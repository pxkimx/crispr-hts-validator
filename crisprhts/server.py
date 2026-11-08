"""Local web server for CRISPR HTS Validator. Standard library HTTP; nothing leaves the computer."""
from __future__ import annotations

import json
import os
import re
import sys
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, quote, urlparse

from . import __version__
from .lifecycle import Lifecycle
from .session import PROJECT, UPLOADS, remove_stale_uploads
from .arrayed_session import ArrayedProject

ARRAYED = ArrayedProject(PROJECT)

WEB = Path(__file__).resolve().parent / "web"
DEFAULT_PORT = 8788
MAX_UPLOAD = 64 * 1024 ** 3        # FASTQ files can be large; they are streamed to disk
MAX_JSON = 20 * 1024 ** 2
TYPES = {".html": "text/html; charset=utf-8", ".js": "text/javascript; charset=utf-8", ".css": "text/css; charset=utf-8",
         ".svg": "image/svg+xml", ".png": "image/png", ".ico": "image/x-icon", ".woff2": "font/woff2"}
LIFE: Lifecycle | None = None


def _clean(obj):
    """JSON-safe: numpy scalars to Python, NaN/inf to None."""
    import math
    import numpy as np
    if isinstance(obj, dict):
        return {str(k): _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _clean(obj.tolist())
    if isinstance(obj, (np.floating, float)):
        f = float(obj)
        return None if math.isnan(f) or math.isinf(f) else f
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


class Handler(BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"
    server_version = "CRISPRHTSValidator"

    def log_message(self, fmt, *args):
        pass

    def _send(self, code: int, body: bytes, ctype: str = "application/json", extra: dict | None = None):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        self.wfile.write(body)

    def _json(self, obj, code=200):
        self._send(code, json.dumps(_clean(obj), allow_nan=False).encode())

    def _guard(self, fn):
        try:
            self._json(fn())
        except (ValueError, FileNotFoundError, RuntimeError, KeyError) as e:
            self._json({"error": str(e).strip("'")}, 400)
        except Exception as e:  # noqa: BLE001
            import traceback
            traceback.print_exc()
            self._json({"error": f"{type(e).__name__}: {e}"}, 500)

    def _local_only(self) -> bool:
        origin = self.headers.get("Origin") or ""
        host = (self.headers.get("Host") or "").split(":")[0]
        if host not in ("localhost", "127.0.0.1"):
            return False
        return not origin or re.match(r"^https?://(localhost|127\.0\.0\.1)(:\d+)?$", origin) is not None

    def _body_json(self) -> dict:
        n = int(self.headers.get("Content-Length") or 0)
        if n > MAX_JSON:
            raise ValueError("Request too large.")
        return json.loads(self.rfile.read(n) or b"{}") if n else {}

    def _upload(self, name: str, sub: str = "") -> Path:
        name = Path(name).name
        if not name or name.startswith(".") or not re.match(r"^[\w .()+,#=-]+$", name):
            raise ValueError("Please use a simple file name (letters, numbers, - _ .).")
        n = int(self.headers.get("Content-Length") or 0)
        if n > MAX_UPLOAD:
            raise ValueError("That file is too large.")
        folder = UPLOADS / sub if sub else UPLOADS
        folder.mkdir(parents=True, exist_ok=True)
        dest = folder / name
        left = n
        with open(dest, "wb") as fh:
            while left > 0:
                chunk = self.rfile.read(min(left, 4 * 1024 * 1024))
                if not chunk:
                    break
                fh.write(chunk)
                left -= len(chunk)
        if left:
            dest.unlink(missing_ok=True)
            raise ValueError("The upload was interrupted.")
        return dest

    # ------------------------------------------------------------------ GET
    def do_GET(self):
        u = urlparse(self.path)
        q = {k: v[0] for k, v in parse_qs(u.query).items()}
        if not self._local_only():
            return self._send(403, b"forbidden", "text/plain")
        P = PROJECT
        A = ARRAYED
        if u.path == "/api/window":
            return self._window()
        if u.path == "/api/health":
            return self._json({"app": "crispr-hts-validator", "version": __version__})
        if u.path == "/api/state":
            return self._guard(P.summary)
        if u.path == "/api/job":
            return self._guard(lambda: P.job(q.get("id", "")))
        if u.path == "/api/qc":
            return self._guard(lambda: P.qc or {"error": "Run the quality checks first."})
        if u.path == "/api/result":
            return self._guard(lambda: P.result(q.get("name", "")))
        if u.path == "/api/gene":
            return self._guard(lambda: P.gene(q.get("name", ""), q.get("gene", "")))
        if u.path == "/api/a/state":
            return self._guard(A.summary)
        if u.path == "/api/a/layout":
            return self._guard(lambda: A.plate_layout(int(q.get("plate", "1"))))
        if u.path == "/api/a/qc":
            return self._guard(lambda: A.qc_view(int(q["plate"]) if q.get("plate") else None))
        if u.path == "/api/a/result":
            return self._guard(lambda: A.result_view(int(q["plate"]) if q.get("plate") else None))
        if u.path == "/api/a/export":
            try:
                body, fn, ctype = A.export(q.get("what", ""), q.get("fmt", "csv"))
            except (ValueError, KeyError) as e:
                return self._json({"error": str(e)}, 400)
            return self._send(200, body, ctype, {"Content-Disposition": f"attachment; filename*=UTF-8''{quote(fn)}"})
        if u.path in ("/api/report", "/api/a/report"):
            try:
                from . import report
                if u.path == "/api/report":
                    if not P.table:
                        raise ValueError("Load a screen first.")
                    body, nm = report.pooled_report(P), P.name or "pooled screen"
                else:
                    if not A.screen:
                        raise ValueError("Load a screen first.")
                    body, nm = report.arrayed_report(A), A.name or "arrayed screen"
            except ValueError as e:
                return self._json({"error": str(e)}, 400)
            fn = re.sub(r"[^\w.+-]+", "_", nm).strip("_") + ".report.pdf"
            disp = "inline" if q.get("view") else "attachment"
            return self._send(200, body, "application/pdf", {"Content-Disposition": f"{disp}; filename*=UTF-8''{quote(fn)}"})
        if u.path == "/api/export":
            try:
                body, fn, ctype = P.export(q.get("what", ""), q.get("name", ""), q.get("fmt", "csv"))
            except (ValueError, KeyError) as e:
                return self._json({"error": str(e)}, 400)
            return self._send(200, body, ctype, {"Content-Disposition": f"attachment; filename*=UTF-8''{quote(fn)}"})
        rel = "index.html" if u.path in ("/", "") else u.path.lstrip("/")
        f = (WEB / rel).resolve()
        if WEB.resolve() not in f.parents or not f.is_file():
            return self._send(404, b"not found", "text/plain")
        self._send(200, f.read_bytes(), TYPES.get(f.suffix, "application/octet-stream"))

    # ------------------------------------------------------------------ POST
    def do_POST(self):
        u = urlparse(self.path)
        q = {k: v[0] for k, v in parse_qs(u.query).items()}
        if not self._local_only():
            return self._send(403, b"forbidden", "text/plain")
        P = PROJECT
        if u.path == "/api/upload":
            kind = q.get("kind", "")

            def up():
                if kind == "counts":
                    P.reset()
                    P.clear_uploads()
                    return P.load_counts(self._upload(q.get("name", "")), q.get("name", ""))
                if kind == "library":
                    if P.table and P.table.fmt != "FASTQ":
                        P.reset()
                    return {"library": P.load_library(self._upload(q.get("name", "")), q.get("name", ""))}
                if kind == "fastq":
                    if P.table and P.table.fmt != "FASTQ":
                        P.reset()
                    return P.add_fastq(self._upload(q.get("name", "")), q.get("name", ""))
                if kind == "positions":
                    return {"positions": P.load_positions(self._upload(q.get("name", "")), q.get("name", ""))}
                if kind == "screeninfo":
                    return {"design": P.import_screen_info(self._upload(q.get("name", "")))}
                raise ValueError("Unknown upload.")
            return self._guard(up)
        A = ARRAYED
        if u.path == "/api/a/upload":
            def aup():
                if A.screen is not None or A.example:
                    A.reset()
                return A.add_file(self._upload(q.get("name", ""), "arrayed"))
            return self._guard(aup)
        try:
            d = self._body_json()
        except (ValueError, json.JSONDecodeError) as e:
            return self._json({"error": str(e)}, 400)
        if u.path == "/api/example":
            return self._guard(lambda: (P.clear_uploads(), P.load_example(d.get("key", "")))[1])
        if u.path == "/api/reset":
            return self._guard(lambda: (P.reset(), P.clear_uploads(), P.summary())[2])
        if u.path == "/api/fastq/remove":
            def rm():
                P.fastq = [f for f in P.fastq if f["name"] != d.get("name")]
                return {"fastq": P.fastq}
            return self._guard(rm)
        if u.path == "/api/count":
            return self._guard(lambda: P.run_count(d))
        if u.path == "/api/design":
            return self._guard(lambda: P.set_design(d))
        if u.path == "/api/thresholds/reset":
            def rt():
                from . import arrayed_stats, prefs, qc
                if d.get("which") == "arrayed":
                    prefs.clear("arrayed_thresholds")
                    return {"thresholds": dict(arrayed_stats.QC_THRESHOLDS)}
                prefs.clear("pooled_thresholds")
                P.thresholds = {}
                return {"thresholds": dict(qc.THRESHOLDS)}
            return self._guard(rt)
        if u.path == "/api/qc":
            return self._guard(lambda: P.run_qc(d.get("thresholds")))
        if u.path == "/api/analyse":
            return self._guard(lambda: P.run_analysis(d))
        if u.path == "/api/controls":
            def ctl():
                m = P.control_mask(d.get("controls", "none"), d.get("control_pattern", ""), d.get("control_list", ""))
                if m is None:
                    return {"guides": 0, "genes": []}
                import numpy as np
                g = sorted(set(np.asarray(P.table.genes, dtype=object)[m].astype(str)))
                return {"guides": int(m.sum()), "genes": g[:30], "n_genes": len(g)}
            return self._guard(ctl)
        if u.path == "/api/a/load":
            return self._guard(A.load)
        if u.path == "/api/a/example":
            return self._guard(lambda: A.load_example(d.get("key", "")))
        if u.path == "/api/a/reset":
            def ar():
                import shutil
                A.reset()
                shutil.rmtree(UPLOADS / "arrayed", ignore_errors=True)
                return A.summary()
            return self._guard(ar)
        if u.path == "/api/a/map":
            return self._guard(lambda: A.set_map(d))
        if u.path == "/api/a/qc":
            return self._guard(lambda: A.run_qc(d.get("thresholds")))
        if u.path == "/api/a/analyse":
            return self._guard(lambda: A.run_analysis(d))
        if u.path == "/api/enrich":
            return self._guard(lambda: P.enrich(d))
        if u.path == "/api/a/enrich":
            return self._guard(lambda: A.enrich(d))
        if u.path == "/api/cancel":
            return self._guard(lambda: (P.cancel(d.get("id", "")), {"ok": True})[1])
        if u.path == "/api/quit":
            self._json({"ok": True})
            if LIFE:
                LIFE.quit_now()
            return
        self._send(404, b"not found", "text/plain")

    def _window(self):
        self.send_response(200)
        self.send_header("Content-Type", "text/event-stream")
        self.send_header("Cache-Control", "no-store")
        self.send_header("Connection", "close")
        self.end_headers()
        if LIFE:
            LIFE.window_opened()
        try:
            while True:
                self.wfile.write(b": hold\n\n")
                self.wfile.flush()
                threading.Event().wait(2.0)
        except (BrokenPipeError, ConnectionResetError, OSError):
            pass
        finally:
            if LIFE:
                LIFE.window_closed()
        self.close_connection = True


def serve(port: int = DEFAULT_PORT, open_browser: bool = True) -> None:
    global LIFE
    remove_stale_uploads()
    PROJECT.clear_uploads()
    httpd = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    httpd.daemon_threads = True

    def quit_():
        sys.stdout.flush()
        import shutil
        shutil.rmtree(UPLOADS, ignore_errors=True)     # uploaded copies are working files only
        os._exit(0)

    LIFE = Lifecycle(on_quit=quit_)
    LIFE.watch()
    url = f"http://localhost:{port}"
    print(f"CRISPR HTS Validator {__version__} at {url}", flush=True)
    if open_browser:
        import webbrowser
        threading.Timer(0.6, lambda: webbrowser.open(url)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        quit_()


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(prog="crispr-hts-validator", description="CRISPR HTS Validator local app")
    ap.add_argument("--port", type=int, default=int(os.environ.get("CH_PORT", DEFAULT_PORT)))
    ap.add_argument("--no-browser", action="store_true")
    a = ap.parse_args()
    serve(a.port, not a.no_browser)


if __name__ == "__main__":
    main()
