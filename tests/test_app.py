"""The app as a whole: API round trip on example data, and quitting when the window closes."""
import json
import os
import socket
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def free_port():
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    p = s.getsockname()[1]
    s.close()
    return p


def start(port, env=None):
    e = dict(os.environ, **(env or {}))
    p = subprocess.Popen([sys.executable, "-m", "crisprhts.server", "--port", str(port), "--no-browser"], cwd=ROOT,
                         env=e, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for _ in range(100):
        try:
            urllib.request.urlopen(f"http://127.0.0.1:{port}/api/health", timeout=0.5)
            return p
        except OSError:
            time.sleep(0.1)
    p.kill()
    raise RuntimeError("server did not start")


def call(port, path, body=None):
    req = urllib.request.Request(f"http://localhost:{port}{path}", data=None if body is None else json.dumps(body).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read())


def wait_job(port, jid):
    for _ in range(600):
        j = call(port, f"/api/job?id={jid}")
        if j["state"] != "running":
            return j
        time.sleep(0.2)
    raise TimeoutError


@pytest.mark.skipif(not (ROOT / "test-data/mageck-demo").exists(), reason="test data not present")
def test_api_round_trip():
    port = free_port()
    p = start(port)
    try:
        s = call(port, "/api/example", {"key": "demo1"})
        assert s["table"]["guides"] == 999 and s["design"]["comparisons"][0]["name"] == "final vs initial"
        j = wait_job(port, call(port, "/api/qc", {})["id"])
        assert j["state"] == "done"
        qc = call(port, "/api/qc")
        assert {c["id"] for c in qc["checks"]} >= {"reads", "zeros", "gini", "coverage", "replicates", "pca", "essentials"}
        j = wait_job(port, call(port, "/api/analyse", {"norm": "median"})["id"])
        assert j["state"] == "done", j
        r = call(port, "/api/result?name=final%20vs%20initial")
        assert len(r["rows"]) == 100
        g = call(port, "/api/gene?name=final%20vs%20initial&gene=" + r["rows"][0][0])
        assert g["guides"]
        with urllib.request.urlopen(f"http://localhost:{port}/api/export?what=genes&fmt=txt&name=final%20vs%20initial") as f:
            head = f.read().decode().splitlines()[0]
        assert head.startswith("id\tnum\tneg|score\tneg|p-value")
        with urllib.request.urlopen(f"http://localhost:{port}/api/export?what=all&fmt=xlsx") as f:
            assert f.read(2) == b"PK"
        # another site cannot drive the app
        req = urllib.request.Request(f"http://localhost:{port}/api/state", headers={"Origin": "https://evil.example"})
        with pytest.raises(urllib.error.HTTPError):
            urllib.request.urlopen(req)
    finally:
        p.kill()


def test_quits_when_window_closes_and_survives_reload():
    port = free_port()
    p = start(port)
    try:
        def open_window():
            s = socket.create_connection(("127.0.0.1", port))
            s.sendall(f"GET /api/window HTTP/1.1\r\nHost: localhost:{port}\r\n\r\n".encode())
            s.recv(100)
            return s
        w = open_window()
        time.sleep(1)
        w.close()
        w = open_window()            # a reload reconnects within the grace period
        time.sleep(7)
        assert p.poll() is None
        w.close()
        p.wait(timeout=12)           # window gone: quits on its own
        assert p.returncode == 0
    finally:
        if p.poll() is None:
            p.kill()


def test_quits_when_launcher_is_gone():
    port = free_port()
    parent = subprocess.Popen(["sleep", "30"])
    p = start(port, {"CH_PARENT_PID": str(parent.pid)})
    try:
        parent.kill()
        parent.wait()
        p.wait(timeout=5)
    finally:
        if p.poll() is None:
            p.kill()


@pytest.mark.skipif(not (ROOT / "test-data/cellhts2/KcViab").exists(), reason="test data not present")
def test_arrayed_round_trip():
    port = free_port()
    p = start(port)
    try:
        s = call(port, "/api/a/example", {"key": "kcviab"})
        assert s["screen"]["plates"] == 57 and s["roles"]["neg"] == "neg"
        lay = call(port, "/api/a/layout?plate=1")
        assert len(lay["wells"]) == 384
        s = call(port, "/api/a/map", {"mark": [{"wells": "P24", "plates": "1", "content": "empty"}], "direction": "lower"})
        assert s["roles"]["empty"] == "empty"
        assert wait_job(port, call(port, "/api/a/qc", {})["id"])["state"] == "done"
        q = call(port, "/api/a/qc?plate=2")
        assert q["plate"] == 2 and {c["id"] for c in q["checks"]} == {"plates", "controls", "edge", "rowcol", "replicates"}
        assert wait_job(port, call(port, "/api/a/analyse", {"method": "bscore"})["id"])["state"] == "done"
        r = call(port, "/api/a/result?plate=3")
        assert r["plate"] == 3 and r["n_hits"] > 0 and r["comparison"]["published_hits"] == 417
        with urllib.request.urlopen(f"http://localhost:{port}/api/a/export?what=hits&fmt=csv") as f:
            assert f.read().decode().splitlines()[0].startswith("plate,well")
    finally:
        p.kill()


@pytest.mark.skipif(not (ROOT / "test-data/hart-bagel").exists(), reason="test data not present")
def test_methods03_round_trip():
    port = free_port()
    p = start(port)
    try:
        s = call(port, "/api/example", {"key": "tkov3"})
        assert s["positions"]["matched"] > 70000 and "bagel" in s["methods"]
        j = wait_job(port, call(port, "/api/analyse", {"method": "bagel", "comparisons": ["BATCH1_HAP1_T18"], "cn_correct": True})["id"])
        assert j["state"] == "done", j
        r = call(port, "/api/result?name=BATCH1_HAP1_T18")
        assert r["method"] == "bagel" and r["summary"]["calls"] > 1500 and r["cn"]["segments"] > 20
        e = call(port, "/api/enrich", {"name": "BATCH1_HAP1_T18", "collection": "Reactome", "test": "ora"})
        assert e["rows"][0]["fdr"] < 1e-50
        j = wait_job(port, call(port, "/api/analyse", {"method": "drugz", "comparisons": ["Starved vs T18 (drugZ)"]})["id"])
        assert j["state"] == "done", j
        g = call(port, "/api/gene?name=Starved%20vs%20T18%20(drugZ)&gene=TSC2")
        assert g["summary"]["normZ"] < -5
        with urllib.request.urlopen(f"http://localhost:{port}/api/report") as f:
            assert f.read(5) == b"%PDF-"
    finally:
        p.kill()
