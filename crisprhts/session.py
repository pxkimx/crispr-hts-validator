"""The open project: loaded data, sample sheet, QC and results, plus background jobs with progress.

One project is open at a time (it is a desktop app). Jobs run in threads of the server process, so closing the
app ends them too. Uploaded files live in ~/Library/Application Support/CRISPR HTS Validator/uploads and are
removed when another project is loaded.
"""
from __future__ import annotations

import io
import json
import os
import re
import shutil
import threading
import time
import traceback
import uuid
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__, cnv, count, enrichment, pooled_io, pooled_methods, pooled_stats, qc

HOME = Path.home() / "Library" / "Application Support" / "CRISPR HTS Validator"
UPLOADS_ROOT = HOME / "uploads"
UPLOADS = UPLOADS_ROOT / str(os.getpid())      # one folder per running copy, so copies never touch each other's files


def remove_stale_uploads():
    """Delete upload folders left by copies of the app that are no longer running."""
    if not UPLOADS_ROOT.exists():
        return
    for d in UPLOADS_ROOT.iterdir():
        if not d.is_dir() or not d.name.isdigit() or int(d.name) == os.getpid():
            continue
        try:
            os.kill(int(d.name), 0)
        except OSError:
            shutil.rmtree(d, ignore_errors=True)
TOOL_DIR = Path(__file__).resolve().parent.parent
TEST_DATA = TOOL_DIR / "test-data"

EXAMPLES = {
    "tkov3": {"title": "HAP1 TKOv3 (Hart lab, BAGEL2 example)",
              "about": "Genome-wide knockout screen in HAP1 cells: T0, T3, T18 and T18 starved, three replicates each. "
                       "70,877 guides, 18,056 genes.",
              "counts": "hart-bagel/tkov3-example/HAP1-TKOv3-EXAMPLE.txt",
              "screeninfo": "hart-bagel/tkov3-example/HAP1-TKOv3-EXAMPLE-SCREENINFO.txt",
              "positions": "hart-bagel/tkov3-example/TKOv3-guide-positions.txt"},
    "hap1": {"title": "HAP1 TKO (Hart lab, BAGEL example)",
             "about": "Genome-wide knockout screen in HAP1 cells: T0 and T18 A/B/C. About 71,000 guides.",
             "counts": "hart-bagel/reads_hap1.txt"},
    "demo1": {"title": "MAGeCK demo 1 (count table)",
              "about": "Small leukaemia screen excerpt: HL60 and KBM7, initial and final. 999 guides.",
              "counts": "mageck-demo/demo1/sample.txt",
              "comparisons": [{"name": "final vs initial", "treatment": ["HL60.final", "KBM7.final"],
                               "control": ["HL60.initial", "KBM7.initial"]}]},
    "demo3": {"title": "MAGeCK demo 3 (count table)",
              "about": "Larger excerpt of the same leukaemia screen. 9,999 guides.",
              "counts": "mageck-demo/demo3/leukemia.new.csv",
              "comparisons": [{"name": "final vs initial", "treatment": ["HL60.final", "KBM7.final"],
                               "control": ["HL60.initial", "KBM7.initial"]}]},
    "demo2": {"title": "MAGeCK demo 2 (raw FASTQ)",
              "about": "Two small FASTQ files (2,500 reads each) and a 2,550-guide library: read counting from raw data.",
              "library": "mageck-demo/demo2/library.txt",
              "fastq": [["mageck-demo/demo2/test1.fastq", "L1"], ["mageck-demo/demo2/test2.fastq", "CTRL"]],
              "comparisons": [{"name": "L1 vs CTRL", "treatment": ["L1"], "control": ["CTRL"]}],
              "reference": ["CTRL"]},
}

CONTROL_PATTERN = r"(?i)^(non[-_ ]?target(ing)?|nontargeting|ntc|neg(ative)?[-_ ]?control|safe[-_ ]?harbou?r|egfp|gfp|lacz|luciferase|luc|rfp|aavs1|ctrl|control)"


class Job:
    def __init__(self, kind: str):
        self.id = uuid.uuid4().hex[:10]
        self.kind = kind
        self.progress = 0.0
        self.message = "Starting…"
        self.state = "running"
        self.error = None
        self.started = time.time()
        self.cancelled = False

    def public(self) -> dict:
        return {"id": self.id, "kind": self.kind, "progress": round(self.progress, 3), "message": self.message,
                "state": self.state, "error": self.error, "elapsed": round(time.time() - self.started, 1)}


class Project:
    def __init__(self):
        self.lock = threading.RLock()
        self.reset()

    def reset(self):
        self.table: pooled_io.CountTable | None = None
        self.library: pooled_io.Library | None = None
        self.fastq: list[dict] = []
        self.count_info: list[dict] = []
        self.design: dict = {"samples": [], "comparisons": []}
        self.qc: dict | None = None
        self.results: dict[str, dict] = {}
        self.settings: dict = {}
        self.name = ""
        self.jobs: dict[str, Job] = {}
        self.thresholds: dict = {}
        self.positions: pd.DataFrame | None = None
        self.positions_name = ""

    # ------------------------------------------------------------------ loading
    def clear_uploads(self):
        shutil.rmtree(UPLOADS, ignore_errors=True)
        UPLOADS.mkdir(parents=True, exist_ok=True)

    def set_table(self, t: pooled_io.CountTable, name: str):
        with self.lock:
            self.table = t
            self.name = name
            self.design = pooled_io.guess_design(t.samples)
            self.qc = None
            self.results = {}

    def load_counts(self, path: Path, name: str) -> dict:
        t = pooled_io.read_count_table(path, name)
        self.fastq, self.count_info = [], []
        self.set_table(t, name)
        return self.summary()

    def load_library(self, path: Path, name: str) -> dict:
        lib = pooled_io.read_library(path, name)
        with self.lock:
            self.library = lib
        return {"guides": len(lib.ids), "genes": len(set(lib.genes)), "notes": lib.notes, "name": name,
                "positions": lib.pos is not None}

    def load_positions(self, path: Path, name: str) -> dict:
        pos = cnv.read_positions(path)
        with self.lock:
            self.positions = pos
            self.positions_name = name
        return self.positions_info()

    def positions_info(self) -> dict | None:
        if self.positions is None:
            return None
        matched = int(pd.Index(self.table.guides).isin(self.positions["guide"]).sum()) if self.table else None
        return {"name": self.positions_name, "guides": int(len(self.positions)), "matched": matched,
                "chromosomes": int(self.positions["chrom"].nunique())}

    def add_fastq(self, path: Path, name: str) -> dict:
        sample = re.sub(r"(\.fastq|\.fq)(\.gz)?$", "", name, flags=re.I)
        with self.lock:
            self.fastq = [f for f in self.fastq if f["name"] != name]
            self.fastq.append({"path": str(path), "name": name, "sample": sample,
                               "size": path.stat().st_size})
        return {"fastq": self.fastq}

    def load_example(self, key: str) -> dict:
        ex = EXAMPLES.get(key)
        if not ex:
            raise ValueError("Unknown example.")
        self.reset()
        if "counts" in ex:
            p = TEST_DATA / ex["counts"]
            if not p.exists():
                raise FileNotFoundError("The example data is not in this copy of the tool (test-data folder missing).")
            self.load_counts(p, Path(ex["counts"]).name)
            if "screeninfo" in ex:
                self.design["comparisons"] = pooled_io.read_screen_info(TEST_DATA / ex["screeninfo"])
                if key == "tkov3":   # a treated-vs-untreated pair for drugZ
                    self.design["comparisons"].append({"name": "Starved vs T18 (drugZ)",
                                                       "treatment": ["T18_A_Starved", "T18_B_Starved", "T18_C_Starved"],
                                                       "control": ["T18_A", "T18_B", "T18_C"]})
            if "positions" in ex:
                self.load_positions(TEST_DATA / ex["positions"], Path(ex["positions"]).name)
        else:
            self.load_library(TEST_DATA / ex["library"], Path(ex["library"]).name)
            for f, s in ex["fastq"]:
                self.fastq.append({"path": str(TEST_DATA / f), "name": Path(f).name, "sample": s,
                                   "size": (TEST_DATA / f).stat().st_size})
        if "comparisons" in ex:
            self._pending_comparisons = ex["comparisons"]
            self._pending_reference = ex.get("reference")
            if self.table:
                self.apply_pending()
        self.name = ex["title"]
        return self.summary()

    def apply_pending(self):
        comps = getattr(self, "_pending_comparisons", None)
        if comps and self.table:
            self.design["comparisons"] = [dict(c) for c in comps]
            refs = set(sum((c["control"] for c in comps), []))
            for r in self.design["samples"]:
                r["reference"] = r["sample"] in refs
        self._pending_comparisons = None

    # ------------------------------------------------------------------ design
    def set_design(self, d: dict) -> dict:
        if not self.table:
            raise ValueError("Load data first.")
        samples = self.table.samples
        rows = []
        for r in d.get("samples", []):
            if r.get("sample") in samples:
                rows.append({"sample": r["sample"], "condition": str(r.get("condition") or r["sample"]).strip(),
                             "reference": bool(r.get("reference")), "replicate": int(r.get("replicate") or 1)})
        comps = []
        for c in d.get("comparisons", []):
            t = [s for s in c.get("treatment", []) if s in samples]
            k = [s for s in c.get("control", []) if s in samples]
            if not t or not k:
                raise ValueError(f"Comparison “{c.get('name', '')}” needs at least one treatment and one control sample.")
            if set(t) & set(k):
                raise ValueError(f"Comparison “{c.get('name', '')}” uses the same sample as treatment and control.")
            comps.append({"name": str(c.get("name") or f"{'+'.join(t)} vs {'+'.join(k)}").strip(), "treatment": t, "control": k})
        names = [c["name"] for c in comps]
        if len(set(names)) != len(names):
            raise ValueError("Two comparisons have the same name.")
        with self.lock:
            self.design = {"samples": rows, "comparisons": comps}
            self.qc = None
        return self.design

    def import_screen_info(self, path: Path) -> dict:
        comps = pooled_io.read_screen_info(path)
        if self.table:
            missing = sorted({s for c in comps for s in c["treatment"] + c["control"]} - set(self.table.samples))
            if missing:
                raise ValueError("These samples in the screen-info file are not in the count table: " + ", ".join(missing))
            refs = set(sum((c["control"] for c in comps), []))
            d = {"samples": [dict(r, reference=r["sample"] in refs) for r in self.design["samples"]], "comparisons": comps}
            return self.set_design(d)
        raise ValueError("Load the count table first, then the screen-info file.")

    # ------------------------------------------------------------------ jobs
    def start(self, kind: str, fn) -> dict:
        with self.lock:
            busy = [j for j in self.jobs.values() if j.state == "running"]
            if busy:
                raise RuntimeError("Another step is still running; wait for it to finish.")
            job = Job(kind)
            self.jobs[job.id] = job

        def run():
            try:
                fn(job)
                job.progress = 1.0
                job.state = "done"
                job.message = "Done."
            except Exception as e:  # noqa: BLE001
                traceback.print_exc()
                job.state = "error"
                job.error = str(e) or type(e).__name__
        threading.Thread(target=run, daemon=True, name=f"job-{kind}").start()
        return job.public()

    def job(self, jid: str) -> dict:
        j = self.jobs.get(jid)
        if not j:
            raise ValueError("No such job.")
        return j.public()

    def cancel(self, jid: str):
        j = self.jobs.get(jid)
        if j:
            j.cancelled = True

    # ---- counting
    def run_count(self, opts: dict) -> dict:
        if not self.library:
            raise ValueError("Load the guide library first.")
        if not self.fastq:
            raise ValueError("Add at least one FASTQ file.")
        samples = opts.get("samples") or {}
        for f in self.fastq:
            if f["name"] in samples and str(samples[f["name"]]).strip():
                f["sample"] = str(samples[f["name"]]).strip()
        names = [f["sample"] for f in self.fastq]
        if len(set(names)) != len(names):
            raise ValueError("Two FASTQ files have the same sample name.")
        lib = self.library

        def work(job: Job):
            res = []
            sizes = [max(1, f["size"]) for f in self.fastq]
            done = 0
            for f, sz in zip(self.fastq, sizes):
                job.message = f"Counting {f['name']}…"
                est = sz / (60 if f["name"].lower().endswith(".gz") is False else 25)  # rough reads estimate

                def prog(n, base=done, sz=sz, est=est):
                    job.progress = min(0.99, (base + min(n / max(est, 1), 1) * sz) / sum(sizes))
                    job.message = f"Counting {f['name']}: {n:,} reads"
                r = count.count_fastq(f["path"], lib, f["sample"], mode=opts.get("mode", "auto"),
                                      trim=opts.get("trim"), mismatch=bool(opts.get("mismatch")),
                                      flank=(opts.get("flank") or count.SCAFFOLD).upper(), progress=prog,
                                      cancel=lambda: job.cancelled)
                res.append(r)
                done += sz
            counts = np.stack([r.counts for r in res], axis=1).astype(float)
            t = pooled_io.CountTable(list(lib.ids), list(lib.genes), [r.sample for r in res], counts, "FASTQ",
                                     ", ".join(f["name"] for f in self.fastq),
                                     notes=[n for r in res for n in (f"{r.sample}: {x}" for x in r.notes)],
                                     reads={r.sample: r.reads for r in res}, mapped={r.sample: r.mapped for r in res})
            info = [{"sample": r.sample, "file": Path(r.file).name, "reads": r.reads, "mapped": r.mapped,
                     "fraction": r.mapped_fraction, "trims": [x + 1 for x in r.trims], "reverse": r.reverse,
                     "mode": r.mode, "mismatch_reads": r.mismatch_reads, "zero": int((r.counts == 0).sum()),
                     "gini": qc.gini(r.counts)} for r in res]
            pend, pref = getattr(self, "_pending_comparisons", None), getattr(self, "_pending_reference", None)
            self.set_table(t, self.name or "FASTQ counts")
            self.count_info = info
            self.settings["count"] = {k: opts.get(k) for k in ("mode", "trim", "mismatch", "flank")}
            if pend:
                self._pending_comparisons = pend
                self.apply_pending()
            if pref:
                for r in self.design["samples"]:
                    r["reference"] = r["sample"] in pref
        return self.start("count", work)

    # ---- QC
    def run_qc(self, thresholds: dict | None = None) -> dict:
        if not self.table:
            raise ValueError("Load data first.")
        th = {k: float(v) for k, v in (thresholds or {}).items() if k in qc.THRESHOLDS}
        self.thresholds.update(th)

        def work(job: Job):
            job.message = "Checking samples…"
            t = self.table
            r = qc.run_qc(t.counts, t.guides, t.genes, t.samples, self.design, t.reads, t.mapped,
                          self.thresholds)
            r["count_info"] = self.count_info
            r["threshold_text"] = qc.threshold_text(self.thresholds)
            self.qc = r
        return self.start("qc", work)

    # ---- analysis
    def control_mask(self, how: str, pattern: str = "", genes_text: str = "") -> np.ndarray | None:
        if not self.table:
            return None
        g = np.asarray(self.table.genes, dtype=object).astype(str)
        if how == "pattern":
            rx = re.compile(pattern or CONTROL_PATTERN)
            return np.array([bool(rx.search(x)) for x in g])
        if how == "list":
            names = {x.strip() for x in re.split(r"[\s,;]+", genes_text or "") if x.strip()}
            ids = np.asarray(self.table.guides, dtype=object).astype(str)
            return np.array([(a in names) or (b in names) for a, b in zip(g, ids)])
        return None

    def run_analysis(self, opts: dict) -> dict:
        if not self.table:
            raise ValueError("Load data first.")
        comps = [c for c in self.design["comparisons"] if not opts.get("comparisons") or c["name"] in opts["comparisons"]]
        if not comps:
            raise ValueError("Set up at least one comparison on the Samples page.")
        norm = opts.get("norm", "median")
        if norm not in pooled_stats.NORM_METHODS:
            raise ValueError("Unknown normalisation.")
        cm = self.control_mask(opts.get("controls", "none"), opts.get("control_pattern", ""), opts.get("control_list", ""))
        if norm == "control" and (cm is None or not cm.any()):
            raise ValueError("Control-guide normalisation needs control guides: none matched.")
        use_null = bool(opts.get("control_null")) and cm is not None and cm.any()
        fdr_thr = float(opts.get("gene_test_fdr_threshold", 0.25))
        remove_zero = opts.get("remove_zero", "both")
        gene_lfc = opts.get("gene_lfc", "median")
        method = opts.get("method", "mageck")
        if method not in pooled_methods.METHODS:
            raise ValueError("Unknown method.")
        cn = bool(opts.get("cn_correct"))
        if cn and self.positions is None:
            raise ValueError("Copy-number correction needs guide positions: load a position file on the Data page.")
        if cn and method == "drugz":
            raise ValueError("Copy-number correction is for comparisons against plasmid/T0; drugZ compares two arms of "
                             "the same cells, where copy number cancels out.")
        paired = opts.get("drugz_paired", True) is not False
        remove = [g for g in re.split(r"[\s,;]+", opts.get("drugz_remove", "")) if g]
        t = self.table

        def work(job: Job):
            out = {}
            for i, c in enumerate(comps):
                job.message = f"{method}: {c['name']} ({i + 1} of {len(comps)})…"
                cols = [t.samples.index(s) for s in c["control"] + c["treatment"]]
                X = t.counts[:, cols]

                def prog(x, i=i):
                    job.progress = min(0.99, (i + x) / len(comps))
                settings = {"method": method, "cn_correct": cn}
                if method == "mageck":
                    cninfo = None
                    if cn:
                        job.message = f"Copy-number correction: {c['name']}…"
                        sh = pooled_methods.cn_shift(X, t.guides, t.genes, len(c["control"]), self.positions)
                        X = pooled_methods.corrected_counts(X, len(c["control"]), sh["shift"])
                        cninfo = {"segments": sh["segments"], "n_positioned": sh["n_positioned"],
                                  "corrected": int(np.sum(np.abs(sh["shift"]) > 0.25))}
                        job.message = f"{method}: {c['name']} ({i + 1} of {len(comps)})…"
                    r = pooled_stats.mageck_test(X, t.genes, t.guides, len(c["control"]), norm=norm, control_mask=cm,
                                                 remove_zero=remove_zero, gene_test_fdr_threshold=fdr_thr,
                                                 gene_lfc=gene_lfc, use_control_null=use_null, progress=prog)
                    r["method"] = "mageck"
                    r["cn"] = cninfo
                    settings.update({"norm": norm, "controls": int(cm.sum()) if cm is not None else 0,
                                     "control_null": use_null, "remove_zero": remove_zero,
                                     "gene_test_fdr_threshold": fdr_thr, "gene_lfc": gene_lfc})
                elif method == "bagel":
                    r = pooled_methods.run_bagel(t, c, self.positions, cn, progress=prog)
                    settings.update({"norm": "10 M reads per sample, pseudocount 5", "folds": 10, "seed": 1,
                                     "cutoff": pooled_methods.BF_CUTOFF})
                else:
                    r = pooled_methods.run_drugz(t, c, paired=paired, remove_genes=remove, progress=prog)
                    settings.update({"norm": "10 M reads per sample, pseudocount 5", "paired": paired,
                                     "removed_genes": remove})
                out[c["name"]] = {"comparison": c, "result": r, "when": time.strftime("%Y-%m-%d %H:%M"),
                                  "settings": settings}
            with self.lock:
                self.results.update(out)
            self.settings["analysis"] = {"norm": norm, "controls": opts.get("controls", "none")}
        return self.start("analysis", work)

    # ------------------------------------------------------------------ views
    def summary(self) -> dict:
        t = self.table
        out = {"version": __version__, "name": self.name, "examples": {k: {"title": v["title"], "about": v["about"]}
                                                                     for k, v in EXAMPLES.items()},
               "examples_available": TEST_DATA.exists(),
               "library": None, "fastq": self.fastq, "count_info": self.count_info, "design": self.design,
               "has_qc": self.qc is not None, "qc_overall": self.qc["overall"] if self.qc else None,
               "results": [{"name": k, "when": v["when"], "settings": v["settings"],
                            "method": v["result"].get("method", "mageck"),
                            "hits_neg": pooled_methods.hit_counts(v)[0], "hits_pos": pooled_methods.hit_counts(v)[1]}
                           for k, v in self.results.items()],
               "methods": pooled_methods.METHODS, "positions": self.positions_info() if self.table else None,
               "enrichment": {k: enrichment.COLLECTIONS[k] for k in enrichment.COLLECTIONS},
               "norm_methods": pooled_stats.NORM_METHODS, "control_pattern": CONTROL_PATTERN,
               "jobs": [j.public() for j in self.jobs.values() if j.state == "running"]}
        if self.library:
            out["library"] = {"guides": len(self.library.ids), "genes": len(set(self.library.genes)),
                              "notes": self.library.notes, "positions": self.library.pos is not None}
        if t:
            out["table"] = {"format": t.fmt, "source": t.source, "guides": len(t.guides), "genes": len(set(t.genes)),
                            "samples": t.samples, "notes": t.notes,
                            "totals": [float(x) for x in t.counts.sum(axis=0)],
                            "preview": [[t.guides[i], t.genes[i]] + [float(x) for x in t.counts[i]] for i in range(min(8, len(t.guides)))]}
        return out

    def result(self, name: str) -> dict:
        r = self.results.get(name)
        if not r:
            raise ValueError("No result with that name.")
        res = r["result"]
        if res.get("method", "mageck") != "mageck":
            return pooled_methods.view(name, r)
        genes = res["genes"]
        return {"name": name, "method": "mageck", "cn": pooled_methods._cn_view(res.get("cn")), "comparison": r["comparison"], "settings": r["settings"], "when": r["when"],
                "alpha_neg": res["alpha_neg"], "alpha_pos": res["alpha_pos"],
                "size_factors": dict(zip(r["comparison"]["control"] + r["comparison"]["treatment"],
                                         [float(x) for x in res["guide"]["factors"]])),
                "norm_warnings": res["guide"]["norm_warnings"],
                "model": [float(x) for x in res["guide"]["model"]],
                "columns": ["gene", "guides", "lfc", "neg_score", "neg_p", "neg_fdr", "neg_rank", "neg_good",
                            "pos_score", "pos_p", "pos_fdr", "pos_rank", "pos_good"],
                "rows": [[g["gene"], g["guides"], round(g["lfc"], 5), g["neg_score"], g["neg_p"], g["neg_fdr"],
                          g["neg_rank"], g["neg_good"], g["pos_score"], g["pos_p"], g["pos_fdr"], g["pos_rank"],
                          g["pos_good"]] for g in genes]}

    def gene(self, name: str, gene: str) -> dict:
        r = self.results.get(name)
        if not r:
            raise ValueError("No result with that name.")
        t = self.table
        if r["result"].get("method", "mageck") != "mageck":
            return pooled_methods.gene_view(t, r, gene)
        g = r["result"]["guide"]
        c = r["comparison"]
        idx = [i for i, x in enumerate(t.genes) if x == gene]
        if not idx:
            raise ValueError("Gene not found.")
        cols = c["control"] + c["treatment"]
        rows = []
        for i in idx:
            rows.append({"guide": t.guides[i], "raw": [float(t.counts[i, t.samples.index(s)]) for s in cols],
                         "norm": [round(float(x), 2) for x in g["normalised"][i]],
                         "control_mean": float(g["control_mean"][i]), "treat_mean": float(g["treat_mean"][i]),
                         "lfc": float(g["lfc"][i]), "p_low": float(g["p_low"][i]), "p_high": float(g["p_high"][i]),
                         "fdr": float(g["fdr"][i]), "score": float(g["score"][i]), "valid": bool(g["valid"][i])})
        row = next((x for x in r["result"]["genes"] if x["gene"] == gene), None)
        return {"gene": gene, "method": "mageck", "samples": cols, "n_control": len(c["control"]), "guides": rows, "summary": row}

    # ------------------------------------------------------------------ enrichment
    def enrich(self, d: dict) -> dict:
        r = self.results.get(d.get("name", ""))
        if not r:
            raise ValueError("Run the analysis first.")
        direction = d.get("direction", "neg")
        cutoff = float(d.get("cutoff", 0.05))
        hits, scores, better = pooled_methods.hits_of(r, direction, cutoff)
        return run_enrichment(d, hits, scores, list(scores), better, self.table.genes)

    # ------------------------------------------------------------------ exports
    def guide_frame(self, name: str) -> pd.DataFrame:
        r = self.results[name]
        if r["result"].get("method", "mageck") != "mageck":
            return pooled_methods.guide_frame(self.table, r)
        g = r["result"]["guide"]
        t = self.table
        c = r["comparison"]
        cols = c["control"] + c["treatment"]
        nc = len(c["control"])
        df = pd.DataFrame({"sgrna": t.guides, "Gene": t.genes,
                           "control_count": ["/".join(f"{x:.5g}" for x in row[:nc]) for row in g["normalised"]],
                           "treatment_count": ["/".join(f"{x:.5g}" for x in row[nc:]) for row in g["normalised"]],
                           "control_mean": g["control_mean"], "treat_mean": g["treat_mean"], "LFC": g["lfc"],
                           "control_var": g["control_var"], "adj_var": g["adj_var"], "score": np.abs(g["score"]),
                           "p.low": g["p_low"], "p.high": g["p_high"], "p.twosided": g["p_two"], "FDR": g["fdr"],
                           "high_in_treatment": g["treat_mean"] > g["control_mean"], "used_in_gene_test": g["valid"]})
        return df.iloc[np.argsort(-np.abs(g["score"]), kind="mergesort")]

    def gene_frame(self, name: str) -> pd.DataFrame:
        if self.results[name]["result"].get("method", "mageck") != "mageck":
            return pooled_methods.gene_frame(self.results[name])
        rows = self.results[name]["result"]["genes"]
        return pd.DataFrame({"id": [r["gene"] for r in rows], "num": [r["guides"] for r in rows],
                             "neg|score": [r["neg_score"] for r in rows], "neg|p-value": [r["neg_p"] for r in rows],
                             "neg|fdr": [r["neg_fdr"] for r in rows], "neg|rank": [r["neg_rank"] for r in rows],
                             "neg|goodsgrna": [r["neg_good"] for r in rows], "neg|lfc": [r["lfc"] for r in rows],
                             "pos|score": [r["pos_score"] for r in rows], "pos|p-value": [r["pos_p"] for r in rows],
                             "pos|fdr": [r["pos_fdr"] for r in rows], "pos|rank": [r["pos_rank"] for r in rows],
                             "pos|goodsgrna": [r["pos_good"] for r in rows], "pos|lfc": [r["lfc"] for r in rows]})

    def counts_frame(self, normalised: bool = False) -> pd.DataFrame:
        t = self.table
        X = qc.median_ratio_normalise(t.counts) if normalised else t.counts
        df = pd.DataFrame(X, columns=t.samples)
        df.insert(0, "Gene", t.genes)
        df.insert(0, "sgRNA", t.guides)
        return df

    def qc_frame(self) -> pd.DataFrame:
        if not self.qc:
            raise ValueError("Run the quality checks first.")
        rows = []
        for c in self.qc["checks"]:
            rows.append({"check": c["name"], "item": "overall", "value": "", "status": c["status"], "why": c["why"]})
            for p in c.get("per_sample", []):
                rows.append({"check": c["name"], "item": p["sample"], "value": p["value"], "status": p["status"], "why": ""})
            for p in c.get("pairs", []):
                rows.append({"check": c["name"], "item": f"{p['a']} vs {p['b']}", "value": p["value"], "status": p["status"], "why": ""})
            for p in c.get("per_condition", []):
                rows.append({"check": c["name"], "item": p["condition"], "value": p["value"], "status": p["status"], "why": ""})
        return pd.DataFrame(rows)

    def export(self, what: str, name: str = "", fmt: str = "csv") -> tuple[bytes, str, str]:
        safe = re.sub(r"[^\w.+-]+", "_", name or self.name or "screen").strip("_") or "screen"
        if what == "genes":
            df, fn = self.gene_frame(name), f"{safe}.gene_summary"
        elif what == "guides":
            df, fn = self.guide_frame(name), f"{safe}.sgrna_summary"
        elif what == "counts":
            df, fn = self.counts_frame(False), f"{safe}.count"
        elif what == "normcounts":
            df, fn = self.counts_frame(True), f"{safe}.count_normalized"
        elif what == "qc":
            df, fn = self.qc_frame(), f"{safe}.qc"
        elif what == "all":
            return self.export_workbook()
        else:
            raise ValueError("Unknown export.")
        if fmt == "xlsx":
            buf = io.BytesIO()
            df.to_excel(buf, index=False, sheet_name=what[:30])
            return buf.getvalue(), fn + ".xlsx", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        sep = "\t" if fmt == "txt" else ","
        return pooled_io.to_csv_bytes(df, sep), fn + (".txt" if fmt == "txt" else ".csv"), "text/plain; charset=utf-8"

    def export_workbook(self) -> tuple[bytes, str, str]:
        buf = io.BytesIO()
        with pd.ExcelWriter(buf) as xw:
            about = pd.DataFrame({"field": ["tool", "version", "project", "exported", "normalisation"],
                                  "value": ["CRISPR HTS Validator", __version__, self.name, time.strftime("%Y-%m-%d %H:%M"),
                                            json.dumps(self.settings.get("analysis", {}))]})
            about.to_excel(xw, index=False, sheet_name="about")
            if self.qc:
                self.qc_frame().to_excel(xw, index=False, sheet_name="qc")
            for k in self.results:
                self.gene_frame(k).to_excel(xw, index=False, sheet_name=re.sub(r"[\[\]:*?/\\]", "_", "genes " + k)[:31])
            self.counts_frame(False).to_excel(xw, index=False, sheet_name="counts")
        safe = re.sub(r"[^\w.+-]+", "_", self.name or "screen").strip("_") or "screen"
        return buf.getvalue(), f"{safe}.results.xlsx", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"


def run_enrichment(d: dict, hits, scores: dict, universe, better: bool, all_genes) -> dict:
    """Shared by pooled and arrayed results. d: collection (GO_BP, GO_MF, GO_CC, Reactome, custom), test (ora, ranked),
    species (auto, human, mouse), gmt (text of a user .gmt when collection == custom)."""
    coll = d.get("collection", "Reactome")
    test = d.get("test", "ora")
    species = d.get("species", "auto")
    if species == "auto":
        species = enrichment.guess_species(list(all_genes))
    if coll == "custom":
        if not d.get("gmt"):
            raise ValueError("Choose a .gmt file of gene sets.")
        sets = enrichment.read_gmt(d["gmt"])
        label = d.get("gmt_name", "your gene sets")
    else:
        sets = enrichment.load(coll, species)
        label = f"{enrichment.COLLECTIONS[coll]} ({species})"
    lo, hi = int(d.get("min_size", 10)), int(d.get("max_size", 500))
    if test == "ranked":
        rows = enrichment.ranked(scores, sets, better, lo, hi)
    else:
        rows = enrichment.ora(hits, universe, sets, lo, hi)
    return {"collection": label, "test": test, "species": species, "n_hits": len(hits), "n_universe": len(set(universe)),
            "n_sets": len(rows), "rows": rows[:300]}


PROJECT = Project()
