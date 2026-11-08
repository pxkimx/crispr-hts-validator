"""The open arrayed screen: loading, plate map edits, QC, analysis and exports (one screen at a time)."""
from __future__ import annotations

import io
import re
import threading
import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__, arrayed_io, arrayed_stats, prefs
from .session import TEST_DATA, UPLOADS, Project

EXAMPLES = {
    "kcviab": {"title": "KcViab — genome-wide RNAi viability screen (cellHTS2)",
               "about": "Drosophila Kc167 cells, 57 × 384-well plates in duplicate, ATP luminescence (Boutros et al., Science 2004). "
                        "Comes with the published 2003 scores for comparison.",
               "dir": "cellhts2/KcViab", "direction": "lower"},
    "dual": {"title": "Dual-channel screen (cellHTS2)",
             "about": "3 plates in duplicate with two reporters: growth/viability and pathway activity. Named control wells.",
             "dir": "cellhts2/DualChannelScreen", "direction": "lower"},
    "envision": {"title": "EnVision plate-reader export (cellHTS2)",
                 "about": "One raw 384-well luminescence export with instrument metadata: shows reading plate-shaped grids.",
                 "dir": "cellhts2/EnVisionExample", "direction": "lower"},
}
ROLES = ("sample", "neg", "pos", "other", "empty")


class ArrayedProject:
    def __init__(self, jobs_host: Project):
        self.host = jobs_host           # reuse the pooled project's job runner (one job at a time overall)
        self.lock = threading.RLock()
        self.reset()

    def reset(self):
        self.screen: arrayed_io.ArrayedScreen | None = None
        self.files: dict[str, Path] = {}
        self.roles: dict[str, str] = {}
        self.direction = "lower"
        self.channel = 1
        self.qc = None
        self.result = None
        self.name = ""
        self.example = None

    # ------------------------------------------------------------------ loading
    def add_file(self, path: Path) -> dict:
        with self.lock:
            self.files[path.name] = path
        return {"files": sorted(self.files)}

    def load(self) -> dict:
        if not self.files:
            raise ValueError("Add the screen's files first.")
        scr = arrayed_io.detect_and_load(self.files, self.name or "")
        self._set(scr)
        return self.summary()

    def load_example(self, key: str) -> dict:
        ex = EXAMPLES.get(key)
        if not ex:
            raise ValueError("Unknown example.")
        d = TEST_DATA / ex["dir"]
        if not d.exists():
            raise FileNotFoundError("The example data is not in this copy of the tool (test-data folder missing).")
        self.reset()
        self.files = {p.name: p for p in d.iterdir() if p.is_file() and p.name != "Analysis2003.txt"
                      and not p.name.startswith("old-")}
        self.name = ex["title"].split(" — ")[0].split(" (")[0]
        self.example = key
        scr = arrayed_io.detect_and_load(self.files, self.name)
        self._set(scr)
        self.direction = ex["direction"]
        return self.summary()

    def _set(self, scr):
        with self.lock:
            self.screen = scr
            self.name = self.name or scr.name or "Arrayed screen"
            kinds = scr.wells.drop_duplicates(["plate", "well"])["content"].value_counts()
            self.roles = {k: (k if k in ROLES else "other") for k in kinds.index}
            self.qc = None
            self.result = None

    # ------------------------------------------------------------------ plate map
    def set_map(self, d: dict) -> dict:
        if not self.screen:
            raise ValueError("Load a screen first.")
        w = self.screen.wells
        for rule in d.get("mark", []):
            wells = [x.strip().upper() for x in re.split(r"[\s,;]+", rule.get("wells", "")) if x.strip()]
            content = str(rule.get("content", "")).strip()
            if not wells or not content:
                continue
            norm = []
            for x in wells:
                rc = arrayed_io.parse_well(x)
                if not rc:
                    raise ValueError(f"‘{x}’ is not a well name (like B02).")
                norm.append(arrayed_io.well_name(*rc))
            plates = rule.get("plates", "all")
            mask = w["well"].isin(norm)
            if plates not in ("all", "", None):
                ps = {int(p) for p in re.split(r"[\s,;]+", str(plates)) if p.strip().isdigit()}
                mask &= w["plate"].isin(ps)
            w.loc[mask, "content"] = content
            self.roles.setdefault(content, content if content in ROLES else "other")
        for k, v in (d.get("roles") or {}).items():
            if v in ROLES:
                self.roles[k] = v
        if d.get("direction") in ("lower", "higher", "both"):
            self.direction = d["direction"]
        if d.get("channel"):
            self.channel = int(d["channel"])
        self.qc = None
        return self.summary()

    def _role_screen(self):
        """A view of the screen where content names are replaced by their role (sample / neg / pos / other)."""
        scr = self.screen
        w = scr.wells.copy()
        w["content_name"] = w["content"]
        w["content"] = w["content"].map(lambda c: self.roles.get(c, "other"))
        return arrayed_io.ArrayedScreen(w, scr.n_wells, scr.name, scr.notes, scr.description)

    # ------------------------------------------------------------------ jobs
    def run_qc(self, thresholds=None) -> dict:
        if not self.screen:
            raise ValueError("Load a screen first.")
        if thresholds is not None:
            th = {k: float(v) for k, v in thresholds.items() if k in arrayed_stats.QC_THRESHOLDS}
            prefs.set("arrayed_thresholds", th, arrayed_stats.QC_THRESHOLDS)
        th = prefs.get("arrayed_thresholds", arrayed_stats.QC_THRESHOLDS)

        def work(job):
            job.message = "Checking plates…"
            self.qc = arrayed_stats.run_qc(self._role_screen(), self.channel, thresholds=th)
            self.qc["threshold_text"] = {}
        return self.host.start("arrayed-qc", work)

    def run_analysis(self, o: dict) -> dict:
        if not self.screen:
            raise ValueError("Load a screen first.")
        method = o.get("method", "robust_z")
        if method not in arrayed_stats.NORM:
            raise ValueError("Unknown normalisation.")
        kw = dict(method=method, channel=self.channel, log2=bool(o.get("log2")), score=o.get("score", "z"),
                  threshold=float(o.get("threshold", 2.0)), replicate_min=float(o.get("replicate_min", 1.0)),
                  summary=o.get("summary", "mean"), direction=o.get("direction") or self.direction,
                  gene_fraction=float(o.get("gene_fraction", 0.5)), use_flagged=bool(o.get("use_flagged")))

        def work(job):
            job.message = f"Normalising ({method}) and scoring…"
            r = arrayed_stats.analyse(self._role_screen(), **kw)
            r["when"] = time.strftime("%Y-%m-%d %H:%M")
            r["comparison"] = self.compare_published(r) if self.example == "kcviab" else None
            self.result = r
        return self.host.start("arrayed-analysis", work)

    def compare_published(self, r) -> dict | None:
        """KcViab only: agreement with the published 2003 scores (Analysis2003.txt; hit = score ≥ 2)."""
        p = TEST_DATA / "cellhts2/KcViab/Analysis2003.txt"
        if not p.exists():
            return None
        pub = pd.read_csv(p, sep="\t")
        w = r["wells"].merge(pub, left_on=["plate", "well"], right_on=["Plate", "Well"])
        ok = w["score"].notna()
        ph = w["Score"] >= 2
        oh = w["hit"].astype(bool)
        return {"wells": int(len(w)), "r": float(np.corrcoef(w.loc[ok, "score"], w.loc[ok, "Score"])[0, 1]),
                "published_hits": int(ph.sum()), "our_hits": int(oh.sum()), "both": int((ph & oh).sum()),
                "recall": float((ph & oh).sum() / max(1, ph.sum())), "precision": float((ph & oh).sum() / max(1, oh.sum())),
                "exact_2dp": float((w.loc[ok, "score"].round(2) == w.loc[ok, "Score"]).mean())}

    # ------------------------------------------------------------------ views
    def summary(self) -> dict:
        out = {"examples": {k: {"title": v["title"], "about": v["about"]} for k, v in EXAMPLES.items()},
               "examples_available": TEST_DATA.exists(), "files": sorted(self.files), "name": self.name,
               "norm_methods": arrayed_stats.NORM, "hit_defaults": arrayed_stats.HIT_DEFAULTS,
               "default_norm": arrayed_stats.DEFAULT_NORM, "hit_defaults": arrayed_stats.HIT_DEFAULTS,
               "qc_thresholds": {**arrayed_stats.QC_THRESHOLDS, **prefs.get("arrayed_thresholds", arrayed_stats.QC_THRESHOLDS)},
               "qc_threshold_defaults": arrayed_stats.QC_THRESHOLDS, "has_qc": self.qc is not None,
               "qc_overall": self.qc["overall"] if self.qc else None, "has_result": self.result is not None,
               "roles": self.roles, "direction": self.direction, "channel": self.channel, "example": self.example}
        if self.screen:
            s = self.screen.summary()
            s["dims"] = list(self.screen.dims)
            out["screen"] = s
        return out

    def plate_layout(self, plate: int) -> dict:
        if not self.screen:
            raise ValueError("Load a screen first.")
        w = self.screen.wells
        g = w[(w["plate"] == plate)].drop_duplicates("well")
        nr, nc = self.screen.dims
        return {"plate": plate, "dims": [nr, nc],
                "wells": [{"r": int(a), "c": int(b), "well": x, "content": ct, "gene": ge, "reagent": rg}
                          for a, b, x, ct, ge, rg in zip(g["row"], g["col"], g["well"], g["content"], g["gene"], g["reagent"])]}

    def qc_view(self, plate: int | None = None) -> dict:
        if not self.qc:
            raise ValueError("Run the plate checks first.")
        q = dict(self.qc)
        maps = q.pop("raw_maps")
        if plate is None:
            plate = q["plates"][0]
        q["maps"] = {k: v for k, v in maps.items() if k.endswith(f":{plate}")}
        q["plate"] = plate
        # replicate scatter (plate-median scaled), sampled to keep the page light
        w = self._role_screen().wells
        w = w[(w["channel"] == self.channel) & (w["content"] == "sample")]
        w = w.assign(v=w["value"].where(~w["flagged"].astype(bool)))
        w["pn"] = w["v"] / w.groupby(["replicate", "plate"])["v"].transform("median")
        piv = w.pivot_table(index=["plate", "well"], columns="replicate", values="pn").dropna()
        if piv.shape[1] >= 2:
            piv = piv.sample(min(len(piv), 12000), random_state=0) if len(piv) > 12000 else piv
            q["rep_scatter"] = {"x": piv.iloc[:, 0].round(4).tolist(), "y": piv.iloc[:, 1].round(4).tolist(),
                                "labels": [f"plate {p} {x}" for p, x in piv.index]}
        # row and column profiles
        nr, nc = self.screen.dims
        q["row_profile"] = [float(x) if np.isfinite(x) else None for x in w.groupby("row")["pn"].median().reindex(range(nr))]
        q["col_profile"] = [float(x) if np.isfinite(x) else None for x in w.groupby("col")["pn"].median().reindex(range(nc))]
        return q

    def result_view(self, plate: int | None = None) -> dict:
        r = self.result
        if not r:
            raise ValueError("Run the analysis first.")
        W = r["wells"]
        reps = r["replicates"]
        plates = sorted(int(x) for x in W["plate"].unique())
        plate = plate or plates[0]
        nr, nc = self.screen.dims
        P = W[W["plate"] == plate]
        M = np.full((nr, nc), np.nan)
        M[P["row"].to_numpy(), P["col"].to_numpy()] = P["score"].to_numpy()
        hits = W[W["hit"]].sort_values("score", ascending=False)
        cols = ["plate", "well", "content", "reagent", "gene", "score"] + [f"score_{x}" for x in reps] + [f"raw_{x}" for x in reps]
        samp = W[W["content"] == "sample"]
        sc = samp[[f"score_{x}" for x in reps]].dropna()
        scat = None
        if len(reps) >= 2:
            sc2 = sc.sample(min(len(sc), 12000), random_state=0) if len(sc) > 12000 else sc
            hit_ix = set(samp.index[samp["hit"]])
            scat = {"x": sc2.iloc[:, 0].round(3).tolist(), "y": sc2.iloc[:, 1].round(3).tolist(),
                    "hit": [i in hit_ix for i in sc2.index],
                    "labels": [f"plate {samp.at[i, 'plate']} {samp.at[i, 'well']} {samp.at[i, 'gene']}" for i in sc2.index]}
        allsc = samp["score"].dropna().to_numpy()
        hist, edges = np.histogram(np.clip(allsc, -8, 12), bins=80)
        g = r["genes"]
        return {"settings": r["settings"], "when": r["when"], "replicates": reps, "plates": plates, "plate": plate,
                "map": M.round(3).tolist(), "dims": [nr, nc],
                "map_info": [{"r": int(a), "c": int(b), "well": x, "content": ct, "gene": ge, "hit": bool(h)}
                             for a, b, x, ct, ge, h in zip(P["row"], P["col"], P["well"], P["content"], P["gene"], P["hit"])],
                "n_hits": int(W["hit"].sum()), "n_samples": int((W["content"] == "sample").sum()),
                "n_gene_hits": int(g["hit"].sum()) if len(g) else 0, "n_genes": int(len(g)),
                "hits_columns": cols, "hits": hits[cols].head(3000).replace({np.nan: None}).values.tolist(),
                "genes_columns": list(g.columns), "genes": g.head(3000).replace({np.nan: None}).values.tolist(),
                "scatter": scat, "hist": {"edges": edges.round(3).tolist(), "counts": hist.tolist()},
                "hits_per_plate": W.groupby("plate")["hit"].sum().astype(int).tolist(),
                "comparison": r.get("comparison")}

    def enrich(self, d: dict) -> dict:
        from .session import run_enrichment
        if not self.result:
            raise ValueError("Run the analysis first.")
        g = self.result["genes"]
        if g.empty:
            raise ValueError("This screen has no gene annotation, so gene sets cannot be tested.")
        hits = list(g.loc[g["hit"], "gene"])
        scores = dict(zip(g["gene"], g["median_score"]))
        out = run_enrichment(d, hits, scores, list(scores), True, list(scores))
        if out["n_sets"] == 0:
            raise ValueError("None of the gene sets overlap this screen's gene names (the bundled sets use human and "
                             "mouse gene symbols). Use a .gmt file with the same gene IDs as the screen.")
        return out

    def export(self, what: str, fmt: str = "csv") -> tuple[bytes, str, str]:
        safe = re.sub(r"[^\w.+-]+", "_", self.name or "arrayed").strip("_") or "arrayed"
        if what == "wells":
            if not self.result:
                raise ValueError("Run the analysis first.")
            df, fn = self.result["wells"], f"{safe}.wells"
        elif what == "hits":
            if not self.result:
                raise ValueError("Run the analysis first.")
            df, fn = self.result["wells"][self.result["wells"]["hit"]], f"{safe}.hits"
        elif what == "genes":
            if not self.result:
                raise ValueError("Run the analysis first.")
            df, fn = self.result["genes"], f"{safe}.genes"
        elif what == "raw":
            df, fn = self.screen.wells, f"{safe}.raw"
        elif what == "qc":
            if not self.qc:
                raise ValueError("Run the plate checks first.")
            df = pd.DataFrame([{"check": c["name"], "status": c["status"], "why": c["why"]} for c in self.qc["checks"]])
            fn = f"{safe}.qc"
        else:
            raise ValueError("Unknown export.")
        if fmt == "xlsx":
            buf = io.BytesIO()
            df.to_excel(buf, index=False, sheet_name=what[:30])
            return buf.getvalue(), fn + ".xlsx", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        return df.to_csv(index=False).encode(), fn + ".csv", "text/plain; charset=utf-8"
