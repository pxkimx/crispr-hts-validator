"""BAGEL2 Bayes factors, drugZ and copy-number correction for the pooled project: running, views and exports."""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import bagel, cnv, drugz, qc

METHODS = {
    "mageck": "MAGeCK-style: negative-binomial guide test and α-RRA gene ranking, both directions (default).",
    "bagel": "BAGEL2 Bayes factors: how strongly each gene behaves like the reference essential genes "
             "(fitness screens against plasmid/T0; needs CEG/NEG reference genes in the library).",
    "drugz": "drugZ: treated vs untreated at the same time point; finds genes whose loss sensitises (synthetic "
             "lethal) or protects (suppressor). Replicates are paired in the order listed.",
}
# Confirmed by the user on 2026-11-07, the value BAGEL2's own tutorial uses. On HAP1 TKOv3 T18 it calls 1,796
# genes at 94.0% recall with no false positives among the reference sets; it also holds 98% precision on the weak
# T3 timepoint, where BF > 1 would let 55 non-essential genes through. The results page can show any cut-off; this
# is the default and the one exports use.
BF_CUTOFF = 5.0


def reference_sets(genes) -> tuple[set, set, str]:
    """CEG/NEG for the species the gene names look like (human symbols upper case, mouse capitalised)."""
    from .enrichment import guess_species
    species = guess_species(list(genes))
    ceg, neg = qc.load_reference_sets(species)
    return ceg, neg, species


# ------------------------------------------------------------------ copy-number correction
def cn_shift(counts: np.ndarray, guides, genes, n_ctrl: int, positions: pd.DataFrame, progress=None) -> dict:
    """Guide fold change (median-normalised mean of treatment vs control) -> CRISPRcleanR-style segment shift."""
    N = qc.median_ratio_normalise(counts)
    lfc = np.log2(N[:, n_ctrl:].mean(axis=1) + 1) - np.log2(N[:, :n_ctrl].mean(axis=1) + 1)
    return cnv.correct(lfc, guides, genes, positions, progress=progress)


def corrected_counts(counts: np.ndarray, n_ctrl: int, shift: np.ndarray) -> np.ndarray:
    """As CRISPRcleanR's ccr.correctCounts: treatment counts scaled so that their fold change loses the segment mean."""
    X = np.asarray(counts, dtype=float).copy()
    X[:, n_ctrl:] = X[:, n_ctrl:] * (2.0 ** (-shift))[:, None]
    return X


# ------------------------------------------------------------------ runners
def run_bagel(t, comp, positions=None, cn=False, seed: int = 1, progress=None) -> dict:
    ctrl = [t.samples.index(s) for s in comp["control"]]
    treat = [t.samples.index(s) for s in comp["treatment"]]
    X = t.counts[:, ctrl + treat]
    fc, keep = bagel.fold_changes(X, list(range(len(ctrl))))
    genes = np.asarray(t.genes, dtype=object)[keep]
    cninfo = None
    if cn:
        sh = cn_shift(X, t.guides, t.genes, len(ctrl), positions)
        fc = fc - sh["shift"][keep][:, None]
        cninfo = {"segments": sh["segments"], "n_positioned": sh["n_positioned"],
                  "corrected": int(np.sum(np.abs(sh["shift"]) > 0.25))}
    ceg, neg, species = reference_sets(t.genes)
    r = bagel.bayes_factors(fc, genes, ceg, neg, seed=seed, progress=progress)
    pr = bagel.precision_recall(r["genes"], r["bf"], ceg, neg)
    n_guides = pd.Series(genes).value_counts()
    return {"method": "bagel", "genes": r["genes"], "bf": r["bf"], "guide_bf": r["guide_bf"], "fc": fc,
            "keep": keep, "fits": r["fits"], "pr": pr, "summary": bagel.summary_at(pr, BF_CUTOFF),
            "guides_per_gene": n_guides.to_dict(), "species": species, "cn": cninfo,
            "ref": {"ceg": ceg, "neg": neg}}


def run_drugz(t, comp, paired=True, remove_genes=(), progress=None) -> dict:
    ctrl = [t.samples.index(s) for s in comp["control"]]
    treat = [t.samples.index(s) for s in comp["treatment"]]
    r = drugz.drugz(t.counts, t.guides, t.genes, ctrl, treat, paired=paired, remove_genes=remove_genes)
    if progress:
        progress(1.0)
    return {"method": "drugz", **r}


# ------------------------------------------------------------------ views
def view(name: str, entry: dict) -> dict:
    r = entry["result"]
    base = {"name": name, "method": r["method"], "comparison": entry["comparison"], "settings": entry["settings"],
            "when": entry["when"]}
    if r["method"] == "bagel":
        ceg, neg = r["ref"]["ceg"], r["ref"]["neg"]
        order = np.argsort(-r["bf"], kind="mergesort")
        rank = np.empty(len(order), dtype=int)
        rank[order] = np.arange(1, len(order) + 1)
        rows = [[g, int(r["guides_per_gene"].get(g, 0)), float(b), int(k),
                 "essential" if g in ceg else ("non-essential" if g in neg else "")]
                for g, b, k in zip(r["genes"], r["bf"], rank)]
        rows.sort(key=lambda x: x[3])
        pr = r["pr"]
        step = max(1, len(pr["bf"]) // 400)
        return {**base, "columns": ["gene", "guides", "bf", "rank", "reference"], "rows": rows,
                "summary": r["summary"], "cutoff": BF_CUTOFF, "species": r["species"],
                "pr": {"recall": pr["recall"][::step].round(4).tolist(), "precision": pr["precision"][::step].round(4).tolist(),
                       "bf": pr["bf"][::step].round(3).tolist()},
                "fits": r["fits"], "cn": _cn_view(r.get("cn"))}
    if r["method"] == "drugz":
        g = r["genes"].reset_index()
        cols = ["GENE", "sumZ", "numObs", "normZ", "pval_synth", "fdr_synth", "rank_synth", "pval_supp", "fdr_supp", "rank_supp"]
        rows = g[cols].to_numpy().tolist()
        return {**base, "columns": ["gene", "sumZ", "numObs", "normZ", "p_synth", "fdr_synth", "rank_synth",
                                    "p_supp", "fdr_supp", "rank_supp"], "rows": rows, "paired": r["paired"]}
    raise ValueError("Unknown method.")


def _cn_view(cn):
    if not cn:
        return None
    seg = sorted(cn["segments"], key=lambda s: s["mean_lfc"])
    return {"n_positioned": cn["n_positioned"], "shifted_guides": cn["corrected"], "segments": len(cn["segments"]),
            "strongest": seg[:6] + seg[-3:]}


def gene_view(t, entry: dict, gene: str) -> dict:
    r = entry["result"]
    c = entry["comparison"]
    cols = c["control"] + c["treatment"]
    idx = [i for i, x in enumerate(t.genes) if x == gene]
    if not idx:
        raise ValueError("Gene not found.")
    rows = []
    if r["method"] == "bagel":
        kept = np.cumsum(r["keep"]) - 1
        for i in idx:
            if not r["keep"][i]:
                continue
            k = kept[i]
            rows.append({"guide": t.guides[i], "raw": [float(t.counts[i, t.samples.index(s)]) for s in cols],
                         "fc": [float(x) for x in r["fc"][k]], "bf": [float(x) for x in r["guide_bf"][k]]})
        j = list(r["genes"]).index(gene) if gene in set(r["genes"]) else None
        summ = {"bf": float(r["bf"][j])} if j is not None else {}
        return {"gene": gene, "method": "bagel", "samples": cols, "n_control": len(c["control"]),
                "replicates": c["treatment"], "guides": rows, "summary": summ}
    if r["method"] == "drugz":
        f = r["guides"]
        sub = f[f["GENE"] == gene]
        zc = r["zcols"]
        fcc = [x for x in f.columns if x.startswith("fc_") or x == "mean_fc"]
        for gid, row in sub.iterrows():
            i = t.guides.index(gid)
            rows.append({"guide": gid, "raw": [float(t.counts[i, t.samples.index(s)]) for s in cols],
                         "fc": [float(row[x]) for x in fcc], "z": [float(row[x]) for x in zc]})
        gs = r["genes"]
        summ = gs.loc[gene].to_dict() if gene in gs.index else {}
        return {"gene": gene, "method": "drugz", "samples": cols, "n_control": len(c["control"]),
                "replicates": [f"{a} → {b}" for a, b in zip(c["control"], c["treatment"])] if r["paired"] else ["mean"],
                "guides": rows, "summary": {k: float(v) for k, v in summ.items()}}
    raise ValueError("Unknown method.")


def gene_frame(entry: dict) -> pd.DataFrame:
    v = view("", entry)
    df = pd.DataFrame(v["rows"], columns=v["columns"])
    if entry["result"]["method"] == "bagel":
        df = df.rename(columns={"gene": "GENE", "bf": "BF"})
    return df


def guide_frame(t, entry: dict) -> pd.DataFrame:
    r = entry["result"]
    if r["method"] == "bagel":
        c = entry["comparison"]
        g = np.asarray(t.guides, dtype=object)[r["keep"]]
        df = pd.DataFrame(r["fc"], columns=[f"FC {s}" for s in c["treatment"]])
        for j, s in enumerate(c["treatment"]):
            df[f"BF {s}"] = r["guide_bf"][:, j]
        df.insert(0, "GENE", np.asarray(t.genes, dtype=object)[r["keep"]])
        df.insert(0, "REAGENT_ID", g)
        return df
    if r["method"] == "drugz":
        return r["guides"].reset_index().rename(columns={"index": "sgRNA"})
    raise ValueError("Unknown method.")


def hits_of(entry: dict, direction: str = "neg", cutoff: float = 0.05) -> tuple[list, dict, bool]:
    """Hit genes, gene scores (for ranked enrichment) and whether higher scores are 'better' for this direction."""
    r = entry["result"]
    m = r.get("method", "mageck")
    if m == "mageck":
        rows = r["genes"]
        key = "neg" if direction == "neg" else "pos"
        hits = [g["gene"] for g in rows if g[f"{key}_fdr"] <= cutoff]
        scores = {g["gene"]: -np.log10(max(g[f"{key}_p"], 1e-300)) for g in rows}
        return hits, scores, True
    if m == "bagel":
        hits = [g for g, b in zip(r["genes"], r["bf"]) if b > BF_CUTOFF]
        return hits, dict(zip(r["genes"], r["bf"])), True
    if m == "drugz":
        g = r["genes"]
        if direction == "neg":
            return list(g.index[g["fdr_synth"] <= cutoff]), (-g["normZ"]).to_dict(), True
        return list(g.index[g["fdr_supp"] <= cutoff]), g["normZ"].to_dict(), True
    raise ValueError("Unknown method.")


def hit_counts(entry: dict) -> tuple[int, int]:
    r = entry["result"]
    m = r.get("method", "mageck")
    if m == "mageck":
        return (sum(1 for g in r["genes"] if g["neg_fdr"] < 0.05), sum(1 for g in r["genes"] if g["pos_fdr"] < 0.05))
    if m == "bagel":
        return int((r["bf"] > BF_CUTOFF).sum()), 0
    if m == "drugz":
        g = r["genes"]
        return int((g["fdr_synth"] < 0.05).sum()), int((g["fdr_supp"] < 0.05).sum())
    return 0, 0
