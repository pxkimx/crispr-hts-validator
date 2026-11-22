"""Compare two pooled screens: which dependencies they share, and which belong to one of them.

The question this answers is "what does this cell line (or this drug) need that the other does not" — for
example an HGPS line against its isogenic corrected control, or treated against untreated.

Two screens cannot be compared by their raw log fold changes. A later timepoint, a faster-growing line or a
deeper run all stretch the scale, so the same gene gets a bigger number in one screen than the other for
reasons that have nothing to do with biology. Each screen is therefore put on a scale of its own making,
using the genes whose behaviour is known in advance: the median core-essential gene (CEGv2) is set to -1 and
the median non-essential gene (NEGv1) to 0 (Hart et al. 2015; the same idea as CERES and Chronos). A scaled
effect of -1 then means "as essential as a ribosomal subunit" in either screen, and the two can be subtracted.

The differential test works on the guides, not the gene summaries. When both screens share a library, each
guide appears in both, and its difference B - A cancels the guide's own efficiency, which is the largest
nuisance in a screen. Those differences are turned into z-scores against an empirical null taken from the
screen itself, exactly as drugZ does it (a standard deviation estimated in a sliding window of guides ordered
by abundance, never allowed to fall as abundance rises), and combined per gene as sumZ / sqrt(n). That is the
same machinery the app already uses for drugZ, and when the two comparisons share one reference sample the
two analyses are asking the same question and give closely matching answers.

Thresholds here are scientific choices and are listed in the app for Paul to confirm.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from . import qc
from .drugz import _empirical_bayes

# Scientific choices, shown in the app. DEFAULTS["dependency"] follows the convention that comes with the
# scaling: -1 is the median core-essential gene, and -0.5 is the usual line for calling a gene a dependency.
DEFAULTS = {
    "cmp_fdr": 0.05,          # FDR below which a difference is called
    "cmp_min_delta": 0.5,     # and the smallest difference in scaled effect worth reporting
    "cmp_min_guides": 3,      # guides a gene needs in both screens
    "cmp_dependency": -0.5,   # scaled effect at or below which a gene is a dependency in a screen
}


def guide_lfc(counts: np.ndarray, samples: list[str], comparison: dict) -> np.ndarray:
    """One comparison's guide log2 fold changes, against its own control samples.

    Each screen keeps its own reference: in a two-arm experiment the HGPS arm is measured against the HGPS
    T0 and the corrected arm against its own, so nothing is borrowed across the pair.
    """
    N = qc.median_ratio_normalise(np.asarray(counts, dtype=float))
    t = [samples.index(s) for s in comparison["treatment"]]
    c = [samples.index(s) for s in comparison["control"]]
    return np.log2(N[:, t].mean(axis=1) + 1.0) - np.log2(N[:, c].mean(axis=1) + 1.0)


def gene_effect(lfc: np.ndarray, genes: np.ndarray) -> pd.Series:
    """A gene's effect is the median of its guides, as elsewhere in the app."""
    return pd.Series(np.asarray(lfc, dtype=float)).groupby(np.asarray(genes, dtype=object).astype(str)).median()


def scale_factors(effects: pd.Series, ceg: set, neg: set) -> dict:
    """Put a screen on the common scale: median core-essential gene at -1, median non-essential at 0."""
    e = effects[[g in ceg for g in effects.index]]
    n = effects[[g in neg for g in effects.index]]
    if len(e) < 10 or len(n) < 10:
        return {"ok": False, "n_ceg": int(len(e)), "n_neg": int(len(n)),
                "median_ceg": float(e.median()) if len(e) else float("nan"),
                "median_neg": float(n.median()) if len(n) else float("nan"), "span": float("nan")}
    me, mn = float(e.median()), float(n.median())
    span = mn - me                       # how far this screen separates the two sets
    return {"ok": span > 0, "n_ceg": int(len(e)), "n_neg": int(len(n)),
            "median_ceg": me, "median_neg": mn, "span": span}


def scaled(effects: pd.Series, f: dict) -> pd.Series:
    if not f.get("ok"):
        return effects.astype(float)
    return (effects - f["median_neg"]) / f["span"]


def compare(counts, guides, genes, samples, comp_a: dict, comp_b: dict, species: str = "auto",
            thresholds: dict | None = None, half_window: int = 500) -> dict:
    """Compare two comparisons taken from one count table (so every guide is in both screens)."""
    th = {**DEFAULTS, **(thresholds or {})}
    counts = np.asarray(counts, dtype=float)
    gene_arr = np.asarray(genes, dtype=object).astype(str)
    if species == "auto":
        from .enrichment import guess_species
        species = guess_species(list(gene_arr))
    try:
        ceg, neg = qc.load_reference_sets(species)
    except ValueError:
        ceg, neg = set(), set()

    la, lb = guide_lfc(counts, samples, comp_a), guide_lfc(counts, samples, comp_b)
    ea, eb = gene_effect(la, gene_arr), gene_effect(lb, gene_arr)
    fa, fb = scale_factors(ea, ceg, neg), scale_factors(eb, ceg, neg)
    sa, sb = scaled(ea, fa), scaled(eb, fb)

    # --- guide-level difference, z-scored against the screen's own spread (drugZ's empirical Bayes)
    d = lb - la
    ctrl_cols = [samples.index(s) for s in comp_a["control"]] + [samples.index(s) for s in comp_b["control"]]
    abundance = counts[:, ctrl_cols].mean(axis=1)
    order = np.argsort(-abundance, kind="mergesort")          # most abundant first, as drugZ sorts
    eb_sd = np.empty(len(d))
    eb_sd[order] = _empirical_bayes(d[order], half_window)
    z = np.divide(d, eb_sd, out=np.zeros_like(d), where=eb_sd > 0)

    df = pd.DataFrame({"gene": gene_arr, "guide": np.asarray(guides, dtype=object).astype(str),
                       "lfc_a": la, "lfc_b": lb, "delta_guide": d, "z": z})
    agg = df.groupby("gene", sort=True)
    sumz = agg["z"].sum()
    nobs = agg["z"].count()
    keep = nobs >= int(th["cmp_min_guides"])
    sumz, nobs = sumz[keep], nobs[keep]
    normz = pd.Series(stats.zscore(sumz / np.sqrt(nobs)), index=sumz.index) if len(sumz) > 1 else sumz * 0.0
    p = 2.0 * stats.norm.sf(np.abs(normz))                     # two-sided: a gene may go either way
    fdr = pd.Series(_bh(p), index=normz.index)
    # how many of a gene's guides move the way the gene does — a gene carried by one guide is not a result
    sign = np.sign(normz.reindex(df["gene"].to_numpy()).to_numpy())
    df["_agree"] = (np.sign(df["z"].to_numpy()) == sign) & (sign != 0)
    agree = df.groupby("gene", sort=True)["_agree"].mean()[keep]

    genes_out = pd.DataFrame({
        "gene": normz.index,
        "guides": nobs.to_numpy(),
        "effect_a": ea.reindex(normz.index).to_numpy(), "effect_b": eb.reindex(normz.index).to_numpy(),
        "scaled_a": sa.reindex(normz.index).to_numpy(), "scaled_b": sb.reindex(normz.index).to_numpy(),
        "sumZ": sumz.to_numpy(), "normZ": normz.to_numpy(), "p": p, "fdr": fdr.to_numpy(),
        "agreement": agree.reindex(normz.index).to_numpy(),
    })
    genes_out["delta"] = genes_out["scaled_b"] - genes_out["scaled_a"]
    dep = float(th["cmp_dependency"])
    da, db = genes_out["scaled_a"] <= dep, genes_out["scaled_b"] <= dep
    # "different" is the test's verdict; "call" is which screen needs the gene. A gene can be a dependency in
    # both and still differ — mitochondrial genes go deeper under starvation — so the two are kept apart
    # rather than folded into one label that would hide it.
    genes_out["different"] = ((genes_out["fdr"] <= float(th["cmp_fdr"]))
                              & (genes_out["delta"].abs() >= float(th["cmp_min_delta"])))
    genes_out["call"] = np.select([da & db, da & ~db, ~da & db], ["shared", "specific to A", "specific to B"],
                                  default="neither")
    genes_out["stronger_in"] = np.where(~genes_out["different"], "",
                                        np.where(genes_out["delta"] < 0, "B", "A"))
    genes_out = genes_out.sort_values("normZ").reset_index(drop=True)

    return {"genes": genes_out, "guide": df.drop(columns=["_agree"]),
            "a": {"name": comp_a["name"], "scale": fa}, "b": {"name": comp_b["name"], "scale": fb},
            "species": species, "thresholds": th,
            "concordance": _concordance(genes_out, th),
            "shared_reference": sorted(comp_a["control"]) == sorted(comp_b["control"])}


def _bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p, kind="mergesort")
    q = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1.0))[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1.0)
    return out


def _concordance(g: pd.DataFrame, th: dict) -> dict:
    a, b = g["scaled_a"].to_numpy(), g["scaled_b"].to_numpy()
    ok = np.isfinite(a) & np.isfinite(b)
    dep = float(th["cmp_dependency"])
    sa, sb = set(g["gene"][(g["scaled_a"] <= dep)]), set(g["gene"][(g["scaled_b"] <= dep)])
    union = len(sa | sb)
    return {
        "n_genes": int(len(g)),
        "pearson": float(np.corrcoef(a[ok], b[ok])[0, 1]) if ok.sum() > 2 else float("nan"),
        "spearman": float(stats.spearmanr(a[ok], b[ok]).statistic) if ok.sum() > 2 else float("nan"),
        "dep_a": len(sa), "dep_b": len(sb), "dep_both": len(sa & sb),
        "jaccard": (len(sa & sb) / union) if union else float("nan"),
        "n_different": int(g["different"].sum()),
        "n_specific_a": int(((g["call"] == "specific to A") & g["different"]).sum()),
        "n_specific_b": int(((g["call"] == "specific to B") & g["different"]).sum()),
        "n_shared_stronger": int(((g["call"] == "shared") & g["different"]).sum()),
    }


def quality_note(qc_result: dict | None, comp_a: dict, comp_b: dict) -> list[dict]:
    """How well each screen separates essential from non-essential genes. A comparison between two screens is
    only as good as the weaker one, so this is shown beside the result rather than buried in the QC page."""
    if not qc_result:
        return []
    ess = next((c for c in qc_result["checks"] if c["id"] == "essentials"), None)
    if not ess or not ess.get("per_condition"):
        return []
    out = []
    for comp in (comp_a, comp_b):
        hit = next((p for p in ess["per_condition"] if set(p["samples"]) & set(comp["treatment"])), None)
        if hit:
            out.append({"comparison": comp["name"], "condition": hit["condition"],
                        "auc": hit["value"], "status": hit["status"]})
    return out
