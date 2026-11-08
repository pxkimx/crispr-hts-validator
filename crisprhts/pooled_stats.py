"""Normalisation, guide-level test and MAGeCK-style gene test for pooled screens.

The guide test follows `mageck test` 0.5.9.5 (crisprFunction.py) closely so results can be checked against it:
  * median-ratio size factors on the samples of one comparison (geometric mean of count+1, upper median);
  * control and treatment means are medians of the normalised replicate counts;
  * a mean-variance model log2(var - mean + 1) = k log2(mean + 1) + b, fitted by least squares weighted by the
    mean, on the control replicates (or all samples of the comparison when there is only one control);
  * guide p-values from a normal approximation truncated at zero (MAGeCK's 'NB' test as shipped);
  * guide scores = normal quantile of the depletion p-value (MAGeCK's qnorm table range), used for ranking;
  * alpha-RRA on the ranked guides for each direction, alpha = fraction of guides with p <= 0.25.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.special import ndtr, ndtri

from .rra import rra

# ----------------------------------------------------------------------------- normalisation
NORM_METHODS = {
    "median": "Median ratio (default): each sample is scaled so that its median guide ratio to the geometric "
              "mean of all samples is 1 (as MAGeCK and DESeq do). Robust to a few guides changing strongly.",
    "total": "Total count: each sample is scaled to the same total number of reads.",
    "control": "Control guides: median ratio computed on non-targeting / safe-harbour control guides only. "
               "Best when many genes drop out (strong screens).",
    "none": "No normalisation: counts are used as given.",
}


def size_factors(counts: np.ndarray, method: str = "median", control_mask=None) -> tuple[np.ndarray, list[str]]:
    """Per-sample factors to multiply counts by (MAGeCK normalizeCounts). Returns (factors, warnings)."""
    X = np.asarray(counts, dtype=float)
    warnings: list[str] = []
    if method == "none":
        return np.ones(X.shape[1]), warnings
    sub = X
    if method == "control":
        if control_mask is None or not np.any(control_mask):
            raise ValueError("Control-guide normalisation needs control guides (none were found).")
        sub = X[np.asarray(control_mask, dtype=bool)]
        method = "median"
    tot = sub.sum(axis=0)
    if np.any(tot <= 0):
        raise ValueError("A sample has no reads, so it cannot be normalised.")
    total_f = tot.mean() / tot
    if method == "total":
        return total_f, warnings
    keep = sub.sum(axis=1) > 0
    gm = np.exp(np.log(sub[keep] + 1.0).mean(axis=1))
    gm[gm <= 0] = 1.0
    ratios = np.sort(sub[keep] / gm[:, None], axis=0)
    med = ratios[len(ratios) // 2, :]
    zero_frac = (sub == 0).sum(axis=0) / len(sub)
    use_total = False
    for j in range(X.shape[1]):
        if med[j] <= 0:
            warnings.append(f"Sample {j + 1} has a median count of zero, so the median ratio cannot be used; "
                            "switched to total-count normalisation.")
            use_total = True
        if zero_frac[j] > 0.45:
            warnings.append(f"Sample {j + 1} has {zero_frac[j]:.0%} zero-count guides; median normalisation is "
                            "unstable, so total-count normalisation was used.")
            use_total = True
        elif zero_frac[j] > 0.3:
            warnings.append(f"Sample {j + 1} has {zero_frac[j]:.0%} zero-count guides; median normalisation may be "
                            "unstable. Consider total-count normalisation.")
    if use_total:
        return total_f, warnings
    f = 1.0 / med
    if np.any((f < 0.2) | (f > 5)):
        warnings.append("Some samples have unusually small or large normalisation factors (<0.2 or >5).")
    return f, warnings


# ----------------------------------------------------------------------------- helpers (MAGeCK maths)
def _median_rows(m: np.ndarray) -> np.ndarray:
    return np.median(m, axis=1)


def _var_rows(m: np.ndarray, centre: np.ndarray) -> np.ndarray:
    return ((m - centre[:, None]) ** 2).sum(axis=1) / max(m.shape[1] - 1, 1.0)


def mean_var_model(m: np.ndarray) -> tuple[float, float]:
    """MAGeCK modelmeanvar(method='linear'): returns (k, b)."""
    mean = _median_rows(m)
    var = _var_rows(m, mean)
    if len(mean) > 10:
        mm = float(np.median(mean))
        vm = float(((mean - mm) ** 2).sum() / max(len(mean) - 1, 1.0))
        high = mm + 4 * math.sqrt(vm)
    else:
        high = float(mean.max()) * 2
    good = (mean < var) & (mean < high)
    mg, vg = mean[good], var[good] - mean[good]
    if len(mg) < 10:
        mg = mean
        vg = np.maximum(var - mean, 0.01)
    x = np.log2(mg + 1)
    y = np.log2(vg + 1)
    w = mg
    nw = w.sum()
    sy, sx = (y * w).sum(), (x * w).sum()
    sx2, sxy = (x * x * w).sum(), (x * y * w).sum()
    den = nw * sx2 - sx * sx
    if den == 0:
        return 1.0, 0.0
    a = (sy * sx2 - sx * sxy) / den
    k = (nw * sxy - sx * sy) / den
    return max(k, 1.0), max(a, 0.0)


_QSTEP = 0.001
_QGRID = None


def _qnorm_log(lp: np.ndarray) -> np.ndarray:
    """MAGeCK QNormConverter: normal quantile of exp(lp) by linear interpolation on a 0.001 grid over
    [-24, -0.001]; NaN outside that range (MAGeCK then keeps the plain z score)."""
    global _QGRID
    if _QGRID is None:
        g = np.round(np.arange(-24000, 0) * _QSTEP, 3)
        _QGRID = (g, ndtri(np.exp(g)))
    g, q = _QGRID
    qmin, qmax = g[0], g[-1]
    step = g[1] - g[0]
    i = np.floor((lp - qmin) / step).astype(np.int64)
    i = np.clip(i, 0, len(g) - 2)
    fr = np.clip((lp - g[i]) / step, 0.0, 1.0)
    val = q[i] * (1 - fr) + q[i + 1] * fr
    return np.where((lp >= qmin) & (lp <= qmax), val, np.nan)


def bh(p: np.ndarray) -> np.ndarray:
    """MAGeCK pFDR(method='fdr'): Benjamini-Hochberg."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(-p, kind="mergesort")
    adj = n / (n - np.arange(n)) * p[order]
    adj = np.minimum.accumulate(adj)
    out = np.empty(n)
    out[order] = adj
    return out


# ----------------------------------------------------------------------------- the test
DEFAULTS = {
    "norm": "median",
    "remove_zero": "both",
    "remove_zero_threshold": 0.0,
    "gene_test_fdr_threshold": 0.25,
    "gene_lfc": "median",
}


def guide_test(counts: np.ndarray, n_ctrl: int, norm: str = "median", control_mask=None,
               remove_zero: str = "both", remove_zero_threshold: float = 0.0, var_cols=None) -> dict:
    """Guide-level test for one comparison. `counts` columns are the control samples followed by the treatment
    samples (raw counts). Returns per-guide arrays."""
    X = np.asarray(counts, dtype=float)
    nt = X.shape[1] - n_ctrl
    if n_ctrl < 1 or nt < 1:
        raise ValueError("A comparison needs at least one control and one treatment sample.")
    f, warns = size_factors(X, norm, control_mask)
    N = X * f[None, :]
    C, T = N[:, :n_ctrl], N[:, n_ctrl:]
    if var_cols is not None and len(var_cols) > 1:
        V = N[:, var_cols]
    elif n_ctrl > 1:
        V = C
    else:
        V = N
    k, b = mean_var_model(V)
    cmean = _median_rows(C)
    valid1 = cmean > remove_zero_threshold
    pos = cmean[cmean > 0]
    if len(pos) == 0:
        raise ValueError("Every guide has zero reads in the control samples.")
    cmean = np.maximum(cmean, pos.min())
    var_c = _var_rows(V, _median_rows(V))
    if (var_cols is not None and len(var_cols) > 1) or n_ctrl > 1 or nt == 1:
        rawvar = var_c
    else:
        var_t = _var_rows(T, _median_rows(T))
        fc = (n_ctrl - 1) / (n_ctrl - 1 + nt - 1)
        rawvar = var_c * fc + var_t * (1 - fc)
    adjvar = cmean ** k * 2 ** b + cmean
    tmean = _median_rows(T)
    valid2 = tmean > remove_zero_threshold
    valid = {"control": valid1, "treatment": valid2, "any": valid1 & valid2,
             "both": valid1 | valid2}.get(remove_zero, np.ones(len(X), dtype=bool))
    if not valid.any():
        raise ValueError("No guides are left after removing zero-count guides.")
    tpos = tmean[tmean > 0]
    tadj = np.maximum(tmean, tpos.min() if len(tpos) else 0.0)
    sd = np.sqrt(adjvar)
    z = (tadj - cmean) / sd
    z0 = (0.0 - cmean) / sd
    p0 = ndtr(z0)
    p_low = (ndtr(z) - p0) / (1 - p0)
    p_high = ndtr(-z) / (1 - p0)
    theta0 = (tmean - cmean) / sd
    nz = p_low[p_low * 0.1 > 0] * 0.1
    floor = nz.min() if len(nz) else 1e-300
    lp = np.log(np.maximum(p_low, floor))
    q = _qnorm_log(lp)
    theta_low = np.where(np.isnan(q), theta0, q)
    theta = np.where(tmean < cmean, theta_low, theta0)
    p_two = np.where(p_low < p_high, 2 * p_low, 2 * p_high)
    fdr = bh(p_two)
    lfc = np.log2(tmean + 1.0) - np.log2(cmean + 1.0)
    return {"factors": f, "norm_warnings": warns, "model": (k, b), "control_mean": cmean, "treat_mean": tmean,
            "control_var": rawvar, "adj_var": adjvar, "score": theta, "p_low": p_low, "p_high": p_high,
            "p_two": p_two, "fdr": fdr, "lfc": lfc, "valid": valid, "normalised": N}


def mageck_test(counts: np.ndarray, genes, guides, n_ctrl: int, norm: str = "median", control_mask=None,
                remove_zero: str = "both", remove_zero_threshold: float = 0.0, gene_test_fdr_threshold: float = 0.25,
                gene_lfc: str = "median", use_control_null: bool = False, progress=None) -> dict:
    """Full MAGeCK-style test of one comparison: guide table plus gene table for both directions."""
    genes = np.asarray(genes, dtype=object)
    guides = np.asarray(guides, dtype=object)
    g = guide_test(counts, n_ctrl, norm, control_mask, remove_zero, remove_zero_threshold)
    valid = g["valid"]
    nv = valid.sum()
    n_low = int(((g["p_low"] <= gene_test_fdr_threshold) & valid).sum())
    a_low = min(max(n_low / nv, 0.05), 0.5)
    n_high_all = int((g["p_high"] <= gene_test_fdr_threshold).sum())
    n_high = int(((g["p_high"] <= gene_test_fdr_threshold) & valid).sum())
    a_high = n_high / nv if n_high_all > 0 else 0.01
    cm = np.asarray(control_mask, dtype=bool)[valid] if (use_control_null and control_mask is not None) else None

    def prog(lo, hi):
        return (lambda x: progress(lo + (hi - lo) * x)) if progress else None

    neg = rra(g["score"][valid], genes[valid], a_low, control_mask=cm, guide_names=guides[valid],
              progress=prog(0.05, 0.5))
    pos = rra(-g["score"][valid], genes[valid], a_high, control_mask=cm, guide_names=guides[valid],
              progress=prog(0.5, 0.95))
    # gene log fold change: median (default) or mean of the valid guides
    lfc_v, gene_v = g["lfc"][valid], genes[valid]
    order = np.argsort(gene_v.astype(str), kind="mergesort")
    gs, ls = gene_v[order].astype(str), lfc_v[order]
    cut = np.flatnonzero(np.r_[True, gs[1:] != gs[:-1], True])
    glfc = {}
    for s, e in zip(cut[:-1], cut[1:]):
        vals = ls[s:e]
        glfc[gs[s]] = float(np.median(vals) if gene_lfc == "median" else np.mean(vals))
    neg_rank = {gname: i + 1 for i, gname in enumerate(neg["gene"][np.argsort(neg["lo"], kind="mergesort")])}
    pos_rank = {gname: i + 1 for i, gname in enumerate(pos["gene"][np.argsort(pos["lo"], kind="mergesort")])}
    pos_ix = {gname: i for i, gname in enumerate(pos["gene"])}
    rows = []
    for i, gname in enumerate(neg["gene"]):
        j = pos_ix[gname]
        rows.append({"gene": gname, "guides": int(neg["n"][i]),
                     "neg_score": float(neg["lo"][i]), "neg_p": float(neg["p"][i]), "neg_fdr": float(neg["fdr"][i]),
                     "neg_rank": neg_rank[gname], "neg_good": int(neg["good"][i]),
                     "pos_score": float(pos["lo"][j]), "pos_p": float(pos["p"][j]), "pos_fdr": float(pos["fdr"][j]),
                     "pos_rank": pos_rank[gname], "pos_good": int(pos["good"][j]),
                     "lfc": glfc.get(gname, 0.0)})
    rows.sort(key=lambda r: r["neg_rank"])
    if progress:
        progress(1.0)
    return {"guide": g, "genes": rows, "alpha_neg": a_low, "alpha_pos": a_high,
            "permutation_rounds": neg["permutation_rounds"]}
