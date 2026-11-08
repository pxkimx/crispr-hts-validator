"""Arrayed (plate) screens: per-plate normalisation, quality checks, scores, replicate summary and hit calls.

Normalisation is done separately for every plate, replicate and channel:
  none        raw values
  poc         percent of control: 100 * x / mean(negative controls on the plate)
  npi         normalised percent inhibition: 100 * (mean(neg) - x) / (mean(neg) - mean(pos))
  median      x / median(sample wells on the plate)                         (cellHTS 'median')
  robust_z    (x - median(samples)) / (1.4826 * MAD(samples))               (robust z-score)
  bscore      two-way median polish of the sample wells removes row and column effects; residual / (1.4826 * MAD)
              (Brideau et al., J Biomol Screen 2003)
  zscore      (x - mean(all wells)) / SD(all wells): the classic plate z-score used for the 2003 KcViab analysis
Scores: z-type methods (robust_z, bscore, zscore) are scores as they are; the others are turned into robust
z-scores over all sample wells of the replicate. 'lower' screens (e.g. viability) flip the sign so that a hit always
has a high score. SSMD: robust SSMD for screens without replicates (Zhang 2011), (x - plate median)/(1.4826 MAD √2).
"""
from __future__ import annotations

import math
import warnings

import numpy as np
import pandas as pd

from .arrayed_io import ArrayedScreen

NORM = {
    "none": "No normalisation: raw values.",
    "poc": "Percent of control: 100 × value ÷ mean of the negative-control wells on the same plate.",
    "npi": "Normalised percent inhibition: 0% = negative controls, 100% = positive controls, on each plate.",
    "median": "Plate median: value ÷ median of the sample wells on the plate (cellHTS 'median').",
    "robust_z": "Robust z-score (default): (value − plate median) ÷ (1.4826 × MAD) of the sample wells on the plate.",
    "bscore": "B-score: median polish removes row and column effects within each plate, then residual ÷ (1.4826 × MAD).",
    "zscore": "Classic plate z-score: (value − plate mean) ÷ plate SD over all wells, as in the 2003 KcViab analysis.",
}
Z_TYPE = {"robust_z", "bscore", "zscore"}

QC_THRESHOLDS = {
    "plate_ratio_warn": 2.0,     # max/min of raw plate medians
    "plate_ratio_fail": 5.0,
    "zprime_pass": 0.5,          # Z' factor (Zhang, Chung & Oldenburg 1999)
    "zprime_fail": 0.0,
    "ssmd_pass": 3.0,            # |SSMD| between positive and negative controls
    "ssmd_fail": 2.0,
    "edge_pass": 0.10,           # |edge / interior - 1| of plate-normalised sample wells
    "edge_fail": 0.20,
    "rowcol_pass": 0.10,         # share of within-plate variance explained by row and column effects
    "rowcol_fail": 0.25,
    "rep_pass": 0.70,            # Pearson r between replicates (normalised sample wells)
    "rep_fail": 0.40,
}
DEFAULT_NORM = "robust_z"

# Confirmed by the user on 2026-11-07: robust z stays the default normalisation (per-plate median and MAD), and a
# well is a hit when its combined score passes 2.0 and every replicate also passes 1.0. On the cellHTS2 KcViab
# screen that gives 668 hits, against 415 for the classic plate z-score the published 2003 analysis used and 1,411
# for B-score, which also removes that screen's row and column effects.
HIT_DEFAULTS = {"score": "z", "threshold": 2.0, "replicate_min": 1.0, "summary": "mean", "direction": "lower",
                "gene_fraction": 0.5}


def _mad(v):
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return float("nan")
    return 1.4826 * float(np.median(np.abs(v - np.median(v))))


def median_polish(m: np.ndarray, iters: int = 10, tol: float = 1e-6) -> tuple[np.ndarray, float, np.ndarray, np.ndarray]:
    """Tukey's median polish (ignoring NaN). Returns residuals, overall, row and column effects."""
    r = m.astype(float).copy()
    overall = 0.0
    re_ = np.zeros(r.shape[0])
    ce = np.zeros(r.shape[1])
    for _ in range(iters):
        rm = np.nanmedian(r, axis=1)
        rm = np.where(np.isfinite(rm), rm, 0.0)
        r -= rm[:, None]
        re_ += rm
        d = np.nanmedian(ce)
        ce -= d
        overall += d
        cm = np.nanmedian(r, axis=0)
        cm = np.where(np.isfinite(cm), cm, 0.0)
        r -= cm[None, :]
        ce += cm
        d = np.nanmedian(re_)
        re_ -= d
        overall += d
        if np.nanmax(np.abs(np.concatenate([rm, cm]))) < tol:
            break
    return r, overall, re_, ce


def _plate_groups(w: pd.DataFrame):
    return w.groupby(["channel", "replicate", "plate"], sort=True).groups


def normalise(scr: ArrayedScreen, method: str = "robust_z", log2: bool = False, neg: str = "neg", pos: str = "pos",
              use_flagged: bool = False) -> pd.Series:
    """Normalised value for every measurement (index aligned with scr.wells). Wells flagged in the screen log become
    NaN unless `use_flagged` (needed to reproduce analyses that ignored the log)."""
    if method not in NORM:
        raise ValueError("Unknown normalisation.")
    w = scr.wells
    x = w["value"].astype(float)
    if not use_flagged:
        x = x.where(~w["flagged"].astype(bool))
    if log2:
        x = np.log2(x.where(x > 0))
    out = pd.Series(np.nan, index=w.index)
    nr, nc = scr.dims
    for key, idx in _plate_groups(w).items():
        v = x.loc[idx].to_numpy()
        cont = w.loc[idx, "content"].to_numpy()
        s = v[cont == "sample"]
        if method == "none":
            y = v
        elif method in ("poc", "npi"):
            n_ = v[cont == neg]
            n_ = n_[np.isfinite(n_)]
            if len(n_) == 0:
                raise ValueError(f"Plate {key[2]} has no negative-control wells (‘{neg}’), which {method.upper()} needs.")
            if method == "poc":
                y = 100.0 * v / n_.mean()
            else:
                p_ = v[cont == pos]
                p_ = p_[np.isfinite(p_)]
                if len(p_) == 0:
                    raise ValueError(f"Plate {key[2]} has no positive-control wells (‘{pos}’), which NPI needs.")
                y = 100.0 * (n_.mean() - v) / (n_.mean() - p_.mean())
        elif method == "median":
            y = v / np.nanmedian(s) if not log2 else v - np.nanmedian(s)
        elif method == "robust_z":
            y = (v - np.nanmedian(s)) / _mad(s)
        elif method == "zscore":
            y = (v - np.nanmean(v)) / np.nanstd(v, ddof=1)
        elif method == "bscore":
            rows = w.loc[idx, "row"].to_numpy()
            cols = w.loc[idx, "col"].to_numpy()
            M = np.full((nr, nc), np.nan)
            samp = cont == "sample"
            M[rows[samp], cols[samp]] = v[samp]
            _, overall, re_, ce = median_polish(M)
            resid = v - (overall + re_[rows] + ce[cols])
            y = resid / _mad(resid[samp])
        out.loc[idx] = y
    return out


def scores(scr: ArrayedScreen, normed: pd.Series, method: str, direction: str = "lower", kind: str = "z") -> pd.Series:
    """Per-measurement score (higher = stronger hit in the chosen direction)."""
    w = scr.wells
    sign = -1.0 if direction == "lower" else 1.0
    out = pd.Series(np.nan, index=w.index)
    if kind == "ssmd":
        for key, idx in _plate_groups(w).items():
            v = normed.loc[idx].to_numpy()
            s = v[w.loc[idx, "content"].to_numpy() == "sample"]
            out.loc[idx] = (v - np.nanmedian(s)) / (_mad(s) * math.sqrt(2))
    elif method in Z_TYPE:
        out = normed.copy()
    else:
        for key, idx in w.groupby(["channel", "replicate"]).groups.items():
            v = normed.loc[idx]
            s = v[w.loc[idx, "content"] == "sample"].to_numpy()
            out.loc[idx] = (v - np.nanmedian(s)) / _mad(s)
    if method == "npi" and kind != "ssmd":
        out = -out            # NPI already counts inhibition upwards; turn it back to 'signal' orientation
    if direction == "both":
        return out
    return sign * out


def well_table(scr: ArrayedScreen, normed: pd.Series, sc: pd.Series, channel: int = 1) -> pd.DataFrame:
    """One row per plate x well: annotation plus normalised value and score for each replicate."""
    w = scr.wells.assign(norm=normed, score=sc)
    w = w[w["channel"] == channel]
    base = w.drop_duplicates(["plate", "well"])[["plate", "well", "row", "col", "content", "reagent", "gene"]]
    reps = sorted(w["replicate"].unique())
    t = base.set_index(["plate", "well"])
    for r in reps:
        sub = w[w["replicate"] == r].set_index(["plate", "well"])
        t[f"raw_{r}"] = sub["value"]
        t[f"norm_{r}"] = sub["norm"]
        t[f"score_{r}"] = sub["score"]
        t[f"flag_{r}"] = sub["flagged"]
    return t.reset_index(), [int(r) for r in reps]


def call_hits(t: pd.DataFrame, reps: list[int], threshold: float, replicate_min: float, summary: str,
              direction: str) -> pd.DataFrame:
    S = t[[f"score_{r}" for r in reps]].to_numpy(float)
    if summary == "median":
        summ = np.nanmedian(S, axis=1)
    elif summary == "min":
        # the replicate closest to zero: a well must score in every replicate
        summ = np.where(np.all(np.isnan(S), axis=1), np.nan,
                        np.take_along_axis(S, np.nanargmin(np.where(np.isnan(S), np.inf, np.abs(S)), axis=1)[:, None], axis=1)[:, 0])
    elif summary == "max":
        summ = np.nanmax(S, axis=1)
    else:
        summ = np.nanmean(S, axis=1)
    t = t.copy()
    t["score"] = summ
    t["n_reps"] = np.isfinite(S).sum(axis=1)
    if direction == "both":
        agree = np.all(np.where(np.isfinite(S), np.abs(S) >= replicate_min, True), axis=1) & \
                (np.all(np.where(np.isfinite(S), np.sign(S) == np.sign(summ)[:, None], True), axis=1))
        hit = (np.abs(summ) >= threshold) & agree
    else:
        agree = np.all(np.where(np.isfinite(S), S >= replicate_min, True), axis=1)
        hit = (summ >= threshold) & agree
    t["agree"] = agree
    t["hit"] = hit & (t["content"] == "sample") & (t["n_reps"] > 0)
    if len(reps) > 1:
        t["rep_sd"] = np.nanstd(S, axis=1, ddof=1)
    return t


def gene_calls(t: pd.DataFrame, fraction: float = 0.5) -> pd.DataFrame:
    s = t[(t["content"] == "sample") & (t["gene"] != "")]
    if s.empty:
        return pd.DataFrame(columns=["gene", "wells", "hit_wells", "median_score", "best_score", "hit"])
    g = s.groupby("gene").agg(wells=("well", "size"), hit_wells=("hit", "sum"), median_score=("score", "median"),
                              best_score=("score", "max")).reset_index()
    g["hit"] = (g["hit_wells"] >= np.maximum(1, np.ceil(g["wells"] * fraction))).astype(bool)
    return g.sort_values(["hit", "median_score"], ascending=[False, False]).reset_index(drop=True)


# ----------------------------------------------------------------------------- quality checks
def _st(v, p, f, higher=True):
    if v is None or not np.isfinite(v):
        return "na"
    if higher:
        return "pass" if v >= p else ("fail" if v < f else "warn")
    return "pass" if v <= p else ("fail" if v > f else "warn")


def _worst(states):
    order = {"fail": 3, "warn": 2, "pass": 1, "na": 0}
    states = list(states)
    return max(states, key=lambda s: order[s]) if states else "na"


def zprime(pos, neg) -> float:
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    pos, neg = pos[np.isfinite(pos)], neg[np.isfinite(neg)]
    if len(pos) < 2 or len(neg) < 2:
        return float("nan")
    d = abs(pos.mean() - neg.mean())
    return float(1 - 3 * (pos.std(ddof=1) + neg.std(ddof=1)) / d) if d > 0 else float("-inf")


def ssmd(pos, neg) -> float:
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    pos, neg = pos[np.isfinite(pos)], neg[np.isfinite(neg)]
    if len(pos) < 2 or len(neg) < 2:
        return float("nan")
    return float((pos.mean() - neg.mean()) / math.sqrt(pos.var(ddof=1) + neg.var(ddof=1)))


def run_qc(scr: ArrayedScreen, channel: int = 1, neg: str = "neg", pos: str = "pos", thresholds: dict | None = None) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return _run_qc(scr, channel, neg, pos, thresholds)


def _run_qc(scr, channel, neg, pos, thresholds):
    th = {**QC_THRESHOLDS, **(thresholds or {})}
    w = scr.wells[scr.wells["channel"] == channel].copy()
    w["v"] = w["value"].astype(float).where(~w["flagged"].astype(bool))
    nr, nc = scr.dims
    checks = []
    samp = w["content"] == "sample"
    reps = sorted(w["replicate"].unique())

    # 1 plate-to-plate medians
    med = w[samp].groupby(["replicate", "plate"])["v"].median()
    ratio = float(med.max() / med.min()) if med.min() > 0 else float("inf")
    st = "pass" if ratio < th["plate_ratio_warn"] else ("fail" if ratio > th["plate_ratio_fail"] else "warn")
    checks.append({"id": "plates", "name": "Plate-to-plate medians", "status": st, "value": ratio,
                   "per_plate": [{"replicate": int(r), "plate": int(p), "value": float(v)} for (r, p), v in med.items()],
                   "why": (f"Raw plate medians vary {ratio:.1f}-fold (from {med.min():,.0f} to {med.max():,.0f})"
                           + (", which per-plate normalisation removes; never compare raw values across plates."
                              if st != "pass" else "; plates are comparable even before normalisation.")),
                   "rule": f"max/min of raw sample-well medians: pass < {th['plate_ratio_warn']:g}×, fail > {th['plate_ratio_fail']:g}×."})

    # 2 controls: Z' and SSMD per plate when possible, otherwise pooled across plates (after plate-median scaling)
    pm = w[samp].groupby(["replicate", "plate"])["v"].median()
    w["pn"] = w["v"] / pd.Series([pm.get((r, p), np.nan) for r, p in zip(w["replicate"], w["plate"])], index=w.index)
    cnt = w.groupby(["replicate", "plate"])["content"].agg(lambda c: ((c == pos).sum(), (c == neg).sum()))
    min_pos = min((a for a, b in cnt), default=0)
    min_neg = min((b for a, b in cnt), default=0)
    contents = sorted(set(w["content"]) - {"sample"})
    if min_pos == 0 and min_neg == 0 and not ((w["content"] == pos).any() or (w["content"] == neg).any()):
        checks.append({"id": "controls", "name": "Control separation (Z′, SSMD)", "status": "na", "per_plate": [],
                       "why": f"No wells are marked ‘{pos}’ and ‘{neg}’. Mark the positive and negative controls on the "
                              f"Plate map page (well types present: {', '.join(contents) or 'samples only'}).",
                       "rule": ""})
    else:
        per = []
        if min_pos >= 2 and min_neg >= 2:
            for (r, p), g in w.groupby(["replicate", "plate"]):
                z = zprime(g.loc[g["content"] == pos, "v"], g.loc[g["content"] == neg, "v"])
                s = ssmd(g.loc[g["content"] == pos, "v"], g.loc[g["content"] == neg, "v"])
                per.append({"replicate": int(r), "plate": int(p), "zprime": z, "ssmd": s,
                            "status": _worst([_st(z, th["zprime_pass"], th["zprime_fail"]),
                                              _st(abs(s), th["ssmd_pass"], th["ssmd_fail"])])})
            n_bad = sum(1 for q in per if q["status"] == "fail")
            zs = [q["zprime"] for q in per]
            why = (f"Z′ per plate from {min(zs):.2f} to {max(zs):.2f} (median {np.median(zs):.2f}); "
                   + (f"{n_bad} plate(s) fail." if n_bad else "controls separate well on every plate."))
            st = _worst(q["status"] for q in per)
            pooled = None
        else:
            pp = w.loc[w["content"] == pos, "pn"]
            nn = w.loc[w["content"] == neg, "pn"]
            z = zprime(pp, nn)
            s = ssmd(pp, nn)
            st = _worst([_st(z, th["zprime_pass"], th["zprime_fail"]), _st(abs(s), th["ssmd_pass"], th["ssmd_fail"])])
            pooled = {"zprime": z, "ssmd": s, "n_pos": int(pp.notna().sum()), "n_neg": int(nn.notna().sum()),
                      "pos_mean": float(pp.mean()), "neg_mean": float(nn.mean())}
            why = (f"Each plate has only {min_pos} positive and {min_neg} negative control well(s), so Z′ and SSMD "
                   f"cannot be computed per plate. Pooled over all plates (after scaling each plate to its median): "
                   f"Z′ = {z:.2f}, SSMD = {s:.1f} from {pooled['n_pos']} positive and {pooled['n_neg']} negative wells; "
                   f"positive controls sit at {pooled['pos_mean']:.0%} and negative at {pooled['neg_mean']:.0%} of the plate median. "
                   "Pooled values include plate-to-plate differences, so they understate how well controls separate within a plate.")
        checks.append({"id": "controls", "name": "Control separation (Z′, SSMD)", "status": st, "per_plate": per,
                       "pooled": pooled, "why": why,
                       "rule": f"Z′: pass ≥ {th['zprime_pass']:g}, fail < {th['zprime_fail']:g}; |SSMD|: pass ≥ {th['ssmd_pass']:g}, fail < {th['ssmd_fail']:g} "
                               "(Zhang et al. 1999; Zhang 2007). Needs ≥ 2 wells of each control per plate, otherwise pooled."})

    # 3 edge effect and 4 row/column effects (plate-median scaled sample wells)
    edge = (w["row"] == 0) | (w["row"] == nr - 1) | (w["col"] == 0) | (w["col"] == nc - 1)
    er = []
    rc = []
    maps = {}
    for (r, p), g in w[samp].groupby(["replicate", "plate"]):
        e = g.loc[edge[g.index], "pn"].median()
        i = g.loc[~edge[g.index], "pn"].median()
        er.append(float(e / i - 1) if i else float("nan"))
        M = np.full((nr, nc), np.nan)
        M[g["row"].to_numpy(), g["col"].to_numpy()] = g["pn"].to_numpy()
        resid, _, _, _ = median_polish(M)
        tot = np.nanvar(M)
        rc.append(float(1 - np.nanvar(resid) / tot) if tot > 0 else 0.0)
    er_med = float(np.nanmedian(er)) if er else float("nan")
    rc_med = float(np.nanmedian(rc)) if rc else float("nan")
    est = _st(abs(er_med), th["edge_pass"], th["edge_fail"], False)
    checks.append({"id": "edge", "name": "Edge effect", "status": est, "value": er_med, "per_plate": er,
                   "why": (f"Edge wells read {abs(er_med):.0%} {'higher' if er_med > 0 else 'lower'} than interior wells "
                           f"(median over plates)" + (" — within the normal range." if est == "pass" else
                                                       "; evaporation or temperature gradients are likely. B-score normalisation or "
                                                       "excluding the outer wells reduces this.")),
                   "rule": f"|edge ÷ interior − 1| of plate-scaled sample wells: pass ≤ {th['edge_pass']:.0%}, fail > {th['edge_fail']:.0%}."})
    rst = _st(rc_med, th["rowcol_pass"], th["rowcol_fail"], False)
    checks.append({"id": "rowcol", "name": "Row and column effects", "status": rst, "value": rc_med, "per_plate": rc,
                   "why": (f"Row and column patterns explain {rc_med:.0%} of the variation within a plate (median over plates)"
                           + (". Little positional bias." if rst == "pass" else
                              ": dispensing or reader patterns. B-score normalisation removes them.")),
                   "rule": f"share of within-plate variance removed by a median polish of rows and columns: pass ≤ {th['rowcol_pass']:.0%}, fail > {th['rowcol_fail']:.0%}."})

    # 5 replicate agreement
    if len(reps) >= 2:
        piv = w[samp].pivot_table(index=["plate", "well"], columns="replicate", values="pn")
        pairs = []
        for i in range(len(reps)):
            for j in range(i + 1, len(reps)):
                a, b = piv[reps[i]], piv[reps[j]]
                ok = a.notna() & b.notna()
                r = float(np.corrcoef(a[ok], b[ok])[0, 1]) if ok.sum() > 2 else float("nan")
                from scipy.stats import spearmanr
                rho = float(spearmanr(a[ok], b[ok]).correlation) if ok.sum() > 2 else float("nan")
                per_plate = piv[ok].groupby(level=0).apply(lambda d: float(np.corrcoef(d[reps[i]], d[reps[j]])[0, 1]) if len(d) > 2 else np.nan)
                pairs.append({"a": int(reps[i]), "b": int(reps[j]), "value": r, "spearman": rho,
                              "status": _st(r, th["rep_pass"], th["rep_fail"]),
                              "worst_plates": [{"plate": int(k), "r": float(v)} for k, v in per_plate.sort_values().head(5).items()]})
        st = _worst(p["status"] for p in pairs)
        checks.append({"id": "replicates", "name": "Replicate agreement", "status": st, "pairs": pairs,
                       "why": "; ".join(f"replicate {p['a']} vs {p['b']}: r = {p['value']:.2f} (Spearman {p['spearman']:.2f}), "
                                        f"weakest plate {p['worst_plates'][0]['plate']} (r = {p['worst_plates'][0]['r']:.2f})" for p in pairs)
                              + ("." if st == "pass" else " — replicates disagree; hits need to be confirmed in both."),
                       "rule": f"Pearson r of plate-scaled sample wells: pass ≥ {th['rep_pass']:g}, fail < {th['rep_fail']:g}."})
    else:
        checks.append({"id": "replicates", "name": "Replicate agreement", "status": "na", "pairs": [],
                       "why": "Only one replicate, so hits cannot be confirmed by a second measurement.", "rule": ""})

    # heatmap data: raw values per plate/replicate
    for (r, p), g in w.groupby(["replicate", "plate"]):
        M = np.full((nr, nc), np.nan)
        M[g["row"].to_numpy(), g["col"].to_numpy()] = g["v"].to_numpy()
        maps[f"{int(r)}:{int(p)}"] = M.round(4).tolist()
    return {"checks": checks, "overall": _worst(c["status"] for c in checks), "thresholds": th, "raw_maps": maps,
            "dims": [nr, nc], "plates": sorted(int(x) for x in w["plate"].unique()), "replicates": [int(x) for x in reps]}


def analyse(scr: ArrayedScreen, method: str = "robust_z", channel: int = 1, log2: bool = False, neg: str = "neg",
            pos: str = "pos", score: str = "z", threshold: float = 2.0, replicate_min: float = 1.0,
            summary: str = "mean", direction: str = "lower", gene_fraction: float = 0.5, use_flagged: bool = False) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return _analyse(scr, method, channel, log2, neg, pos, score, threshold, replicate_min, summary, direction,
                        gene_fraction, use_flagged)


def _analyse(scr, method, channel, log2, neg, pos, score, threshold, replicate_min, summary, direction, gene_fraction,
             use_flagged):
    n = normalise(scr, method, log2, neg, pos, use_flagged)
    sc = scores(scr, n, method, direction, score)
    t, reps = well_table(scr, n, sc, channel)
    t = call_hits(t, reps, threshold, replicate_min, summary, direction)
    genes = gene_calls(t, gene_fraction)
    return {"wells": t, "replicates": reps, "genes": genes,
            "settings": {"method": method, "channel": channel, "log2": log2, "neg": neg, "pos": pos, "score": score,
                         "threshold": threshold, "replicate_min": replicate_min, "summary": summary,
                         "direction": direction, "gene_fraction": gene_fraction, "use_flagged": use_flagged}}
