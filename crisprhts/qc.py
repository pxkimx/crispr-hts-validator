"""Quality checks for pooled screens. Each check returns pass / warn / fail with one plain sentence why.

Thresholds are scientific choices: the defaults below come from published guidance (MAGeCK/MAGeCK-VISPR QC
recommendations, Hart et al. 2015/2017 for essential-gene separation, Joung et al. 2017 Nat Protoc for library
skew and representation) and are shown in the app, where they can be changed and are then remembered.

Confirmed by the user on 2026-11-07: read depth 300/100 reads per guide; mapped reads 80%/65%; Gini 0.10/0.20
at baseline; zero-count guides 1%/5% in the reference; replicate r 0.60/0.40 (Pearson, with Spearman shown
beside it); representation 500x/200x cells per guide; MOI 0.4/0.6; 4 guides per gene; 100 non-targeting
controls. The arrayed defaults they confirmed are in arrayed_stats.py.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parent / "data"

THRESHOLDS = {
    "coverage_pass": 300.0,    # reads per guide in a sample
    "coverage_fail": 100.0,
    "mapped_pass": 0.80,       # fraction of reads matching a library guide (common screen guidance: 80-90%)
    "mapped_fail": 0.65,
    "zero_ref_pass": 0.01,     # fraction of zero-count guides in reference (plasmid / T0) samples
    "zero_ref_fail": 0.05,
    "zero_late_pass": 0.05,    # ... in later samples, where essential-gene guides are expected to drop out
    "zero_late_fail": 0.15,
    "gini_ref_pass": 0.10,     # MAGeCK Gini index (on log counts), reference samples
    "gini_ref_fail": 0.20,
    "gini_late_pass": 0.20,
    "gini_late_fail": 0.30,
    "skew_pass": 10.0,         # 90th / 10th percentile of reference counts
    "skew_fail": 20.0,
    "rep_r_pass": 0.60,        # Pearson r of guide log fold changes between replicates
    "rep_r_fail": 0.40,
    "auc_pass": 0.90,          # essential (CEGv2) vs non-essential (NEGv1) separation
    "auc_fail": 0.80,
    "cells_pass": 500.0,       # cells carrying each guide at the tightest point of the experiment
    "cells_fail": 200.0,
    "moi_pass_hi": 0.4,        # multiplicity of infection: one guide per cell
    "moi_fail_hi": 0.6,
    "moi_low": 0.05,           # below this, most cells are untransduced (wasteful, not wrong)
    "guides_gene_pass": 4.0,   # independent guides per gene
    "guides_gene_fail": 2.0,
    "ntc_pass": 100.0,         # non-targeting control guides in the library
    "ntc_fail": 20.0,
    "ntc_shift": 0.25,         # |median LFC| of the non-targeting guides: they should sit at zero
}

THRESHOLD_TEXT = {
    "coverage": "Reads per guide: pass ≥ {coverage_pass:g}×, warn {coverage_fail:g}–{coverage_pass:g}×, fail < {coverage_fail:g}×.",
    "mapped": "Mapped reads: pass ≥ {mapped_pass:.0%}, warn {mapped_fail:.0%}–{mapped_pass:.0%}, fail < {mapped_fail:.0%} (MAGeCK QC guideline).",
    "zero": "Zero-count guides: reference samples pass ≤ {zero_ref_pass:.0%}, fail > {zero_ref_fail:.0%}; later samples pass ≤ {zero_late_pass:.0%}, fail > {zero_late_fail:.0%}.",
    "gini": "Gini index (MAGeCK, on log counts): reference pass ≤ {gini_ref_pass:g}, fail > {gini_ref_fail:g}; later samples pass ≤ {gini_late_pass:g}, fail > {gini_late_fail:g}.",
    "skew": "Library evenness (90th/10th percentile of reference counts): pass ≤ {skew_pass:g}, fail > {skew_fail:g} (Joung et al. 2017).",
    "rep": "Replicate agreement (Pearson r of guide log fold changes): pass ≥ {rep_r_pass:g}, fail < {rep_r_fail:g}.",
    "auc": "Essential vs non-essential separation (ROC AUC, CEGv2 vs NEGv1): pass ≥ {auc_pass:g}, warn {auc_fail:g}–{auc_pass:g}, fail < {auc_fail:g}.",
    "cells": "Library representation: pass ≥ {cells_pass:g} cells per guide at every step, fail < {cells_fail:g}×. MOI: pass ≤ {moi_pass_hi:g}, fail > {moi_fail_hi:g} (one guide per cell).",
    "guides_gene": "Guides per gene: pass ≥ {guides_gene_pass:g}, fail < {guides_gene_fail:g} (gene-level tests need several independent guides).",
    "ntc": "Non-targeting controls: pass ≥ {ntc_pass:g} guides, fail < {ntc_fail:g}; their median log₂ fold change should be within ±{ntc_shift:g} of zero.",
}


def threshold_text(th: dict | None = None) -> dict:
    t = {**THRESHOLDS, **(th or {})}
    return {k: v.format(**t) for k, v in THRESHOLD_TEXT.items()}


# ------------------------------------------------------------------------------- metrics
def gini(x) -> float:
    """MAGeCK's Gini index (mageckcount_gini) of log(count + 1). 0 = perfectly even, 1 = all reads in one guide."""
    v = np.sort(np.log(np.asarray(x, dtype=float) + 1.0))
    n = len(v)
    if n < 2:
        return 0.0
    s = v.sum() or 1.0
    return float(1.0 - 2.0 * (n - ((np.arange(n) + 1.0) * v).sum() / s) / (n - 1))


def gini_raw(x) -> float:
    """Classical Gini index on raw counts (reported alongside, e.g. for comparison with BAGEL's QC)."""
    v = np.sort(np.asarray(x, dtype=float))
    n = len(v)
    if n < 2 or v.sum() == 0:
        return 0.0
    return float((2 * ((np.arange(n) + 1) * v).sum()) / (n * v.sum()) - (n + 1) / n)


def skew_ratio(x) -> float:
    v = np.asarray(x, dtype=float)
    p10, p90 = np.percentile(v, [10, 90])
    return float("inf") if p10 <= 0 else float(p90 / p10)


def auc(pos_scores, neg_scores) -> float:
    """Probability that a positive (essential) gene scores lower than a negative one (Mann-Whitney, ties = 1/2)."""
    a = np.asarray(pos_scores, dtype=float)
    b = np.asarray(neg_scores, dtype=float)
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    allv = np.concatenate([a, b])
    from scipy.stats import rankdata
    r = rankdata(-allv)              # lower LFC = higher rank value
    ra = r[:len(a)].sum()
    return float((ra - len(a) * (len(a) + 1) / 2) / (len(a) * len(b)))


def roc_pr(pos_scores, neg_scores, points: int = 200) -> dict:
    """ROC and precision-recall curves, ranking genes from most depleted (lowest score) upward."""
    a = np.asarray(pos_scores, dtype=float)
    b = np.asarray(neg_scores, dtype=float)
    s = np.concatenate([a, b])
    y = np.concatenate([np.ones(len(a)), np.zeros(len(b))])
    o = np.argsort(s, kind="mergesort")
    y = y[o]
    tp = np.cumsum(y)
    fp = np.cumsum(1 - y)
    tpr = tp / max(len(a), 1)
    fpr = fp / max(len(b), 1)
    prec = tp / np.maximum(tp + fp, 1)
    idx = np.unique(np.linspace(0, len(s) - 1, min(points, len(s))).astype(int))
    # recall at 5% FDR (precision >= 0.95)
    ok = np.where(prec >= 0.95)[0]
    rec95 = float(tpr[ok].max()) if len(ok) else 0.0
    return {"fpr": [0.0] + fpr[idx].round(4).tolist(), "tpr": [0.0] + tpr[idx].round(4).tolist(),
            "recall": tpr[idx].round(4).tolist(), "precision": prec[idx].round(4).tolist(), "recall_at_fdr5": rec95}


def load_reference_sets(species: str = "human") -> tuple[set[str], set[str]]:
    """Hart lab reference genes: human CEGv2 / NEGv1, mouse CEG / NEG (BAGEL repository, MIT)."""
    files = {"human": ("CEGv2.txt", "NEGv1.txt"), "mouse": ("CEG_mouse.txt", "NEG_mouse.txt")}
    if species not in files:
        raise ValueError("Reference essential genes are bundled for human and mouse only.")
    def read(p):
        out = set()
        for i, ln in enumerate(p.read_text().splitlines()):
            f = ln.split("\t")[0].strip()
            if i == 0 or not f:
                continue
            out.add(f)
        return out
    return read(DATA / files[species][0]), read(DATA / files[species][1])


def pearson(x, y) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3 or np.std(x[ok]) == 0 or np.std(y[ok]) == 0:
        return float("nan")
    return float(np.corrcoef(x[ok], y[ok])[0, 1])


def _status(value, pass_at, fail_at, higher_is_better=True) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "na"
    if higher_is_better:
        return "pass" if value >= pass_at else ("fail" if value < fail_at else "warn")
    return "pass" if value <= pass_at else ("fail" if value > fail_at else "warn")


def _worst(states) -> str:
    order = {"fail": 3, "warn": 2, "pass": 1, "na": 0}
    states = list(states)
    return max(states, key=lambda s: order[s]) if states else "na"


# ------------------------------------------------------------------------------- the checks
def median_ratio_normalise(counts: np.ndarray) -> np.ndarray:
    from .pooled_stats import size_factors
    f, _ = size_factors(counts, "median")
    return counts * f[None, :]


def guide_lfc_matrix(counts: np.ndarray, samples: list[str], ref_samples: list[str]) -> np.ndarray:
    """log2 fold change of every sample against the mean of the reference samples (median-ratio normalised,
    pseudocount 1)."""
    N = median_ratio_normalise(counts)
    ri = [samples.index(s) for s in ref_samples]
    ref = N[:, ri].mean(axis=1)
    return np.log2(N + 1.0) - np.log2(ref + 1.0)[:, None]


def run_qc(counts: np.ndarray, guides: list[str], genes: list[str], samples: list[str], design: dict,
           reads: dict | None = None, mapped: dict | None = None, thresholds: dict | None = None,
           species: str = "auto", screen: dict | None = None, control_mask=None) -> dict:
    """All pooled-screen checks. `design` = {"samples": [{sample, condition, reference, replicate}], ...}."""
    th = {**THRESHOLDS, **(thresholds or {})}
    X = np.asarray(counts, dtype=float)
    ng = X.shape[0]
    info = {r["sample"]: r for r in design.get("samples", [])}
    is_ref = {s: bool(info.get(s, {}).get("reference")) for s in samples}
    refs = [s for s in samples if is_ref[s]]
    checks = []

    # 0 library representation (cells per guide) and MOI -- from the numbers the user enters, not from the reads
    checks.append(_representation_check(screen or {}, ng, th))

    # 1 reads per sample
    tot = {s: float(reads[s]) if reads and s in reads else float(X[:, j].sum()) for j, s in enumerate(samples)}
    per = []
    for j, s in enumerate(samples):
        cov = tot[s] / ng
        st = _status(cov, th["coverage_pass"], th["coverage_fail"])
        per.append({"sample": s, "value": tot[s], "coverage": cov, "status": st})
    low = [p["sample"] for p in per if p["status"] != "pass"]
    checks.append({
        "id": "reads", "name": "Reads per sample", "status": _worst(p["status"] for p in per), "per_sample": per,
        "why": ("Every sample has at least {:g} reads per guide.".format(th["coverage_pass"]) if not low else
                f"{len(low)} sample(s) have fewer than {th['coverage_pass']:g} reads per guide ({', '.join(low[:4])}"
                f"{'…' if len(low) > 4 else ''}); low depth makes guide counts noisy."),
        "rule": threshold_text(th)["coverage"]})

    # 2 mapped fraction
    if reads and mapped:
        per = []
        for s in samples:
            fr = mapped.get(s, 0) / reads[s] if reads.get(s) else float("nan")
            per.append({"sample": s, "value": fr, "status": _status(fr, th["mapped_pass"], th["mapped_fail"])})
        bad = [p for p in per if p["status"] != "pass"]
        checks.append({"id": "mapped", "name": "Mapped reads", "status": _worst(p["status"] for p in per),
                       "per_sample": per,
                       "why": ("At least {:.0%} of reads match a library guide in every sample.".format(th["mapped_pass"])
                               if not bad else f"{len(bad)} sample(s) map poorly (lowest {min(p['value'] for p in bad):.0%}): "
                               "check the guide position in the read, the library file, or contamination."),
                       "rule": threshold_text(th)["mapped"]})
    else:
        checks.append({"id": "mapped", "name": "Mapped reads", "status": "na", "per_sample": [],
                       "why": "Not available: the data came in as a count table, so the unmapped reads are not known.",
                       "rule": threshold_text(th)["mapped"]})

    # 3 zero-count guides, 4 Gini
    zper, gper = [], []
    for j, s in enumerate(samples):
        z = float((X[:, j] == 0).mean())
        g = gini(X[:, j])
        if is_ref[s]:
            zs = _status(z, th["zero_ref_pass"], th["zero_ref_fail"], False)
            gs = _status(g, th["gini_ref_pass"], th["gini_ref_fail"], False)
        else:
            zs = _status(z, th["zero_late_pass"], th["zero_late_fail"], False)
            gs = _status(g, th["gini_late_pass"], th["gini_late_fail"], False)
        zper.append({"sample": s, "value": z, "count": int((X[:, j] == 0).sum()), "status": zs, "reference": is_ref[s]})
        gper.append({"sample": s, "value": g, "raw": gini_raw(X[:, j]), "status": gs, "reference": is_ref[s]})
    zbad = [p for p in zper if p["status"] != "pass"]
    checks.append({"id": "zeros", "name": "Zero-count guides", "status": _worst(p["status"] for p in zper),
                   "per_sample": zper,
                   "why": ("Few guides are missing in any sample." if not zbad else
                           "; ".join(f"{p['sample']}: {p['count']:,} guides ({p['value']:.1%}) have no reads" for p in zbad[:3])
                           + (" — in the reference this means the library was not fully represented." if any(p["reference"] for p in zbad)
                              else " — some drop-out is expected late in a screen, but this much suggests a bottleneck.")),
                   "rule": threshold_text(th)["zero"]})
    gbad = [p for p in gper if p["status"] != "pass"]
    checks.append({"id": "gini", "name": "Read evenness (Gini index)", "status": _worst(p["status"] for p in gper),
                   "per_sample": gper,
                   "why": ("Reads are spread evenly over the guides." if not gbad else
                           f"{', '.join(p['sample'] for p in gbad[:4])}: reads are concentrated on fewer guides "
                           "(Gini {}), which happens with strong selection or a population bottleneck.".format(
                               ", ".join(f"{p['value']:.2f}" for p in gbad[:4]))),
                   "rule": threshold_text(th)["gini"]})

    # 5 guide coverage / evenness in the reference
    if refs:
        per = []
        for s in refs:
            v = X[:, samples.index(s)]
            sk = skew_ratio(v)
            det = float((v > 0).mean())
            per.append({"sample": s, "value": sk, "detected": det,
                        "status": _status(sk, th["skew_pass"], th["skew_fail"], False)})
        checks.append({"id": "coverage", "name": "Guide coverage in the reference", "status": _worst(p["status"] for p in per),
                       "per_sample": per,
                       "why": "; ".join(f"{p['sample']}: {p['detected']:.1%} of guides detected, 90/10 skew {p['value']:.1f}"
                                         for p in per) + (". The library is evenly represented." if all(p["status"] == "pass" for p in per)
                                                          else ". An uneven starting library lowers power for the guides at the bottom."),
                       "rule": threshold_text(th)["skew"]})
    else:
        checks.append({"id": "coverage", "name": "Guide coverage in the reference", "status": "na", "per_sample": [],
                       "why": "No reference (plasmid or T0) sample is marked in the sample sheet.",
                       "rule": threshold_text(th)["skew"]})

    # 5b guides per gene: gene-level tests need several independent guides
    gene_arr0 = np.asarray(genes, dtype=object).astype(str)
    ctrl = np.zeros(ng, dtype=bool) if control_mask is None else np.asarray(control_mask, dtype=bool)
    tgt = ~ctrl
    if tgt.any():
        uniq, cnt = np.unique(gene_arr0[tgt], return_counts=True)
        med = float(np.median(cnt))
        thin = uniq[cnt < th["guides_gene_pass"]]
        st = _status(med, th["guides_gene_pass"], th["guides_gene_fail"])
        frac_thin = float(len(thin)) / max(len(uniq), 1)
        if st == "pass" and frac_thin > 0.10:
            st = "warn"
        hist = [{"n": int(k), "genes": int(v)} for k, v in zip(*np.unique(cnt, return_counts=True))][:24]
        head = "Median {:g} guides per gene (range {}–{}).".format(med, int(cnt.min()), int(cnt.max()))
        if not len(thin):
            why = head + " Every gene is covered by several independent guides."
        elif st == "pass":
            why = (head + " {:,} of {:,} genes ({:.1%}) have fewer than {:g} — read those gene calls with care, "
                   "since one bad or off-target guide carries more weight there.".format(
                       len(thin), len(uniq), frac_thin, th["guides_gene_pass"]))
        else:
            why = (head + " {:,} of {:,} genes ({:.1%}) have fewer than {:g} (e.g. {}). A gene called from one or "
                   "two guides cannot be told apart from a single bad or off-target guide.".format(
                       len(thin), len(uniq), frac_thin, th["guides_gene_pass"], ", ".join(sorted(thin)[:4])))
        checks.append({"id": "guides_gene", "name": "Guides per gene", "status": st, "median": med,
                       "n_genes": int(len(uniq)), "n_thin": int(len(thin)), "hist": hist,
                       "examples": sorted(thin)[:20], "min": int(cnt.min()), "max": int(cnt.max()),
                       "why": why, "rule": threshold_text(th)["guides_gene"]})
    else:
        checks.append({"id": "guides_gene", "name": "Guides per gene", "status": "na", "hist": [],
                       "why": "Every guide in the library is marked as a control.",
                       "rule": threshold_text(th)["guides_gene"]})

    # 6 replicate correlation (guide LFC vs reference)
    conds: dict[str, list[str]] = {}
    for s in samples:
        conds.setdefault(info.get(s, {}).get("condition", s), []).append(s)
    lfc = guide_lfc_matrix(X, samples, refs) if refs else None
    logn = np.log2(median_ratio_normalise(X) + 1.0)
    pairs = []
    for c, ss in conds.items():
        if len(ss) < 2 or all(is_ref[s] for s in ss):
            continue
        for i in range(len(ss)):
            for j in range(i + 1, len(ss)):
                a, b = samples.index(ss[i]), samples.index(ss[j])
                r_lfc = pearson(lfc[:, a], lfc[:, b]) if lfc is not None else float("nan")
                r_cnt = pearson(logn[:, a], logn[:, b])
                rho = spearman(lfc[:, a], lfc[:, b]) if lfc is not None else float("nan")
                pairs.append({"condition": c, "a": ss[i], "b": ss[j], "value": r_lfc, "r_counts": r_cnt,
                              "spearman": rho,
                              "status": _status(r_lfc, th["rep_r_pass"], th["rep_r_fail"])})
    if pairs:
        bad = [p for p in pairs if p["status"] != "pass"]
        checks.append({"id": "replicates", "name": "Replicate agreement", "status": _worst(p["status"] for p in pairs),
                       "pairs": pairs,
                       "why": ("Replicates agree (guide fold-change r from {:.2f} to {:.2f}).".format(
                           min(p["value"] for p in pairs), max(p["value"] for p in pairs)) if not bad else
                           "Some replicates agree poorly ({}): their fold changes are dominated by noise.".format(
                               ", ".join(f"{p['a']}–{p['b']} r={p['value']:.2f}" for p in bad[:3]))),
                       "rule": threshold_text(th)["rep"]})
    else:
        checks.append({"id": "replicates", "name": "Replicate agreement", "status": "na", "pairs": [],
                       "why": "No condition has two or more replicates" + ("" if refs else " (or no reference is marked)") + ".",
                       "rule": threshold_text(th)["rep"]})

    # 7 sample PCA on log normalised counts
    pca = sample_pca(logn, samples)
    cond_of = {s: info.get(s, {}).get("condition", s) for s in samples}
    if len(samples) >= 3 and any(len(v) > 1 for v in conds.values()):
        D = np.corrcoef(logn.T)
        lonely, alike = [], set()
        for j, s in enumerate(samples):
            own = [samples.index(o) for o in conds[cond_of[s]] if o != s]
            if not own:
                continue
            m_own = D[j, own].mean()
            for c2, ss2 in conds.items():
                if c2 == cond_of[s]:
                    continue
                m_other = D[j, [samples.index(o) for o in ss2]].mean()
                if m_other > m_own:
                    lonely.append((s, c2))
        by_cond: dict[str, set] = {}
        for s, c2 in lonely:
            by_cond.setdefault(cond_of[s], set()).add(s)
        swaps = []
        for c, ss in by_cond.items():
            if len(ss) == len(conds[c]):      # the whole condition overlaps another: conditions look alike
                partner = sorted({c2 for s, c2 in lonely if cond_of[s] == c})
                alike.add(tuple(sorted([c] + partner[:1])))
            else:
                swaps.extend(sorted(ss))
        parts = []
        if swaps:
            parts.append(f"{', '.join(swaps)} {'is' if len(swaps) == 1 else 'are'} more similar to another condition than "
                         f"to {'its' if len(swaps) == 1 else 'their'} own replicates — check for a sample swap or a failed replicate.")
        if alike:
            parts.append(" ".join(f"{a} and {b} cannot be told apart by their counts: if they should differ, the "
                                  "treatment had little effect on most guides." for a, b in sorted(alike)))
        st = "pass" if not parts else "warn"
        why = "Replicates group together and conditions are distinct." if not parts else " ".join(parts)
    else:
        st, why = "na", "Needs at least three samples and a condition with replicates."
    checks.append({"id": "pca", "name": "Sample PCA", "status": st, "pca": pca, "why": why,
                   "rule": "Each sample should correlate better (log counts) with its own replicates than with any other condition."})

    # 7b non-targeting controls: the null distribution every gene is judged against
    n_ntc = int(ctrl.sum())
    if n_ntc == 0:
        checks.append({"id": "ntc", "name": "Non-targeting controls", "status": "na", "per_condition": [], "n": 0,
                       "why": "No control guides were found. Name them on the Analyse page (or in the library) so they "
                              "can set the null distribution — by default any gene matching non-targeting, NTC, "
                              "safe-harbour, EGFP, LacZ, luciferase or AAVS1.",
                       "rule": threshold_text(th)["ntc"]})
    else:
        st = _status(float(n_ntc), th["ntc_pass"], th["ntc_fail"])
        per = []
        if lfc is not None:
            for c, ss in conds.items():
                if all(is_ref[s] for s in ss):
                    continue
                cols = [samples.index(s) for s in ss]
                v = lfc[:, cols].mean(axis=1)
                cv, tv = v[ctrl], v[tgt]
                med = float(np.median(cv))
                sd = float(np.subtract(*np.percentile(cv, [75, 25])) / 1.349) if len(cv) > 3 else float("nan")
                shifted = abs(med) > th["ntc_shift"]
                per.append({"condition": c, "median": med, "sd": sd,
                            "median_targeting": float(np.median(tv)),
                            "floor": float(2 * sd) if sd == sd else float("nan"),
                            "status": "warn" if shifted else "pass"})
            if any(p["status"] != "pass" for p in per):
                st = _worst([st, "warn"])
        bad = [p for p in per if p["status"] != "pass"]
        if st == "fail":
            why = (f"Only {n_ntc:,} control guides: too few to measure the spread of a gene with no effect, so "
                   "p-values rest entirely on the targeting guides.")
        elif n_ntc < th["ntc_pass"]:
            why = (f"{n_ntc:,} control guides — usable, but {th['ntc_pass']:.0f} or more give a steadier null "
                   "distribution for the gene test.")
        elif bad:
            why = ("Control guides should sit at a log₂ fold change of zero; "
                   + "; ".join(f"in {p['condition']} their median is {p['median']:+.2f}" for p in bad[:3])
                   + ". That is a normalisation shift, not biology — the whole library moved.")
        elif per:
            why = ("{:,} control guides sit at zero as they should (median {}). Their spread sets the noise floor: "
                   "a gene effect smaller than about {} is not distinguishable from no effect.".format(
                       n_ntc, ", ".join(f"{p['median']:+.2f}" for p in per[:3]),
                       "±{:.2f}".format(np.nanmean([p["floor"] for p in per])) if per else "n/a"))
        else:
            why = f"{n_ntc:,} control guides are in the library (needs a reference sample to check how they behave)."
        checks.append({"id": "ntc", "name": "Non-targeting controls", "status": st, "n": n_ntc,
                       "per_condition": per, "why": why, "rule": threshold_text(th)["ntc"]})

    # 8 essential vs non-essential separation
    ess = None
    if species == "auto":
        from .enrichment import guess_species
        species = guess_species(list(genes))
    try:
        ceg, neg = load_reference_sets(species)
    except ValueError:
        ceg, neg = set(), set()
    gene_arr = np.asarray(genes, dtype=object)
    in_ceg = np.array([g in ceg for g in gene_arr])
    in_neg = np.array([g in neg for g in gene_arr])
    n_ceg = len(set(gene_arr[in_ceg]))
    n_neg = len(set(gene_arr[in_neg]))
    if lfc is None:
        checks.append({"id": "essentials", "name": "Essential vs non-essential genes", "status": "na", "per_condition": [],
                       "why": "Needs a reference (plasmid or T0) sample to compute fold changes.",
                       "rule": threshold_text(th)["auc"]})
    elif n_ceg < 50 or n_neg < 50:
        checks.append({"id": "essentials", "name": "Essential vs non-essential genes", "status": "na", "per_condition": [],
                       "why": f"Only {n_ceg} core-essential and {n_neg} non-essential reference genes are in this library "
                              "(needs ≥ 50 of each; the bundled lists are human and mouse gene symbols).",
                       "rule": threshold_text(th)["auc"]})
    else:
        per = []
        for c, ss in conds.items():
            if all(is_ref[s] for s in ss):
                continue
            cols = [samples.index(s) for s in ss]
            glfc = gene_means(lfc[:, cols].mean(axis=1), gene_arr)
            e = [v for g, v in glfc.items() if g in ceg]
            n = [v for g, v in glfc.items() if g in neg]
            a = auc(e, n)
            curves = roc_pr(e, n)
            per.append({"condition": c, "samples": ss, "value": a, "status": _status(a, th["auc_pass"], th["auc_fail"]),
                        "n_ess": len(e), "n_non": len(n), "median_ess": float(np.median(e)), "median_non": float(np.median(n)),
                        "recall_at_fdr5": curves["recall_at_fdr5"], "curves": curves,
                        "hist": _hist(e, n)})
        if per:
            best = max(p["value"] for p in per)
            notes = []
            for p in per:
                if p["status"] != "pass" and best >= th["auc_pass"]:
                    notes.append(f"{p['condition']} (AUC {p['value']:.3f}) separates essential genes much less than the "
                                 f"best condition (AUC {best:.3f}); if it is an earlier timepoint, it is probably too early "
                                 "for fitness effects to show.")
                elif p["status"] != "pass":
                    notes.append(f"{p['condition']}: AUC {p['value']:.3f} — essential genes are not clearly depleted, so "
                                 "the screen may lack the power to find weaker hits.")
            why = ("Core-essential genes drop out clearly in every condition (AUC {}).".format(
                ", ".join(f"{p['condition']} {p['value']:.3f}" for p in per)) if not notes else " ".join(notes))
            ess = {"id": "essentials", "name": "Essential vs non-essential genes", "status": _worst(p["status"] for p in per),
                   "per_condition": per, "why": why, "n_ceg": n_ceg, "n_neg": n_neg, "rule": threshold_text(th)["auc"]}
        else:
            ess = {"id": "essentials", "name": "Essential vs non-essential genes", "status": "na", "per_condition": [],
                   "why": "Every condition is marked as a reference.", "rule": threshold_text(th)["auc"]}
        checks.append(ess)

    overall = _worst(c["status"] for c in checks)
    return {"checks": checks, "overall": overall, "thresholds": th, "distributions": _distributions(X, samples)}


def spearman(a, b) -> float:
    """Rank correlation; robust to the few huge fold changes that dominate a Pearson r."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return float("nan")
    return pearson(_rank(a[ok]), _rank(b[ok]))


def _rank(x) -> np.ndarray:
    order = np.argsort(x, kind="mergesort")
    r = np.empty(len(x), dtype=float)
    r[order] = np.arange(len(x), dtype=float)
    # average ties, as Spearman requires
    v = x[order]
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[j + 1] == v[i]:
            j += 1
        if j > i:
            r[order[i:j + 1]] = (i + j) / 2.0
        i = j + 1
    return r


def _representation_check(screen: dict, n_guides: int, th: dict) -> dict:
    """Library representation (cells per guide) and MOI. These are experiment numbers, not read counts:
    sequencing deeply cannot recover guides that were never in enough cells."""
    def num(k):
        try:
            v = float(screen.get(k) or 0)
            return v if v > 0 else None
        except (TypeError, ValueError):
            return None
    lib = int(num("library_guides") or n_guides)
    moi = num("moi")
    steps = [("at transduction", num("cells_transduced")), ("at the tightest passage", num("cells_bottleneck")),
             ("at harvest", num("cells_harvested"))]
    given = [(name, c / lib) for name, c in steps if c]
    rows = [{"step": name, "cells": c, "per_guide": c / lib,
             "status": _status(c / lib, th["cells_pass"], th["cells_fail"])} for name, c in steps if c]
    moi_row = None
    if moi is not None:
        if moi > th["moi_fail_hi"]:
            mst = "fail"
        elif moi > th["moi_pass_hi"] or moi < th["moi_low"]:
            mst = "warn"
        else:
            mst = "pass"
        moi_row = {"moi": moi, "status": mst,
                   "one_guide": 100.0 * _single_fraction(moi)}
    if not rows and moi_row is None:
        return {"id": "representation", "name": "Library representation and MOI", "status": "na",
                "rows": [], "moi": None, "library_guides": lib,
                "why": "Not checked: fill in the screen design on the Samples page (cells transduced, cells carried "
                       "at the tightest passage, cells harvested, and the MOI). Sequencing cannot recover a guide "
                       "that was lost because too few cells carried it.",
                "rule": threshold_text(th)["cells"]}
    st = _worst([r["status"] for r in rows] + ([moi_row["status"]] if moi_row else []))
    parts = []
    worst = [r for r in rows if r["status"] != "pass"]
    if rows and not worst:
        parts.append("At least {:g} cells carry each guide throughout ({}).".format(
            th["cells_pass"], ", ".join("{:,.0f}× {}".format(r["per_guide"], r["step"]) for r in rows)))
    elif worst:
        parts.append("Only {} — below {:g} cells per guide, guides are lost by chance rather than by selection, "
                     "which shows up later as false dropouts.".format(
                         "; ".join("{:,.0f}× {}".format(r["per_guide"], r["step"]) for r in worst),
                         th["cells_pass"]))
    if moi_row:
        if moi_row["status"] == "pass":
            parts.append("At MOI {:g}, about {:.0f}% of the cells that took up virus carry a single guide.".format(
                moi, moi_row["one_guide"]))
        elif moi < th["moi_low"]:
            parts.append("MOI {:g} is very low: almost all cells are untransduced, so most of the culture is wasted "
                         "and reaching the needed representation takes far more cells.".format(moi))
        else:
            parts.append("At MOI {:g}, only about {:.0f}% of transduced cells carry one guide; the rest carry two or "
                         "more, so a knockout's effect is mixed with its neighbours' and real hits are diluted.".format(
                             moi, moi_row["one_guide"]))
    sel = str(screen.get("selection") or "").strip()
    if sel:
        parts.append(f"Selection: {sel}.")
    return {"id": "representation", "name": "Library representation and MOI", "status": st, "rows": rows,
            "moi": moi_row, "library_guides": lib, "selection": sel, "why": " ".join(parts),
            "rule": threshold_text(th)["cells"]}


def _single_fraction(moi: float) -> float:
    """Poisson: of the cells that were infected at all, the fraction that took exactly one virus."""
    p1 = moi * math.exp(-moi)
    any_ = 1.0 - math.exp(-moi)
    return p1 / any_ if any_ > 0 else 0.0


def gene_means(values, genes) -> dict:
    genes = np.asarray(genes, dtype=object).astype(str)
    order = np.argsort(genes, kind="mergesort")
    g, v = genes[order], np.asarray(values, dtype=float)[order]
    cut = np.flatnonzero(np.r_[True, g[1:] != g[:-1], True])
    return {g[s]: float(v[s:e].mean()) for s, e in zip(cut[:-1], cut[1:])}


def sample_pca(logn: np.ndarray, samples: list[str]) -> dict:
    Y = logn - logn.mean(axis=1, keepdims=True)
    var = Y.var(axis=1)
    top = np.argsort(-var)[: min(5000, len(var))]
    M = Y[top].T
    M = M - M.mean(axis=0)
    if M.shape[0] < 2:
        return {"samples": samples, "pc1": [0.0] * len(samples), "pc2": [0.0] * len(samples), "explained": [1.0, 0.0]}
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    pcs = U * S
    ev = (S ** 2) / max((S ** 2).sum(), 1e-12)
    pc2 = pcs[:, 1] if pcs.shape[1] > 1 else np.zeros(len(samples))
    return {"samples": samples, "pc1": pcs[:, 0].round(4).tolist(), "pc2": pc2.round(4).tolist(),
            "explained": ev[:2].round(4).tolist() + ([0.0] if len(ev) < 2 else [])}


def _hist(e, n, bins: int = 40) -> dict:
    allv = np.concatenate([e, n])
    lo, hi = np.percentile(allv, [0.5, 99.5])
    edges = np.linspace(lo, hi, bins + 1)
    he, _ = np.histogram(np.clip(e, lo, hi), edges, density=True)
    hn, _ = np.histogram(np.clip(n, lo, hi), edges, density=True)
    return {"edges": edges.round(4).tolist(), "ess": he.round(4).tolist(), "non": hn.round(4).tolist()}


def _distributions(X: np.ndarray, samples: list[str]) -> dict:
    """log2(count + 1) quantiles per sample, for box plots."""
    L = np.log2(X + 1.0)
    q = np.percentile(L, [1, 10, 25, 50, 75, 90, 99], axis=0)
    return {"samples": samples, "quantiles": q.T.round(3).tolist()}
