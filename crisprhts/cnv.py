"""Copy-number correction from guide positions, after CRISPRcleanR (Iorio et al., BMC Genomics 2018).

In cancer cell lines, guides in amplified regions deplete because Cas9 cuts many copies, whatever the gene does.
CRISPRcleanR finds such regions without any copy-number data: guide fold changes are ordered along each
chromosome and segmented; in a segment whose guides target at least `min_genes` different genes, the segment's mean
fold change is subtracted from every guide. A single essential gene (a few adjacent guides) never forms such a
segment, so true hits are kept.

Segmentation here is recursive binary segmentation: the split with the largest two-sample t statistic is kept when
its permutation p-value (shuffled fold changes, fixed seed) is below `alpha`, and both halves are searched again.
CRISPRcleanR uses circular binary segmentation (DNAcopy); the two find the same large segments, binary segmentation
can split a short internal amplicon in two steps instead of one.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

from .pooled_io import _text


def read_positions(src) -> pd.DataFrame:
    """Guide positions: a guide-ID column plus chromosome and start columns (any order; names like chr / chrm /
    chromosome and start / startpos / pos / position). CRISPRcleanR library files work as they are."""
    lines = [ln for ln in _text(src).split("\n") if ln.strip()]
    sep = "\t" if lines[0].count("\t") >= lines[0].count(",") else ","
    head = [h.strip().strip('"') for h in lines[0].split(sep)]
    hl = [h.lower() for h in head]
    ci = next((i for i, h in enumerate(hl) if h in ("chr", "chrm", "chrom", "chromosome")), None)
    pi = next((i for i, h in enumerate(hl) if h in ("start", "startpos", "pos", "position", "cut", "cut_position")), None)
    gi = next((i for i, h in enumerate(hl) if h in ("", "id", "sgrna", "guide", "guide_id", "sgrna_id", "code") and i not in (ci, pi)), 0)
    if ci is None or pi is None:
        raise ValueError("The position file needs a chromosome column (chr / chromosome) and a start column (start / pos).")
    rows = []
    bad = 0
    for ln in lines[1:]:
        f = [x.strip().strip('"') for x in ln.split(sep)]
        if len(f) <= max(ci, pi, gi):
            bad += 1
            continue
        try:
            rows.append((f[gi], re.sub(r"^chr", "", f[ci], flags=re.I), float(f[pi])))
        except ValueError:
            bad += 1
    df = pd.DataFrame(rows, columns=["guide", "chrom", "pos"]).drop_duplicates("guide")
    df.attrs["skipped"] = bad
    return df


def _max_t(x: np.ndarray, min_size: int) -> tuple[float, int]:
    """Largest |two-sample t| over split points (left = x[:k]); vectorised over rows of a 2-D array as well."""
    X = np.atleast_2d(x)
    n = X.shape[1]
    k = np.arange(min_size, n - min_size + 1)
    if len(k) == 0:
        return (np.zeros(X.shape[0]), np.zeros(X.shape[0], dtype=int))
    cs = np.cumsum(X, axis=1)
    tot = cs[:, -1:]
    left = cs[:, k - 1] / k
    right = (tot - cs[:, k - 1]) / (n - k)
    # robust noise estimate from successive differences (insensitive to the step itself)
    sd = 1.4826 * np.median(np.abs(np.diff(X, axis=1) - np.median(np.diff(X, axis=1), axis=1, keepdims=True)), axis=1) / np.sqrt(2)
    sd = np.where(sd > 0, sd, X.std(axis=1) + 1e-12)
    t = np.abs(left - right) / (sd[:, None] * np.sqrt(1.0 / k + 1.0 / (n - k)))
    j = np.argmax(t, axis=1)
    return t[np.arange(X.shape[0]), j], k[j]


def segment(x: np.ndarray, alpha: float = 0.01, nperm: int = 200, min_size: int = 10, seed: int = 0) -> list[tuple[int, int]]:
    """Recursive binary segmentation of a 1-D signal. Returns [(start, end)] half-open index ranges."""
    rng = np.random.default_rng(seed)
    out = []
    stack = [(0, len(x))]
    while stack:
        a, b = stack.pop()
        seg = x[a:b]
        if b - a < 2 * min_size:
            out.append((a, b))
            continue
        tmax, k = _max_t(seg, min_size)
        perms = np.stack([rng.permutation(seg) for _ in range(nperm)])
        tp, _ = _max_t(perms, min_size)
        p = (1 + np.sum(tp >= tmax[0])) / (1 + nperm)
        if p < alpha:
            stack.append((a + int(k[0]), b))
            stack.append((a, a + int(k[0])))
        else:
            out.append((a, b))
    return sorted(out)


def correct(lfc: np.ndarray, guides, genes, positions: pd.DataFrame, min_genes: int = 3, alpha: float = 0.01,
            min_size: int = 10, progress=None) -> dict:
    """CRISPRcleanR-style correction of guide fold changes. Guides without a position are left as they are.
    Returns corrected LFC, the per-guide correction (segment mean subtracted) and the segments."""
    lfc = np.asarray(lfc, dtype=float)
    pos = positions.set_index("guide")
    idx = pd.Index(list(guides))
    have = idx.isin(pos.index)
    if have.sum() < 100:
        raise ValueError(f"Only {int(have.sum())} guides have a genomic position; copy-number correction needs the "
                         "guide IDs in the position file to match the count table.")
    chrom = pd.Series(pos.reindex(idx)["chrom"].to_numpy(), dtype=object)
    p = pos.reindex(idx)["pos"].to_numpy()
    genes = np.asarray(genes, dtype=object)
    shift = np.zeros(len(lfc))
    segments = []
    chroms = [c for c in pd.unique(chrom[have]) if isinstance(c, str)]
    chroms.sort(key=lambda c: (0, int(c)) if c.isdigit() else (1, c))
    for n, c in enumerate(chroms):
        ii = np.where((chrom == c).to_numpy() & have & np.isfinite(lfc))[0]
        ii = ii[np.argsort(p[ii], kind="mergesort")]
        for a, b in segment(lfc[ii], alpha=alpha, min_size=min_size, seed=n):
            g = ii[a:b]
            ng = len(set(genes[g]))
            m = float(np.mean(lfc[g]))
            corrected = ng >= min_genes
            if corrected:
                shift[g] = m
            segments.append({"chrom": c, "start": float(p[g[0]]), "end": float(p[g[-1]]), "guides": int(b - a),
                             "genes": ng, "mean_lfc": m, "corrected": corrected})
        if progress:
            progress((n + 1) / len(chroms))
    return {"lfc": lfc - shift, "shift": shift, "segments": segments, "n_positioned": int(have.sum())}
