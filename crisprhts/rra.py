"""Robust rank aggregation (alpha-RRA), written to follow MAGeCK 0.5.9.5's RRA program step by step.

MAGeCK (Li et al., Genome Biology 2014; BSD licence) ranks guides, turns each guide's rank into a percentile, and
scores a gene by the smallest beta-distribution probability of its top-ranked guides that fall under the alpha
cut-off. A gene's p-value comes from scoring random percentiles (same number of guides) many times. To make the
p-values reproducible and comparable with `mageck test`, the random numbers come from the same generator MAGeCK
uses (Park & Miller's Lehmer generator, stream 0, seed 123456) and are drawn in the same order.
"""
from __future__ import annotations

import numpy as np
from scipy.special import betainc

MODULUS = 2147483647
MULTIPLIER = 48271
_BLOCK = 1 << 16
_POW: np.ndarray | None = None


def _powers() -> np.ndarray:
    global _POW
    if _POW is None:
        p = np.empty(_BLOCK, dtype=np.int64)
        x = 1
        for i in range(_BLOCK):
            x = (x * MULTIPLIER) % MODULUS
            p[i] = x
        _POW = p
    return _POW


class Lehmer:
    """x <- 48271 x mod (2^31 - 1); Random() = x / (2^31 - 1). Same sequence as rngs.c stream 0."""

    def __init__(self, seed: int = 123456):
        self.x = int(seed) % MODULUS

    def ints(self, n: int) -> np.ndarray:
        pw = _powers()
        out = np.empty(n, dtype=np.int64)
        x = self.x
        for start in range(0, n, _BLOCK):
            m = min(_BLOCK, n - start)
            out[start:start + m] = (x * pw[:m]) % MODULUS
            x = int(out[start + m - 1])
        if n:
            self.x = x
        return out

    def uniform(self, n: int) -> np.ndarray:
        return self.ints(n) / MODULUS


def btree_index(a: np.ndarray, v) -> np.ndarray:
    """Vector version of RRA's bTreeSearchingF: index of the element of sorted `a` nearest to v
    (ties in distance go to the upper neighbour; values outside the range go to the ends)."""
    v = np.asarray(v, dtype=float)
    n = len(a)
    if n == 1:
        return np.zeros(v.shape, dtype=np.int64)
    j = np.searchsorted(a, v, side="right") - 1
    j = np.clip(j, 0, n - 2)
    lo, hi = j, j + 1
    pick = np.where(np.abs(a[hi] - v) > np.abs(a[lo] - v), lo, hi)
    pick = np.where(v <= a[0], 0, pick)
    pick = np.where(v >= a[-1], n - 1, pick)
    return pick


def lo_values(pct: np.ndarray, alpha: float, min_pct: float = -1.0) -> tuple[np.ndarray, np.ndarray]:
    """RRA ComputeLoValue for each row of a (genes x guides) percentile matrix.

    Sorted percentiles x_1..x_k; rho = min over i of BetaCDF(x_i; i, k - i + 1), using x_1 always and the next ones
    while they are within the alpha cut-off. Returns (lo-value, number of 'good' guides under alpha)."""
    pct = np.asarray(pct, dtype=float)
    if pct.ndim == 1:
        pct = pct[:, None]
    r, k = pct.shape
    if k == 0:
        return np.ones(r), np.zeros(r, dtype=int)
    s = np.sort(pct, axis=1)
    out_of = (s > alpha) | (s < min_pct)
    use = ~out_of
    use[:, 0] = True
    # once a guide is out of range (after the first) RRA stops; sorted order makes this a prefix
    use = np.logical_and.accumulate(use, axis=1) | (np.arange(k) == 0)
    rows, cols = np.nonzero(use)          # only evaluate the beta CDF where RRA does
    f = np.full(s.shape, np.inf)
    f[rows, cols] = betainc((cols + 1).astype(float), (k - cols).astype(float), s[rows, cols])
    lo = np.minimum(f.min(axis=1), 1.0)
    good = np.logical_and.accumulate(~out_of, axis=1).sum(axis=1)
    return lo, good


def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg, as RRA's calculateFDRVector."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p, kind="mergesort")
    q = p[order] / (np.arange(n) + 1.0) * n
    q[-1] = min(q[-1], 1.0)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = q
    return out


def rra(scores, genes, alpha: float, control_mask=None, guide_names=None, skip_genes=("NA", "na"),
        permutations: int | None = None, max_guides_permuted: int = 100, seed: int = 123456,
        progress=None) -> dict:
    """Gene-level alpha-RRA on guide scores (smaller score = stronger selection).

    `scores`, `genes` cover only the guides that enter the test (MAGeCK's 'valid' guides).
    `control_mask` marks negative-control guides: their percentiles then form the null instead of uniform
    random numbers (MAGeCK --control-sgrna).
    Returns arrays per gene (in order of each gene's best guide): gene, n, lo, p, fdr, good.
    """
    scores = np.asarray(scores, dtype=float)
    genes = np.asarray(genes, dtype=object)
    n = len(scores)
    if n == 0:
        raise ValueError("No guides left for the gene test.")
    order = np.argsort(scores, kind="mergesort")       # the RRA input file is written in this order
    sv = np.sort(scores)
    lo_i = btree_index(sv, scores - 1e-9)
    hi_i = btree_index(sv, scores + 1e-9)
    pct = (lo_i + hi_i + 1.0) / (2.0 * n)

    # groups in order of first appearance in the (sorted) input file
    gseq = genes[order]
    names, first = np.unique(gseq.astype(str), return_index=True)
    names = names[np.argsort(first)]
    idx_of = {g: i for i, g in enumerate(names)}
    gid = np.array([idx_of[str(g)] for g in genes])
    G = len(names)
    sizes = np.bincount(gid, minlength=G)

    lo = np.ones(G)
    good = np.zeros(G, dtype=int)
    for k in np.unique(sizes):
        members = np.where(sizes == k)[0]
        mat = np.empty((len(members), k))
        pos = {g: r for r, g in enumerate(members)}
        fill = np.zeros(len(members), dtype=int)
        for gi_, p_ in zip(gid[order], pct[order]):
            r = pos.get(gi_)
            if r is not None:
                mat[r, fill[r]] = p_
                fill[r] += 1
        lo[members], good[members] = lo_values(mat, alpha)

    skip = np.array([str(g) in set(skip_genes) for g in names])
    null = None
    if control_mask is not None and np.any(control_mask):
        cm = np.asarray(control_mask, dtype=bool)
        cnames = np.asarray(guide_names if guide_names is not None else np.arange(n)).astype(str)
        # std::map order: control guides sorted by name
        ci = np.argsort(cnames[cm], kind="mergesort")
        null = pct[cm][ci]
        only_ctrl = np.bincount(gid[cm], minlength=G) == sizes
        skip |= only_ctrl
    else:
        only_ctrl = np.zeros(G, dtype=bool)

    rand_pass = 100 if permutations is None else int(permutations)
    if permutations is None and rand_pass * G < 100000:
        rand_pass = 100000 // G + 1
    scan = (rand_pass * G) // G + 1
    nrand = G * scan
    rng = Lehmer(seed)
    pval = np.ones(G)
    last = None
    ks = sorted(set(int(x) for x in sizes[~skip]))
    for step, k in enumerate(ks):
        if k <= max_guides_permuted:
            rl = np.empty(nrand)
            chunk = max(1, 2_000_000 // k)
            for s0 in range(0, nrand, chunk):
                m = min(chunk, nrand - s0)
                u = rng.uniform(m * k).reshape(m, k)
                if null is not None:
                    ix = np.minimum((len(null) * u).astype(np.int64), len(null) - 1)
                    u = null[ix]
                rl[s0:s0 + m] = lo_values(u, alpha)[0]
            rl.sort()
            last = rl
        if last is None:
            continue
        sel = np.where((sizes == k) & ~skip)[0]
        a = btree_index(last, lo[sel] - 1e-9)
        b = btree_index(last, lo[sel] + 1e-9)
        pval[sel] = (a + b + 1.0) / 2.0 / len(last)
        if progress:
            progress((step + 1) / len(ks))
    fdr = bh_fdr(pval)
    keep = ~only_ctrl
    return {"gene": names[keep], "n": sizes[keep], "lo": lo[keep], "p": pval[keep], "fdr": fdr[keep],
            "good": good[keep], "alpha": alpha, "permutation_rounds": scan, "null_size": nrand}
