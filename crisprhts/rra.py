"""Robust rank aggregation (alpha-RRA), written to follow MAGeCK 0.5.9.5's RRA program step by step.

MAGeCK (Li et al., Genome Biology 2014; BSD licence) ranks guides, turns each guide's rank into a percentile, and
scores a gene by the smallest beta-distribution probability of its top-ranked guides that fall under the alpha
cut-off. A gene's p-value comes from scoring random percentiles (same number of guides) many times. To make the
p-values reproducible and comparable with `mageck test`, the random numbers come from the same generator MAGeCK
uses (Park & Miller's Lehmer generator, stream 0, seed 123456) and are drawn in the same order.
"""
from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy.special import betainc

# Permutation blocks run on the performance cores, leaving the machine usable. CH_RRA_WORKERS=1 forces the
# single-threaded path; the answers are the same either way.
def _workers() -> int:
    try:
        n = int(os.environ.get("CH_RRA_WORKERS", "0"))
    except ValueError:
        n = 0
    if n > 0:
        return n
    return max(1, min(6, (os.cpu_count() or 2) - 1))


_WORKERS = _workers()

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

    def jumped(self, n: int) -> "Lehmer":
        """A generator standing where this one would stand after n draws, reached without making them.
        x_n = x_0 * 48271^n mod (2^31 - 1), so a block of the stream can be produced on its own thread and
        the numbers are the very ones a single run would have used."""
        g = Lehmer.__new__(Lehmer)
        g.x = (self.x * pow(MULTIPLIER, int(n), MODULUS)) % MODULUS
        return g

    def jump(self, n: int) -> None:
        self.x = self.jumped(n).x


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
    s = _sort_rows(pct, k)
    # once a guide is out of range (after the first) RRA stops; sorted order makes this a prefix
    ok = np.logical_and.accumulate((s <= alpha) & (s >= min_pct), axis=1)
    good = ok.sum(axis=1)
    # Column by column on contiguous memory. The earlier version evaluated the beta CDF through a
    # (row, column) index pair over the whole matrix; gathering and scattering those scattered cells cost
    # several times the arithmetic itself. Results are unchanged, bit for bit.
    lo = s[:, 0].copy() if k == 1 else betainc(1.0, float(k), s[:, 0])   # betainc(1,1,x) is exactly x
    for i in range(1, k):
        m = ok[:, i]
        if not m.any():                   # prefix mask: every later column is empty too
            break
        v = betainc(i + 1.0, float(k - i), s[:, i][m])
        cur = lo[m]
        lo[m] = np.where(v < cur, v, cur)
    return np.minimum(lo, 1.0), good


def _sort_rows(a: np.ndarray, k: int) -> np.ndarray:
    """Sort each row. Nearly every gene has 1-6 guides, and for those a sorting network of min/max beats
    np.sort, which pays a per-row overhead it cannot amortise over so few elements."""
    if k == 1:
        return a
    if k == 2:
        x, y = a[:, 0], a[:, 1]
        return np.stack([np.minimum(x, y), np.maximum(x, y)], axis=1)
    if k == 3:
        x, y, z = a[:, 0], a[:, 1], a[:, 2]
        lo1, hi1 = np.minimum(x, y), np.maximum(x, y)
        return np.stack([np.minimum(lo1, z), np.minimum(hi1, np.maximum(lo1, z)), np.maximum(hi1, z)], axis=1)
    if k == 4:
        x, y, z, w = a[:, 0], a[:, 1], a[:, 2], a[:, 3]
        l1, h1 = np.minimum(x, y), np.maximum(x, y)
        l2, h2 = np.minimum(z, w), np.maximum(z, w)
        m1, m2 = np.maximum(l1, l2), np.minimum(h1, h2)
        return np.stack([np.minimum(l1, l2), np.minimum(m1, m2), np.maximum(m1, m2), np.maximum(h1, h2)], axis=1)
    return np.sort(a, axis=1)


def _null_lo(k: int, alpha: float, nrand: int, rng: "Lehmer", null) -> np.ndarray:
    """The lo-values of `nrand` random genes with k guides each — the null a gene's score is ranked against.

    This is where nearly all the time goes on a genome-wide screen, and it is dominated by the few genes with
    many guides: one gene with 94 guides needs as much work as all 17,000 genes with four. The blocks are
    independent, so they run on several threads; each starts from the generator state it would have reached by
    drawing in turn (Lehmer.jumped), so the numbers, and the p-values, are exactly those of a single-threaded
    run. NumPy releases the interpreter lock inside the sort, the comparisons and the beta CDF, which is where
    the time is spent.
    """
    chunk = max(1, 2_000_000 // k)
    spans = [(s0, min(chunk, nrand - s0)) for s0 in range(0, nrand, chunk)]
    rl = np.empty(nrand)
    base = Lehmer.__new__(Lehmer)
    base.x = rng.x

    def block(span):
        s0, m = span
        u = base.jumped(s0 * k).uniform(m * k).reshape(m, k)
        if null is not None:
            ix = np.minimum((len(null) * u).astype(np.int64), len(null) - 1)
            u = null[ix]
        rl[s0:s0 + m] = lo_values(u, alpha)[0]      # disjoint slices: no two blocks touch the same cells

    workers = min(len(spans), _WORKERS)
    if workers > 1:
        with ThreadPoolExecutor(workers) as pool:
            for _ in pool.map(block, spans):
                pass
    else:
        for span in spans:
            block(span)
    rng.jump(nrand * k)      # leave the stream where drawing one after another would have left it
    return rl


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
            rl = _null_lo(k, alpha, nrand, rng, null)
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
