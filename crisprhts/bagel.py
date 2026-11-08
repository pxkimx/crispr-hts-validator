"""BAGEL2-style Bayes factors for gene essentiality (Kim & Hart, Genome Medicine 2021; BAGEL2 is MIT-licensed).

Ported step by step from BAGEL.py (github.com/hart-lab/bagel, 2021) so results can be checked against it:
  fc  counts + pseudocount 5, every sample scaled to 10 million reads, control columns summed,
      log2(sample / control), written with three decimals (the bf step reads those rounded values);
  bf  10-fold cross-validation over genes (numpy's legacy random generator, seeded): in each fold, kernel
      densities of guide fold changes for training essential (CEGv2) and non-essential (NEGv1) genes give a
      log-likelihood ratio, which is fitted by a straight line over the informative range; the line scores the
      held-out genes' guides in every replicate; a gene's Bayes factor is the sum over its guides and replicates;
  pr  precision and recall against the reference sets, genes ranked by BF.
A BF above 5 is commonly read as 'essential' (BAGEL2 tutorial); the user chooses the cut-off.
"""
from __future__ import annotations

import numpy as np
from scipy import stats


def fold_changes(counts: np.ndarray, control_cols: list[int], pseudo: float = 5.0, min_reads: float = 0.0,
                 round3: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """BAGEL `fc` for every non-control column. Returns (log2 FC guides x non-control samples, kept-guide mask)."""
    X = np.asarray(counts, dtype=float)
    ctrl = X[:, control_cols].sum(axis=1)
    others = [j for j in range(X.shape[1]) if j not in control_cols]
    M = np.column_stack([X[:, others], ctrl]) + pseudo
    N = M / M.sum(axis=0) * 1e7
    keep = M[:, -1] >= min_reads
    fc = np.log2(N[:, :-1] / N[:, [-1]])
    if round3:
        fc = np.array([[float(f"{v:.3f}") for v in row] for row in fc])   # BAGEL writes '%4.3f' and reads it back
    return fc[keep], keep


class _Training:
    """BAGEL's Training class (cross-validation draws), same calls to the random generator."""

    def __init__(self, X, rng, cvnum=10):
        self._n = len(X)
        self._cvnum = cvnum
        self._bid = int(self._n / cvnum)
        self._bucket = np.arange(len(X))
        self._X = X
        self._step = 0
        self.rng = rng

    def cross_validation(self):
        if self._bid < 1:
            raise ValueError(f"Too few genes for {self._cvnum}-fold cross-validation.")
        mask = np.array([True] * self._n)
        for _ in range(self._bid):
            select = self.rng.randint(len(self._bucket))
            mask[self._bucket[select]] = False
            self._bucket = np.delete(self._bucket, select)
        if self._step < self._n % self._cvnum:
            select = self.rng.randint(len(self._bucket))
            mask[self._bucket[select]] = False
            self._bucket = np.delete(self._bucket, select)
        self._step += 1
        return self._X[mask], self._X[~mask]


def _round_hundredth(x):
    return np.around(x * 100) / 100.0


def bayes_factors(fc: np.ndarray, genes, ess: set[str], non: set[str], seed: int = 1, folds: int = 10,
                  progress=None) -> dict:
    """BAGEL `bf` (default settings: 10-fold cross-validation, no multi-targeting filter, no network boost).

    `fc`: guides x replicates log2 fold changes (as from `fold_changes`), `genes`: gene of each guide.
    Returns genes (in order of first appearance), BF per gene, and per-guide mean BF per replicate."""
    genes = np.asarray(genes, dtype=object).astype(str)
    order_genes = list(dict.fromkeys(genes))          # BAGEL keeps genes in file order
    genes_array = np.array(order_genes)
    gidx = {g: i for i, g in enumerate(order_genes)}
    guide_gene = np.array([gidx[g] for g in genes])
    guides_of = [[] for _ in order_genes]
    for i, g in enumerate(guide_gene):
        guides_of[g].append(i)
    ess_arr, non_arr = np.array(sorted(ess)), np.array(sorted(non))
    is_ess_gene = np.isin(genes_array, ess_arr)
    is_non_gene = np.isin(genes_array, non_arr)
    if is_ess_gene.sum() < 20 or is_non_gene.sum() < 20:
        raise ValueError(f"BAGEL needs reference genes in the library: found {int(is_ess_gene.sum())} essential and "
                         f"{int(is_non_gene.sum())} non-essential genes.")
    fc_thresh = 2 ** (-1.1535 * np.log(len(np.intersect1d(genes_array, non_arr)) + 13.324) + 0.7728)
    rng = np.random.RandomState(seed)
    training = _Training(np.arange(len(order_genes)), rng, cvnum=folds)
    nrep = fc.shape[1]
    guide_bf = np.full((len(genes), nrep), np.nan)
    fits = []
    for loop in range(folds):
        train_idx, test_idx = training.cross_validation()
        tr_ess = train_idx[is_ess_gene[train_idx]]
        tr_non = train_idx[is_non_gene[train_idx]]
        ess_vals = np.concatenate([fc[guides_of[g]].ravel() for g in tr_ess])
        non_vals = np.concatenate([fc[guides_of[g]].ravel() for g in tr_non])
        kess = stats.gaussian_kde(ess_vals)
        knon = stats.gaussian_kde(non_vals)
        x = np.arange(-10, 2, 0.01)
        nonfitx = knon.evaluate(x)
        f = np.where(nonfitx > fc_thresh)
        xmin = _round_hundredth(min(x[f]))
        subx = np.arange(xmin, max(x[f]), 0.01)
        logratio_sample = np.log2(kess.evaluate(subx) / knon.evaluate(subx))
        f = np.where(logratio_sample == logratio_sample.min())
        xmax = _round_hundredth(subx[f])[0]
        grid = np.arange(xmin, xmax + 0.01, 0.01)
        keys = np.around(grid * 100)
        lookup = {k: float(np.log2(kess.evaluate(i) / knon.evaluate(i))[0]) for k, i in zip(keys, grid)}
        tr_guides = np.concatenate([guides_of[g] for g in train_idx])
        vals = fc[tr_guides].ravel()
        sel = vals[(vals >= xmin) & (vals <= xmax)]
        k = np.around(sel * 100)
        testx = k / 100
        testy = np.array([lookup[kk] for kk in k])
        slope, intercept, *_ = stats.linregress(testx, testy)
        fits.append({"fold": loop, "slope": float(slope), "intercept": float(intercept), "xmin": float(xmin),
                     "xmax": float(xmax), "train_ess": int(len(tr_ess)), "train_non": int(len(tr_non))})
        te_guides = np.concatenate([guides_of[g] for g in test_idx])
        guide_bf[te_guides] = slope * fc[te_guides] + intercept
        if progress:
            progress((loop + 1) / folds)
    gene_bf = np.zeros(len(order_genes))
    np.add.at(gene_bf, guide_gene, guide_bf.sum(axis=1))
    return {"genes": genes_array, "bf": gene_bf, "guide_bf": guide_bf, "fits": fits, "fc_thresh": float(fc_thresh),
            "n_ess": int(is_ess_gene.sum()), "n_non": int(is_non_gene.sum())}


def precision_recall(genes, bf, ess: set[str], non: set[str]) -> dict:
    """BAGEL `pr`: genes sorted by BF (high first); recall over reference essentials in the screen."""
    genes = np.asarray(genes, dtype=object).astype(str)
    bf = np.asarray(bf, dtype=float)
    order = np.argsort(-bf, kind="mergesort")
    g, b = genes[order], bf[order]
    is_e = np.array([x in ess for x in g])
    is_n = np.array([(x in non) and (x not in ess) for x in g])
    tp, fp = np.cumsum(is_e), np.cumsum(is_n)
    tot = max(1, int(is_e.sum()))
    recall = tp / tot
    precision = np.where(tp + fp > 0, tp / np.maximum(tp + fp, 1), 1.0)
    return {"gene": g, "bf": b, "recall": recall, "precision": precision, "fdr": 1 - precision}


def summary_at(pr: dict, cutoff: float = 5.0) -> dict:
    """Essential calls at BF > cutoff and how they score against the reference sets."""
    sel = pr["bf"] > cutoff
    n = int(sel.sum())
    if n == 0:
        return {"cutoff": cutoff, "calls": 0, "recall": 0.0, "precision": None}
    i = np.where(sel)[0][-1]
    return {"cutoff": cutoff, "calls": n, "recall": float(pr["recall"][i]), "precision": float(pr["precision"][i])}
