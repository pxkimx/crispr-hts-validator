"""drugZ: gene-level chemogenetic interaction scores, treated vs untreated (Colic et al., Genome Medicine 2019;
drugZ is MIT-licensed). Ported step by step from drugz.py (github.com/hart-lab/drugz, 2021):

  * counts of the control and treated samples are scaled to 10 million reads each;
  * paired (default): for each replicate pair, guide log2 fold change = log2((treated + 5) / (control + 5));
    unpaired: the same with the mean of the treated and of the control samples;
  * guides are sorted by control read count (high to low) and the fold-change SD is estimated in a sliding window
    of 1,000 guides moved in steps of 25, never allowed to decrease ('empirical Bayes' smoothing); z = FC / SD;
  * per gene: sum of z over guides and replicates, divided by √(number of observations), then z-scored across
    genes = normZ. Negative normZ = synthetic sick/lethal (the drug hurts more without the gene), positive =
    suppressor. p-values from the normal distribution, FDR by Benjamini-Hochberg in each direction.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

NORM_VALUE = 1e7


def _empirical_bayes(fc_sorted: np.ndarray, half_window: int = 500) -> np.ndarray:
    n = len(fc_sorted)
    s = pd.Series(fc_sorted)
    eb = np.zeros(n)
    eb[0:half_window] = s.iloc[0: half_window * 2].std()
    for i in range(half_window, n - half_window + 25, 25):
        sd = s.iloc[i - half_window: i + half_window].std()
        if sd >= eb[i - 1]:
            eb[i: i + 25] = sd
        else:
            eb[i: i + 25] = eb[i - 1]
    eb[n - half_window:] = eb[n - (half_window + 1)]
    return eb


def _gene_scores(z: pd.DataFrame, genes: pd.Series, min_obs: int = 1) -> pd.DataFrame:
    vals = z.to_numpy()
    df = pd.DataFrame({"gene": genes.to_numpy()})
    rows = {}
    for g, idx in df.groupby("gene", sort=True).indices.items():
        v = vals[idx]
        rows[g] = (np.nansum(v), np.count_nonzero(v))
    res = pd.DataFrame.from_dict(rows, orient="index", columns=["sumZ", "numObs"])
    res = res.loc[res["numObs"] >= min_obs].copy()
    n = len(res)
    res["normZ"] = stats.zscore(res["sumZ"] / np.sqrt(res["numObs"]))
    res = res.sort_values("normZ", ascending=True)
    res["pval_synth"] = stats.norm.sf(res["normZ"] * -1)
    res["rank_synth"] = np.arange(1, n + 1)
    fdr = res["pval_synth"] / (res["rank_synth"] / float(n))
    res["fdr_synth"] = np.minimum.accumulate(fdr[::-1])[::-1]
    res = res.sort_values("normZ", ascending=False)
    res["pval_supp"] = stats.norm.sf(res["normZ"])
    res["rank_supp"] = np.arange(1, n + 1)
    fdr = res["pval_supp"] / (res["rank_supp"] / float(n))
    res["fdr_supp"] = np.minimum.accumulate(fdr[::-1])[::-1]
    res = res.sort_values("normZ", ascending=True)
    res.index.name = "GENE"
    return res


def drugz(counts: np.ndarray, guides, genes, control_cols: list[int], treated_cols: list[int], paired: bool = True,
          pseudo: float = 5.0, half_window: int = 500, remove_genes=(), min_obs: int = 1) -> dict:
    """Run drugZ. `counts` guides x samples; columns given by index. Returns gene table and guide z-scores."""
    if paired and len(control_cols) != len(treated_cols):
        raise ValueError("Paired drugZ needs the same number of control and treated samples (pair them by replicate), "
                         "or choose the unpaired form.")
    reads = pd.DataFrame(np.asarray(counts, dtype=float), index=pd.Index(list(guides), name="guide"))
    reads.insert(0, "GENE", list(genes))
    if remove_genes:
        reads = reads.loc[~reads["GENE"].isin(set(remove_genes))]
    cols = list(control_cols) + list(treated_cols)
    sub = reads[cols]
    norm = NORM_VALUE * sub / sub.sum().values
    n = len(reads)
    if n < 2 * half_window:
        raise ValueError(f"drugZ's variance window needs at least {2 * half_window} guides (this screen has {n}).")
    zcols = []
    if paired:
        frames = []
        for i, (c, t) in enumerate(zip(control_cols, treated_cols)):
            fc = pd.DataFrame(index=reads.index)
            fc["GENE"] = reads["GENE"]
            fc["ctrl_reads"] = reads[c]
            fc = fc.sort_values("ctrl_reads", ascending=False)            # same (unstable) sort as drugz.py
            f = np.log2((norm[t] + pseudo) / (norm[c] + pseudo))
            fcv = f.loc[fc.index].to_numpy()
            eb = _empirical_bayes(fcv, half_window)
            fc[f"fc_{i}"] = fcv
            fc[f"eb_std_{i}"] = eb
            fc[f"zscore_fc_{i}"] = fcv / eb
            frames.append(fc.drop(columns=["ctrl_reads"]))
            zcols.append(f"zscore_fc_{i}")
        allfc = pd.concat(frames, axis=1, sort=False)
        allfc = allfc.loc[:, ~allfc.columns.duplicated()]
    else:
        fc = pd.DataFrame(index=reads.index)
        fc["GENE"] = reads["GENE"]
        fc["ctrl_mean_reads"] = reads[list(control_cols)].mean(axis=1)
        fc = fc.sort_values("ctrl_mean_reads", ascending=False)
        f = np.log2((norm[list(treated_cols)].mean(axis=1) + pseudo) / (norm[list(control_cols)].mean(axis=1) + pseudo))
        fcv = f.loc[fc.index].to_numpy()
        eb = _empirical_bayes(fcv, half_window)
        fc["mean_fc"] = fcv
        fc["eb_std"] = eb
        fc["zscore"] = fcv / eb
        allfc = fc
        zcols = ["zscore"]
    genes_res = _gene_scores(allfc[zcols], allfc["GENE"], min_obs)
    return {"genes": genes_res, "guides": allfc, "zcols": zcols, "paired": paired}
