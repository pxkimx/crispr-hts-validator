"""Pathway and gene-set enrichment, fully offline, with bundled GO and Reactome sets (CC BY 4.0) or a user .gmt.

Two tests:
  * over-representation (ORA): are the hit genes over-represented in a set, compared with all genes tested in the
    screen? One-sided hypergeometric test; the background is the screen, not the genome.
  * ranked: do a set's genes score better than the other genes? Mann-Whitney U test on the gene scores (two-sided,
    with the direction reported), which needs no hit cut-off.
Benjamini-Hochberg FDR over the sets tested (sets with 10-500 genes present in the screen).
"""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
from scipy import stats

from .rra import bh_fdr

DATA = Path(__file__).resolve().parent / "data" / "genesets"
COLLECTIONS = {
    "GO_BP": "GO biological process", "GO_MF": "GO molecular function", "GO_CC": "GO cellular component",
    "Reactome": "Reactome pathways",
}


def available(species: str) -> list[str]:
    return [k for k in COLLECTIONS if (DATA / f"{k}_{species}.gmt.gz").exists()]


def read_gmt(path_or_text) -> dict[str, set[str]]:
    if isinstance(path_or_text, (str, Path)) and Path(path_or_text).exists():
        p = Path(path_or_text)
        opener = gzip.open if p.suffix == ".gz" else open
        with opener(p, "rt", encoding="utf-8", errors="replace") as fh:
            text = fh.read()
    else:
        text = path_or_text.decode("utf-8", "replace") if isinstance(path_or_text, bytes) else str(path_or_text)
    sets = {}
    for ln in text.splitlines():
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 3:
            sets[f[0]] = {g.strip() for g in f[2:] if g.strip()}
    if not sets:
        raise ValueError("No gene sets found (expected .gmt: name, description, genes…, tab separated).")
    return sets


def load(collection: str, species: str = "human") -> dict[str, set[str]]:
    p = DATA / f"{collection}_{species}.gmt.gz"
    if not p.exists():
        raise ValueError(f"No {COLLECTIONS.get(collection, collection)} gene sets for {species}.")
    return read_gmt(p)


def guess_species(genes) -> str:
    g = [x for x in genes if isinstance(x, str) and x]
    if not g:
        return "human"
    upper = sum(1 for x in g if x.isupper() or not any(c.isalpha() for c in x))
    return "human" if upper / len(g) > 0.6 else "mouse"


def _filter(sets, universe, min_size, max_size):
    out = {}
    for name, genes in sets.items():
        inter = genes & universe
        if min_size <= len(inter) <= max_size:
            out[name] = inter
    return out


def ora(hits, universe, sets: dict[str, set[str]], min_size: int = 10, max_size: int = 500) -> list[dict]:
    U = set(universe)
    H = set(hits) & U
    if not H:
        raise ValueError("No hit genes to test.")
    S = _filter(sets, U, min_size, max_size)
    N, n = len(U), len(H)
    rows = []
    for name, genes in S.items():
        K = len(genes)
        k = len(genes & H)
        p = float(stats.hypergeom.sf(k - 1, N, K, n)) if k > 0 else 1.0
        rows.append({"set": name, "size": K, "hits": k, "expected": K * n / N, "ratio": (k / (K * n / N)) if K else 0.0,
                     "p": p, "genes": sorted(genes & H)})
    return _finish(rows)


def ranked(scores: dict[str, float], sets: dict[str, set[str]], higher_is_better: bool = True, min_size: int = 10,
           max_size: int = 500) -> list[dict]:
    names = np.array(list(scores))
    vals = np.array([scores[g] for g in names], dtype=float)
    ok = np.isfinite(vals)
    names, vals = names[ok], vals[ok]
    U = set(names)
    S = _filter(sets, U, min_size, max_size)
    ranks = stats.rankdata(vals if higher_is_better else -vals)    # high rank = better
    pos = {g: i for i, g in enumerate(names)}
    n = len(vals)
    rows = []
    for name, genes in S.items():
        idx = np.fromiter((pos[g] for g in genes), dtype=int)
        k = len(idx)
        r_in = ranks[idx].sum()
        u = r_in - k * (k + 1) / 2
        mu = k * (n - k) / 2
        sd = np.sqrt(k * (n - k) * (n + 1) / 12)
        z = (u - mu) / sd
        p = float(2 * stats.norm.sf(abs(z)))
        top = sorted(genes, key=lambda g: -ranks[pos[g]])[:12]
        rows.append({"set": name, "size": k, "auc": float(u / (k * (n - k))), "z": float(z),
                     "direction": "better" if z > 0 else "worse", "p": p, "genes": top})
    return _finish(rows)


def _finish(rows):
    if not rows:
        return []
    fdr = bh_fdr(np.array([r["p"] for r in rows]))
    for r, q in zip(rows, fdr):
        r["fdr"] = float(q)
    rows.sort(key=lambda r: r["p"])
    return rows
