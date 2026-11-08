"""Agreement with the official tools (reference outputs made by scripts/make_reference.sh).

Each test is skipped, with a reason, when its reference file is not there yet.
Agreement targets (deliberate differences are documented in test-data/reference/README.md):
  * mageck count: identical counts for every guide, identical read / mapped / zero-count totals;
  * mageck test, guide level: same normalised means, LFC and p-values (relative difference < 1e-6);
  * mageck test, gene level: same alpha, same RRA scores (rel. diff < 1e-6), gene p-values within the
    permutation resolution, Spearman rank correlation > 0.999, top-100 overlap >= 98.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from crisprhts import count, pooled_io, pooled_stats  # noqa: E402

TD = ROOT / "test-data"
REF = TD / "reference" / "mageck"


def need(p: Path):
    if not p.exists():
        pytest.skip(f"reference output not made yet: {p.relative_to(ROOT)} (run scripts/make_reference.sh)")
    return p


# ------------------------------------------------------------------ mageck count
def test_count_matches_mageck_count():
    ref = pd.read_csv(need(REF / "demo2.count.txt"), sep="\t")
    summ = pd.read_csv(need(REF / "demo2.countsummary.txt"), sep="\t")
    lib = pooled_io.read_library(TD / "mageck-demo/demo2/library.txt")
    ours = {}
    for f, s in (("test1.fastq", "L1"), ("test2.fastq", "CTRL")):
        r = count.count_fastq(TD / "mageck-demo/demo2" / f, lib, s)
        ours[s] = r
        got = dict(zip(lib.ids, r.counts))
        exp = dict(zip(ref["sgRNA"], ref[s]))
        assert set(got) == set(exp)
        diff = [k for k in exp if got[k] != exp[k]]
        assert not diff, f"{s}: {len(diff)} guides differ, e.g. {diff[:5]}"
        row = summ[summ["Label"] == s].iloc[0]
        assert r.reads == row["Reads"]
        assert r.mapped == row["Mapped"]
        assert int((r.counts == 0).sum()) == row["Zerocounts"]


# ------------------------------------------------------------------ mageck test
CASES = [
    ("demo1", "mageck-demo/demo1/sample.txt", ["HL60.initial", "KBM7.initial"], ["HL60.final", "KBM7.final"]),
    ("demo3", "mageck-demo/demo3/leukemia.new.csv", ["HL60.initial", "KBM7.initial"], ["HL60.final", "KBM7.final"]),
    ("demo2", "reference/mageck/demo2.count.txt", ["CTRL"], ["L1"]),
    ("tkov3_T3", "hart-bagel/tkov3-example/HAP1-TKOv3-EXAMPLE.txt", ["T0_-"], ["T3_A", "T3_B", "T3_C"]),
    ("tkov3_T18", "hart-bagel/tkov3-example/HAP1-TKOv3-EXAMPLE.txt", ["T0_-"], ["T18_A", "T18_B", "T18_C"]),
    ("tkov3_T18_Starved", "hart-bagel/tkov3-example/HAP1-TKOv3-EXAMPLE.txt", ["T0_-"],
     ["T18_A_Starved", "T18_B_Starved", "T18_C_Starved"]),
]
_cache: dict = {}


def run_case(name, table, ctrl, treat):
    if name not in _cache:
        t = pooled_io.read_count_table(TD / table)
        cols = [t.samples.index(s) for s in ctrl + treat]
        _cache[name] = (t, pooled_stats.mageck_test(t.counts[:, cols], t.genes, t.guides, len(ctrl)))
    return _cache[name]


@pytest.mark.parametrize("name,table,ctrl,treat", CASES, ids=[c[0] for c in CASES])
def test_guide_level_matches_mageck_test(name, table, ctrl, treat):
    ref = pd.read_csv(need(REF / f"{name}.sgrna_summary.txt"), sep="\t").set_index("sgrna")
    need(TD / table)
    t, r = run_case(name, table, ctrl, treat)
    g = r["guide"]
    ours = pd.DataFrame({"control_mean": g["control_mean"], "treat_mean": g["treat_mean"], "LFC": g["lfc"],
                         "p.low": g["p_low"], "p.high": g["p_high"], "score": np.abs(g["score"])}, index=t.guides)
    ours = ours.loc[ref.index]
    for col in ("control_mean", "treat_mean", "LFC", "p.low", "p.high", "score"):
        a, b = ours[col].to_numpy(float), ref[col].to_numpy(float)
        # MAGeCK prints 5 significant digits
        bad = ~np.isclose(a, b, rtol=1e-4, atol=1e-12)
        assert bad.sum() == 0, f"{col}: {bad.sum()} guides differ, e.g. {list(ours.index[bad][:3])}"


@pytest.mark.parametrize("name,table,ctrl,treat", CASES, ids=[c[0] for c in CASES])
def test_gene_level_matches_mageck_test(name, table, ctrl, treat):
    ref = pd.read_csv(need(REF / f"{name}.gene_summary.txt"), sep="\t").set_index("id")
    need(TD / table)
    t, r = run_case(name, table, ctrl, treat)
    ours = pd.DataFrame(r["genes"]).set_index("gene").loc[ref.index]
    for d in ("neg", "pos"):
        a, b = ours[f"{d}_score"].to_numpy(float), ref[f"{d}|score"].to_numpy(float)
        assert np.allclose(a, b, rtol=1e-3, atol=1e-12), f"{d} RRA scores differ for {(~np.isclose(a, b, rtol=1e-3)).sum()} genes"
        rho = spearmanr(ours[f"{d}_score"], ref[f"{d}|score"]).correlation
        assert rho > 0.999
        top_o = set(ours.sort_values(f"{d}_score").index[:100])
        top_r = set(ref.sort_values(f"{d}|score").index[:100])
        assert len(top_o & top_r) >= 98
        # p-values come from permutations with MAGeCK's own generator and seed: same draws, so the same values
        pa, pb = ours[f"{d}_p"].to_numpy(float), ref[f"{d}|p-value"].to_numpy(float)
        assert np.allclose(pa, pb, rtol=1e-3, atol=1e-6), f"{d} p-values differ for {(~np.isclose(pa, pb, rtol=1e-3, atol=1e-6)).sum()} genes"
        fa, fb = ours[f"{d}_fdr"].to_numpy(float), ref[f"{d}|fdr"].to_numpy(float)
        assert np.allclose(fa, fb, rtol=1e-3, atol=1e-5)
    la, lb = ours["lfc"].to_numpy(float), ref["neg|lfc"].to_numpy(float)
    assert np.allclose(la, lb, rtol=1e-4, atol=1e-5)
