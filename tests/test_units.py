"""Unit tests on small synthetic data with known answers (every QC metric, normalisation, the RRA pieces)."""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from crisprhts import count, pooled_io, pooled_stats, qc, rra  # noqa: E402


# ------------------------------------------------------------------ QC metrics
def test_gini_even_and_uneven():
    assert qc.gini([100] * 50) == pytest.approx(0.0, abs=1e-12)
    # MAGeCK formula on log(x+1): one guide with all reads, others zero
    v = [0] * 9 + [1000]
    lv = sorted(math.log(x + 1) for x in v)
    n = len(lv)
    expect = 1 - 2 * (n - sum((i + 1) * x for i, x in enumerate(lv)) / sum(lv)) / (n - 1)
    assert qc.gini(v) == pytest.approx(expect)
    assert qc.gini(v) == pytest.approx(1.0)


def test_gini_raw():
    assert qc.gini_raw([5, 5, 5, 5]) == pytest.approx(0.0)
    assert qc.gini_raw([0, 0, 0, 10]) == pytest.approx(0.75)    # (n-1)/n for one non-zero value


def test_skew_ratio():
    v = np.arange(1, 101, dtype=float)
    assert qc.skew_ratio(v) == pytest.approx(np.percentile(v, 90) / np.percentile(v, 10))
    assert qc.skew_ratio([0] * 20 + [5] * 80) == float("inf")


def test_auc_perfect_random_and_ties():
    assert qc.auc([-3, -2, -1], [1, 2, 3]) == 1.0
    assert qc.auc([1, 2, 3], [-3, -2, -1]) == 0.0
    assert qc.auc([0, 0], [0, 0]) == 0.5
    assert qc.auc([-1, 1], [0]) == 0.5


def test_roc_pr_recall_at_fdr():
    r = qc.roc_pr([-5, -4, -3, 0.5], [-0.1, 0, 1, 2, 3])
    # three essentials come first at 100% precision, so recall at 5% FDR = 3/4
    assert r["recall_at_fdr5"] == pytest.approx(0.75)
    assert r["tpr"][-1] == 1.0 and r["fpr"][-1] == 1.0


def test_pearson():
    assert qc.pearson([1, 2, 3, 4], [2, 4, 6, 8]) == pytest.approx(1.0)
    assert math.isnan(qc.pearson([1, 1, 1], [1, 2, 3]))


def _synthetic_screen(seed=0, n_genes=400, reps=3):
    """Reference + one condition; the first 100 'genes' are essential (drop 8-fold), named from CEGv2."""
    rng = np.random.default_rng(seed)
    ceg, neg = qc.load_reference_sets()
    ess = sorted(ceg)[:100]
    non = sorted(neg)[:100]
    other = [f"G{i}" for i in range(n_genes - 200)]
    genes = [g for g in ess + non + other for _ in range(4)]
    base = rng.lognormal(6, 0.4, len(genes))
    effect = np.array([0.125 if g in set(ess) else 1.0 for g in genes])
    cols = [rng.poisson(base)]
    for _ in range(reps):
        cols.append(rng.poisson(base * effect * rng.lognormal(0, 0.1, len(genes))))
    X = np.stack(cols, axis=1).astype(float)
    samples = ["T0"] + [f"T14_{c}" for c in "ABC"[:reps]]
    design = {"samples": [{"sample": "T0", "condition": "T0", "reference": True}] +
              [{"sample": s, "condition": "T14", "reference": False} for s in samples[1:]]}
    guides = [f"g{i}" for i in range(len(genes))]
    return X, guides, genes, samples, design


def test_run_qc_known_answers():
    X, guides, genes, samples, design = _synthetic_screen()
    r = qc.run_qc(X, guides, genes, samples, design)
    by = {c["id"]: c for c in r["checks"]}
    # depth: about 400 reads per guide in every sample -> pass
    assert by["reads"]["status"] == "pass"
    assert by["mapped"]["status"] == "na"              # count table: unmapped reads unknown
    assert by["zeros"]["per_sample"][0]["count"] == int((X[:, 0] == 0).sum())
    assert by["essentials"]["status"] == "pass"
    assert by["essentials"]["per_condition"][0]["value"] > 0.99
    assert by["replicates"]["status"] in ("pass", "warn", "fail")
    assert len(by["replicates"]["pairs"]) == 3


def test_run_qc_flags_problems():
    X, guides, genes, samples, design = _synthetic_screen()
    X = X.copy()
    X[:, 0] = X[:, 0] / 10          # shallow reference: ~40 reads per guide
    X[: len(X) // 10, 0] = 0        # 10% zero guides in the reference
    r = qc.run_qc(X, guides, genes, samples, design)
    by = {c["id"]: c for c in r["checks"]}
    assert by["reads"]["per_sample"][0]["status"] == "fail"
    assert by["zeros"]["per_sample"][0]["status"] == "fail"
    assert r["overall"] == "fail"


def test_replicate_disagreement_fails():
    X, guides, genes, samples, design = _synthetic_screen()
    rng = np.random.default_rng(3)
    X = X.copy()
    X[:, 3] = rng.permutation(X[:, 3])      # scramble one replicate
    r = qc.run_qc(X, guides, genes, samples, design)
    rep = next(c for c in r["checks"] if c["id"] == "replicates")
    bad = [p for p in rep["pairs"] if "T14_C" in (p["a"], p["b"])]
    assert all(p["status"] == "fail" for p in bad)


def test_essentials_not_applicable_for_non_human_names():
    X, guides, genes, samples, design = _synthetic_screen()
    genes = [f"m{g}" for g in genes]
    r = qc.run_qc(X, guides, genes, samples, design)
    assert next(c for c in r["checks"] if c["id"] == "essentials")["status"] == "na"


def test_pca_detects_swap():
    X, guides, genes, samples, design = _synthetic_screen(reps=3)
    X = np.concatenate([X, X[:, :1] * 1.01, X[:, :1] * 0.99], axis=1)     # two more reference replicates
    samples = samples + ["T0b", "T0c"]
    design["samples"][0]["condition"] = "T0"
    design["samples"] += [{"sample": "T0b", "condition": "T0", "reference": True},
                          {"sample": "T0c", "condition": "T0", "reference": True}]
    ok = qc.run_qc(X, guides, genes, samples, design)
    assert next(c for c in ok["checks"] if c["id"] == "pca")["status"] == "pass"
    # label one late sample as a reference replicate: it no longer groups with its condition
    design["samples"][1]["condition"] = "T0"
    bad = qc.run_qc(X, guides, genes, samples, design)
    assert next(c for c in bad["checks"] if c["id"] == "pca")["status"] == "warn"


# ------------------------------------------------------------------ normalisation
def test_median_ratio_size_factors():
    base = np.array([10, 20, 30, 40, 50.0])
    X = np.stack([base, base * 2, base * 4], axis=1)
    f, w = pooled_stats.size_factors(X, "median")
    N = X * f
    assert np.allclose(N[:, 0] / N[:, 1], 1.0, rtol=0.02)
    assert np.allclose(N[:, 0] / N[:, 2], 1.0, rtol=0.02)
    t, _ = pooled_stats.size_factors(X, "total")
    assert np.allclose((X * t).sum(axis=0), X.sum(axis=0).mean())


def test_control_normalisation_uses_controls_only():
    X = np.array([[100, 100], [100, 100], [100, 10], [100, 10.0]])
    ctl = np.array([True, True, False, False])
    f, _ = pooled_stats.size_factors(X, "control", ctl)
    assert f[0] == pytest.approx(f[1])


def test_too_many_zeros_falls_back_to_total():
    X = np.zeros((100, 2))
    X[:, 0] = 50
    X[:40, 1] = 50          # 60% zeros in sample 2
    f, w = pooled_stats.size_factors(X, "median")
    assert any("total" in x for x in w)


# ------------------------------------------------------------------ RRA pieces
def test_lehmer_matches_park_miller_check_value():
    g = rra.Lehmer(1)
    g.ints(10000)
    assert g.x == 399268537          # rngs.c TestRandom check value


def test_lehmer_streams_are_contiguous():
    a = rra.Lehmer(123456).uniform(200000)
    g = rra.Lehmer(123456)
    b = np.concatenate([g.uniform(70001), g.uniform(129999)])
    assert np.array_equal(a, b)


def test_lo_values_match_brute_force():
    from scipy.stats import beta
    rng = np.random.default_rng(1)
    P = rng.random((50, 5))
    lo, good = rra.lo_values(P, 0.3)
    for row, l, g in zip(P, lo, good):
        s = np.sort(row)
        best = 1.0
        n_good = 0
        for i, x in enumerate(s):
            if x > 0.3 and i > 0:
                break
            if x <= 0.3:
                n_good += 1
            best = min(best, beta.cdf(x, i + 1, len(s) - i))
        assert l == pytest.approx(best, rel=1e-9)
        assert g == n_good


def test_btree_index_semantics():
    a = np.array([0.1, 0.2, 0.2, 0.4])
    assert rra.btree_index(a, 0.2 - 1e-9) == 1        # first of the tied block
    assert rra.btree_index(a, 0.2 + 1e-9) == 2        # last of the tied block
    assert rra.btree_index(a, -1) == 0 and rra.btree_index(a, 9) == 3


def test_rra_finds_planted_gene():
    rng = np.random.default_rng(0)
    genes = np.repeat([f"G{i}" for i in range(300)], 4)
    scores = rng.normal(size=len(genes))
    scores[:4] -= 6                     # G0 strongly selected
    r = rra.rra(scores, genes, 0.1)
    i = list(r["gene"]).index("G0")
    assert r["p"][i] == min(r["p"])
    assert r["fdr"][i] < 0.05


def test_bh():
    p = np.array([0.01, 0.04, 0.03, 0.5])
    from statsmodels.stats.multitest import multipletests
    assert np.allclose(rra.bh_fdr(p), multipletests(p, method="fdr_bh")[1])
    assert np.allclose(pooled_stats.bh(p), multipletests(p, method="fdr_bh")[1])


# ------------------------------------------------------------------ inputs
def test_count_table_reader(tmp_path):
    f = tmp_path / "c.txt"
    f.write_text("sgRNA\tGene\tA\tB\ng1\tX\t1\t2\ng1\tX\t5\t5\ng2\tY\t3\nbad\tZ\tq\t1\ng3\tY\t4\t0\n")
    t = pooled_io.read_count_table(f)
    assert t.fmt == "MAGeCK" and t.guides == ["g1", "g3"] and t.counts[0].tolist() == [1, 2]
    assert len(t.notes) == 3


def test_bagel_reader_and_csv(tmp_path):
    f = tmp_path / "c.csv"
    f.write_text("SEQUENCE,GENE,T0,T18\nACGT,X,1.0,2.5\n")
    t = pooled_io.read_count_table(f)
    assert t.fmt == "BAGEL" and t.samples == ["T0", "T18"]


def test_library_and_design_guess(tmp_path):
    f = tmp_path / "lib.csv"
    f.write_text("id,seq,gene\na,ACGTACGTAC,G1\nb,ACGTACGTAC,G1\nc,TTTTACGTAC,G2\n")
    lib = pooled_io.read_library(f, "lib.csv")
    assert lib.ids == ["a", "c"]
    d = pooled_io.guess_design(["T0_-", "T18_A", "T18_B", "T18_A_Starved", "T18_B_Starved"])
    assert {c["name"] for c in d["comparisons"]} == {"T18 vs T0_-", "T18_Starved vs T0_-"}


def test_fastq_counting_offsets_mismatch_and_reverse(tmp_path):
    lib_f = tmp_path / "lib.txt"
    guides = ["ACCGTAGGTACCATGCAGTT", "TTTTGGGGCCCCAAAATTTT", "GATTACAGATTACAGATTAC"]
    lib_f.write_text("".join(f"g{i}\t{s}\tG{i}\n" for i, s in enumerate(guides)))
    lib = pooled_io.read_library(lib_f)
    reads = ["CCC" + guides[0] + "GTTTTAGAGC"] * 5 + ["CCC" + guides[1] + "GTTTTAGAGC"] * 3 + \
            ["CCC" + "A" + guides[2][1:] + "GTTTTAGAGC"] * 2       # one mismatch at the first base
    fq = tmp_path / "r.fastq"
    fq.write_text("".join(f"@r{i}\n{r}\n+\n{'I' * len(r)}\n" for i, r in enumerate(reads)))
    r = count.count_fastq(fq, lib, "s")
    assert r.trims == [3] and r.counts.tolist() == [5, 3, 0] and r.reads == 10
    r = count.count_fastq(fq, lib, "s", mismatch=True)
    assert r.counts.tolist() == [5, 3, 2]
    rc = tmp_path / "rc.fastq"
    rc.write_text("".join(f"@r{i}\n{count.revcomp(x)}\n+\n{'I' * len(x)}\n" for i, x in enumerate(reads)))
    r = count.count_fastq(rc, lib, "s")
    assert r.reverse and r.counts.tolist() == [5, 3, 0]
