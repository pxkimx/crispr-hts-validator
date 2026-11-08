"""Version 0.3: BAGEL2 Bayes factors, drugZ, copy-number correction, gene-set enrichment and the PDF report.
Unit tests use small synthetic data with known answers; validation tests compare with the official BAGEL2 and drugZ
outputs in test-data/reference/ (made by scripts/make_reference_hart.sh) and skip until those exist."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from crisprhts import bagel, cnv, drugz, enrichment, pooled_io, qc  # noqa: E402

TD = ROOT / "test-data"
TK = TD / "hart-bagel/tkov3-example/HAP1-TKOv3-EXAMPLE.txt"
REF_B = TD / "reference/bagel"
REF_D = TD / "reference/drugz"


def need(p):
    if not p.exists():
        pytest.skip(f"reference output not made yet: {p.relative_to(ROOT)} (run scripts/make_reference_hart.sh)")
    return p


# ------------------------------------------------------------------ BAGEL2
def test_bagel_fold_changes_formula():
    X = np.array([[95.0, 45.0, 5.0], [195.0, 95.0, 195.0], [0.0, 0.0, 0.0]])
    fc, keep = bagel.fold_changes(X, [0], round3=False)
    M = X + 5
    N = M / M.sum(axis=0) * 1e7
    assert np.allclose(fc, np.log2(N[:, 1:] / N[:, [0]]))
    assert keep.all()
    fc3, _ = bagel.fold_changes(X, [0])
    assert np.allclose(fc3, np.round(fc, 3), atol=5e-4)


def test_bagel_cross_validation_covers_every_gene_once():
    rng = np.random.RandomState(7)
    tr = bagel._Training(np.arange(103), rng, cvnum=10)
    seen = []
    for _ in range(10):
        train, test = tr.cross_validation()
        assert len(set(train) & set(test)) == 0 and len(train) + len(test) == 103
        seen.extend(test)
    assert sorted(seen) == list(range(103))


def _synthetic_fitness(seed=0):
    rng = np.random.default_rng(seed)
    ceg, neg = qc.load_reference_sets()
    ess, non = sorted(ceg)[:300], sorted(neg)[:300]
    other = [f"G{i}" for i in range(600)]
    genes, fc = [], []
    for g in ess + non + other:
        mu = -2.5 if g in set(ess) else 0.0
        for _ in range(4):
            genes.append(g)
            fc.append(rng.normal(mu, 0.6, 3))
    return np.array(fc), np.array(genes), ceg, neg, set(ess)


def test_bagel_separates_planted_essentials():
    fc, genes, ceg, neg, ess = _synthetic_fitness()
    r = bagel.bayes_factors(np.round(fc, 3), genes, ceg, neg, seed=1)
    bf = dict(zip(r["genes"], r["bf"]))
    assert np.median([bf[g] for g in ess]) > 20
    assert np.median([v for g, v in bf.items() if g in neg]) < -20
    pr = bagel.precision_recall(r["genes"], r["bf"], ceg, neg)
    s = bagel.summary_at(pr, 5)
    assert s["recall"] > 0.95 and s["precision"] > 0.99
    # reproducible with the same seed
    assert np.array_equal(r["bf"], bagel.bayes_factors(np.round(fc, 3), genes, ceg, neg, seed=1)["bf"])


@pytest.mark.parametrize("cond", ["T3", "T18", "T18_Starved"])
def test_bagel_matches_official(cond):
    ref = pd.read_csv(need(REF_B / f"tkov3_{cond}.bf"), sep="\t", index_col=0)
    t = pooled_io.read_count_table(TK)
    fc, keep = bagel.fold_changes(t.counts, [0])
    others = t.samples[1:]
    cols = {"T3": ["T3_A", "T3_B", "T3_C"], "T18": ["T18_A", "T18_B", "T18_C"],
            "T18_Starved": ["T18_A_Starved", "T18_B_Starved", "T18_C_Starved"]}[cond]
    ceg, neg = qc.load_reference_sets()
    r = bagel.bayes_factors(fc[:, [others.index(c) for c in cols]], np.asarray(t.genes)[keep], ceg, neg, seed=1)
    ours = pd.Series(r["bf"], index=r["genes"]).loc[ref.index]
    # same seed, same draws: equal to BAGEL's three printed decimals
    assert np.allclose(ours.round(3).to_numpy(), ref["BF"].to_numpy(), atol=5e-4)
    pr_ref = pd.read_csv(need(REF_B / f"tkov3_{cond}.pr"), sep="\t", index_col=0)
    sel = pr_ref["BF"] > 5
    pr = bagel.precision_recall(r["genes"], r["bf"], ceg, neg)
    assert bagel.summary_at(pr, 5)["calls"] == int(sel.sum())
    assert set(ours.index[ours > 5]) == set(ref.index[ref["BF"] > 5])


# ------------------------------------------------------------------ drugZ
def test_drugz_window_sd_never_decreases_and_scores():
    rng = np.random.default_rng(0)
    n = 3000
    fc = rng.normal(0, np.linspace(1.0, 0.3, n))           # noise shrinks with read count
    eb = drugz._empirical_bayes(fc, 500)
    assert np.all(np.diff(eb[500:n - 500]) >= -1e-12)
    assert eb[0] == pytest.approx(np.std(fc[:1000], ddof=1))


def test_drugz_finds_planted_interactions():
    rng = np.random.default_rng(1)
    genes = np.repeat([f"G{i}" for i in range(1000)], 4)
    base = rng.lognormal(6, 0.5, len(genes))
    ctrl = np.stack([rng.poisson(base) for _ in range(3)], axis=1)
    eff = np.ones(len(genes))
    eff[genes == "G1"] = 0.15          # synthetic lethal
    eff[genes == "G2"] = 5.0           # suppressor
    trt = np.stack([rng.poisson(base * eff) for _ in range(3)], axis=1)
    X = np.concatenate([ctrl, trt], axis=1).astype(float)
    r = drugz.drugz(X, [f"g{i}" for i in range(len(genes))], genes, [0, 1, 2], [3, 4, 5])
    g = r["genes"]
    assert g.index[0] == "G1" and g.loc["G1", "fdr_synth"] < 0.01
    assert g.index[-1] == "G2" and g.loc["G2", "fdr_supp"] < 0.01
    # normZ is a z-score across genes
    assert g["normZ"].mean() == pytest.approx(0, abs=1e-9) and g["normZ"].std(ddof=0) == pytest.approx(1, abs=1e-9)


@pytest.mark.parametrize("paired,fname", [(True, "tkov3_starved.drugz.txt"), (False, "tkov3_starved_unpaired.drugz.txt")])
def test_drugz_matches_official(paired, fname):
    ref = pd.read_csv(need(REF_D / fname), sep="\t", index_col=0)
    t = pooled_io.read_count_table(TK)
    s = t.samples
    r = drugz.drugz(t.counts, t.guides, t.genes, [s.index(x) for x in ("T18_A", "T18_B", "T18_C")],
                    [s.index(x) for x in ("T18_A_Starved", "T18_B_Starved", "T18_C_Starved")], paired=paired)
    ours = r["genes"].loc[ref.index]
    assert np.allclose(ours["normZ"].round(2), ref["normZ"], atol=5e-3)    # printed with 2 decimals
    assert np.allclose(ours["sumZ"].round(2), ref["sumZ"], atol=5e-3)
    assert (ours["rank_synth"].to_numpy() == ref["rank_synth"].to_numpy()).all()
    assert np.allclose(ours["pval_synth"], ref["pval_synth"], rtol=6e-3)     # printed with 3 significant figures
    for k in ("fdr_synth", "fdr_supp"):
        assert ((ours[k] < 0.05) == (ref[k] < 0.05)).all()


# ------------------------------------------------------------------ copy-number correction
def test_segmentation_finds_a_step():
    rng = np.random.default_rng(0)
    x = np.concatenate([rng.normal(0, 0.3, 200), rng.normal(-1.0, 0.3, 60), rng.normal(0, 0.3, 200)])
    segs = cnv.segment(x, min_size=10)
    starts = [a for a, b in segs]
    assert any(abs(a - 200) <= 5 for a in starts) and any(abs(b - 260) <= 5 for a, b in segs)


def test_cn_correction_removes_amplicon_but_keeps_single_essential_gene():
    rng = np.random.default_rng(2)
    n_genes = 400
    genes = np.repeat([f"G{i}" for i in range(n_genes)], 4)
    pos = pd.DataFrame({"guide": [f"g{i}" for i in range(len(genes))], "chrom": "1",
                        "pos": np.arange(len(genes)) * 1000.0})
    lfc = rng.normal(0, 0.3, len(genes))
    amp = (genes >= "G100") & (genes < "G112") & (np.char.str_len(genes.astype(str)) == 4)   # 12 neighbouring genes
    lfc[amp] -= 1.5
    ess = genes == "G300"
    lfc[ess] -= 3.0
    r = cnv.correct(lfc, pos["guide"], genes, pos)
    assert abs(r["lfc"][amp].mean()) < 0.25                 # amplicon effect gone
    assert r["lfc"][ess].mean() < -2.5                      # the single essential gene is kept


def test_read_positions_crisprcleanr_format():
    p = TD / "hart-bagel/tkov3-example/TKOv3-guide-positions.txt"
    if not p.exists():
        pytest.skip("positions file not present")
    pos = cnv.read_positions(p)
    assert len(pos) > 70000 and set(pos["chrom"]) >= {"1", "X"}


# ------------------------------------------------------------------ enrichment
def test_ora_matches_hypergeometric():
    universe = [f"G{i}" for i in range(200)]
    sets = {"A": set(universe[:20]), "B": set(universe[100:130])}
    hits = universe[:10] + universe[150:160]
    rows = {r["set"]: r for r in enrichment.ora(hits, universe, sets)}
    assert rows["A"]["hits"] == 10
    assert rows["A"]["p"] == pytest.approx(stats.hypergeom.sf(9, 200, 20, 20))
    assert rows["B"]["p"] == pytest.approx(1.0)


def test_ranked_matches_mann_whitney():
    rng = np.random.default_rng(0)
    scores = {f"G{i}": float(v) for i, v in enumerate(rng.normal(size=300))}
    for i in range(15):
        scores[f"G{i}"] += 2
    sets = {"top": {f"G{i}" for i in range(15)} | {f"G{i}" for i in range(200, 205)}}
    r = enrichment.ranked(scores, sets)[0]
    inside = [scores[g] for g in sets["top"]]
    outside = [v for g, v in scores.items() if g not in sets["top"]]
    p = stats.mannwhitneyu(inside, outside, alternative="two-sided", method="asymptotic", use_continuity=False).pvalue
    assert r["p"] == pytest.approx(p, rel=1e-6) and r["direction"] == "better"


def test_bundled_gene_sets():
    for coll in ("GO_BP", "GO_MF", "GO_CC", "Reactome"):
        sets = enrichment.load(coll, "human")
        assert len(sets) > 500 and all(10 <= len(v) <= 500 for v in list(sets.values())[:50])
    assert len(enrichment.load("GO_BP", "mouse")) > 500
    assert (ROOT / "crisprhts/data/genesets/SOURCES.txt").exists()
    assert enrichment.guess_species(["Trp53", "Actb", "Myc"]) == "mouse"
    assert enrichment.guess_species(["TP53", "ACTB", "MYC"]) == "human"


# ------------------------------------------------------------------ report
def test_pdf_reports():
    if not TK.exists():
        pytest.skip("test data not present")
    import time
    from crisprhts import report
    from crisprhts.session import Project
    P = Project()
    P.load_example("demo1")

    def wait(j):
        while P.jobs[j["id"]].state == "running":
            time.sleep(0.1)
        assert P.jobs[j["id"]].state == "done", P.jobs[j["id"]].error
    wait(P.run_qc())
    wait(P.run_analysis({"method": "mageck"}))
    pdf = report.pooled_report(P)
    assert pdf[:5] == b"%PDF-" and len(pdf) > 20000
