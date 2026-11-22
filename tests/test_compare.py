"""Version 1.5: comparing two pooled screens — common scaling, the differential test, and the calls."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from crisprhts import compare, drugz, pooled_io, qc  # noqa: E402

TK = ROOT / "test-data/hart-bagel/tkov3-example/HAP1-TKOv3-EXAMPLE.txt"
A = {"name": "T18", "treatment": ["T18_A", "T18_B", "T18_C"], "control": ["T0_-"]}
B = {"name": "T18 starved", "treatment": ["T18_A_Starved", "T18_B_Starved", "T18_C_Starved"], "control": ["T0_-"]}


# ------------------------------------------------------------------ synthetic pair with planted differences
def _pair(seed=0, n_genes=900, guides=4, a_specific=(), b_specific=(), essential=(), depth=1.0, noise=0.25):
    """Two screens from one library: T0, then arm A and arm B, three replicates each."""
    rng = np.random.default_rng(seed)
    ceg, neg = qc.load_reference_sets()
    ess_names, non_names = sorted(ceg)[:200], sorted(neg)[:200]
    others = [f"G{i}" for i in range(n_genes)]
    names = ess_names + non_names + others
    genes = np.repeat(names, guides)
    gl = [f"{g}_{k}" for g in names for k in range(guides)]
    base = rng.lognormal(7.0, 0.35, len(genes))
    t0 = np.stack([rng.poisson(base) for _ in range(1)], axis=1).astype(float)

    def arm(extra, scale):
        eff = np.ones(len(genes))
        eff[np.isin(genes, ess_names)] = 0.12                      # core-essential genes drop out
        eff[np.isin(genes, list(essential))] = 0.12
        eff[np.isin(genes, list(extra))] = 0.10                    # this arm's own dependencies
        eff = eff ** scale                                         # a screen with more or less dynamic range
        return np.stack([rng.poisson(base * eff * np.exp(rng.normal(0, noise, len(genes))))
                         for _ in range(3)], axis=1).astype(float)

    X = np.concatenate([t0, arm(a_specific, 1.0), arm(b_specific, depth)], axis=1)
    samples = ["T0"] + [f"A{i}" for i in range(3)] + [f"B{i}" for i in range(3)]
    ca = {"name": "arm A", "treatment": ["A0", "A1", "A2"], "control": ["T0"]}
    cb = {"name": "arm B", "treatment": ["B0", "B1", "B2"], "control": ["T0"]}
    return X.astype(float), gl, list(genes), samples, ca, cb


def test_planted_specific_genes_are_recovered():
    aonly, bonly = ["G1", "G2", "G3"], ["G10", "G11", "G12", "G13"]
    X, gl, genes, s, ca, cb = _pair(a_specific=aonly, b_specific=bonly)
    r = compare.compare(X, gl, genes, s, ca, cb)
    g = r["genes"].set_index("gene")
    for name in aonly:
        assert g.loc[name, "call"] == "specific to A", (name, g.loc[name].to_dict())
        assert g.loc[name, "different"] and g.loc[name, "delta"] > 0
    for name in bonly:
        assert g.loc[name, "call"] == "specific to B", (name, g.loc[name].to_dict())
        assert g.loc[name, "different"] and g.loc[name, "delta"] < 0
    # and the genes that were not planted are overwhelmingly left alone
    planted = set(aonly) | set(bonly)
    other = g[~g.index.isin(planted)]
    assert other["different"].mean() < 0.02


def test_two_identical_screens_give_no_differences():
    """The same biology twice, differing only by sampling noise: nothing should be called."""
    X, gl, genes, s, ca, cb = _pair(seed=4)
    r = compare.compare(X, gl, genes, s, ca, cb)
    g = r["genes"]
    assert g["different"].sum() == 0, g[g["different"]][["gene", "delta", "fdr"]].head().to_dict()
    assert r["concordance"]["n_specific_a"] == 0 and r["concordance"]["n_specific_b"] == 0
    assert r["concordance"]["pearson"] > 0.9


def test_scaling_puts_essential_at_minus_one_and_non_essential_at_zero():
    X, gl, genes, s, ca, cb = _pair(seed=1, depth=1.6)        # arm B deliberately has more dynamic range
    r = compare.compare(X, gl, genes, s, ca, cb)
    ceg, neg = qc.load_reference_sets()
    g = r["genes"].set_index("gene")
    for col in ("scaled_a", "scaled_b"):
        e = g[col][[x in ceg for x in g.index]]
        n = g[col][[x in neg for x in g.index]]
        assert e.median() == pytest.approx(-1.0, abs=0.02), (col, e.median())
        assert n.median() == pytest.approx(0.0, abs=0.02), (col, n.median())
    # the raw effects differ in range between the two arms; the scaled ones agree
    assert abs(r["a"]["scale"]["span"] - r["b"]["scale"]["span"]) > 0.3
    assert abs(g["scaled_a"].std() - g["scaled_b"].std()) < 0.25


def test_a_screen_without_the_reference_gene_sets_is_not_scaled_silently():
    X, gl, genes, s, ca, cb = _pair(seed=2)
    genes = [f"X{i}" for i in range(len(genes))]              # no CEG or NEG genes at all
    r = compare.compare(X, gl, genes, s, ca, cb)
    assert r["a"]["scale"]["ok"] is False and r["b"]["scale"]["ok"] is False


# ------------------------------------------------------------------ against what the app already does
@pytest.mark.skipif(not TK.exists(), reason="test data not present")
def test_shared_reference_comparison_agrees_with_drugz():
    """When both screens are measured against the same T0, comparing them asks drugZ's question by another
    route — guide differences instead of treated-vs-untreated counts. The two should rank genes alike."""
    t = pooled_io.read_count_table(TK)
    r = compare.compare(t.counts, t.guides, t.genes, t.samples, A, B)
    assert r["shared_reference"] is True
    s = t.samples
    dz = drugz.drugz(t.counts, t.guides, t.genes,
                     [s.index(x) for x in A["treatment"]], [s.index(x) for x in B["treatment"]])["genes"]
    m = r["genes"].set_index("gene").join(dz[["normZ"]].rename(columns={"normZ": "dz"}), how="inner")
    assert len(m) > 17000
    assert pd.Series(m["normZ"]).corr(pd.Series(m["dz"]), method="spearman") > 0.9
    top = set(r["genes"].head(6)["gene"])
    assert {"TSC1", "TSC2", "KEAP1"} <= top, top        # starvation makes mTORC1 restraint essential


@pytest.mark.skipif(not TK.exists(), reason="test data not present")
def test_tkov3_starvation_comparison_reads_sensibly():
    t = pooled_io.read_count_table(TK)
    r = compare.compare(t.counts, t.guides, t.genes, t.samples, A, B)
    c = r["concordance"]
    assert c["pearson"] > 0.8 and 0.5 < c["jaccard"] < 0.95
    assert c["dep_both"] > 1000 and c["n_specific_b"] > c["n_specific_a"]
    g = r["genes"]
    # agreement is reported, not used to filter: on TKOv3 four of the 88 called genes rest on a 2-of-4 guide
    # split, which the table has to show rather than hide
    agree = g["agreement"][g["different"]]
    assert (agree >= 0.5).all() and agree.median() >= 0.75
    assert g["guides"].min() >= compare.DEFAULTS["cmp_min_guides"]


def test_minimum_guides_is_respected():
    X, gl, genes, s, ca, cb = _pair(seed=3, guides=2)
    r = compare.compare(X, gl, genes, s, ca, cb, thresholds={"cmp_min_guides": 3})
    assert len(r["genes"]) == 0                              # every gene has two guides, so none qualifies
    r2 = compare.compare(X, gl, genes, s, ca, cb, thresholds={"cmp_min_guides": 2})
    assert len(r2["genes"]) > 100


# ------------------------------------------------------------------ through the app
def _project_with_comparison():
    import time
    from crisprhts import session
    P = session.Project()
    P.load_example("tkov3")              # 174 non-targeting controls: the case that broke the QC export
    def run(start):                      # one at a time: the app refuses to start a second job
        j = start()
        while P.jobs[j["id"]].state == "running":
            time.sleep(0.05)
        assert P.jobs[j["id"]].state == "done", P.jobs[j["id"]].error
    run(P.run_qc)
    run(lambda: P.run_compare({"a": "BATCH1_HAP1_T18", "b": "BATCH1_HAP1_T18_STARVED"}))
    return P


def test_comparison_exports_on_a_screen_with_non_targeting_controls():
    P = _project_with_comparison()
    body, fname, _ = P.export("comparison", fmt="csv")
    assert b"specific to B" in body and b"TSC2" in body and fname.endswith(".csv")
    book, _, _ = P.export("all", fmt="xlsx")          # the workbook gains a comparison sheet
    assert len(book) > 10000
    import io

    import openpyxl
    assert "comparison" in openpyxl.load_workbook(io.BytesIO(book)).sheetnames


def test_comparison_view_and_gene_panel():
    P = _project_with_comparison()
    v = P.comparison_view()
    assert v["total"] > 17000 and len(v["points"]) == v["total"]
    assert v["concordance"]["n_specific_b"] > 0 and len(v["quality"]) == 2
    assert all(isinstance(p[4], float) for p in v["points"][:50])     # FDR keeps full precision
    assert min(p[4] for p in v["points"]) > 0                         # ... so nothing rounds to zero
    hit = P.comparison_view(search="TSC2")
    assert any(r[0] == "TSC2" for r in hit["rows"])
    g = P.comparison_gene("TSC2")
    assert len(g["guides"]) == 4 and all(x["lfc_b"] < x["lfc_a"] for x in g["guides"])


def test_the_report_carries_the_comparison():
    from crisprhts import report
    P = _project_with_comparison()
    pdf = report.pooled_report(P)
    assert pdf[:5] == b"%PDF-" and len(pdf) > 50000


def test_comparing_a_screen_with_itself_is_refused():
    import pytest as _pt
    P = _project_with_comparison()
    with _pt.raises(ValueError):
        P.run_compare({"a": "BATCH1_HAP1_T18", "b": "BATCH1_HAP1_T18"})
    with _pt.raises(ValueError):
        P.run_compare({"a": "BATCH1_HAP1_T18", "b": "not a comparison"})
