"""Version 1.1: the screen-design checks — library representation (cells per guide), MOI, guides per gene,
non-targeting controls, and Spearman replicate agreement. Synthetic data with known answers."""
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from crisprhts import qc, session  # noqa: E402

TH = qc.THRESHOLDS


def by_id(res, key):
    return next(c for c in res["checks"] if c["id"] == key)


# ------------------------------------------------------------------ representation and MOI
def test_representation_not_checked_until_the_numbers_are_given():
    c = qc._representation_check({}, 70000, TH)
    assert c["status"] == "na" and "Samples page" in c["why"]


def test_cells_per_guide_pass_warn_fail():
    lib = 70000
    for cells, want in [(lib * 1000, "pass"), (lib * 500, "pass"), (lib * 300, "warn"), (lib * 100, "fail")]:
        c = qc._representation_check({"cells_bottleneck": cells, "library_guides": lib}, lib, TH)
        assert c["status"] == want, (cells / lib, c["status"])
    # the tightest step decides, not the largest
    c = qc._representation_check({"cells_transduced": lib * 1000, "cells_bottleneck": lib * 50,
                                  "library_guides": lib}, lib, TH)
    assert c["status"] == "fail" and "tightest passage" in c["why"]


def test_moi_single_guide_fraction_is_poisson():
    # of the cells infected at all, the share that took exactly one virus
    for moi in (0.1, 0.3, 1.0, 2.0):
        p1 = stats.poisson.pmf(1, moi)
        any_ = 1 - stats.poisson.pmf(0, moi)
        assert qc._single_fraction(moi) == pytest.approx(p1 / any_, rel=1e-12)
    assert qc._single_fraction(0.3) > 0.85           # the usual target
    assert qc._single_fraction(1.0) < 0.6


def test_moi_status_and_wording():
    lib = 70000
    assert qc._representation_check({"moi": 0.3}, lib, TH)["status"] == "pass"
    assert qc._representation_check({"moi": 0.5}, lib, TH)["status"] == "warn"
    high = qc._representation_check({"moi": 1.2}, lib, TH)
    assert high["status"] == "fail" and "two or more" in high["why"]
    low = qc._representation_check({"moi": 0.01}, lib, TH)
    assert low["status"] == "warn" and "untransduced" in low["why"]


# ------------------------------------------------------------------ Spearman
def test_spearman_matches_scipy_including_ties():
    rng = np.random.default_rng(0)
    a = rng.normal(size=400)
    b = a + rng.normal(size=400)
    assert qc.spearman(a, b) == pytest.approx(stats.spearmanr(a, b).statistic, rel=1e-10)
    at, bt = np.round(a, 1), np.round(b, 1)                      # plenty of ties
    assert qc.spearman(at, bt) == pytest.approx(stats.spearmanr(at, bt).statistic, rel=1e-10)


def test_spearman_is_not_dragged_by_a_few_huge_values():
    rng = np.random.default_rng(1)
    a = rng.normal(size=300)
    b = rng.normal(size=300)                                      # unrelated
    a[0], b[0] = 400.0, 400.0                                     # one shared outlier
    assert qc.pearson(a, b) > 0.9                                 # Pearson is fooled
    assert abs(qc.spearman(a, b)) < 0.2                           # rank correlation is not


# ------------------------------------------------------------------ the checks on a whole screen
def _screen(guides_per_gene=4, n_ntc=200, n_essential=30, seed=0):
    """T0 + two treatment replicates; essential genes drop out, controls do not."""
    rng = np.random.default_rng(seed)
    genes, names = [], []
    for i in range(300):
        for k in range(guides_per_gene):
            genes.append(f"G{i}")
            names.append(f"G{i}_{k}")
    for k in range(n_ntc):
        genes.append("NonTargetingControl")
        names.append(f"NTC_{k}")
    genes = np.array(genes)
    base = rng.lognormal(6.0, 0.4, len(genes))
    t0 = rng.poisson(base).astype(float)
    eff = np.ones(len(genes))
    eff[np.isin(genes, [f"G{i}" for i in range(n_essential)])] = 0.1
    later = [rng.poisson(base * eff).astype(float) for _ in range(2)]
    X = np.column_stack([t0] + later)
    design = {"samples": [{"sample": "T0", "condition": "T0", "reference": True, "replicate": 1},
                          {"sample": "A", "condition": "treated", "reference": False, "replicate": 1},
                          {"sample": "B", "condition": "treated", "reference": False, "replicate": 2}],
              "comparisons": []}
    return X, list(names), list(genes), ["T0", "A", "B"], design


def _run(X, names, genes, samples, design, **kw):
    ctrl = np.array([g.startswith("NonTargeting") for g in genes])
    return qc.run_qc(X, names, genes, samples, design, control_mask=ctrl, **kw)


def test_guides_per_gene_check():
    r = _run(*_screen(guides_per_gene=4))
    c = by_id(r, "guides_gene")
    assert c["status"] == "pass" and c["median"] == 4 and c["n_thin"] == 0
    assert c["n_genes"] == 300                                     # controls are not counted as a gene
    c2 = by_id(_run(*_screen(guides_per_gene=2)), "guides_gene")   # thin, but some libraries do use 2
    assert c2["status"] == "warn" and c2["median"] == 2
    assert "off-target guide" in c2["why"]
    c3 = by_id(_run(*_screen(guides_per_gene=1)), "guides_gene")   # one guide cannot be a gene call
    assert c3["status"] == "fail"


def test_non_targeting_control_check():
    r = _run(*_screen(n_ntc=200))
    c = by_id(r, "ntc")
    assert c["status"] == "pass" and c["n"] == 200
    p = c["per_condition"][0]
    assert abs(p["median"]) < 0.25                                 # controls sit at zero
    # the noise floor is a real number, and well inside the planted 10-fold dropout it has to detect
    assert 0 < p["floor"] < abs(np.log2(0.1))
    assert "noise floor" in c["why"]


def test_too_few_controls_fails():
    c = by_id(_run(*_screen(n_ntc=10)), "ntc")
    assert c["status"] == "fail" and "too few" in c["why"]


def test_no_controls_is_not_applicable_with_advice():
    X, names, genes, samples, design = _screen(n_ntc=0)
    r = qc.run_qc(X, names, genes, samples, design, control_mask=np.zeros(len(genes), dtype=bool))
    c = by_id(r, "ntc")
    assert c["status"] == "na" and "Analyse page" in c["why"]


def test_shifted_controls_are_flagged_as_normalisation_not_biology():
    """When most of the library drops out, median-ratio normalisation pins the typical guide at zero and so
    lifts the non-targeting controls. Their fold change is then a property of the normalisation, not of the
    cells, and the app has to say so rather than report the controls as enriched."""
    X, names, genes, samples, design = _screen(n_essential=250)    # 250 of 300 genes essential
    c = by_id(_run(X, names, genes, samples, design), "ntc")
    assert c["status"] == "warn" and "normalisation" in c["why"]
    assert c["per_condition"][0]["median"] > qc.THRESHOLDS["ntc_shift"]


def test_representation_reaches_the_qc_result():
    X, names, genes, samples, design = _screen()
    r = _run(X, names, genes, samples, design, screen={"cells_bottleneck": 1000.0 * len(genes), "moi": 0.3,
                                                       "library_guides": float(len(genes))})
    c = by_id(r, "representation")
    assert c["status"] == "pass" and c["rows"][0]["per_guide"] == pytest.approx(1000)
    assert r["checks"][0]["id"] == "representation"                # shown first: it is the experiment, not the reads


def test_replicate_check_reports_spearman():
    r = _run(*_screen())
    p = by_id(r, "replicates")["pairs"][0]
    assert 0 < p["spearman"] <= 1 and p["spearman"] == p["spearman"]


# ------------------------------------------------------------------ the numbers the user types
def test_screen_fields_are_validated():
    P = session.Project()
    P.load_example("demo1")
    P.set_design({**P.design, "screen": {"moi": "0.3", "cells_bottleneck": "20,000,000",
                                         "selection": "2 ug/ml puromycin"}})
    assert P.screen["moi"] == 0.3 and P.screen["cells_bottleneck"] == 2e7
    assert P.screen["selection"].startswith("2 ug")
    for bad in ({"moi": "lots"}, {"moi": "-1"}, {"moi": "40"}):
        with pytest.raises(ValueError):
            P.set_design({**P.design, "screen": bad})
    P.set_design({**P.design, "screen": {}})
    assert P.screen == {}
