"""Version 1.2: thresholds the user changes are remembered between sessions, and can be put back."""
import json
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from crisprhts import arrayed_stats, prefs, qc, session  # noqa: E402


@pytest.fixture(autouse=True)
def clean():
    prefs.clear("pooled_thresholds")
    prefs.clear("arrayed_thresholds")
    yield
    prefs.clear("pooled_thresholds")
    prefs.clear("arrayed_thresholds")


def test_only_real_changes_to_known_thresholds_are_kept():
    saved = prefs.set("pooled_thresholds",
                      {"coverage_pass": 500, "auc_pass": qc.THRESHOLDS["auc_pass"], "nonsense": 3, "gini_ref_pass": "x"},
                      qc.THRESHOLDS)
    assert saved == {"coverage_pass": 500.0}          # default value, unknown name and non-number all dropped
    assert prefs.get("pooled_thresholds", qc.THRESHOLDS) == {"coverage_pass": 500.0}


def test_a_saved_value_that_is_no_longer_a_threshold_is_ignored():
    """A later version may drop or rename a threshold; the saved file must not resurrect it."""
    prefs.set("pooled_thresholds", {"coverage_pass": 500}, qc.THRESHOLDS)
    with open(prefs.path()) as fh:
        d = json.load(fh)
    d["pooled_thresholds"]["removed_in_a_later_version"] = 7
    with open(prefs.path(), "w") as fh:
        json.dump(d, fh)
    assert prefs.get("pooled_thresholds", qc.THRESHOLDS) == {"coverage_pass": 500.0}


def test_a_damaged_file_does_not_stop_the_app():
    with open(prefs.path(), "w") as fh:
        fh.write("{not json")
    assert prefs.get("pooled_thresholds", qc.THRESHOLDS) == {}
    prefs.set("pooled_thresholds", {"coverage_pass": 400}, qc.THRESHOLDS)    # and writing repairs it
    assert prefs.get("pooled_thresholds", qc.THRESHOLDS) == {"coverage_pass": 400.0}


def test_a_new_project_starts_from_the_saved_thresholds():
    prefs.set("pooled_thresholds", {"coverage_pass": 777}, qc.THRESHOLDS)
    P = session.Project()
    assert P.thresholds == {"coverage_pass": 777.0}
    P.reset()                                        # loading other data must not forget them
    assert P.thresholds == {"coverage_pass": 777.0}


def test_running_the_checks_saves_what_was_typed_and_reset_puts_it_back():
    P = session.Project()
    P.load_example("demo1")
    j = P.run_qc({"coverage_pass": 42})
    _wait(P, j)
    assert P.qc["thresholds"]["coverage_pass"] == 42
    assert prefs.get("pooled_thresholds", qc.THRESHOLDS) == {"coverage_pass": 42.0}
    assert P.qc["threshold_defaults"]["coverage_pass"] == qc.THRESHOLDS["coverage_pass"]
    assert "42" in P.qc["threshold_text"]["coverage"]
    # a second project in a later session picks them up
    assert session.Project().thresholds == {"coverage_pass": 42.0}
    prefs.clear("pooled_thresholds")
    assert session.Project().thresholds == {}


def test_arrayed_thresholds_are_remembered_too():
    from crisprhts.arrayed_session import ArrayedProject
    A = ArrayedProject(session.Project())
    A.load_example("kcviab")
    _wait(A.host, A.run_qc({"zprime_pass": 0.8, "rep_pass": 0.9}))
    assert A.qc["thresholds"]["zprime_pass"] == 0.8
    assert prefs.get("arrayed_thresholds", arrayed_stats.QC_THRESHOLDS) == {"zprime_pass": 0.8, "rep_pass": 0.9}
    # a later run with nothing posted still uses them, and the state hands the app both sets
    _wait(A.host, A.run_qc())
    assert A.qc["thresholds"]["zprime_pass"] == 0.8
    st = A.summary()
    assert st["qc_thresholds"]["zprime_pass"] == 0.8
    assert st["qc_threshold_defaults"]["zprime_pass"] == arrayed_stats.QC_THRESHOLDS["zprime_pass"]
    assert st["default_norm"] == "robust_z"


def test_clearing_one_section_leaves_the_other():
    prefs.set("pooled_thresholds", {"coverage_pass": 500}, qc.THRESHOLDS)
    prefs.set("arrayed_thresholds", {"zprime_pass": 0.6}, arrayed_stats.QC_THRESHOLDS)
    prefs.clear("pooled_thresholds")
    assert prefs.get("arrayed_thresholds", arrayed_stats.QC_THRESHOLDS) == {"zprime_pass": 0.6}


def test_the_confirmed_defaults():
    """The values the user confirmed on 2026-11-07."""
    assert qc.THRESHOLDS["coverage_pass"] == 300 and qc.THRESHOLDS["coverage_fail"] == 100
    assert qc.THRESHOLDS["mapped_pass"] == 0.80 and qc.THRESHOLDS["mapped_fail"] == 0.65
    assert qc.THRESHOLDS["rep_r_pass"] == 0.60 and qc.THRESHOLDS["rep_r_fail"] == 0.40
    assert qc.THRESHOLDS["cells_pass"] == 500 and qc.THRESHOLDS["guides_gene_pass"] == 4
    assert arrayed_stats.HIT_DEFAULTS["threshold"] == 2.0
    assert arrayed_stats.HIT_DEFAULTS["replicate_min"] == 1.0
    assert arrayed_stats.DEFAULT_NORM == "robust_z" and "robust_z" in arrayed_stats.NORM


def test_preferences_are_not_written_to_the_real_home_during_tests():
    assert os.environ.get("CH_PREFS")                      # conftest points it away
    assert "Application Support/CRISPR HTS Validator/preferences.json" not in str(prefs.path())


def _wait(P, j):
    import time
    while P.jobs[j["id"]].state == "running":
        time.sleep(0.05)
    assert P.jobs[j["id"]].state == "done", P.jobs[j["id"]].error
