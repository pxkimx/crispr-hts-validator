"""Arrayed screens: readers, plate statistics and QC on synthetic plates with known answers, and agreement with the
published KcViab analysis (cellHTS2, Analysis2003.txt)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from crisprhts import arrayed_io as aio, arrayed_stats as ast  # noqa: E402

KC = ROOT / "test-data/cellhts2/KcViab"


def synthetic(n_plates=4, reps=2, seed=0, edge=0.0, colfx=None, hits=(), scale=None, hit_level=0.4):
    """384-well plates, sample wells N(1000, 50), neg controls (B01, C01) at 1000, pos (B02, C02) at 200.
    `hits`: wells (plate, well) reduced to 40%. `scale`: per-plate multipliers."""
    rng = np.random.default_rng(seed)
    rows = []
    for p in range(1, n_plates + 1):
        for r in range(1, reps + 1):
            for i in range(16):
                for j in range(24):
                    w = aio.well_name(i, j)
                    content = {"B01": "neg", "C01": "neg", "B02": "pos", "C02": "pos"}.get(w, "sample")
                    v = 200.0 if content == "pos" else rng.normal(1000, 50)
                    if content == "neg":
                        v = rng.normal(1000, 50)
                    if edge and (i in (0, 15) or j in (0, 23)):
                        v *= 1 + edge
                    if colfx is not None:
                        v *= colfx[j]
                    if (p, w) in hits:
                        v *= hit_level
                    if scale:
                        v *= scale[p - 1]
                    rows.append((p, r, 1, w, v, content, f"G{p}_{w}"))
    df = pd.DataFrame(rows, columns=["plate", "replicate", "channel", "well", "value", "content", "gene"])
    return aio._finish(df, 384, "synthetic", [])


# ------------------------------------------------------------------ readers
def test_well_names_and_formats():
    assert aio.parse_well("B03") == (1, 2) and aio.parse_well("p24") == (15, 23)
    assert aio.well_name(0, 0) == "A01" and aio.well_name(31, 47) == "AF48"
    assert aio.plate_format(8, 12) == 96 and aio.plate_format(16, 24) == 384 and aio.plate_format(32, 48) == 1536


def test_plate_conf_wildcards_later_rules_win(tmp_path):
    f = tmp_path / "Plateconf.txt"
    f.write_text("Wells: 384\nPlates: 3\nPlate\tWell\tContent\n*\t*\tsample\n*\tA0[1-2]\tother\n*\tB01\tneg\n2\tB01\tpos\n")
    content_of, n, _ = aio.read_plate_conf(f)
    assert n == 384
    assert content_of(1, "A02") == "other" and content_of(1, "B01") == "neg" and content_of(2, "B01") == "pos"
    assert content_of(3, "C05") == "sample"


def test_grid_reader_picks_results_block(tmp_path):
    hdr = "," + ",".join(f"{c:02d}" for c in range(1, 25)) + ","
    body = "\n".join(f"{chr(65 + i)}," + ",".join(str(i * 100 + j) for j in range(24)) + "," for i in range(16))
    f = tmp_path / "plate.csv"
    f.write_text(f"Calculated results: x\n{hdr}\n{body}\n\nResults for LUM\n{hdr}\n" + body.replace(",1", ",9") + "\n")
    blocks = aio.read_grids(f)
    assert len(blocks) == 2 and blocks[0]["values"].shape == (16, 24)
    assert aio._default_block([b["title"] for b in blocks]) == 1


def test_well_table_with_replicate_columns(tmp_path):
    f = tmp_path / "wells.csv"
    lines = ["Plate,Well,Gene,Type,rep1,rep2"]
    for i in range(16):
        for j in range(24):
            w = aio.well_name(i, j)
            lines.append(f"P1,{w},G{i}_{j},{'negative' if w == 'B01' else 'sample'},{100 + i},{101 + i}")
    lines.append("P1,ZZ9,x,sample,1,1")
    f.write_text("\n".join(lines))
    s = aio.load_well_table(f)
    assert s.summary()["replicates"] == [1, 2] and s.n_wells == 384
    assert (s.wells["content"] == "neg").sum() == 2          # B01 in both replicates
    assert any("skipped" in n for n in s.notes)


# ------------------------------------------------------------------ statistics
def test_median_polish_removes_additive_effects():
    rows, cols = np.arange(16.0)[:, None], np.arange(24.0)[None, :] * 2
    M = 5 + rows + cols
    resid, overall, re_, ce = ast.median_polish(M)
    assert np.allclose(resid, 0, atol=1e-9)
    assert np.allclose(overall + re_[:, None] + ce[None, :], M)


def test_zprime_and_ssmd_known_values():
    pos = [10.0, 12.0]
    neg = [100.0, 104.0]
    sp, sn = np.std(pos, ddof=1), np.std(neg, ddof=1)
    assert ast.zprime(pos, neg) == pytest.approx(1 - 3 * (sp + sn) / 91)
    assert ast.ssmd(pos, neg) == pytest.approx((11 - 102) / np.sqrt(sp ** 2 + sn ** 2))
    assert np.isnan(ast.zprime([1.0], [5.0, 6.0]))


def test_normalisations_on_synthetic_plate():
    s = synthetic(n_plates=2, scale=[1.0, 3.0])
    w = s.wells
    med = ast.normalise(s, "median")
    samp = w["content"] == "sample"
    assert np.allclose(med[samp].groupby([w["plate"], w["replicate"]]).median(), 1.0)
    poc = ast.normalise(s, "poc")
    negs = w["content"] == "neg"
    assert np.allclose(poc[negs].groupby([w["plate"], w["replicate"]]).mean(), 100.0)
    npi = ast.normalise(s, "npi")
    assert np.allclose(npi[negs].groupby([w["plate"], w["replicate"]]).mean(), 0.0)
    assert np.allclose(npi[w["content"] == "pos"].groupby([w["plate"], w["replicate"]]).mean(), 100.0)
    rz = ast.normalise(s, "robust_z")
    assert abs(float(np.median(rz[samp]))) < 1e-9


def test_bscore_removes_column_gradient_and_finds_hits():
    """A ±30% column gradient inflates the plate's spread, so moderate hits (−20%) disappear in a plain robust z;
    the B-score removes the gradient first and finds them."""
    colfx = np.linspace(0.7, 1.3, 24)
    hits = {(1, "H12"), (1, "D14"), (2, "K11"), (2, "F13")}
    s = synthetic(n_plates=2, colfx=colfx, hits=hits, hit_level=0.8)
    found = {}
    for m in ("bscore", "robust_z"):
        r = ast.analyse(s, m)
        found[m] = {(p, w) for p, w, h in zip(r["wells"]["plate"], r["wells"]["well"], r["wells"]["hit"]) if h}
    assert hits <= found["bscore"]
    assert not (hits & found["robust_z"])
    assert len(found["bscore"] - hits) <= 10        # noise-level false calls only


def test_hit_calls_need_replicate_agreement():
    s = synthetic(n_plates=1, hits={(1, "E05")})
    w = s.wells
    # a well that drops in replicate 1 only
    w.loc[(w["plate"] == 1) & (w["well"] == "F06") & (w["replicate"] == 1), "value"] *= 0.2
    r = ast.analyse(s, "robust_z", threshold=2, replicate_min=1)
    t = r["wells"].set_index("well")
    assert bool(t.at["E05", "hit"]) and not bool(t.at["F06", "hit"])
    assert t.at["F06", "score"] >= 2          # strong on average, but not confirmed by replicate 2


def test_gene_calls():
    t = pd.DataFrame({"content": ["sample"] * 5, "gene": ["A", "A", "B", "B", "B"], "well": list("abcde"),
                      "hit": [True, False, True, False, False], "score": [3, 1, 2.5, 0, 0]})
    g = ast.gene_calls(t, 0.5).set_index("gene")
    assert bool(g.at["A", "hit"]) and not bool(g.at["B", "hit"])


def test_qc_known_answers():
    clean = ast.run_qc(synthetic(n_plates=3))
    by = {c["id"]: c for c in clean["checks"]}
    assert by["plates"]["status"] == "pass"
    assert all(p["zprime"] > 0.5 for p in by["controls"]["per_plate"])      # 2 wells of each control per plate
    assert by["edge"]["status"] == "pass" and by["rowcol"]["status"] == "pass"
    assert by["replicates"]["status"] == "fail"                               # pure noise: replicates disagree
    bad = ast.run_qc(synthetic(n_plates=3, edge=0.3, scale=[1, 2.5, 6]))
    by = {c["id"]: c for c in bad["checks"]}
    assert by["plates"]["status"] == "fail" and by["plates"]["value"] == pytest.approx(6, rel=0.05)
    assert by["edge"]["status"] == "fail" and by["edge"]["value"] == pytest.approx(0.3, abs=0.03)


# ------------------------------------------------------------------ KcViab against the 2003 analysis
@pytest.fixture(scope="module")
def kc():
    if not KC.exists():
        pytest.skip("KcViab test data not present")
    files = {p.name: p for p in KC.iterdir() if p.name != "Analysis2003.txt" and not p.name.startswith("old-")}
    return aio.detect_and_load(files, "KcViab"), pd.read_csv(KC / "Analysis2003.txt", sep="\t")


def test_kcviab_reading(kc):
    s, _ = kc
    sm = s.summary()
    assert sm["plates"] == 57 and sm["replicates"] == [1, 2] and len(s.wells) == 57 * 384 * 2
    assert sm["contents"] == {"sample": 21660, "other": 114, "neg": 57, "pos": 57}
    assert sm["flagged"] == 48
    assert any("10 row(s)" in n for n in sm["notes"])           # the ragged annotation rows from line 21121


def test_kcviab_classic_zscore_reproduces_published_scores(kc):
    s, pub = kc
    r = ast.analyse(s, "zscore", direction="lower", summary="mean", threshold=2, replicate_min=-1e9, use_flagged=True)
    m = r["wells"].merge(pub, left_on=["plate", "well"], right_on=["Plate", "Well"])
    assert len(m) == len(pub) == 21060
    assert np.corrcoef(m["score"], m["Score"])[0, 1] > 0.9999
    assert (m["score"].round(2) == m["Score"]).mean() > 0.999
    pub_hit, our_hit = m["Score"] >= 2, m["hit"].astype(bool)
    assert (pub_hit & our_hit).sum() >= 0.995 * pub_hit.sum()
    assert (our_hit & ~pub_hit).sum() <= 3                       # only rounding at the threshold


def test_kcviab_robust_default_recovers_published_hits(kc):
    s, pub = kc
    r = ast.analyse(s)                                          # robust z, lower signal, mean ≥ 2, each replicate ≥ 1
    m = r["wells"].merge(pub, left_on=["plate", "well"], right_on=["Plate", "Well"])
    pub_hit, our_hit = m["Score"] >= 2, m["hit"].astype(bool)
    assert np.corrcoef(m["score"].fillna(0), m["Score"])[0, 1] > 0.97
    assert (pub_hit & our_hit).sum() / pub_hit.sum() > 0.98
