"""Version 1.6: comparing screens held in two projects, and so in two different libraries.

The hard part is not the plumbing, it is that two libraries target a gene with different sequences, so the
guide-by-guide cancellation that makes the within-project comparison sensitive is gone. These tests hold the
method to that: it has to recover planted differences when the data can support them, refuse when it cannot,
and say which case it is in.
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from crisprhts import compare, crosscompare, pooled_io, qc  # noqa: E402

TK = ROOT / "test-data/hart-bagel/tkov3-example/HAP1-TKOv3-EXAMPLE.txt"
HAP1 = ROOT / "test-data/hart-bagel/reads_hap1.txt"
FED = {"name": "T18 fed", "treatment": ["T18_A", "T18_B", "T18_C"], "control": ["T0_-"]}
STARVED = {"name": "T18 starved", "treatment": ["T18_A_Starved", "T18_B_Starved", "T18_C_Starved"],
           "control": ["T0_-"]}


# ------------------------------------------------------------------ two libraries of known truth
def two_libraries(n_guides=4, seed=0, n_other=1200, shared=-2.5, specific=-3.5, n_spec=30, depth_b=1.0):
    """Two screens with no guide in common. Guide efficiency is modelled explicitly, because it is the one
    thing a cross-library comparison cannot cancel: each guide realises only part of its gene's effect, drawn
    independently in each library."""
    rng = np.random.default_rng(seed)
    ceg, neg = qc.load_reference_sets()
    ess, non = sorted(ceg)[:200], sorted(neg)[:200]
    others = [f"G{i}" for i in range(n_other)]
    a_spec, b_spec = others[:n_spec], others[n_spec:2 * n_spec]
    names = ess + non + others
    truth = {}
    for g in names:
        base = shared if g in set(ess) else 0.0
        truth[g] = (base + (specific if g in set(a_spec) else 0.0),
                    (base + (specific if g in set(b_spec) else 0.0)) * depth_b)

    def lib(tag, which):
        genes = np.repeat(names, n_guides)
        guides = [f"{tag}_{g}_{k}" for g in names for k in range(n_guides)]
        eff = rng.beta(2.2, 1.0, len(genes))                    # most guides cut well, some are duds
        t0 = rng.lognormal(7.0, 0.35, len(genes))
        lfc = np.array([truth[g][which] for g in genes]) * eff + rng.normal(0, 0.30, len(genes))
        ct = np.stack([rng.poisson(t0 * np.exp2(lfc) * np.exp(rng.normal(0, .18, len(genes))))
                       for _ in range(3)], 1)
        X = np.concatenate([rng.poisson(t0)[:, None], ct], 1).astype(float)
        samples = [f"{tag}_T0"] + [f"{tag}_R{i}" for i in range(3)]
        return {"counts": X, "guides": guides, "genes": list(genes), "samples": samples, "project": tag,
                "comparison": {"name": f"{tag} screen", "treatment": samples[1:], "control": [samples[0]]}}
    return lib("libA", 0), lib("libB", 1), set(a_spec), set(b_spec), set(ess), set(non)


def test_planted_library_specific_genes_are_recovered():
    A, B, aset, bset, _, _ = two_libraries(n_guides=6, seed=1)
    r = crosscompare.cross_compare(A, B, thresholds={"xc_min_delta": 0.5})
    g = r["genes"].set_index("gene")
    called = set(g.index[g["different"]])
    assert len(called & bset) >= 0.7 * len(bset), (len(called & bset), r["power"])
    assert len(called & aset) >= 0.7 * len(aset), (len(called & aset), r["power"])
    for name in sorted(bset)[:5]:
        assert g.loc[name, "delta"] < 0 and g.loc[name, "call"] == "specific to B"
    for name in sorted(aset)[:5]:
        assert g.loc[name, "delta"] > 0 and g.loc[name, "call"] == "specific to A"
    # and the genes that were not planted are left alone
    other = g[~g.index.isin(aset | bset)]
    assert other["different"].mean() < 0.01, other[other["different"]].head().to_dict()


def test_power_improves_with_guides_per_gene_and_is_reported_honestly():
    """The resolution the app prints has to track what it can really find, or it is worse than no number."""
    seen = []
    for ng in (4, 10):
        A, B, _, bset, _, _ = two_libraries(n_guides=ng, seed=20 + ng)
        r = crosscompare.cross_compare(A, B, thresholds={"xc_min_delta": 0.0})
        g = r["genes"].set_index("gene")
        seen.append((r["power"]["resolution"], len(set(g.index[g["different"]]) & bset),
                     float(g.loc[sorted(bset), "delta"].abs().median())))
    (res4, found4, d4), (res10, found10, d10) = seen
    assert res10 < res4, seen                      # deeper library, finer resolution
    assert found10 >= found4, seen                 # and more of the planted genes found
    # the printed resolution is in the right place: the planted difference clears it, and most are found
    assert d10 > res10 and found10 >= 0.8 * 30, seen


def test_a_difference_in_dynamic_range_alone_is_absorbed_by_the_scale():
    """One screen with far more dynamic range than the other is the commonest false-positive trap."""
    A, B, _, _, ess, non = two_libraries(n_guides=6, seed=7, n_spec=0, depth_b=1.7)
    r = crosscompare.cross_compare(A, B)
    g = r["genes"].set_index("gene")
    assert abs(r["a"]["scale"]["span"] - r["b"]["scale"]["span"]) > 0.5     # the raw screens really do differ
    for col in ("scaled_a", "scaled_b"):
        assert g[col][[x in ess for x in g.index]].median() == pytest.approx(-1.0, abs=0.05)
        assert g[col][[x in non for x in g.index]].median() == pytest.approx(0.0, abs=0.05)
    assert int(g["different"].sum()) == 0, g[g["different"]][["delta", "fdr"]].to_dict()


def test_the_scale_is_required_and_a_thin_reference_set_is_refused():
    A, B, _, _, _, _ = two_libraries(n_guides=4, seed=3)
    thin = dict(B)
    thin["genes"] = [g if i % 97 == 0 else f"X{i}" for i, g in enumerate(B["genes"])]
    with pytest.raises(ValueError, match="reference genes"):
        crosscompare.cross_compare(A, thin)
    none = dict(B)
    none["genes"] = [f"X{i}" for i in range(len(B["genes"]))]
    with pytest.raises(ValueError, match="common essential-gene scale|reference genes"):
        crosscompare.cross_compare(A, none)


def test_minimum_guides_applies_to_both_libraries():
    A, B, _, _, _, _ = two_libraries(n_guides=4, seed=5)
    thin = dict(B)
    keep = [i for i, _ in enumerate(B["genes"]) if i % 4 < 2]          # two guides a gene on one side only
    thin["counts"] = B["counts"][keep]
    thin["guides"] = [B["guides"][i] for i in keep]
    thin["genes"] = [B["genes"][i] for i in keep]
    with pytest.raises(ValueError, match="at least 3 guides"):
        crosscompare.cross_compare(A, thin, thresholds={"xc_min_guides": 3})
    r = crosscompare.cross_compare(A, thin, thresholds={"xc_min_guides": 2})
    assert len(r["genes"]) > 500 and int(r["genes"]["guides_b"].min()) == 2


# ------------------------------------------------------------------ matching gene names
def test_gene_names_are_matched_without_guessing():
    m = crosscompare.match_genes(["TP53", " brca1 ", "ENSG00000141510.12"], ["tp53", "BRCA1", "ENSG00000141510"])
    assert set(m["shared"]) == {"TP53", "BRCA1", "ENSG00000141510"}
    assert m["kind_a"] == "gene symbols"
    m2 = crosscompare.match_genes(["TP53", "BRCA1"], ["ENSG00000141510", "ENSG00000012048"])
    assert not m2["shared"] and "Ensembl" in m2["note"] and "mapping file" in m2["note"]


def test_a_mapping_file_is_read_either_way_round():
    a, b = ["TP53", "BRCA1", "MYC"], ["ENSG00000141510", "ENSG00000012048", "ENSG00000136997"]
    fwd = {"TP53": "ENSG00000141510", "BRCA1": "ENSG00000012048"}
    m = crosscompare.match_genes(a, b, fwd)
    assert set(m["shared"]) == {"TP53", "BRCA1"} and m["n_mapped"] == 2
    rev = {v: k for k, v in fwd.items()}
    assert set(crosscompare.match_genes(a, b, rev)["shared"]) == {"TP53", "BRCA1"}


def test_the_mapping_file_reader_takes_csv_or_tab_with_or_without_a_header():
    for text in (b"gene,other\nTP53,ENSG1\nBRCA1,ENSG2\n", b"TP53\tENSG1\nBRCA1\tENSG2\n"):
        assert crosscompare.read_mapping(text) == {"TP53": "ENSG1", "BRCA1": "ENSG2"}
    with pytest.raises(ValueError):
        crosscompare.read_mapping(b"onlyonecolumn\nagain\n")


@pytest.mark.skipif(not (TK.exists() and HAP1.exists()), reason="test data not present")
def test_the_same_library_is_recognised_through_two_id_conventions():
    """reads_hap1 keys guides on the bare 20-mer and the TKOv3 table on GENE_20mer, so no ID string is in
    both — but they are one library, and the user has to be told so rather than run a weaker analysis."""
    t1, t2 = pooled_io.read_count_table(HAP1), pooled_io.read_count_table(TK)
    assert not (set(map(str, t1.guides)) & set(map(str, t2.guides)))       # nothing matches on the ID
    sg = crosscompare.shared_guides(t1.guides, t2.guides)
    assert sg["fraction"] > 0.99 and sg["by_sequence"]
    assert crosscompare.guide_key("EEF2_TGCAGACGGAGACAGTGCTG") == "TGCAGACGGAGACAGTGCTG"
    assert crosscompare.guide_key("AAACAAACCACCGAAACCCT") == "AAACAAACCACCGAAACCCT"
    assert crosscompare.guide_key("BRDN0001") == "BRDN0001"                # no protospacer to find


# ------------------------------------------------------------------ against real data
@pytest.mark.skipif(not TK.exists(), reason="test data not present")
def test_splitting_a_real_library_in_half_shows_the_cost_and_does_not_manufacture_results():
    """The hardest case there is: TKOv3's four guides per gene split two and two, so the fed screen is
    measured with guides the starved screen never sees. The biology is the one the within-project comparison
    already answered, so the two can be held against each other — and the honest answer at two guides a side
    is that almost nothing can be called, which is what the resolution figure has to say."""
    t = pooled_io.read_count_table(TK)
    genes = np.asarray(t.genes, dtype=object).astype(str)
    rank = pd.Series(np.arange(len(genes))).groupby(genes).cumcount().to_numpy()
    h1, h2 = rank % 2 == 0, rank % 2 == 1
    guides = np.asarray(t.guides, dtype=object)
    assert not (set(guides[h1].astype(str)) & set(guides[h2].astype(str)))

    def side(mask, comp, name):
        return {"counts": t.counts[mask], "guides": list(guides[mask]), "genes": list(genes[mask]),
                "samples": t.samples, "comparison": comp, "project": name}
    r = crosscompare.cross_compare(side(h1, FED, "half 1"), side(h2, STARVED, "half 2"),
                                   thresholds={"xc_min_guides": 2})
    g = r["genes"].set_index("gene")
    assert r["shared_guides"]["n"] == 0
    assert r["a"]["scale"]["ok"] and r["b"]["scale"]["ok"]
    # the ranking still carries the biology: starvation makes mTORC1 restraint essential
    top = list(g.index[:15])
    assert "TSC2" in top and "TSC1" in top, top
    # but two guides a side cannot clear a genome-wide bar, and the app says so rather than pretending
    assert r["power"]["resolution"] > 2.0, r["power"]
    assert int(g["different"].sum()) == 0, g[g["different"]][["delta", "fdr"]].to_dict()
    # the within-project route, which pairs each guide with itself, does find them
    within = compare.compare(t.counts, t.guides, t.genes, t.samples, FED, STARVED)["genes"].set_index("gene")
    assert {"TSC1", "TSC2"} <= set(within.index[within["different"]])


@pytest.mark.skipif(not (TK.exists() and HAP1.exists()), reason="test data not present")
def test_one_screen_in_two_tables_yields_no_differences():
    """The strongest negative control available: the same HAP1 screen read from two count tables that share
    no guide ID. Anything called here would be an artefact of the method."""
    t1, t2 = pooled_io.read_count_table(TK), pooled_io.read_count_table(HAP1)
    A = {"counts": t1.counts, "guides": t1.guides, "genes": t1.genes, "samples": t1.samples,
         "project": "TKOv3 table", "comparison": FED}
    B = {"counts": t2.counts, "guides": t2.guides, "genes": t2.genes, "samples": t2.samples,
         "project": "BAGEL table",
         "comparison": {"name": "HAP1 T18", "treatment": ["HAP1_T18A", "HAP1_T18B", "HAP1_T18C"],
                        "control": ["HAP1_T0"]}}
    r = crosscompare.cross_compare(A, B)
    c = r["concordance"]
    assert c["pearson"] > 0.95 and c["jaccard"] > 0.9
    assert c["n_different"] == 0
    notes = {n["id"]: n for n in crosscompare.library_notes(r)}
    assert notes["same_library"]["status"] == "warn"           # and it says to use the paired route instead
    assert notes["over_agree"]["status"] == "warn"             # the two are not independent measurements


# ------------------------------------------------------------------ through the app
def _project_with_cross():
    from crisprhts import session
    P = session.Project()
    P.load_example("cross")

    def run(start):
        j = start()
        while P.jobs[j["id"]].state == "running":
            time.sleep(0.05)
        assert P.jobs[j["id"]].state == "done", P.jobs[j["id"]].error
    run(P.run_qc)
    run(lambda: P.run_cross({"a": "BATCH1_HAP1_T18", "b_name": "HAP1 T18 (BAGEL table)",
                             "b_treatment": ["HAP1_T18A", "HAP1_T18B", "HAP1_T18C"],
                             "b_control": ["HAP1_T0"]}))
    return P


@pytest.mark.skipif(not (TK.exists() and HAP1.exists()), reason="test data not present")
def test_the_example_loads_both_projects_and_reports_the_overlap():
    from crisprhts import session
    P = session.Project()
    P.load_example("cross")
    info = P.partner_info()
    assert info["guides"] > 70000 and info["shared"]["n"] > 17000
    assert info["comparisons"], "the other project's design should be guessed for the user"
    s = P.summary()
    assert s["partner"]["name"] and s["has_cross"] is False
    assert set(s["cross_defaults"]) == set(crosscompare.DEFAULTS)


@pytest.mark.skipif(not (TK.exists() and HAP1.exists()), reason="test data not present")
def test_view_gene_panel_and_exports():
    P = _project_with_cross()
    v = P.cross_view()
    assert v["total"] > 17000 and len(v["points"]) == v["total"]
    assert all(isinstance(p[4], float) for p in v["points"][:50])          # FDR keeps full precision
    assert {n["id"] for n in v["notes"]} >= {"overlap", "guides", "resolution", "null", "scale"}
    assert v["power"]["resolution"] > 0 and v["power"]["median_df"] > 0
    hit = P.cross_view(search="TSC2")
    assert any(r[0] == "TSC2" for r in hit["rows"])
    g = P.cross_gene("TSC2")
    assert g["guides_a"] and g["guides_b"]
    assert all("scaled" in x for x in g["guides_a"])

    body, fname, _ = P.export("cross", fmt="csv")
    assert b"TSC2" in body and b"screen_a" in body and fname.endswith(".csv")
    book, _, _ = P.export("all", fmt="xlsx")
    import io

    import openpyxl
    assert "cross-project" in openpyxl.load_workbook(io.BytesIO(book)).sheetnames


@pytest.mark.skipif(not (TK.exists() and HAP1.exists()), reason="test data not present")
def test_the_report_carries_the_cross_project_section():
    from crisprhts import report
    P = _project_with_cross()
    pdf = report.pooled_report(P)
    assert pdf[:5] == b"%PDF-" and len(pdf) > 50000


@pytest.mark.skipif(not (TK.exists() and HAP1.exists()), reason="test data not present")
def test_bad_requests_are_refused_clearly():
    from crisprhts import session
    P = session.Project()
    P.load_example("cross")
    with pytest.raises(ValueError, match="this project's comparisons"):
        P.run_cross({"a": "nope", "b_treatment": ["HAP1_T18A"], "b_control": ["HAP1_T0"]})
    with pytest.raises(ValueError, match="treated and reference"):
        P.run_cross({"a": "BATCH1_HAP1_T18", "b_treatment": [], "b_control": ["HAP1_T0"]})
    with pytest.raises(ValueError, match="cannot be both"):
        P.run_cross({"a": "BATCH1_HAP1_T18", "b_treatment": ["HAP1_T0"], "b_control": ["HAP1_T0"]})
    with pytest.raises(ValueError, match="Run a cross-project comparison first"):
        P.cross_view()
    P2 = session.Project()
    P2.load_example("tkov3")
    with pytest.raises(ValueError, match="other project's count table"):
        P2.run_cross({"a": "BATCH1_HAP1_T18", "b_treatment": ["x"], "b_control": ["y"]})
