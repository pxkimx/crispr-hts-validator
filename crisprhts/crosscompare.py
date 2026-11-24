"""Compare two screens that come from different projects — and, in general, from different libraries.

Why this needs its own machinery
--------------------------------
`compare.py` puts two screens from one count table side by side. There every guide appears in both screens,
so the difference B - A is taken guide by guide and each guide's own cutting efficiency cancels. That
cancellation is the single largest reason the within-project comparison is sensitive, and it is exactly what
is lost here: two libraries target the same gene with different sequences, so nothing pairs up.

What is left to build on is the common scale. Hart's reference sets were defined for this purpose: a screen
is scaled so that its median core-essential gene (CEGv2) sits at -1 and its median non-essential gene (NEGv1)
at 0, and a scaled effect of -1 then means "as essential as a ribosomal subunit" whichever library, timepoint
or cell line produced it (Hart et al. 2015; the same device as CERES and Chronos). Here that scale is not a
convenience, it is the only bridge between the two screens — so if either screen cannot be placed on it, this
module refuses to compare rather than return numbers that have no shared meaning.

The test
--------
Each gene's effect is the median of its guides, as elsewhere in the app, and the difference is

    delta = scaled(B) - scaled(A)

Its uncertainty now has to come from inside each screen, because there is no paired difference to z-score.
And it has to be the uncertainty about the gene's *true* effect given the handful of guides this library
happens to carry — which is dominated by guide efficiency, not by read noise. Splitting a real library in two
shows how brutal that is: EEF2's guides in one half of TKOv3 read -0.6 and -0.2 (they barely cut) while its
guides in the other half read -5.2 and -7.3. Measured effect -0.4 against -5.3, from one cell line at one
timepoint. Any test that treats a gene's measured effect as precise will call that a biological difference.

The spread of a gene's guides is what carries this, so for each gene in each screen:

  * the spread of its own guides (ordinary SD, which an outlying guide inflates — conservative on purpose);
  * the spread a gene of that effect size usually shows, read off a trend fitted across the screen. Guide
    efficiency only shows itself where there is something to measure: every guide of a harmless gene reads
    near zero whether it cuts or not, while against an essential gene a dead guide reads 0 and a good one -6.
    So the spread climbs steeply with the effect, and shrinking every gene toward one global number would make
    exactly the discordant essential genes look far more certain than they are. This is the mean-variance
    trend that MAGeCK and limma-voom fit, here on the gene effect rather than the count mean;
  * a floor from the typical guide spread at that read abundance, estimated in a sliding window of guides
    ordered by abundance and never allowed to fall as abundance rises — drugZ's empirical-Bayes step, reused —
    so a gene carried only by faint guides cannot look precise either;
  * the three combined as sd^2 = w*own^2 + (1-w)*prior^2 with w = (n-1)/((n-1) + d0), limma's weighting by
    degrees of freedom, so a two-guide gene leans mostly on the trend and a twenty-guide gene speaks for
    itself. How strong that prior is, d0, is not a number anyone should have to guess: it is estimated from
    the data by Smyth's (2004) moment method — how far the observed guide spreads scatter about the trend
    tells you how much the trend can be trusted — and the figure it arrives at is reported;
  * the standard error of the gene's effect = sqrt(pi/2) * sd / sqrt(n), the median's SE for normal data.

Then t = delta / sqrt(se_A^2 + se_B^2), referred to a t distribution on (n_A - 1) + (n_B - 1) + d0
degrees of freedom — limma's moderated t. A normal approximation would be wrong here and wrong in the
dangerous direction: with four guides a side the standard error is itself barely known, and treating it as
exact manufactures significance. The t distribution is what stops a two-guide gene from ever looking certain.

That z is not standard normal: the scale factor is itself estimated, the median's SE is asymptotic, and guides
within a gene are not independent of the gene's effect. So the null is calibrated from the data in Efron's
sense — centred on the median of z and widened by its MAD/0.6745 — and the widening is **floored at 1**, so
the calibration may only ever make the test more conservative, never more liberal. How wide the calibration
had to be is reported, and so is the spread of t among non-essential genes, which should land near 1.

The scale is not free either. Both medians are estimates, so each screen's scale carries a standard error of
its own, and that error is propagated into every gene: a relative error in the span becomes an absolute error
proportional to how deep the gene is, which is why a thin library's deepest genes are the least certain rather
than the most. Below MIN_REFERENCE genes in either reference set the comparison is refused outright — with
seventeen core-essential genes the error in the scale alone exceeds the difference the test could resolve, and
for two different libraries the scale is the whole of the bridge, not a cosmetic step.

Two honest limits, both shown in the app
----------------------------------------
  * This is a blunter instrument than the within-project comparison. With four guides a side, the smallest
    difference it can actually detect is around one whole unit of the scale — the distance from a harmless
    gene to halfway down to a core-essential one. The app computes that number from the data in front of it
    and prints it, so a null result can be read as "no difference this big" rather than "no difference".
  * Gene names are matched exactly, after normalising case, whitespace and Ensembl version suffixes, plus an
    optional two-column mapping file. No alias table is bundled and nothing is downloaded, so names that do
    not line up are reported and diagnosed, never guessed at.

Thresholds here are scientific choices and are listed in the app for Paul to confirm.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import special, stats

from . import qc
from .compare import _bh, _concordance, gene_effect, guide_lfc, scale_factors, scaled
from .drugz import _empirical_bayes

DATA = Path(__file__).resolve().parent / "data"

# Scientific choices, confirmed by the user on 2026-11-24 and stated in Methods and on the page.
#   xc_min_delta 1.0  deliberately twice the within-project 0.5: a cross-library difference smaller than that
#                     is inside the systematic error between two libraries, and below the size this test can
#                     resolve with four guides a side anyway — on the real HAP1 pair, 1.83 (see "resolution").
#                     Considered and rejected: 0.5 for one rule across both routes, which would leave the floor
#                     doing nothing while the FDR bar did all the work; and a floor that rises to each pair's
#                     own resolution, which would stop two runs being comparable.
DEFAULTS = {
    "xc_fdr": 0.05,            # FDR below which a difference is called
    "xc_min_delta": 1.0,       # and the smallest difference in scaled effect worth reporting
    "xc_min_guides": 3,        # guides a gene needs in each library
    "xc_dependency": -0.5,     # scaled effect at or below which a gene is a dependency in a screen
}

MEDIAN_SE = float(np.sqrt(np.pi / 2.0))     # SE of a median against a mean's, for normal data
POWER_Z = 2.80                              # 80% power, two-sided 5%, for the nominal statement
# Confirmed 2026-11-24: refuse rather than warn. The bundled leukaemia table carries 17 core-essential and 19
# non-essential genes, and at that size the scale's own error is larger than any difference the test could find,
# so a table of numbers would be worse than none. The refusal names the counts and points at the paired route.
MIN_REFERENCE = 30                          # reference genes each screen needs before its scale is load-bearing

_SPACER = re.compile(r"^[ACGT]{19,25}$")
_ENSEMBL = re.compile(r"^ENS[A-Z]*G\d{6,}")
_ENTREZ = re.compile(r"^\d+$")
_SPACE = re.compile(r"\s+")


# ---------------------------------------------------------------- gene names
def normalise_name(x) -> str:
    """One gene name reduced to a form two libraries can agree on: no surrounding space or quotes, upper
    case, and an Ensembl version suffix dropped (ENSG00000141510.12 and ENSG00000141510 are one gene)."""
    s = _SPACE.sub("", str(x)).strip().strip("'\"").upper()
    if _ENSEMBL.match(s) and "." in s:
        s = s.split(".", 1)[0]
    return s


def id_kind(names) -> str:
    """What kind of identifier a library is using, so a failure to match can be explained."""
    n = [normalise_name(x) for x in names if str(x).strip()]
    if not n:
        return "nothing"
    if sum(1 for x in n if _ENSEMBL.match(x)) / len(n) > 0.6:
        return "Ensembl gene IDs"
    if sum(1 for x in n if _ENTREZ.match(x)) / len(n) > 0.6:
        return "Entrez gene IDs"
    return "gene symbols"


def guide_key(x) -> str:
    """The protospacer inside a guide ID, when there is one to find.

    Guide IDs are not standardised: the Hart lab's BAGEL tables key on the bare 20-mer, their TKOv3 table on
    GENE_20mer, and other libraries use a catalogue number. Comparing ID strings therefore says nothing about
    whether two tables hold the same guides, which is worth knowing: if they do, the paired route on the
    Compare page is far more sensitive and the user should be told to use it.
    """
    s = str(x).strip().upper()
    if _SPACER.match(s):
        return s
    for part in re.split(r"[_\-:|.]", s):
        if _SPACER.match(part):
            return part
    return s


def shared_guides(guides_a, guides_b) -> dict:
    """How much of the two libraries is really the same guides, judged on protospacers where possible."""
    ka = {guide_key(g) for g in guides_a}
    kb = {guide_key(g) for g in guides_b}
    both = ka & kb
    by_seq = sum(1 for k in both if _SPACER.match(k))
    return {"n": len(both), "fraction": len(both) / max(1, min(len(ka), len(kb))),
            "by_sequence": by_seq == len(both) and by_seq > 0}


def read_mapping(src) -> dict:
    """A two-column gene mapping: a name in one library, the matching name in the other. Header optional."""
    from .pooled_io import _delim, _text
    lines = [ln for ln in _text(src).split("\n") if ln.strip() != ""]
    if not lines:
        raise ValueError("The mapping file is empty.")
    sep = _delim(lines) or "\t"
    rows = list(csv.reader(lines, delimiter=sep))
    out = {}
    for i, r in enumerate(rows):
        if len(r) < 2 or not r[0].strip() or not r[1].strip():
            continue
        a, b = r[0].strip(), r[1].strip()
        if i == 0 and a.lower() in ("gene", "gene_a", "a", "from", "gene1"):
            continue                                     # a header line
        out[a] = b
    if not out:
        raise ValueError("The mapping file needs two columns: a gene name in one library, then the other.")
    return out


def match_genes(genes_a, genes_b, mapping: dict | None = None) -> dict:
    """Which names in the two libraries are the same gene, and a diagnosis when they are not.

    Exact match after normalising, then the mapping file if one was given. The mapping is tried both ways
    round, because a user cannot be expected to remember which library is A, and whichever direction resolves
    more names is used and reported.
    """
    na = {normalise_name(g) for g in genes_a if str(g).strip()}
    nb = {normalise_name(g) for g in genes_b if str(g).strip()}
    direct = na & nb
    xlate: dict[str, str] = {}                           # a name in B -> the matching name in A
    direction = ""
    if mapping:
        fwd = {normalise_name(v): normalise_name(k) for k, v in mapping.items()}       # columns are A, B
        rev = {normalise_name(k): normalise_name(v) for k, v in mapping.items()}       # columns are B, A
        gain_f = len({b for b in nb - direct if fwd.get(b) in na})
        gain_r = len({b for b in nb - direct if rev.get(b) in na})
        xlate = fwd if gain_f >= gain_r else rev
        direction = "columns read as (A, B)" if gain_f >= gain_r else "columns read as (B, A)"
        xlate = {b: a for b, a in xlate.items() if b in nb and a in na and b not in direct}
    shared = sorted(direct | {xlate[b] for b in xlate})
    kind_a, kind_b = id_kind(genes_a), id_kind(genes_b)
    note = ""
    if not shared:
        note = (f"No gene name is in both libraries. A uses {kind_a} and B uses {kind_b}"
                + ("; supply a mapping file with the two names side by side."
                   if kind_a != kind_b else ", so the names should match — check for a prefix or suffix "
                   "added by one of the libraries."))
    elif len(shared) < 0.5 * min(len(na), len(nb)):
        note = (f"Only {len(shared):,} of {min(len(na), len(nb)):,} genes line up. A uses {kind_a} and B uses "
                f"{kind_b}." + (" A mapping file would recover the rest."
                                if kind_a != kind_b else " The unmatched names below may be aliases; a mapping "
                                                         "file would recover them."))
    return {"shared": shared, "n_a": len(na), "n_b": len(nb), "n_direct": len(direct),
            "n_mapped": len(xlate), "xlate": xlate, "kind_a": kind_a, "kind_b": kind_b,
            "unmatched_a": sorted(na - set(shared))[:25], "unmatched_b": sorted(nb - set(shared))[:25],
            "mapping_direction": direction, "note": note}


# ---------------------------------------------------------------- one screen
def _window(n: int, want: int) -> int:
    """drugZ's sliding window, narrowed for a small library so it is a window and not the whole table."""
    return int(max(25, min(want, max(25, n // 4))))


def _trend(effect: pd.Series, own_sd: pd.Series, window: int = 400) -> pd.Series:
    """The guide-to-guide spread a gene of this effect size usually shows, as a smooth trend.

    Fitted as a running median of the observed spread over genes ordered by their effect, then read back for
    every gene. Signed rather than absolute, because an enriched gene and a depleted one of the same size do
    not behave alike. A gene's own effect is used to look up its own prior, which is the same mild circularity
    every mean-variance trend carries (voom does it on the count mean) and is the price of not pretending that
    a harmless gene and a core-essential one are measured equally well.
    """
    ok = own_sd.notna() & effect.notna()
    if ok.sum() < 20:                                    # too few genes to see a trend; one number will do
        v = float(own_sd[ok].median()) if ok.any() else 0.0
        return pd.Series(v, index=effect.index)
    e = effect[ok].to_numpy()
    s = own_sd[ok].to_numpy()
    o = np.argsort(e, kind="mergesort")
    w = int(max(25, min(window, max(25, ok.sum() // 10))))
    run = pd.Series(s[o]).rolling(2 * w + 1, center=True, min_periods=max(5, w // 2)).median()
    run = run.ffill().bfill().to_numpy()                 # the window is short of genes at either end
    return pd.Series(np.interp(effect.to_numpy(), e[o], run), index=effect.index)


def _trigamma_inverse(v: np.ndarray) -> np.ndarray:
    """Solve trigamma(x) = v for x, by Newton on 1/x as limma's trigammaInverse does."""
    v = np.asarray(v, dtype=float)
    x = 0.5 + 1.0 / v
    for _ in range(50):
        tri = special.polygamma(1, x)
        step = tri * (1.0 - tri / v) / special.polygamma(2, x)
        x = x + step
        if np.all(np.abs(step) < 1e-10 * np.maximum(x, 1.0)):
            break
    return x


def _prior_df(own: pd.Series, trend: pd.Series, n: pd.Series, cap: float = 200.0) -> float:
    """How many degrees of freedom the guide-spread trend is worth (Smyth 2004, the moment estimator limma
    uses). The log of a gene's observed spread scatters about the trend by a known amount if the trend is the
    truth; any scatter beyond that is the trend being imperfect, and the shortfall gives d0 directly. A large
    d0 means the trend is trustworthy and a two-guide gene may lean on it; a small one means genes genuinely
    differ in how much their guides disagree, and each gene has to speak for itself."""
    ok = own.notna() & (own > 0) & trend.notna() & (trend > 0) & (n >= 2)
    if ok.sum() < 50:
        return 3.0                                    # too few genes to estimate it; a weak prior
    d = (n[ok] - 1.0).to_numpy(dtype=float)
    e = np.log(np.square(own[ok].to_numpy())) - np.log(np.square(trend[ok].to_numpy()))
    # e = [psi(d/2) - log(d/2)] - [psi(d0/2) - log(d0/2)] + noise, var(e) = trigamma(d/2) + trigamma(d0/2)
    excess = float(np.var(e, ddof=1) - np.mean(special.polygamma(1, d / 2.0)))
    if excess <= 0:
        return cap                                    # the trend explains the spreads as well as it could
    return float(min(cap, max(1.0, 2.0 * _trigamma_inverse(np.array([excess]))[0])))


def _side(counts, guides, genes, samples, comparison, ceg: set, neg: set, prior: float,
          half_window: int) -> dict:
    """One screen reduced to what the comparison needs: a scaled effect per gene, and its standard error."""
    counts = np.asarray(counts, dtype=float)
    lfc = guide_lfc(counts, samples, comparison)
    gn = np.array([normalise_name(x) for x in genes], dtype=object)
    eff = gene_effect(lfc, gn)
    f = scale_factors(eff, ceg, neg)
    s_eff = scaled(eff, f)
    # guide values on the same scale, so a spread measured on them is a spread in scaled units
    s_lfc = (lfc - f["median_neg"]) / f["span"] if f.get("ok") else lfc.astype(float)

    # the scale's own standard error: both medians are estimated from a finite reference set
    def _med_se(sel: set) -> float:
        v = eff[[g in sel for g in eff.index]].to_numpy(dtype=float)
        v = v[np.isfinite(v)]
        if len(v) < 3:
            return float("nan")
        sigma = float(np.median(np.abs(v - np.median(v)))) / 0.6744897501960817
        return MEDIAN_SE * sigma / np.sqrt(len(v))
    se_ceg, se_neg = _med_se(ceg), _med_se(neg)
    span = f.get("span", float("nan"))
    scale_se = {"ceg": se_ceg, "neg": se_neg,
                "span": float(np.hypot(se_ceg, se_neg)),
                "rel_span": float(np.hypot(se_ceg, se_neg) / span) if span and span > 0 else float("nan"),
                "offset": float(se_neg / span) if span and span > 0 else float("nan")}

    ctrl = [samples.index(s) for s in comparison["control"]]
    abundance = counts[:, ctrl].mean(axis=1)
    df = pd.DataFrame({"gene": gn, "guide": np.asarray(guides, dtype=object).astype(str),
                       "lfc": lfc, "scaled": s_lfc, "abundance": abundance})

    # the typical guide spread at a given abundance: residuals from each gene's own median, smoothed in a
    # sliding window of guides ordered from most to least abundant (drugZ's empirical-Bayes step)
    resid = (df["scaled"] - df.groupby("gene")["scaled"].transform("median")).to_numpy()
    order = np.argsort(-abundance, kind="mergesort")
    hw = _window(len(resid), half_window)
    sd_g = np.empty(len(resid))
    sd_g[order] = _empirical_bayes(resid[order], hw)
    df["sd_global2"] = np.square(sd_g)

    agg = df.groupby("gene", sort=True)
    n = agg["scaled"].count()
    own = agg["scaled"].std(ddof=1)                      # NaN for a single guide
    s_gene = s_eff.reindex(n.index)
    trend = _trend(s_gene, own)                          # the spread a gene this strong usually shows
    floor = np.sqrt(agg["sd_global2"].mean())            # and the spread its guides' abundance implies
    prior_sd = np.maximum(trend, floor)
    d0 = _prior_df(own, prior_sd, n) if prior <= 0 else float(prior)
    dfree = (n - 1).clip(lower=0).astype(float)
    w = dfree / (dfree + d0)
    sd = np.sqrt(w * own.fillna(0.0).pow(2) + (1.0 - w) * prior_sd.pow(2))
    se = MEDIAN_SE * sd / np.sqrt(n)

    genes_out = pd.DataFrame({"n": n, "effect": eff.reindex(n.index), "scaled": s_gene,
                              "sd_guide": sd, "sd_own": own, "sd_prior": prior_sd, "se": se})
    return {"genes": genes_out, "guide": df.drop(columns=["sd_global2"]), "scale": f,
            "name": comparison["name"], "window": hw, "prior_df": d0, "scale_se": scale_se,
            "guides_per_gene": float(n.median()), "n_guides": int(len(df))}


# ---------------------------------------------------------------- the comparison
def cross_compare(a: dict, b: dict, thresholds: dict | None = None, mapping: dict | None = None,
                  species: str = "auto", half_window: int = 500) -> dict:
    """Compare two screens held in separate count tables, with no guide in common assumed.

    `a` and `b` each carry counts, guides, genes, samples, a comparison (treatment/control sample names) and
    a project name. Both are measured against their own reference samples, then placed on the common
    essential-gene scale, which is the only thing the two libraries share.
    """
    th = {**DEFAULTS, **(thresholds or {})}
    if species == "auto":
        from .enrichment import guess_species
        species = guess_species(list(a["genes"]) + list(b["genes"]))
    try:
        ceg_raw, neg_raw = qc.load_reference_sets(species)
    except ValueError:
        ceg_raw, neg_raw = set(), set()
    # the reference sets are matched in the same normalised space as the libraries (mouse symbols are mixed
    # case, so this has to happen after the species is settled, never before)
    ceg = {normalise_name(g) for g in ceg_raw}
    neg = {normalise_name(g) for g in neg_raw}

    prior, hw = float(th.get("xc_prior_df", 0.0)), int(half_window)    # 0 = estimate it from the data
    A = _side(a["counts"], a["guides"], a["genes"], a["samples"], a["comparison"], ceg, neg, prior, hw)
    B = _side(b["counts"], b["guides"], b["genes"], b["samples"], b["comparison"], ceg, neg, prior, hw)

    # Without the common scale there is no link at all between two libraries, so this is refused outright
    # rather than reported as a caveat under a table of numbers.
    thin = [(s["name"], s["scale"]) for s in (A, B)
            if s["scale"].get("ok") and min(s["scale"]["n_ceg"], s["scale"]["n_neg"]) < MIN_REFERENCE]
    if thin:
        who = " and ".join(n for n, _ in thin)
        f = thin[0][1]
        raise ValueError(
            f"{who} carries only {f['n_ceg']} core-essential and {f['n_neg']} non-essential reference genes, and "
            f"at least {MIN_REFERENCE} of each are needed here. Two different libraries have nothing in common "
            "except the essential-gene scale, so a scale estimated from a handful of genes is not good enough — "
            "its own error would be larger than any difference the test could resolve. A focused or sub-genome "
            "library usually fails this; compare it against a screen of the same library instead, on the Compare "
            "page.")
    bad = [(s["name"], s["scale"]) for s in (A, B) if not s["scale"].get("ok")]
    if bad:
        who = " and ".join(n for n, _ in bad)
        f = bad[0][1]
        raise ValueError(
            f"{who} cannot be placed on the common essential-gene scale, and that scale is the only thing two "
            f"different libraries share — so there is nothing to compare. Found {f['n_ceg']} core-essential and "
            f"{f['n_neg']} non-essential reference genes"
            + (" (at least 10 of each are needed)." if min(f['n_ceg'], f['n_neg']) < 10 else
               "; the screen does not separate them (the median essential gene is not below the median "
               "non-essential one), which an early timepoint or a failed screen looks like."))

    m = match_genes(a["genes"], b["genes"], mapping)
    ga, gb = A["genes"], B["genes"].rename(index=lambda g: m["xlate"].get(g, g))
    gb = gb[~gb.index.duplicated(keep="first")]
    keep = [g for g in m["shared"] if g in ga.index and g in gb.index]
    ga, gb = ga.loc[keep], gb.loc[keep]
    enough = (ga["n"] >= int(th["xc_min_guides"])) & (gb["n"] >= int(th["xc_min_guides"]))
    ga, gb = ga[enough], gb[enough]
    if len(ga) < 20:
        raise ValueError(
            f"Only {len(ga)} gene(s) are in both libraries with at least {int(th['xc_min_guides'])} guides each"
            + (f". {m['note']}" if m["note"] else ". Lower the guides-per-gene threshold, or check that the two "
                                                 "libraries name genes the same way."))

    delta = (gb["scaled"] - ga["scaled"]).to_numpy()
    # each gene's own measurement error, plus the error each screen's scale contributes: a constant from the
    # non-essential offset, and a term proportional to the gene's depth from the span
    def _var(side, g):
        e = side["scale_se"]
        off = 0.0 if not np.isfinite(e["offset"]) else e["offset"]
        rel = 0.0 if not np.isfinite(e["rel_span"]) else e["rel_span"]
        return np.square(g["se"].to_numpy()) + off ** 2 + np.square(g["scaled"].to_numpy() * rel)
    se = np.sqrt(_var(A, ga) + _var(B, gb))
    t_raw = np.divide(delta, se, out=np.zeros_like(delta), where=se > 0)

    # Efron's empirical null: centre on the bulk, widen by its MAD — and only ever widen, so the calibration
    # cannot make the test more liberal than the nominal model. Confirmed 2026-11-24: the floor at 1 stays even
    # where the data say the test is over-cautious (on two tables of one HAP1 screen the statistic scatters at
    # 0.09 where 1 is expected). That is reported in words instead; narrowing is how a tool invents hits.
    centre = float(np.median(t_raw))
    mad = float(np.median(np.abs(t_raw - centre)))
    spread = max(1.0, mad / 0.6744897501960817) if mad > 0 else 1.0
    t = (t_raw - centre) / spread
    d0 = min(A["prior_df"], B["prior_df"])               # the weaker of the two priors, to be conservative
    dfree = (ga["n"].to_numpy() - 1) + (gb["n"].to_numpy() - 1) + d0
    p = 2.0 * stats.t.sf(np.abs(t), dfree)
    fdr = _bh(p)

    out = pd.DataFrame({
        "gene": list(ga.index),
        "guides_a": ga["n"].to_numpy(), "guides_b": gb["n"].to_numpy(),
        "effect_a": ga["effect"].to_numpy(), "effect_b": gb["effect"].to_numpy(),
        "scaled_a": ga["scaled"].to_numpy(), "scaled_b": gb["scaled"].to_numpy(),
        "se_a": ga["se"].to_numpy(), "se_b": gb["se"].to_numpy(),
        "delta": delta, "se_delta": se, "t": t, "df": dfree, "p": p, "fdr": fdr,
    })
    out["normZ"] = out["t"]                              # the column the shared UI and exports sort on
    dep = float(th["xc_dependency"])
    da, db = out["scaled_a"] <= dep, out["scaled_b"] <= dep
    out["different"] = (out["fdr"] <= float(th["xc_fdr"])) & (out["delta"].abs() >= float(th["xc_min_delta"]))
    out["call"] = np.select([da & db, da & ~db, ~da & db], ["shared", "specific to A", "specific to B"],
                            default="neither")
    out["stronger_in"] = np.where(~out["different"], "", np.where(out["delta"] < 0, "B", "A"))
    out = out.sort_values("t").reset_index(drop=True)

    # How far off nominal the null was, and what this pair of screens can actually resolve. The bar a gene
    # has to clear is not a nominal 5%: it is 5% after correcting for every gene tested, so the resolution is
    # quoted against that bar. Quoting the nominal number would flatter the method by a factor of two.
    neg_t = out["t"][[g in neg for g in out["gene"]]]
    med_df = float(np.median(dfree))
    bar = float(stats.t.isf(float(th["xc_fdr"]) / (2.0 * max(1, len(out))), med_df))
    power = {
        "median_se": float(np.median(se)),
        "resolution": float(np.median(se) * bar * spread),
        "min_detectable": float(np.median(se) * POWER_Z * spread),
        "bar": bar, "median_df": med_df, "prior_df": float(d0),
        "scale_rel_a": float(A["scale_se"]["rel_span"]), "scale_rel_b": float(B["scale_se"]["rel_span"]),
        "prior_df_a": float(A["prior_df"]), "prior_df_b": float(B["prior_df"]),
        "calibration": float(spread),
        "centre": centre,
        "null_sd_neg": float(neg_t.std(ddof=1)) if len(neg_t) > 2 else float("nan"),
        "n_neg": int(len(neg_t)),
    }
    shared = shared_guides(a["guides"], b["guides"])
    return {
        "genes": out,
        "a": {"name": A["name"], "project": a.get("project", ""), "scale": A["scale"],
              "guides_per_gene": A["guides_per_gene"], "n_guides": A["n_guides"], "window": A["window"],
              "prior_df": A["prior_df"]},
        "b": {"name": B["name"], "project": b.get("project", ""), "scale": B["scale"],
              "guides_per_gene": B["guides_per_gene"], "n_guides": B["n_guides"], "window": B["window"],
              "prior_df": B["prior_df"]},
        "guide_a": A["guide"], "guide_b": B["guide"],
        "species": species, "thresholds": th, "matching": m, "power": power,
        "shared_guides": shared,
        "concordance": _concordance(out, {"cmp_dependency": dep}),
    }


def library_notes(r: dict) -> list[dict]:
    """What the user has to know about this pair before reading the table, as pass/warn/fail lines."""
    out = []
    m, p, th = r["matching"], r["power"], r["thresholds"]
    n_shared = int(len(r["genes"]))
    frac = n_shared / max(1, min(m["n_a"], m["n_b"]))
    out.append({"id": "overlap", "label": "Genes in both libraries",
                "value": f"{n_shared:,}",
                "status": "pass" if frac >= 0.6 else "warn" if frac >= 0.25 else "fail",
                "why": (f"{n_shared:,} of {min(m['n_a'], m['n_b']):,} genes, each with at least "
                        f"{int(th['xc_min_guides'])} guides in both. " + m["note"]).strip()})
    gmin = min(r["a"]["guides_per_gene"], r["b"]["guides_per_gene"])
    out.append({"id": "guides", "label": "Guides per gene (the thinner library)",
                "value": f"{gmin:.0f}",
                "status": "pass" if gmin >= 4 else "warn" if gmin >= 3 else "fail",
                "why": f"{r['a']['guides_per_gene']:.0f} in A and {r['b']['guides_per_gene']:.0f} in B. "
                       "Fewer guides means a larger standard error, and a larger difference needed to see it."})
    out.append({"id": "resolution", "label": "Smallest difference this pair can resolve",
                "value": f"{p['resolution']:.2f}",
                "status": "pass" if p["resolution"] <= 0.75 else "warn" if p["resolution"] <= 1.5 else "fail",
                "why": (f"A gene with the typical guide spread needs a difference of {p['resolution']:.2f} to clear "
                        f"the FDR bar across all {n_shared:,} genes tested (t {p['bar']:.2f} on {p['median_df']:.0f} "
                        "degrees of freedom). On this scale 1.0 is the distance from a harmless gene halfway down "
                        "to a core-essential one, so anything smaller than that will not be called here whatever "
                        "the biology — read a null result as 'no difference this large', not 'no difference'. "
                        "Guides per gene is the lever: the figure falls roughly as its square root.")})
    out.append({"id": "null", "label": "Null width", "value": f"{p['calibration']:.2f}",
                "status": "pass" if p["calibration"] <= 1.5 else "warn" if p["calibration"] <= 2.5 else "fail",
                "why": (f"The spread of t had to be widened {p['calibration']:.2f}× to behave like a null "
                        "distribution; 1 means the guide-level noise already explains how far apart the two "
                        "screens are, and a large number means something systematic separates them — different "
                        "timepoints, cell lines or library chemistry."
                        + (f" Among {p['n_neg']} non-essential reference genes the spread is "
                           f"{p['null_sd_neg']:.2f}, where 1 is right." if p["n_neg"] > 2 else ""))})
    sg = r["shared_guides"]
    rel = max(p["scale_rel_a"], p["scale_rel_b"])
    out.append({"id": "scale", "label": "Precision of the common scale",
                "value": f"\u00b1{rel * 100:.1f}%",
                "status": "pass" if rel <= 0.05 else "warn" if rel <= 0.12 else "fail",
                "why": f"The scale rests on {r['a']['scale']['n_ceg']} core-essential and "
                       f"{r['a']['scale']['n_neg']} non-essential reference genes in A, and "
                       f"{r['b']['scale']['n_ceg']} and {r['b']['scale']['n_neg']} in B. Its error is carried "
                       "into every gene in proportion to how deep that gene is, so the deepest genes are the "
                       "least certain, not the most."})
    if np.isfinite(p["null_sd_neg"]) and p["null_sd_neg"] < 0.5:
        out.append({"id": "over_agree", "label": "The screens agree more than their noise allows",
                    "value": f"{p['null_sd_neg']:.2f}",
                    "status": "warn",
                    "why": "Among non-essential genes the test statistic scatters far less than the 1 it should, "
                           "so the two screens are not independent measurements — the same data on both sides, "
                           "or shared samples, will do this. Nothing is wrong with the result, but the "
                           "resolution quoted above is pessimistic for this pair, and a difference it does call "
                           "is on firmer ground than the numbers suggest."})
    if sg["fraction"] > 0.5:
        out.append({"id": "same_library", "label": "The two screens share their guides",
                    "value": f"{sg['n']:,} ({sg['fraction'] * 100:.0f}%)",
                    "status": "warn",
                    "why": ("Matched on the protospacer sequence inside each guide ID, so a difference in ID "
                            "style does not hide it. " if sg["by_sequence"] else "")
                           + "These two screens are mostly the same library, so this is a cross-project "
                             "comparison but not a cross-library one. Put both screens in one count table and "
                             "use the Compare page instead: pairing each guide with itself removes guide "
                             "efficiency, the largest source of error here, and resolves differences several "
                             "times smaller than this page can."})
    return out
