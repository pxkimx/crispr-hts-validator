"""PDF report of a screen: data, every quality check with its verdict and reason, results, methods and settings.
Print-friendly (white background) whatever the app's theme. Figures are drawn with matplotlib, pages with reportlab."""
from __future__ import annotations

import io
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from reportlab.lib import colors  # noqa: E402
from reportlab.lib.pagesizes import A4  # noqa: E402
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet  # noqa: E402
from reportlab.lib.units import mm  # noqa: E402
from reportlab.platypus import (Image, KeepTogether, PageBreak, Paragraph, SimpleDocTemplate, Spacer,  # noqa: E402
                                Table, TableStyle)

from . import __version__  # noqa: E402

INK, MUTED, ACCENT = colors.HexColor("#0b0b0b"), colors.HexColor("#5f5e5a"), colors.HexColor("#0f7b74")
STATUS = {"pass": ("#0ca30c", "PASS"), "warn": ("#b77d00", "WARN"), "fail": ("#d03b3b", "FAIL"), "na": ("#898781", "N/A")}
NEG, POS, NS = "#2a78d6", "#e34948", "#c3c2b7"
W = 174 * mm


def _styles():
    ss = getSampleStyleSheet()
    base = ParagraphStyle("b", parent=ss["Normal"], fontName="Helvetica", fontSize=9, leading=12.5, textColor=INK)
    return {
        "h1": ParagraphStyle("h1", parent=base, fontName="Times-Bold", fontSize=22, leading=26, spaceAfter=4),
        "h2": ParagraphStyle("h2", parent=base, fontName="Times-Bold", fontSize=15, leading=19, spaceBefore=10, spaceAfter=5),
        "h3": ParagraphStyle("h3", parent=base, fontName="Helvetica-Bold", fontSize=10, leading=13, spaceBefore=6, spaceAfter=3),
        "kick": ParagraphStyle("k", parent=base, fontName="Courier-Bold", fontSize=8, textColor=ACCENT, spaceAfter=2),
        "p": base,
        "small": ParagraphStyle("s", parent=base, fontSize=7.5, leading=10, textColor=MUTED),
        "cell": ParagraphStyle("c", parent=base, fontSize=7.5, leading=9.5),
    }


def _fig(fig, width=W, dpi=200):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    buf.seek(0)
    w, h = fig.get_size_inches()
    return Image(buf, width=width, height=width * h / w)


def _ax_style(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.spines["left"].set_color("#c3c2b7")
    ax.spines["bottom"].set_color("#c3c2b7")
    ax.tick_params(colors="#52514e", labelsize=7)
    ax.grid(axis="y", color="#e1e0d9", linewidth=0.6)
    ax.set_axisbelow(True)


def _table(rows, widths, st, header=True, zebra=True):
    data = [[Paragraph(str(c), st["cell"]) if not isinstance(c, Paragraph) else c for c in r] for r in rows]
    t = Table(data, colWidths=widths, repeatRows=1 if header else 0)
    style = [("VALIGN", (0, 0), (-1, -1), "TOP"), ("LINEBELOW", (0, 0), (-1, 0), 0.6, colors.HexColor("#c3c2b7")),
             ("TOPPADDING", (0, 0), (-1, -1), 2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2)]
    if header:
        style.append(("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"))
    if zebra:
        for i in range(2, len(rows), 2):
            style.append(("BACKGROUND", (0, i), (-1, i), colors.HexColor("#f6f6f3")))
    t.setStyle(TableStyle(style))
    return t


def _status(s, st):
    col, lab = STATUS.get(s, STATUS["na"])
    return Paragraph(f'<font color="{col}"><b>{lab}</b></font>', st["cell"])


def _fmt(v):
    if v is None:
        return "–"
    if isinstance(v, (int, np.integer)):
        return f"{int(v):,}"
    try:
        v = float(v)
    except (TypeError, ValueError):
        return str(v)
    if not np.isfinite(v):
        return "–"
    a = abs(v)
    if a != 0 and (a < 1e-3 or a >= 1e6):
        return f"{v:.2e}"
    return f"{v:.3g}"


def _footer(name):
    def draw(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 7)
        canvas.setFillColor(MUTED)
        canvas.drawString(18 * mm, 10 * mm, f"CRISPR HTS Validator {__version__} · {name}")
        canvas.drawRightString(192 * mm, 10 * mm, f"page {doc.page}")
        canvas.restoreState()
    return draw


def _qc_section(story, qcres, st, title="Quality checks"):
    story.append(Paragraph(title, st["h2"]))
    if not qcres:
        story.append(Paragraph("The quality checks were not run.", st["p"]))
        return
    rows = [["", "Check", "Why"]]
    for c in qcres["checks"]:
        rows.append([_status(c["status"], st), Paragraph(f"<b>{c['name']}</b>", st["cell"]),
                     Paragraph(c["why"] + (f'<br/><font color="#898781">{c.get("rule", "")}</font>' if c.get("rule") else ""), st["cell"])])
    story.append(_table(rows, [14 * mm, 38 * mm, W - 52 * mm], st))


# ----------------------------------------------------------------------------- pooled
def pooled_report(P) -> bytes:
    from . import pooled_methods
    st = _styles()
    buf = io.BytesIO()
    name = P.name or "Pooled screen"
    doc = SimpleDocTemplate(buf, pagesize=A4, leftMargin=18 * mm, rightMargin=18 * mm, topMargin=16 * mm,
                            bottomMargin=16 * mm, title=f"{name} — CRISPR HTS Validator report", author="CRISPR HTS Validator")
    s = []
    t = P.table
    s.append(Paragraph("POOLED CRISPR SCREEN · QUALITY AND HITS", st["kick"]))
    s.append(Paragraph(name, st["h1"]))
    s.append(Paragraph(f"Report made {time.strftime('%Y-%m-%d %H:%M')} with CRISPR HTS Validator {__version__}. "
                       "Everything was computed on this computer.", st["small"]))
    s.append(Spacer(1, 6))
    s.append(Paragraph("Data", st["h2"]))
    rows = [["Input", "Guides", "Genes", "Samples"],
            [f"{t.fmt} ({t.source})", f"{len(t.guides):,}", f"{len(set(t.genes)):,}", ", ".join(t.samples)]]
    s.append(_table(rows, [50 * mm, 22 * mm, 22 * mm, W - 94 * mm], st))
    if P.count_info:
        rows = [["Sample", "Reads", "Mapped", "%", "Zero-count guides"]] + [
            [c["sample"], f"{c['reads']:,}", f"{c['mapped']:,}", f"{100 * c['fraction']:.1f}", f"{c['zero']:,}"] for c in P.count_info]
        s.append(Spacer(1, 4))
        s.append(_table(rows, [40 * mm, 30 * mm, 30 * mm, 20 * mm, 30 * mm], st))
    refs = [r["sample"] for r in P.design["samples"] if r["reference"]]
    s.append(Paragraph("Samples and comparisons", st["h3"]))
    s.append(Paragraph("Reference (plasmid / T0): " + (", ".join(refs) or "none marked") + "<br/>" +
                       "<br/>".join(f"<b>{c['name']}</b>: {', '.join(c['treatment'])} vs {', '.join(c['control'])}"
                                    for c in P.design["comparisons"]), st["p"]))
    _qc_section(s, P.qc, st)
    if P.qc:
        figs = _pooled_qc_figs(P.qc)
        for f in figs:
            s.append(Spacer(1, 4))
            s.append(f)
    for nm, entry in P.results.items():
        s.append(PageBreak())
        _pooled_result(s, nm, entry, P, st)
    s.append(PageBreak())
    s.append(Paragraph("Methods", st["h2"]))
    for txt in METHODS_POOLED:
        s.append(Paragraph(txt, st["p"]))
        s.append(Spacer(1, 3))
    doc.build(s, onFirstPage=_footer(name), onLaterPages=_footer(name))
    return buf.getvalue()


def _pooled_qc_figs(q):
    out = []
    by = {c["id"]: c for c in q["checks"]}
    c = by.get("reads")
    e = by.get("essentials")
    n = 2 if (e and e.get("per_condition")) else 1
    fig, axes = plt.subplots(1, n, figsize=(7.2, 2.6))
    axes = np.atleast_1d(axes)
    ax = axes[0]
    ps = c["per_sample"]
    ax.barh([p["sample"] for p in ps][::-1], [p["coverage"] for p in ps][::-1],
            color=[STATUS[p["status"]][0] for p in ps][::-1], height=0.6)
    ax.axvline(q["thresholds"]["coverage_pass"], color="#52514e", ls="--", lw=0.8)
    ax.set_xlabel("reads per guide", fontsize=7.5)
    _ax_style(ax)
    ax.grid(axis="x", color="#e1e0d9", linewidth=0.6)
    ax.grid(axis="y", visible=False)
    if n == 2:
        ax = axes[1]
        pal = ["#0f7b74", "#eb6834", "#2a78d6", "#e87ba4", "#008300", "#6250d6"]
        for i, p in enumerate(e["per_condition"]):
            ax.plot(p["curves"]["fpr"], p["curves"]["tpr"], color=pal[i % 6], lw=1.6, label=f"{p['condition']} (AUC {p['value']:.3f})")
        ax.plot([0, 1], [0, 1], color="#c3c2b7", ls="--", lw=0.8)
        ax.set_xlabel("false positive rate (non-essential)", fontsize=7.5)
        ax.set_ylabel("true positive rate (essential)", fontsize=7.5)
        ax.legend(fontsize=6.5, frameon=False, loc="lower right")
        _ax_style(ax)
    fig.tight_layout()
    out.append(_fig(fig))
    return out


def _pooled_result(s, nm, entry, P, st):
    r = entry["result"]
    m = r.get("method", "mageck")
    c = entry["comparison"]
    s.append(Paragraph(f"RESULT · {m.upper()}", st["kick"]))
    s.append(Paragraph(nm, st["h2"]))
    s.append(Paragraph(f"Treatment: {', '.join(c['treatment'])}. Control: {', '.join(c['control'])}. Run {entry['when']}. "
                       f"Settings: " + "; ".join(f"{k} = {v}" for k, v in entry["settings"].items() if v not in (None, [], "")), st["small"]))
    cn = r.get("cn")
    if cn:
        s.append(Paragraph(f"Copy-number correction: {cn['n_positioned']:,} guides with positions, "
                           f"{len(cn['segments'])} segments; {cn['corrected']:,} guides shifted by more than 0.25 log2 units.", st["small"]))
    if m == "mageck":
        rows = r["genes"]
        lfc = np.array([g["lfc"] for g in rows])
        p = np.array([g["neg_p"] if g["lfc"] < 0 else g["pos_p"] for g in rows])
        q = np.array([g["neg_fdr"] if g["lfc"] < 0 else g["pos_fdr"] for g in rows])
        fig, ax = plt.subplots(figsize=(7.2, 3.0))
        col = np.where(q <= 0.05, np.where(lfc < 0, NEG, POS), NS)
        o = np.argsort(q <= 0.05)
        ax.scatter(lfc[o], -np.log10(np.maximum(p[o], 1e-300)), s=4, c=col[o], linewidths=0)
        ax.axvline(0, color="#c3c2b7", lw=0.8)
        ax.set_xlabel("gene log2 fold change", fontsize=7.5)
        ax.set_ylabel("−log10 p (RRA)", fontsize=7.5)
        _ax_style(ax)
        s.append(_fig(fig))
        neg = sum(1 for g in rows if g["neg_fdr"] < 0.05)
        pos = sum(1 for g in rows if g["pos_fdr"] < 0.05)
        s.append(Paragraph(f"{neg:,} genes depleted and {pos:,} enriched at FDR &lt; 0.05 (of {len(rows):,}). "
                           f"α = {r['alpha_neg']:.3f} / {r['alpha_pos']:.3f}.", st["p"]))
        for key, title in (("neg", "Top depleted genes"), ("pos", "Top enriched genes")):
            top = sorted(rows, key=lambda g: g[f"{key}_rank"])[:20]
            tr = [["Rank", "Gene", "Guides", "LFC", "RRA score", "p", "FDR"]] + [
                [g[f"{key}_rank"], g["gene"], g["guides"], _fmt(g["lfc"]), _fmt(g[f"{key}_score"]), _fmt(g[f"{key}_p"]), _fmt(g[f"{key}_fdr"])] for g in top]
            s.append(KeepTogether([Paragraph(title, st["h3"]), _table(tr, [14 * mm, 34 * mm, 16 * mm, 22 * mm, 28 * mm, 28 * mm, 28 * mm], st)]))
    elif m == "bagel":
        pr = r["pr"]
        bf = np.sort(r["bf"])[::-1]
        fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8))
        ax = axes[0]
        ax.plot(np.arange(1, len(bf) + 1), bf, color="#0f7b74", lw=1.4)
        ax.axhline(5, color="#52514e", ls="--", lw=0.8)
        ax.set_xlabel("gene rank", fontsize=7.5)
        ax.set_ylabel("Bayes factor", fontsize=7.5)
        _ax_style(ax)
        ax = axes[1]
        ax.plot(pr["recall"], pr["precision"], color="#eb6834", lw=1.4)
        ax.set_xlabel("recall (CEG)", fontsize=7.5)
        ax.set_ylabel("precision", fontsize=7.5)
        ax.set_ylim(0, 1.02)
        _ax_style(ax)
        fig.tight_layout()
        s.append(_fig(fig))
        sm = r["summary"]
        s.append(Paragraph(f"{sm['calls']:,} genes with BF &gt; {sm['cutoff']:g}; at this cut-off recall of reference essential "
                           f"genes is {100 * sm['recall']:.1f}% and precision {100 * (sm['precision'] or 0):.1f}%.", st["p"]))
        order = np.argsort(-r["bf"])[:30]
        tr = [["Rank", "Gene", "BF", "Reference set"]] + [
            [i + 1, r["genes"][j], _fmt(r["bf"][j]), "essential" if r["genes"][j] in r["ref"]["ceg"] else ("non-essential" if r["genes"][j] in r["ref"]["neg"] else "")]
            for i, j in enumerate(order)]
        s.append(KeepTogether([Paragraph("Top genes by Bayes factor", st["h3"]), _table(tr, [16 * mm, 44 * mm, 30 * mm, 40 * mm], st)]))
    elif m == "drugz":
        g = r["genes"]
        z = g["normZ"].to_numpy()
        fig, ax = plt.subplots(figsize=(7.2, 2.6))
        col = np.where(g["fdr_synth"] < 0.05, NEG, np.where(g["fdr_supp"] < 0.05, POS, NS))
        ax.scatter(np.arange(1, len(z) + 1), z, s=4, c=col, linewidths=0)
        ax.axhline(0, color="#c3c2b7", lw=0.8)
        ax.set_xlabel("gene rank", fontsize=7.5)
        ax.set_ylabel("drugZ normZ", fontsize=7.5)
        _ax_style(ax)
        s.append(_fig(fig))
        s.append(Paragraph(f"{int((g['fdr_synth'] < 0.05).sum()):,} synthetic sick/lethal and {int((g['fdr_supp'] < 0.05).sum()):,} "
                           f"suppressor genes at FDR &lt; 0.05 ({'paired' if r['paired'] else 'unpaired'} analysis).", st["p"]))
        for asc, title, pk, fk in ((True, "Synthetic sick / lethal (most negative normZ)", "pval_synth", "fdr_synth"),
                                   (False, "Suppressors (most positive normZ)", "pval_supp", "fdr_supp")):
            top = g.sort_values("normZ", ascending=asc).head(20)
            tr = [["Gene", "normZ", "observations", "p", "FDR"]] + [[i, _fmt(x.normZ), int(x.numObs), _fmt(x[pk]), _fmt(x[fk])] for i, x in top.iterrows()]
            s.append(KeepTogether([Paragraph(title, st["h3"]), _table(tr, [44 * mm, 26 * mm, 26 * mm, 30 * mm, 30 * mm], st)]))


METHODS_POOLED = [
    "<b>Normalisation and guide test (MAGeCK-style).</b> Median-ratio size factors; control and treatment means are "
    "medians of normalised replicates; variance from a mean–variance model fitted on control replicates; guide "
    "p-values from a normal distribution truncated at zero; α-RRA gene ranking with permutation p-values using "
    "MAGeCK 0.5.9.5's random-number generator and seed; Benjamini–Hochberg FDR (Li et al., Genome Biology 2014). "
    "Validated against mageck test 0.5.9.5 to printing precision.",
    "<b>BAGEL2 Bayes factors.</b> Counts + 5, scaled to 10 million reads; log2 fold change against the summed "
    "controls; 10-fold cross-validated kernel-density log-likelihood ratio between reference essential (CEGv2) and "
    "non-essential (NEGv1) genes; a gene's BF is the sum over its guides and replicates (Kim &amp; Hart, Genome "
    "Medicine 2021). A gene is called essential at BF &gt; 5, the value BAGEL2's tutorial uses, confirmed for this "
    "tool on 7 November 2026.",
    "<b>drugZ.</b> Paired replicate fold changes of treated vs untreated, guide z-scores from an empirical-Bayes "
    "standard deviation in a sliding window of 1,000 guides, gene normZ, normal p-values and FDR in both directions "
    "(Colic et al., Genome Medicine 2019).",
    "<b>Copy-number correction.</b> CRISPRcleanR-style (Iorio et al., BMC Genomics 2018): guide fold changes ordered "
    "along each chromosome are segmented (binary segmentation, permutation test); segments covering at least three "
    "genes have their mean removed.",
    "<b>Quality checks.</b> Each check reports pass, warn or fail against the thresholds shown in the table, which "
    "follow published guidance (MAGeCK QC; Hart et al. 2015, 2017; Joung et al. 2017) and were confirmed for this "
    "tool on 7 November 2026. Any threshold can be changed in the app, and a changed one is marked in this table.",
    "<b>Library representation and MOI.</b> Cells per guide is the number of cells entered for each step divided by "
    "the guides in the library; at least 500\u00d7 is expected at every step, because below that guides are lost by "
    "sampling rather than by selection. The fraction of transduced cells carrying a single guide is the Poisson "
    "p(1)/p(\u22651) at the MOI entered (Joung et al., Nat Protoc 2017; Doench, Nat Rev Genet 2018).",
    "<b>Guides per gene.</b> Counted over targeting guides only; genes with fewer than four independent guides are "
    "listed, since a gene-level score built on one or two guides cannot separate a real effect from a bad guide.",
    "<b>Non-targeting controls.</b> Control guides are found by name (non-targeting, NTC, safe-harbour, EGFP, LacZ, "
    "luciferase, AAVS1, or a pattern you set). Their median log2 fold change should be zero; their robust standard "
    "deviation (IQR/1.349) gives the noise floor quoted as the smallest distinguishable effect.",
]


# ----------------------------------------------------------------------------- arrayed
def arrayed_report(A) -> bytes:
    st = _styles()
    buf = io.BytesIO()
    name = A.name or "Arrayed screen"
    doc = SimpleDocTemplate(buf, pagesize=A4, leftMargin=18 * mm, rightMargin=18 * mm, topMargin=16 * mm,
                            bottomMargin=16 * mm, title=f"{name} — CRISPR HTS Validator report", author="CRISPR HTS Validator")
    s = []
    sm = A.screen.summary()
    s.append(Paragraph("ARRAYED SCREEN · PLATE QUALITY AND HITS", st["kick"]))
    s.append(Paragraph(name, st["h1"]))
    s.append(Paragraph(f"Report made {time.strftime('%Y-%m-%d %H:%M')} with CRISPR HTS Validator {__version__}.", st["small"]))
    s.append(Paragraph("Data", st["h2"]))
    rows = [["Plates", "Format", "Replicates", "Channels", "Wells by type", "Flagged"],
            [sm["plates"], f"{sm['wells']}-well", len(sm["replicates"]), len(sm["channels"]),
             ", ".join(f"{k} {v:,}" for k, v in sm["contents"].items()), sm["flagged"]]]
    s.append(_table(rows, [16 * mm, 20 * mm, 20 * mm, 18 * mm, W - 94 * mm, 20 * mm], st))
    if sm["notes"]:
        s.append(Paragraph("<br/>".join(sm["notes"]), st["small"]))
    d = sm.get("description") or {}
    if d:
        s.append(Paragraph(" · ".join(f"<b>{k}</b>: {d[k]}" for k in ("Title", "Organism", "Celltype", "Assay", "Reference") if d.get(k)), st["small"]))
    _qc_section(s, A.qc, st, "Plate checks")
    if A.qc:
        try:
            q = A.qc_view()
            fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.4))
            for ax, prof, lab in ((axes[0], q["row_profile"], "plate row"), (axes[1], q["col_profile"], "plate column")):
                v = [np.nan if x is None else x for x in prof]
                ax.plot(np.arange(1, len(v) + 1), v, color="#0f7b74", marker="o", ms=3, lw=1.4)
                ax.axhline(1, color="#c3c2b7", ls="--", lw=0.8)
                ax.set_xlabel(lab, fontsize=7.5)
                ax.set_ylabel("÷ plate median", fontsize=7.5)
                _ax_style(ax)
            fig.tight_layout()
            s.append(Spacer(1, 4))
            s.append(_fig(fig))
        except Exception:  # noqa: BLE001
            pass
    if A.result:
        r = A.result
        W_ = r["wells"]
        st_ = r["settings"]
        s.append(PageBreak())
        s.append(Paragraph("Hits", st["h2"]))
        s.append(Paragraph(f"Normalisation: {st_['method']}{' (log2)' if st_['log2'] else ''}; score: {st_['score']}; "
                           f"{st_['direction']} signal is a hit; combined score ({st_['summary']} of replicates) ≥ {st_['threshold']} "
                           f"and every replicate ≥ {st_['replicate_min']}. Run {r['when']}.", st["small"]))
        nh = int(W_["hit"].sum())
        ns = int((W_["content"] == "sample").sum())
        g = r["genes"]
        s.append(Paragraph(f"<b>{nh:,}</b> hit wells of {ns:,} sample wells ({100 * nh / max(1, ns):.1f}%)"
                           + (f"; <b>{int(g['hit'].sum()):,}</b> hit genes of {len(g):,}." if len(g) else "."), st["p"]))
        reps = r["replicates"]
        if len(reps) >= 2:
            samp = W_[W_["content"] == "sample"]
            fig, ax = plt.subplots(figsize=(4.2, 3.4))
            hit = samp["hit"].to_numpy()
            ax.scatter(samp[f"score_{reps[0]}"][~hit], samp[f"score_{reps[1]}"][~hit], s=2, c=NS, linewidths=0)
            ax.scatter(samp[f"score_{reps[0]}"][hit], samp[f"score_{reps[1]}"][hit], s=4, c=POS, linewidths=0)
            ax.set_xlabel(f"score, replicate {reps[0]}", fontsize=7.5)
            ax.set_ylabel(f"score, replicate {reps[1]}", fontsize=7.5)
            _ax_style(ax)
            s.append(_fig(fig, width=100 * mm))
        top = W_[W_["hit"]].sort_values("score", ascending=False).head(40)
        tr = [["Plate", "Well", "Reagent", "Gene", "Score"] + [f"rep {x}" for x in reps]] + [
            [int(x.plate), x.well, x.reagent, x.gene, _fmt(x.score)] + [_fmt(x[f"score_{k}"]) for k in reps] for _, x in top.iterrows()]
        s.append(Paragraph("Strongest hit wells", st["h3"]))
        s.append(_table(tr, [14 * mm, 14 * mm, 30 * mm, 30 * mm, 20 * mm] + [20 * mm] * len(reps), st))
        cmp = r.get("comparison")
        if cmp:
            s.append(Paragraph("Compared with the published analysis", st["h3"]))
            s.append(Paragraph(f"r = {cmp['r']:.4f} with the published scores over {cmp['wells']:,} wells; published hits "
                               f"(score ≥ 2) {cmp['published_hits']}, this run {cmp['our_hits']}, in both {cmp['both']} "
                               f"(recall {100 * cmp['recall']:.1f}%).", st["p"]))
    doc.build(s, onFirstPage=_footer(name), onLaterPages=_footer(name))
    return buf.getvalue()
