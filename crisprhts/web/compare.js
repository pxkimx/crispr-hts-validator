/* Compare two pooled screens: shared dependencies, and the ones that belong to only one of them.
   Adds the Compare page. Extends app.js (go, renderResults, api, charts). */
(() => {
  let CMP = null, selected = null, sortCol = "normZ", sortAsc = true;
  const CALL_KEY = { shared: "ess", "specific to A": "neg", "specific to B": "pos", neither: "ns" };

  SECTION.pooled.compare = "p-compare";
  SECTION.arrayed.compare = "p-aresults";      // no arrayed comparison yet; the tab is hidden there

  const _go = go;
  go = function (p) {
    _go(p);
    $("#tabs button[data-p=compare]").classList.toggle("hidden", MODE !== "pooled");
    if (p === "compare" && MODE === "pooled" && (window.XCMODE || "in") === "in") renderCompare();
  };
  const _setMode = setMode;
  setMode = function (m) {
    _setMode(m);
    $("#tabs button[data-p=compare]").classList.toggle("hidden", m !== "pooled");
  };
  $("#tabs button[data-p=compare]").classList.toggle("hidden", MODE !== "pooled");

  /* ---------------- set-up */
  function renderCompare() {
    const comps = (S && S.design && S.design.comparisons) || [];
    $("#cmpEmpty").classList.toggle("hidden", comps.length >= 2);
    $("#cmpBody").classList.toggle("hidden", comps.length < 2);
    if (comps.length < 2) return;
    for (const [sel, pref] of [["#cmpA", 1], ["#cmpB", 2]]) {
      const el = $(sel);
      const keep = el.value;
      el.innerHTML = comps.map((c) => `<option>${esc(c.name)}</option>`).join("");
      const shown = CMP && (pref === 1 ? CMP.a.name : CMP.b.name);
      el.value = shown && comps.some((c) => c.name === shown) ? shown
        : keep && comps.some((c) => c.name === keep) ? keep
          : (comps[Math.min(pref, comps.length) - 1] || comps[0]).name;
    }
    const th = (S && S.compare_thresholds) || {};
    const def = (S && S.compare_defaults) || {};
    $("#cmpTh").innerHTML = [
      ["cmp_fdr", "FDR at most"], ["cmp_min_delta", "smallest |Δ| worth calling"],
      ["cmp_min_guides", "guides per gene, at least"], ["cmp_dependency", "a dependency scores at or below"],
    ].map(([k, l]) => `<label class="small">${l}<input type="text" data-cth="${k}" value="${th[k] ?? def[k] ?? ""}" size="5"></label>`).join("");
    if (!CMP && S.has_comparison) load();
    else if (CMP) draw();
  }
  $("#cmpRun").onclick = async () => {
    const thresholds = {};
    $$("[data-cth]").forEach((i) => { if (i.value.trim() !== "" && isFinite(+i.value)) thresholds[i.dataset.cth] = +i.value; });
    try {
      const j = await api("/api/compare", { a: $("#cmpA").value, b: $("#cmpB").value, thresholds });
      watchJob(j.id, async () => { S = await api("/api/state"); await load(); });
    } catch (e) { toast(e.message, true); }
  };
  async function load(search) {
    CMP = await api("/api/compare" + (search ? "?search=" + encodeURIComponent(search) : ""));
    selected = null;
    for (const [sel, nm] of [["#cmpA", CMP.a.name], ["#cmpB", CMP.b.name]]) {
      if ([...$(sel).options].some((o) => o.value === nm)) $(sel).value = nm;
    }
    draw();
  }
  $("#cmpSearch").oninput = debounce(() => load($("#cmpSearch").value.trim()), 250);

  /* ---------------- the result */
  function draw() {
    if (!CMP) { $("#cmpResult").classList.add("hidden"); return; }
    $("#cmpResult").classList.remove("hidden");
    const c = CMP.concordance, th = CMP.thresholds;
    const tile = ([k, v, s]) => `<div class="card tile"><div class="k">${k}</div><div class="v">${esc(v)}</div><div class="s">${esc(s)}</div></div>`;
    $("#cmpTiles").innerHTML = [
      ["Shared dependencies", c.dep_both.toLocaleString(), `of ${c.dep_a.toLocaleString()} in A and ${c.dep_b.toLocaleString()} in B`],
      ["Specific to A", c.n_specific_a.toLocaleString(), `a dependency in ${esc(CMP.a.name)} only, and the difference passes the test`],
      ["Specific to B", c.n_specific_b.toLocaleString(), `a dependency in ${esc(CMP.b.name)} only, and the difference passes the test`],
      ["Shared but deeper", (c.n_shared_stronger ?? 0).toLocaleString(), "needed by both, more by one"],
    ].map(tile).join("");

    // the comparison is only as good as the weaker screen
    const q = CMP.quality || [];
    const weak = q.filter((x) => x.status !== "pass");
    $("#cmpQuality").innerHTML = !q.length
      ? `<span class="muted">Run the quality checks to see how well each screen separates essential from non-essential genes.</span>`
      : q.map((x) => `${statusTag(x.status)} <b>${esc(x.comparison)}</b> <span class="muted">essential-gene AUC ${x.auc.toFixed(3)}</span>`).join('<span class="faint"> · </span>')
        + (weak.length ? `<div class="small warnline mt">${weak.map((x) => esc(x.comparison)).join(" and ")} separates essential genes poorly, so differences involving it are unreliable — an early timepoint usually looks like this.</div>` : "");

    const scaleNote = [CMP.a, CMP.b].filter((x) => !x.scale.ok);
    $("#cmpNotes").innerHTML = `<b>${esc(CMP.a.name)}</b> (A) against <b>${esc(CMP.b.name)}</b> (B) · ${CMP.total.toLocaleString()} genes with ≥ ${th.cmp_min_guides} guides · run ${esc(CMP.when)}`
      + (CMP.shared_reference ? "<br>Both screens share one reference sample, so each guide's difference is a direct treated-vs-untreated contrast." : "<br>The screens use different reference samples; each is measured against its own.")
      + (scaleNote.length ? `<br><b>Not put on the common scale:</b> ${scaleNote.map((x) => esc(x.name)).join(" and ")} has too few reference essential or non-essential genes, so the values below are raw log₂ fold changes.`
        : `<br>Both screens scaled so the median core-essential gene is −1 and the median non-essential gene is 0 (${esc(CMP.species)} reference sets).`);

    drawScatter(); drawVolcano(); drawConcordance(); drawTable();
    if (selected) showGeneCmp(selected);
  }

  function colours() { return { ess: C.css("--ess"), neg: C.css("--neg"), pos: C.css("--pos"), ns: C.css("--ns") }; }

  function drawScatter() {
    const pts = CMP.points.map(([g, a, b, d, fdr, call, diff]) => ({ x: a, y: b, k: CALL_KEY[call] || "ns", label: g, d, fdr, call, diff }));
    const top = pts.filter((p) => p.diff).sort((a, b) => Math.abs(b.d) - Math.abs(a.d)).slice(0, 10);
    C.scatter($("#cmpScatter"), pts, {
      colors: colours(), xl: `${CMP.a.name} — scaled effect`, yl: `${CMP.b.name} — scaled effect`, height: 360,
      diag: true, labels: top, selected, z: (p) => (p.k === "ns" ? 0 : p.diff ? 2 : 1),
      tip: (p) => `<b>${esc(p.label)}</b><br>A ${fmt(p.x)} · B ${fmt(p.y)}<br>Δ ${fmt(p.d)} · FDR ${fmt(p.fdr)}<br>${esc(p.call)}`,
      click: (p) => showGeneCmp(p.label),
    });
    $("#cmpScatter").insertAdjacentHTML("beforeend",
      `<div class="legend"><span><i style="background:${colours().ess}"></i>dependency in both</span><span><i style="background:${colours().neg}"></i>in ${esc(CMP.a.name)} only</span><span><i style="background:${colours().pos}"></i>in ${esc(CMP.b.name)} only</span><span><i style="background:${colours().ns}"></i>neither</span><span class="faint">colour shows which screen needs the gene; the tiles above count only those that also pass the test</span></div>`);
  }

  function drawVolcano() {
    const pts = CMP.points.map(([g, a, b, d, fdr, call, diff]) => ({
      x: d, y: -Math.log10(Math.max(fdr, 1e-300)), k: diff ? (d < 0 ? "pos" : "neg") : "ns", label: g, call, fdr,
    }));
    const top = pts.filter((p) => p.k !== "ns").sort((a, b) => b.y - a.y).slice(0, 10);
    C.scatter($("#cmpVolcano"), pts, {
      colors: colours(), xl: `Δ scaled effect (B − A)`, yl: "−log₁₀ FDR", height: 360, vline: 0,
      labels: top, selected, z: (p) => (p.k === "ns" ? 0 : 1),
      tip: (p) => `<b>${esc(p.label)}</b><br>Δ ${fmt(p.x)} · FDR ${fmt(p.fdr)}<br>${esc(p.call)}`,
      click: (p) => showGeneCmp(p.label),
    });
    $("#cmpVolcano").insertAdjacentHTML("beforeend",
      `<div class="legend"><span><i style="background:${colours().pos}"></i>needed more by ${esc(CMP.b.name)}</span><span><i style="background:${colours().neg}"></i>needed more by ${esc(CMP.a.name)}</span><span><i style="background:${colours().ns}"></i>no difference</span></div>`);
  }

  function drawConcordance() {
    const c = CMP.concordance, th = CMP.thresholds;
    $("#cmpConc").innerHTML = `<table class="tbl">
      <tr><th>How alike are the two screens?</th><th class="num"></th></tr>
      <tr><td>Pearson r of scaled gene effects</td><td class="num">${fmt(c.pearson)}</td></tr>
      <tr><td>Spearman ρ</td><td class="num">${fmt(c.spearman)}</td></tr>
      <tr><td>Dependencies shared (Jaccard)</td><td class="num">${fmt(c.jaccard)}</td></tr>
      <tr><td>Genes called different</td><td class="num">${c.n_different.toLocaleString()}</td></tr>
      <tr><td class="muted small" colspan="2">A gene counts as a dependency at a scaled effect of ${th.cmp_dependency} or below, where −1 is the median core-essential gene. A difference is called at FDR ≤ ${th.cmp_fdr} with |Δ| ≥ ${th.cmp_min_delta}.</td></tr>
      </table>`;
  }

  function drawTable() {
    const I = Object.fromEntries(CMP.columns.map((c, i) => [c, i]));
    const rows = [...CMP.rows].sort((a, b) => {
      const x = a[I[sortCol]], y = b[I[sortCol]];
      const v = typeof x === "number" ? x - y : String(x).localeCompare(String(y));
      return sortAsc ? v : -v;
    });
    const head = [["gene", "Gene"], ["guides", "Guides"], ["scaled_a", "A"], ["scaled_b", "B"], ["delta", "Δ"],
      ["normZ", "normZ"], ["fdr", "FDR"], ["agreement", "Guides agreeing"], ["call", "Needed by"]];
    $("#cmpTable").innerHTML = `<tr>${head.map(([k, l]) => `<th class="${k === "gene" || k === "call" ? "" : "num"} click" data-cs="${k}">${l}${sortCol === k ? (sortAsc ? " ▲" : " ▼") : ""}</th>`).join("")}</tr>`
      + rows.map((r) => `<tr class="click ${r[I.gene] === selected ? "sel" : ""}" data-cg="${esc(r[I.gene])}">
          <td><b>${esc(r[I.gene])}</b></td><td class="num">${r[I.guides]}</td>
          <td class="num">${fmt(r[I.scaled_a])}</td><td class="num">${fmt(r[I.scaled_b])}</td>
          <td class="num"><b>${fmt(r[I.delta])}</b></td><td class="num">${fmt(r[I.normZ])}</td>
          <td class="num">${fmt(r[I.fdr])}</td>
          <td class="num ${r[I.agreement] <= 0.5 ? "weak" : ""}">${Math.round(r[I.agreement] * r[I.guides])} of ${r[I.guides]}</td>
          <td>${esc(r[I.call])}${r[I.stronger_in] ? ` <span class="faint">(more in ${esc(r[I.stronger_in])})</span>` : ""}</td></tr>`).join("");
    $("#cmpCount").textContent = `${CMP.rows.length.toLocaleString()} shown of ${CMP.total.toLocaleString()}`;
    $$("#cmpTable [data-cs]").forEach((h) => (h.onclick = () => {
      if (sortCol === h.dataset.cs) sortAsc = !sortAsc; else { sortCol = h.dataset.cs; sortAsc = true; }
      drawTable();
    }));
    $$("#cmpTable [data-cg]").forEach((tr) => (tr.onclick = () => showGeneCmp(tr.dataset.cg)));
  }

  async function showGeneCmp(gene) {
    selected = gene;
    let g;
    try { g = await api("/api/compare/gene?gene=" + encodeURIComponent(gene)); }
    catch (e) { return toast(e.message, true); }
    const s = g.summary;
    $("#cmpGene").innerHTML = `<div class="row"><h3>${esc(gene)}</h3><span class="muted small">${s.guides} guides · ${esc(s.call)}${s.stronger_in ? `, more in ${esc(s.stronger_in)}` : ""}</span></div>
      <div class="small muted">scaled effect ${fmt(s.scaled_a)} in ${esc(g.a)} and ${fmt(s.scaled_b)} in ${esc(g.b)} · Δ ${fmt(s.delta)} · normZ ${fmt(s.normZ)} · FDR ${fmt(s.fdr)}</div>
      <div class="plot mt" id="cmpGuides"></div>`;
    pairBars($("#cmpGuides"), g.guides.map((x) => ({ label: x.guide, a: x.lfc_a, b: x.lfc_b, d: x.delta_guide })), g.a, g.b);
    drawTable();
  }

  /* a paired bar per guide: its log fold change in each screen, side by side.
     Named apart from charts.js's guideBars, which the Results page uses for one screen. */
  function pairBars(el, rows, aName, bName) {
    if (!rows.length) { el.innerHTML = ""; return; }
    const all = rows.flatMap((r) => [r.a, r.b]).filter(isFinite);
    const lo = Math.min(-0.5, ...all), hi = Math.max(0.5, ...all), span = hi - lo || 1;
    const x = (v) => ((v - lo) / span) * 100;
    el.innerHTML = `<div class="gbars">${rows.map((r) => `
      <div class="gb"><span class="nm">${esc(r.label)}</span>
        <div class="tr"><i class="z" style="left:${x(0)}%"></i>
          <i class="ba" style="left:${Math.min(x(r.a), x(0))}%;width:${Math.abs(x(r.a) - x(0))}%"></i>
          <i class="bb" style="left:${Math.min(x(r.b), x(0))}%;width:${Math.abs(x(r.b) - x(0))}%"></i></div>
        <span class="vv">${fmt(r.a)} / ${fmt(r.b)}</span></div>`).join("")}
      <div class="legend"><span><i style="background:${C.css("--neg")}"></i>${esc(aName)}</span><span><i style="background:${C.css("--pos")}"></i>${esc(bName)}</span><span class="faint">guide log\u2082 fold change</span></div></div>`;
  }

  function debounce(fn, ms) { let t; return () => { clearTimeout(t); t = setTimeout(fn, ms); }; }
  window.CMP_RENDER = renderCompare;          // cross.js calls this when the mode switches back
})();
