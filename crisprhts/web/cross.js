/* Compare a screen here with one in another project — and so, in general, in another library.
   Owns the mode switch on the Compare page. Extends app.js (go, api, charts) and compare.js (CMP_RENDER). */
(() => {
  let XC = null, selected = null, sortCol = "t", sortAsc = true;
  let bTreat = new Set(), bCtrl = new Set();
  const CALL_KEY = { shared: "ess", "specific to A": "neg", "specific to B": "pos", neither: "ns" };
  window.XCMODE = "in";

  /* ---------------- the mode switch */
  $$("#cmpPill button").forEach((b) => (b.onclick = () => {
    $$("#cmpPill button").forEach((x) => x.classList.toggle("on", x === b));
    window.XCMODE = b.dataset.m;
    $("#cmpIn").classList.toggle("hidden", window.XCMODE !== "in");
    $("#cmpAcross").classList.toggle("hidden", window.XCMODE !== "across");
    if (window.XCMODE === "across") renderCross(); else window.CMP_RENDER && window.CMP_RENDER();
  }));

  const _go = go;
  go = function (p) {
    _go(p);
    if (p === "compare" && MODE === "pooled" && window.XCMODE === "across") renderCross();
  };

  /* ---------------- set-up */
  function renderCross() {
    const pt = S && S.partner;
    $("#xcClear").classList.toggle("hidden", !pt);
    $("#xcSetup").classList.toggle("hidden", !pt);
    if (!pt) {
      $("#xcPartner").innerHTML = `<span class="muted">No second project loaded yet.</span>`;
      $("#xcMatch").innerHTML = `<span class="muted">Load the other project&rsquo;s count table to see how many genes the two libraries have in common.</span>`;
      $("#xcResult").classList.add("hidden");
      return;
    }
    $("#xcPartner").innerHTML = `<b>${esc(pt.name)}</b><br><span class="muted">${pt.guides.toLocaleString()} guides · ${pt.genes.toLocaleString()} genes · ${pt.samples.length} samples</span>`
      + (pt.notes && pt.notes.length ? `<div class="faint small mt">${pt.notes.map(esc).join("<br>")}</div>` : "");
    drawMatch(pt);

    // this project's comparisons
    const comps = (S.design && S.design.comparisons) || [];
    const el = $("#xcA"), keep = el.value;
    el.innerHTML = comps.map((c) => `<option>${esc(c.name)}</option>`).join("");
    const want = XC && XC.a && XC.a.name;
    el.value = want && comps.some((c) => c.name === want) ? want
      : keep && comps.some((c) => c.name === keep) ? keep : (comps[0] || {}).name || "";
    if (!comps.length) $("#xcA").innerHTML = `<option disabled>Set up a comparison on the Samples page first</option>`;

    // the other project's samples, pre-filled from its guessed design
    if (!bTreat.size && !bCtrl.size) {
      const g = (pt.comparisons || [])[0];
      if (g) { bTreat = new Set(g.treatment); bCtrl = new Set(g.control); if (!$("#xcName").value) $("#xcName").value = `${pt.name} — ${g.name}`; }
    }
    drawChips(pt);

    const th = (S.cross_thresholds) || {}, def = (S.cross_defaults) || {};
    $("#xcTh").innerHTML = [
      ["xc_fdr", "FDR at most"], ["xc_min_delta", "smallest |Δ| worth calling"],
      ["xc_min_guides", "guides per gene in each library, at least"], ["xc_dependency", "a dependency scores at or below"],
    ].map(([k, l]) => `<label class="small">${l}<input type="text" data-xth="${k}" value="${th[k] ?? def[k] ?? ""}" size="5"></label>`).join("");

    if (!XC && S.has_cross) load(); else if (XC) draw();
  }

  function drawMatch(pt) {
    const sh = pt.shared;
    if (!sh) { $("#xcMatch").innerHTML = `<span class="muted">Load this project&rsquo;s count table too.</span>`; return; }
    const frac = sh.n / Math.max(1, Math.min(sh.n_a, sh.n_b));
    const st = frac >= 0.6 ? "pass" : frac >= 0.25 ? "warn" : "fail";
    $("#xcMatch").innerHTML = `${statusTag(st)} <b>${sh.n.toLocaleString()}</b> genes in both`
      + `<div class="muted small mt">${sh.n_a.toLocaleString()} here (${esc(sh.kind_a)}) · ${sh.n_b.toLocaleString()} there (${esc(sh.kind_b)})`
      + (sh.n_mapped ? ` · ${sh.n_mapped.toLocaleString()} joined by the mapping file` : "") + `</div>`
      + (sh.note ? `<div class="small warnline mt">${esc(sh.note)}</div>` : "")
      + (sh.unmatched_b && sh.unmatched_b.length
        ? `<div class="faint small mt">Not matched, there: ${sh.unmatched_b.slice(0, 6).map(esc).join(", ")}${sh.n_b - sh.n > 6 ? " …" : ""}</div>` : "");
    $("#xcMapInfo").textContent = pt.mapping
      ? `${esc(pt.mapping)} — ${(pt.mapping_rows || 0).toLocaleString()} name pairs`
      : "two columns: a name in one library, the matching name in the other";
  }

  function drawChips(pt) {
    $("#xcChips").innerHTML = pt.samples.map((s) =>
      `<button class="chip ${bTreat.has(s) ? "t" : bCtrl.has(s) ? "c" : ""}" data-xs="${esc(s)}">${esc(s)}${bTreat.has(s) ? " · treated" : bCtrl.has(s) ? " · reference" : ""}</button>`).join("");
    $$("#xcChips [data-xs]").forEach((b) => (b.onclick = () => {
      const s = b.dataset.xs;
      if (bTreat.has(s)) { bTreat.delete(s); bCtrl.add(s); }
      else if (bCtrl.has(s)) bCtrl.delete(s);
      else bTreat.add(s);
      drawChips(pt);
    }));
  }

  $("#xcClear").onclick = async () => {
    try { S = await api("/api/cross/clear", {}); XC = null; bTreat = new Set(); bCtrl = new Set(); renderCross(); }
    catch (e) { toast(e.message, true); }
  };

  $("#xcRun").onclick = async () => {
    const thresholds = {};
    $$("[data-xth]").forEach((i) => { if (i.value.trim() !== "" && isFinite(+i.value)) thresholds[i.dataset.xth] = +i.value; });
    try {
      const j = await api("/api/cross", {
        a: $("#xcA").value, b_name: $("#xcName").value.trim(),
        b_treatment: [...bTreat], b_control: [...bCtrl], thresholds,
      });
      watchJob(j.id, async () => { S = await api("/api/state"); await load(); });
    } catch (e) { toast(e.message, true); }
  };

  async function load(search) {
    XC = await api("/api/cross" + (search ? "?search=" + encodeURIComponent(search) : ""));
    selected = null;
    draw();
  }
  $("#xcSearch").oninput = debounce(() => load($("#xcSearch").value.trim()), 250);

  /* ---------------- the result */
  function draw() {
    if (!XC) { $("#xcResult").classList.add("hidden"); return; }
    $("#xcResult").classList.remove("hidden");
    const c = XC.concordance, p = XC.power, th = XC.thresholds;

    // what you have to know before reading the table
    $("#xcChecks").innerHTML = `<div class="row"><h3>Can these two screens be compared?</h3><span class="spacer"></span><span class="faint small">read this first</span></div>`
      + `<table class="tbl">` + XC.notes.map((n) =>
        `<tr><td>${statusTag(n.status)} <b>${esc(n.label)}</b></td><td class="num"><b>${esc(n.value)}</b></td><td class="muted small">${esc(n.why)}</td></tr>`).join("") + `</table>`;

    const tile = ([k, v, s]) => `<div class="card tile"><div class="k">${k}</div><div class="v">${esc(v)}</div><div class="s">${esc(s)}</div></div>`;
    $("#xcTiles").innerHTML = [
      ["Shared dependencies", c.dep_both.toLocaleString(), `of ${c.dep_a.toLocaleString()} in A and ${c.dep_b.toLocaleString()} in B`],
      ["Specific to A", c.n_specific_a.toLocaleString(), `a dependency in ${XC.a.name} only, and the difference passes the test`],
      ["Specific to B", c.n_specific_b.toLocaleString(), `a dependency in ${XC.b.name} only, and the difference passes the test`],
      ["Shared but deeper", (c.n_shared_stronger ?? 0).toLocaleString(), "needed by both, more by one"],
    ].map(tile).join("");

    drawScatter(); drawVolcano(); drawConc(); drawTable();
    $("#xcFoot").innerHTML = `<b>${esc(XC.a.project || "this project")} — ${esc(XC.a.name)}</b> (A) against <b>${esc(XC.b.project || "the other project")} — ${esc(XC.b.name)}</b> (B) · ${XC.total.toLocaleString()} genes with ≥ ${th.xc_min_guides} guides in each library · run ${esc(XC.when)}`
      + `<br>Each screen is measured against its own reference samples, then scaled so that the median core-essential gene is −1 and the median non-essential gene is 0 (${esc(XC.species)} reference sets) — the only thing two different libraries share. A gene&rsquo;s effect is the median of its guides, and the difference B − A is judged against the spread of those guides, widened by the trend that genes of that strength usually show and by the error in each screen&rsquo;s own scale, on a t distribution with ${p.median_df.toFixed(0)} degrees of freedom (prior strength d₀ = ${p.prior_df.toFixed(1)}, estimated from the data).`
      + (p.calibration > 1.001 ? `<br>The null had to be widened ${p.calibration.toFixed(2)}× to behave, so something systematic separates these screens beyond guide noise.` : "");
    if (selected) showGene(selected);
  }

  function colours() { return { ess: C.css("--ess"), neg: C.css("--neg"), pos: C.css("--pos"), ns: C.css("--ns") }; }

  function drawScatter() {
    const pts = XC.points.map(([g, a, b, d, fdr, call, diff]) => ({ x: a, y: b, k: CALL_KEY[call] || "ns", label: g, d, fdr, call, diff }));
    const top = pts.filter((p) => p.diff).sort((a, b) => Math.abs(b.d) - Math.abs(a.d)).slice(0, 10);
    C.scatter($("#xcScatter"), pts, {
      colors: colours(), xl: `${XC.a.name} — scaled effect`, yl: `${XC.b.name} — scaled effect`, height: 360,
      diag: true, labels: top, selected, z: (p) => (p.k === "ns" ? 0 : p.diff ? 2 : 1),
      tip: (p) => `<b>${esc(p.label)}</b><br>A ${fmt(p.x)} · B ${fmt(p.y)}<br>Δ ${fmt(p.d)} · FDR ${fmt(p.fdr)}<br>${esc(p.call)}`,
      click: (p) => showGene(p.label),
    });
    $("#xcScatter").insertAdjacentHTML("beforeend",
      `<div class="legend"><span><i style="background:${colours().ess}"></i>dependency in both</span><span><i style="background:${colours().neg}"></i>in A only</span><span><i style="background:${colours().pos}"></i>in B only</span><span><i style="background:${colours().ns}"></i>neither</span><span class="faint">colour shows which screen needs the gene; the tiles above count only those that also pass the test</span></div>`);
  }

  function drawVolcano() {
    const res = XC.power.resolution;
    const pts = XC.points.map(([g, a, b, d, fdr, call, diff]) => ({
      x: d, y: -Math.log10(Math.max(fdr, 1e-300)), k: diff ? (d < 0 ? "pos" : "neg") : "ns", label: g, call, fdr,
    }));
    const top = pts.filter((p) => p.k !== "ns").sort((a, b) => b.y - a.y).slice(0, 10);
    // when nothing is significant every FDR is 1 and the axis would span floating-point dust, so hold it
    // open as far as the threshold line at least
    const bar = -Math.log10(XC.thresholds.xc_fdr);
    const hi = Math.max(bar * 1.3, ...pts.map((p) => p.y).filter(isFinite));
    C.scatter($("#xcVolcano"), pts, {
      colors: colours(), xl: `Δ scaled effect (B − A)`, yl: "−log₁₀ FDR", height: 360, vline: 0,
      yr: [0, hi], hline: bar,
      labels: top, selected, z: (p) => (p.k === "ns" ? 0 : 1),
      tip: (p) => `<b>${esc(p.label)}</b><br>Δ ${fmt(p.x)} · FDR ${fmt(p.fdr)}<br>${esc(p.call)}`,
      click: (p) => showGene(p.label),
    });
    $("#xcVolcano").insertAdjacentHTML("beforeend",
      `<div class="legend"><span><i style="background:${colours().pos}"></i>needed more by B</span><span><i style="background:${colours().neg}"></i>needed more by A</span><span><i style="background:${colours().ns}"></i>no difference</span><span class="faint">nothing inside ±${fmt(res)} can be called with these two libraries</span></div>`);
  }

  function drawConc() {
    const c = XC.concordance, p = XC.power, th = XC.thresholds, m = XC.matching;
    $("#xcConc").innerHTML = `<table class="tbl">
      <tr><th>How alike are the two screens?</th><th class="num"></th></tr>
      <tr><td>Pearson r of scaled gene effects</td><td class="num">${fmt(c.pearson)}</td></tr>
      <tr><td>Spearman ρ</td><td class="num">${fmt(c.spearman)}</td></tr>
      <tr><td>Dependencies shared (Jaccard)</td><td class="num">${fmt(c.jaccard)}</td></tr>
      <tr><td>Genes called different</td><td class="num">${c.n_different.toLocaleString()}</td></tr>
      <tr><th>What this pair can resolve</th><th class="num"></th></tr>
      <tr><td>Smallest difference callable</td><td class="num"><b>${fmt(p.resolution)}</b></td></tr>
      <tr><td>Typical standard error of Δ</td><td class="num">${fmt(p.median_se)}</td></tr>
      <tr><td>Scale precision (A / B)</td><td class="num">±${(p.scale_rel_a * 100).toFixed(1)}% / ±${(p.scale_rel_b * 100).toFixed(1)}%</td></tr>
      <tr><td>Spread of t among non-essential genes</td><td class="num">${fmt(p.null_sd_neg)}</td></tr>
      <tr><td class="muted small" colspan="2">A gene counts as a dependency at a scaled effect of ${th.xc_dependency} or below, where −1 is the median core-essential gene. A difference is called at FDR ≤ ${th.xc_fdr} with |Δ| ≥ ${th.xc_min_delta}, over the ${m.n_direct.toLocaleString()} genes matched by name${m.n_mapped ? ` and ${m.n_mapped.toLocaleString()} by the mapping file` : ""}. The spread of t among non-essential genes should sit near 1: well below means the two screens are not independent measurements, and the resolution above is then pessimistic.</td></tr>
      </table>`;
  }

  function drawTable() {
    const I = Object.fromEntries(XC.columns.map((c, i) => [c, i]));
    const rows = [...XC.rows].sort((a, b) => {
      const x = a[I[sortCol]], y = b[I[sortCol]];
      const v = typeof x === "number" ? x - y : String(x).localeCompare(String(y));
      return sortAsc ? v : -v;
    });
    const head = [["gene", "Gene"], ["guides_a", "Guides A/B"], ["scaled_a", "A"], ["scaled_b", "B"],
      ["delta", "Δ"], ["se_delta", "± SE"], ["t", "t"], ["fdr", "FDR"], ["call", "Needed by"]];
    $("#xcTable").innerHTML = `<tr>${head.map(([k, l]) => `<th class="${k === "gene" || k === "call" ? "" : "num"} click" data-xs2="${k}">${l}${sortCol === k ? (sortAsc ? " ▲" : " ▼") : ""}</th>`).join("")}</tr>`
      + rows.map((r) => `<tr class="click ${r[I.gene] === selected ? "sel" : ""}" data-xg="${esc(r[I.gene])}">
          <td><b>${esc(r[I.gene])}</b></td><td class="num">${r[I.guides_a]} / ${r[I.guides_b]}</td>
          <td class="num">${fmt(r[I.scaled_a])}</td><td class="num">${fmt(r[I.scaled_b])}</td>
          <td class="num"><b>${fmt(r[I.delta])}</b></td><td class="num faint">${fmt(r[I.se_delta])}</td>
          <td class="num">${fmt(r[I.t])}</td><td class="num">${fmt(r[I.fdr])}</td>
          <td>${esc(r[I.call])}${r[I.stronger_in] ? ` <span class="faint">(more in ${esc(r[I.stronger_in])})</span>` : ""}</td></tr>`).join("");
    $("#xcCount").textContent = `${XC.rows.length.toLocaleString()} shown of ${XC.total.toLocaleString()}`;
    $$("#xcTable [data-xs2]").forEach((h) => (h.onclick = () => {
      if (sortCol === h.dataset.xs2) sortAsc = !sortAsc; else { sortCol = h.dataset.xs2; sortAsc = true; }
      drawTable();
    }));
    $$("#xcTable [data-xg]").forEach((tr) => (tr.onclick = () => showGene(tr.dataset.xg)));
  }

  async function showGene(gene) {
    selected = gene;
    let g;
    try { g = await api("/api/cross/gene?gene=" + encodeURIComponent(gene)); }
    catch (e) { return toast(e.message, true); }
    const s = g.summary;
    $("#xcGene").innerHTML = `<div class="row"><h3>${esc(gene)}</h3><span class="muted small">${s.guides_a} guides in A, ${s.guides_b} in B · ${esc(s.call)}${s.stronger_in ? `, more in ${esc(s.stronger_in)}` : ""}</span></div>
      <div class="small muted">scaled effect ${fmt(s.scaled_a)} in ${esc(g.a)} and ${fmt(s.scaled_b)} in ${esc(g.b)} · Δ ${fmt(s.delta)} ± ${fmt(s.se_delta)} · t ${fmt(s.t)} on ${Math.round(s.df)} df · FDR ${fmt(s.fdr)}</div>
      <div class="faint small mt">Different sequences on each side, so look at how tightly each library&rsquo;s own guides agree: guides that disagree inside one library are what a cross-library difference is most often made of.</div>
      <div class="plot mt" id="xcGuidesA"></div><div class="plot" id="xcGuidesB"></div>`;
    libBars($("#xcGuidesA"), g.guides_a, g.a, "--neg");
    libBars($("#xcGuidesB"), g.guides_b, g.b, "--pos");
    drawTable();
  }

  /* one library's guides for a gene, on the common scale. Kept apart from compare.js's paired bars: here
     there is nothing to pair, and pretending otherwise would be the whole mistake. */
  function libBars(el, rows, name, tok) {
    if (!rows || !rows.length) { el.innerHTML = `<span class="faint small">No guides for this gene in ${esc(name)}.</span>`; return; }
    const vals = rows.map((r) => r.scaled).filter(isFinite);
    const lo = Math.min(-0.5, ...vals), hi = Math.max(0.5, ...vals), span = hi - lo || 1;
    const x = (v) => ((v - lo) / span) * 100;
    el.innerHTML = `<div class="gbars"><div class="small"><b>${esc(name)}</b> <span class="faint">${rows.length} guides</span></div>`
      + rows.map((r) => `<div class="gb"><span class="nm">${esc(r.guide)}</span>
        <div class="tr"><i class="z" style="left:${x(0)}%"></i>
          <i class="ba" style="left:${Math.min(x(r.scaled), x(0))}%;width:${Math.abs(x(r.scaled) - x(0))}%;background:${C.css(tok)}"></i></div>
        <span class="vv">${fmt(r.scaled)}</span></div>`).join("")
      + `<div class="legend"><span class="faint">scaled effect per guide, −1 = median core-essential gene</span></div></div>`;
  }

  function debounce(fn, ms) { let t; return () => { clearTimeout(t); t = setTimeout(fn, ms); }; }
})();
