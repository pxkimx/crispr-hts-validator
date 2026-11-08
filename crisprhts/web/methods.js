/* Version 0.3 additions to the pooled pages: method choice (MAGeCK, BAGEL2, drugZ), copy-number correction,
   method-specific results, pathway enrichment for pooled and arrayed results. Extends app.js and arrayed.js. */
(() => {
  let METHOD = "mageck";
  cut = function () { const e = $("#fdrCut"); const v = e ? +e.value : 0.05; return v > 0 && v < 1 ? v : 0.05; };

  /* ---------------- data page: guide positions */
  const _renderData = renderData;
  renderData = function () {
    _renderData();
    const p = S && S.positions;
    $("#posInfo").textContent = p ? `${p.name}: ${p.guides.toLocaleString()} guides on ${p.chromosomes} chromosomes` + (p.matched !== null ? `, ${p.matched.toLocaleString()} match this screen` : "") : "none loaded";
  };

  /* ---------------- analyse page */
  const _renderAnalyse = renderAnalyse;
  renderAnalyse = function () {
    _renderAnalyse();
    if (!S) return;
    $("#methodPick").innerHTML = Object.entries(S.methods).map(([k, v]) => {
      const [title, ...rest] = v.split(":");
      return `<label class="${k === METHOD ? "on" : ""}"><input type="radio" name="method" value="${k}" ${k === METHOD ? "checked" : ""}><b>${esc(title)}</b><span>${esc(rest.join(":").trim())}</span></label>`;
    }).join("");
    $$("input[name=method]").forEach((i) => (i.onchange = () => { METHOD = i.value; renderAnalyse(); }));
    const pos = S.positions;
    $("#cnCorrect").disabled = !pos || METHOD === "drugz";
    if (!pos || METHOD === "drugz") $("#cnCorrect").checked = false;
    $("#cnNote").textContent = METHOD === "drugz" ? "(not used for drugZ: both arms share the same copy number)"
      : pos ? `(${pos.matched.toLocaleString()} guides with positions; off by default — it can also remove clusters of neighbouring essential genes)`
        : "(load guide positions on the Data page)";
    $("#drugzOpts").classList.toggle("hidden", METHOD !== "drugz");
    $("#mageckOpts").querySelector(".card").classList.toggle("hidden", METHOD !== "mageck");
    $$("#mageckOpts .form").forEach((f) => f.classList.toggle("hidden", METHOD !== "mageck"));
    $("#methodHint").textContent = METHOD === "bagel" ? "BAGEL2 compares each comparison's treatment samples with its control (plasmid / T0); a gene's Bayes factor sums its guides and replicates."
      : METHOD === "drugz" ? "drugZ: the comparison's control = untreated, treatment = treated (same time point). Replicates are paired in the order listed on the Samples page."
        : "";
  };
  $("#runBtn").onclick = async () => {
    const opts = {
      method: METHOD, cn_correct: $("#cnCorrect").checked, drugz_paired: $("#dzPaired").checked, drugz_remove: $("#dzRemove").value,
      norm: ($("input[name=norm]:checked") || {}).value || "median", controls: ($("input[name=ctl]:checked") || {}).value,
      control_pattern: $("#ctlPattern").value, control_list: $("#ctlList").value, control_null: $("#ctlNull").checked,
      remove_zero: $("#aZero").value, gene_lfc: $("#aLfc").value, gene_test_fdr_threshold: +$("#aAlpha").value || 0.25,
      comparisons: $$("[data-ac]:checked").map((i) => i.dataset.ac),
    };
    if (!opts.comparisons.length) return toast("Choose at least one comparison.", true);
    try {
      const j = await api("/api/analyse", opts);
      watchJob(j.id, async () => { await refresh(); R = null; toast("Analysis finished."); go("results"); });
    } catch (e) { toast(e.message, true); }
  };

  /* ---------------- results */
  const _renderResults = renderResults;
  renderResults = function () {
    $$("#rComp option").forEach((o) => { const r = (S.results || []).find((x) => x.name === o.value); if (r) o.textContent = `${r.name} · ${r.method}`; });
    if (R.method === "mageck") {
      $("#volcano").closest(".card").querySelector("h3").textContent = "Volcano";
      if (!$("#fdrCut")) { $("#volcano").closest(".card").querySelector(".row .faint").innerHTML = `gene LFC vs −log₁₀ p · FDR ≤ <input type="text" id="fdrCut" value="0.05" size="4">`; $("#fdrCut").onchange = () => renderResults(); }
      $("#rankplot").closest(".card").querySelector("h3").textContent = "Gene rank";
      $("#rankDir").classList.remove("hidden");
      _renderResults(); cnNote(); pathways("#pwCard", "/api/enrich", true); return;
    }
    renderOther();
    pathways("#pwCard", "/api/enrich", true);
  };
  function cnNote() {
    if (R.cn) $("#resNotes").insertAdjacentHTML("beforeend", `<br>Copy-number correction: ${R.cn.n_positioned.toLocaleString()} guides with positions, ${R.cn.segments} segments, ${R.cn.shifted_guides.toLocaleString()} guides shifted by more than 0.25.`);
  }
  const colorsOf = () => ({ neg: C.css("--neg"), pos: C.css("--pos"), ns: C.css("--ns"), ess: C.css("--ess"), non: C.css("--non") });
  function renderOther() {
    const I = R.idx = Object.fromEntries(R.columns.map((c, i) => [c, i]));
    const col = colorsOf();
    if (R.method === "bagel") {
      const sm = R.summary;
      $("#resTiles").innerHTML = [["Essential genes", sm.calls.toLocaleString(), `Bayes factor > ${sm.cutoff}`],
        ["Recall at BF > " + sm.cutoff, (sm.recall * 100).toFixed(1) + "%", "of reference essential genes (CEG)"],
        ["Precision at BF > " + sm.cutoff, sm.precision === null ? "–" : (sm.precision * 100).toFixed(1) + "%", "essential vs non-essential reference genes"],
        ["Genes scored", R.rows.length.toLocaleString(), `${R.species} reference sets · 10-fold cross-validation`]].map(tile).join("");
      const pts = R.rows.map((r) => ({ x: r[I.rank], y: r[I.bf], k: r[I.reference] === "essential" ? "ess" : r[I.reference] === "non-essential" ? "non" : "ns", label: r[0], r }));
      const bfs = pts.map((p) => p.y).sort((a, b) => a - b), lo = bfs[Math.floor(bfs.length * 0.005)];
      const nOut = bfs.filter((v) => v < lo).length;
      $("#volcano").closest(".card").querySelector(".row .faint").textContent = nOut ? `${nOut} genes below ${fmt(lo)} not shown` : "";
      C.scatter($("#volcano"), pts, { colors: col, xl: "gene rank", yl: "Bayes factor", hline: R.cutoff, height: 340, yr: [lo, bfs[bfs.length - 1]], z: (p) => (p.k === "ns" ? 0 : 1), selected: selGene,
        labels: pts.slice(0, 8), tip: (p) => `<b>${esc(p.label)}</b><br>BF ${fmt(p.y)} · rank ${p.x}${p.r[I.reference] ? "<br>reference: " + p.r[I.reference] : ""}`, click: (p) => showGene(p.label) });
      $("#volcano").insertAdjacentHTML("beforeend", `<div class="legend"><span><i style="background:${col.ess}"></i>reference essential</span><span><i style="background:${col.non}"></i>reference non-essential</span><span><i style="background:${col.ns}"></i>other genes</span></div>`);
      $("#volcano").closest(".card").querySelector("h3").textContent = "Bayes factors";
      $("#rankplot").closest(".card").querySelector("h3").textContent = "Precision–recall";
      $("#rankDir").classList.add("hidden");
      C.lines($("#rankplot"), [{ name: "BAGEL2", x: R.pr.recall, y: R.pr.precision, color: C.css("--accent") }], { xl: "recall", yl: "precision" });
      $("#resNotes").innerHTML = `<b>${esc(R.name)}</b> · BAGEL2 Bayes factors · treatment ${R.comparison.treatment.map(esc).join(", ")} vs ${R.comparison.control.map(esc).join(", ")} · ${esc(R.settings.norm)} · run ${esc(R.when)}`;
      if (R.cn) $("#resNotes").insertAdjacentHTML("beforeend", `<br>Copy-number correction: ${R.cn.n_positioned.toLocaleString()} guides with positions, ${R.cn.segments} segments, ${R.cn.shifted_guides.toLocaleString()} guides shifted by more than 0.25.`);
      sortCol = "rank"; sortAsc = true;
    } else {
      const f = cut();
      const ns = R.rows.filter((r) => r[I.fdr_synth] <= f).length, np_ = R.rows.filter((r) => r[I.fdr_supp] <= f).length;
      $("#resTiles").innerHTML = [["Synthetic sick / lethal", ns, `FDR ≤ ${f} (gene loss + treatment hurts more)`], ["Suppressors", np_, `FDR ≤ ${f} (gene loss protects)`],
        ["Genes scored", R.rows.length.toLocaleString(), R.paired ? "paired replicates" : "unpaired (means)"],
        ["Comparison", `${R.comparison.treatment.length} vs ${R.comparison.control.length}`, "treated vs untreated"]].map(tile).join("");
      const pts = R.rows.map((r) => { const syn = r[I.normZ] < 0; const p = syn ? r[I.p_synth] : r[I.p_supp], q = syn ? r[I.fdr_synth] : r[I.fdr_supp];
        return { x: r[I.normZ], y: -Math.log10(Math.max(p, 1e-300)), k: q <= f ? (syn ? "neg" : "pos") : "ns", label: r[0], r }; });
      const top = [...pts].filter((p) => p.k !== "ns").sort((a, b) => b.y - a.y).slice(0, 10);
      C.scatter($("#volcano"), pts, { colors: col, xl: "drugZ normZ", yl: "−log₁₀ p", vline: 0, labels: top, selected: selGene, z: (p) => (p.k === "ns" ? 0 : 1),
        tip: (p) => `<b>${esc(p.label)}</b><br>normZ ${fmt(p.x)}<br>synthetic: p ${fmt(p.r[I.p_synth])}, FDR ${fmt(p.r[I.fdr_synth])}<br>suppressor: p ${fmt(p.r[I.p_supp])}, FDR ${fmt(p.r[I.fdr_supp])}`, click: (p) => showGene(p.label) });
      $("#volcano").insertAdjacentHTML("beforeend", `<div class="legend"><span><i style="background:${col.neg}"></i>synthetic sick/lethal</span><span><i style="background:${col.pos}"></i>suppressor</span><span><i style="background:${col.ns}"></i>not significant</span></div>`);
      $("#volcano").closest(".card").querySelector("h3").textContent = "drugZ volcano";
      $("#volcano").closest(".card").querySelector(".row .faint").innerHTML = `normZ vs −log₁₀ p · FDR ≤ <input type="text" id="fdrCut" value="${f}" size="4">`;
      $("#fdrCut").onchange = () => renderResults();
      $("#rankplot").closest(".card").querySelector("h3").textContent = "Gene rank";
      $("#rankDir").classList.add("hidden");
      const rk = R.rows.map((r, i) => ({ x: i + 1, y: r[I.normZ], k: r[I.fdr_synth] <= f ? "neg" : r[I.fdr_supp] <= f ? "pos" : "ns", label: r[0], r }));
      C.scatter($("#rankplot"), rk, { colors: col, xl: "gene rank", yl: "normZ", height: 340, labels: rk.slice(0, 6), selected: selGene,
        tip: (p) => `<b>${esc(p.label)}</b><br>normZ ${fmt(p.y)}`, click: (p) => showGene(p.label) });
      $("#resNotes").innerHTML = `<b>${esc(R.name)}</b> · drugZ (${R.paired ? "paired" : "unpaired"}) · treated ${R.comparison.treatment.map(esc).join(", ")} vs untreated ${R.comparison.control.map(esc).join(", ")} · run ${esc(R.when)}`;
      sortCol = "normZ"; sortAsc = true;
    }
    drawTable();
  }
  const tile = ([k, v, s]) => `<div class="card tile"><div class="k">${k}</div><div class="v">${esc(v)}</div><div class="s">${esc(s)}</div></div>`;

  const _drawTable = drawTable;
  drawTable = function () {
    if (R.method === "mageck") { $("#rankDir").classList.remove("hidden"); return _drawTable(); }
    const I = R.idx, q = $("#gSearch").value.trim().toUpperCase();
    let rows = R.rows;
    if (q) rows = rows.filter((r) => String(r[0]).toUpperCase().includes(q));
    const ci = I[sortCol] ?? 0;
    rows = [...rows].sort((a, b) => (a[ci] < b[ci] ? -1 : a[ci] > b[ci] ? 1 : 0) * (sortAsc ? 1 : -1));
    const per = 100, pages = Math.max(1, Math.ceil(rows.length / per));
    gpage = Math.min(gpage, pages - 1);
    const show = rows.slice(gpage * per, gpage * per + per);
    const f = cut();
    $("#geneTable").innerHTML = `<tr>${R.columns.map((k) => `<th class="${k === "gene" || k === "reference" ? "" : "num"} ${k === sortCol ? "sorted" + (sortAsc ? " asc" : "") : ""}" data-sort="${k}">${esc(k.replace("_", " "))}</th>`).join("")}</tr>` +
      show.map((r) => `<tr class="click ${r[0] === selGene ? "sel" : ""}" data-g="${esc(r[0])}">${r.map((v, i) => {
        const k = R.columns[i];
        if (i === 0) return `<td><b>${esc(v)}</b></td>`;
        const hit = (k === "bf" && v > (R.cutoff || 5)) || ((k === "fdr_synth" || k === "fdr_supp") && v <= f);
        return `<td class="${typeof v === "number" ? "num" : ""}" ${hit ? `style="color:${k === "fdr_supp" ? "var(--pos)" : "var(--neg)"};font-weight:600"` : ""}>${typeof v === "number" ? fmt(v) : esc(v ?? "")}</td>`;
      }).join("")}</tr>`).join("");
    $("#gPager").innerHTML = `${rows.length.toLocaleString()} genes · page ${gpage + 1} of ${pages} <span class="spacer"></span><button class="btn sm" id="gPrev" ${gpage ? "" : "disabled"}>‹</button><button class="btn sm" id="gNext" ${gpage < pages - 1 ? "" : "disabled"}>›</button>`;
    $("#gPrev").onclick = () => { gpage--; drawTable(); };
    $("#gNext").onclick = () => { gpage++; drawTable(); };
    $$("#geneTable th").forEach((th) => (th.onclick = () => { const k = th.dataset.sort; if (k === sortCol) sortAsc = !sortAsc; else { sortCol = k; sortAsc = !(k === "bf" || k === "sumZ"); } drawTable(); }));
    $$("#geneTable tr.click").forEach((tr) => (tr.onclick = () => showGene(tr.dataset.g)));
  };

  const _showGene = showGene;
  showGene = async function (g) {
    if (R.method === "mageck") return _showGene(g);
    selGene = g;
    drawTable(); renderOtherPlotsOnly();
    try {
      const d = await api(`/api/gene?name=${encodeURIComponent(R.name)}&gene=${encodeURIComponent(g)}`);
      const s = d.summary || {};
      const head = d.method === "bagel" ? `Bayes factor ${fmt(s.bf)}` : `normZ ${fmt(s.normZ)} · synthetic FDR ${fmt(s.fdr_synth)} · suppressor FDR ${fmt(s.fdr_supp)}`;
      const valCols = d.method === "bagel" ? d.replicates.map((x) => "FC " + x).concat(d.replicates.map((x) => "BF " + x)) : d.replicates.map((x) => "FC " + x).concat(d.replicates.map((x) => "z " + x));
      $("#geneView").innerHTML = `<div class="gv-head"><h2>${esc(d.gene)}</h2><span class="muted small">${d.guides.length} guides · ${head}</span></div>
        <div class="plot" id="gvBars"></div>
        <h3 class="mt">Guides</h3>
        <div class="scroll"><table class="tbl"><tr><th>Guide</th>${d.samples.map((x, i) => `<th class="num" style="color:${i < d.n_control ? "var(--neg)" : "var(--pos)"}">${esc(x)}</th>`).join("")}${valCols.map((c) => `<th class="num">${esc(c)}</th>`).join("")}</tr>
        ${d.guides.map((r) => `<tr><td class="mono small">${esc(r.guide)}</td>${r.raw.map((x) => `<td class="num">${fmt(x)}</td>`).join("")}${r.fc.concat(d.method === "bagel" ? r.bf : r.z).map((x) => `<td class="num">${fmt(x)}</td>`).join("")}</tr>`).join("")}</table></div>
        <p class="faint small">Raw counts: blue = ${d.method === "drugz" ? "untreated" : "control"}, red = ${d.method === "drugz" ? "treated" : "treatment"}. ${d.method === "bagel" ? "FC: log₂ fold change against the control; BF: the guide's Bayes factor in each replicate." : "FC: log₂ fold change treated vs untreated; z: the fold change divided by the smoothed standard deviation."}</p>`;
      C.guideBars($("#gvBars"), d.guides.map((r) => ({ guide: r.guide, lfc: r.fc.reduce((a, b) => a + b, 0) / r.fc.length, valid: true, p_low: NaN, p_high: NaN, control_mean: NaN, treat_mean: NaN })));
    } catch (e) { toast(e.message, true); }
  };
  function renderOtherPlotsOnly() { const keep = $("#geneTable").innerHTML; renderOther(); }

  /* ---------------- pathways (pooled and arrayed) */
  function pathways(sel, url, pooled) {
    const box = $(sel);
    if (!box || box.dataset.ready === (pooled ? R && R.name : "arrayed")) return;
    box.dataset.ready = pooled ? R && R.name : "arrayed";
    const coll = (S && S.enrichment) || { GO_BP: "GO biological process", GO_MF: "GO molecular function", GO_CC: "GO cellular component", Reactome: "Reactome pathways" };
    const dirs = !pooled ? [] : R.method === "bagel" ? [["neg", "essential genes (BF > 5)"]] : R.method === "drugz" ? [["neg", "synthetic sick/lethal"], ["pos", "suppressors"]] : [["neg", "depleted"], ["pos", "enriched"]];
    box.innerHTML = `<div class="row"><h3>Pathways and gene sets</h3><span class="spacer"></span><span class="faint small">offline · GO and Reactome (CC BY 4.0) or your own .gmt</span></div>
      <div class="form">
        <label>Gene sets <select class="pwColl">${Object.entries(coll).map(([k, v]) => `<option value="${k}" ${k === "Reactome" ? "selected" : ""}>${esc(v)}</option>`).join("")}<option value="custom">My .gmt file…</option></select></label>
        <label>Test <select class="pwTest"><option value="ora">over-representation of hits</option><option value="ranked">ranked (no cut-off)</option></select></label>
        ${dirs.length ? `<label>Genes <select class="pwDir">${dirs.map(([k, v]) => `<option value="${k}">${esc(v)}</option>`).join("")}</select></label>` : ""}
        <label>Species <select class="pwSp"><option value="auto">from gene names</option><option value="human">human</option><option value="mouse">mouse</option></select></label>
        <input type="file" class="pwFile" accept=".gmt,.txt" hidden>
        <button class="btn primary sm pwRun">Test</button><span class="muted small pwMsg"></span>
      </div><div class="scroll tall mt"><table class="tbl pwTable"></table></div>`;
    let gmt = null, gmtName = "";
    const q = (c) => box.querySelector(c);
    q(".pwColl").onchange = () => { if (q(".pwColl").value === "custom") q(".pwFile").click(); };
    q(".pwFile").onchange = async (e) => { const f = e.target.files[0]; if (!f) return; gmt = await f.text(); gmtName = f.name; q(".pwMsg").textContent = `${f.name} loaded`; };
    q(".pwRun").onclick = async () => {
      q(".pwRun").disabled = true; q(".pwMsg").textContent = "Testing…";
      try {
        const body = { name: pooled ? R.name : "", collection: q(".pwColl").value, test: q(".pwTest").value, direction: q(".pwDir") ? q(".pwDir").value : "neg",
          species: q(".pwSp").value, gmt, gmt_name: gmtName, cutoff: pooled ? cut() : 0.05 };
        const r = await api(url, body);
        q(".pwMsg").textContent = `${r.collection}: ${r.n_sets.toLocaleString()} sets tested · ${r.test === "ora" ? `${r.n_hits.toLocaleString()} hits among ${r.n_universe.toLocaleString()} genes` : `${r.n_universe.toLocaleString()} ranked genes`}`;
        const ora = r.test === "ora";
        q(".pwTable").innerHTML = `<tr><th>Gene set</th><th class="num">Size</th>${ora ? `<th class="num">Hits</th><th class="num">Expected</th><th class="num">Fold</th>` : `<th class="num">AUC</th><th>Direction</th>`}<th class="num">p</th><th class="num">FDR</th><th>Genes</th></tr>` +
          r.rows.map((x) => `<tr><td>${esc(x.set)}</td><td class="num">${x.size}</td>${ora ? `<td class="num">${x.hits}</td><td class="num">${fmt(x.expected)}</td><td class="num">${fmt(x.ratio)}</td>` : `<td class="num">${fmt(x.auc)}</td><td>${x.direction === "better" ? "toward hits" : "away from hits"}</td>`}<td class="num">${fmt(x.p)}</td><td class="num" ${x.fdr <= 0.05 ? 'style="font-weight:600;color:var(--accent)"' : ""}>${fmt(x.fdr)}</td><td class="genes">${x.genes.slice(0, 12).map(esc).join(", ")}${x.genes.length > 12 ? "…" : ""}</td></tr>`).join("");
      } catch (e) { q(".pwMsg").textContent = ""; toast(e.message, true); }
      q(".pwRun").disabled = false;
    };
  }
  window.pathways = pathways;
  // arrayed results get the same panel
  const _aRender = window.A.render;
  window.A.render = function (p) {
    const out = _aRender(p);
    if (p === "results") setTimeout(() => { const b = $("#apwCard"); if (b && !$("#arBody").classList.contains("hidden")) { b.dataset.ready = ""; pathways("#apwCard", "/api/a/enrich", false); } }, 400);
    return out;
  };
})();
