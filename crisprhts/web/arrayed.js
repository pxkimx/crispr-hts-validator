/* CRISPR HTS Validator — arrayed (plate) screens. Uses helpers from app.js and charts.js. */
window.A = (() => {
  let AS = null;          // arrayed state
  let AQ = null;          // plate QC view
  let AR = null;          // result view
  let mapPlate = null, qcPlate = null, resPlate = null;

  async function refresh() {
    AS = await api("/api/a/state");
    readout();
    renderData();
  }
  function readout() {
    if (!AS) return;
    const s = AS.screen;
    $("#readout").innerHTML = s ? `<span>${s.plates}</span> plates · <span>${s.wells}</span>-well · <span>${s.replicates.length}</span> rep.` : "no screen";
    $("#qcDot").className = "dot " + (AS.qc_overall || "");
  }

  /* ---------------- data */
  function renderData() {
    $("#aExamples").innerHTML = Object.entries(AS.examples).map(([k, v]) => `<div class="card"><div class="kicker">public data</div><h3>${esc(v.title)}</h3><p class="muted small">${esc(v.about)}</p><button class="btn sm" data-aex="${k}" ${AS.examples_available ? "" : "disabled"}>Load</button></div>`).join("");
    $$("[data-aex]").forEach((b) => (b.onclick = async () => {
      b.disabled = true;
      try { await api("/api/a/example", { key: b.dataset.aex }); AQ = AR = null; await refresh(); toast("Example loaded."); window.scrollTo({ top: $("#aLoaded").offsetTop - 70, behavior: "smooth" }); }
      catch (e) { toast(e.message, true); }
      b.disabled = false;
    }));
    $("#aFiles").textContent = AS.files.length ? `${AS.files.length} file(s) ready` : "";
    $("#aLoadBtn").disabled = !AS.files.length || !!AS.screen;
    const s = AS.screen;
    $("#aLoaded").classList.toggle("hidden", !s);
    if (!s) return;
    const ctl = Object.entries(s.contents).filter(([k]) => k !== "sample").map(([k, v]) => `${v} ${k}`).join(", ");
    $("#aTiles").innerHTML = [["Plates", s.plates, `${s.wells}-well format`], ["Replicates", s.replicates.length, s.channels.length > 1 ? `${s.channels.length} channels` : "1 channel"],
      ["Sample wells", (s.contents.sample || 0).toLocaleString(), s.genes ? `${s.genes.toLocaleString()} genes annotated` : "no gene annotation"],
      ["Control / other wells", Object.entries(s.contents).filter(([k]) => k !== "sample").reduce((a, [, v]) => a + v, 0), ctl || "none marked"]]
      .map(([k, v, t]) => `<div class="card tile"><div class="k">${k}</div><div class="v">${esc(v)}</div><div class="s">${esc(t)}</div></div>`).join("");
    $("#aName").textContent = AS.name;
    $("#aNotes").innerHTML = s.notes.length || s.flagged ? `<div class="notes">${s.notes.map(esc).join("<br>")}</div>` : "";
    const d = s.description || {};
    $("#aDesc").innerHTML = ["Title", "Organism", "Celltype", "Assay", "Reference"].filter((k) => d[k]).map((k) => `<b>${k}</b>: ${esc(d[k])}`).join(" · ");
  }
  async function uploadFiles(files) {
    for (const f of files) {
      showJob({ message: `Copying ${f.name}…`, progress: 0.1 });
      try {
        await new Promise((res, rej) => {
          const x = new XMLHttpRequest();
          x.open("POST", `/api/a/upload?name=${encodeURIComponent(f.name)}`);
          x.onload = () => { let j = {}; try { j = JSON.parse(x.responseText); } catch (e) {} x.status === 200 && !j.error ? res(j) : rej(new Error(j.error || "Upload failed.")); };
          x.onerror = () => rej(new Error("Upload failed."));
          x.send(f);
        });
      } catch (e) { toast(`${f.name}: ${e.message}`, true); }
    }
    hideJob();
    await refresh();
    await loadScreen();
  }
  async function loadScreen() {
    try { await api("/api/a/load", {}); AQ = AR = null; await refresh(); toast("Screen read."); }
    catch (e) { toast(e.message, true); await refresh(); }
  }
  $("#aLoadBtn").onclick = loadScreen;

  /* ---------------- plate map */
  async function renderMap() {
    const has = AS && AS.screen;
    $("#amapEmpty").classList.toggle("hidden", !!has);
    $("#amapBody").classList.toggle("hidden", !has);
    if (!has) return;
    const n = AS.screen.plates;
    if (!mapPlate) mapPlate = 1;
    const plates = await plateNumbers();
    $("#amPlate").innerHTML = plates.map((p) => `<option ${p === mapPlate ? "selected" : ""}>${p}</option>`).join("");
    $("#amPlate").onchange = () => { mapPlate = +$("#amPlate").value; drawLayout(); };
    drawLayout();
    $("#amRoles").innerHTML = `<tr><th>Well type</th><th class="num">Wells per screen</th><th>Role</th></tr>` + Object.entries(AS.roles).map(([k, v]) =>
      `<tr><td>${esc(k)}</td><td class="num">${(AS.screen.contents[k] || 0).toLocaleString()}</td><td><select data-role="${esc(k)}">${["sample", "neg", "pos", "other", "empty"].map((r) => `<option value="${r}" ${r === v ? "selected" : ""}>${{ sample: "sample", neg: "negative control", pos: "positive control", other: "other (not scored)", empty: "empty" }[r]}</option>`).join("")}</select></td></tr>`).join("");
    $$(`input[name=adir]`).forEach((i) => (i.checked = i.value === AS.direction));
    $("#amChannel").innerHTML = AS.screen.channels.map((c) => `<option ${c === AS.channel ? "selected" : ""}>${c}</option>`).join("");
    $("#amChanL").classList.toggle("hidden", AS.screen.channels.length < 2);
  }
  let _plates = null;
  async function plateNumbers() {
    if (_plates && _plates.n === AS.screen.plates && _plates.name === AS.name) return _plates.list;
    const list = [];
    // plates are numbered 1..n in every reader
    for (let i = 1; i <= AS.screen.plates; i++) list.push(i);
    _plates = { n: AS.screen.plates, name: AS.name, list };
    return list;
  }
  async function drawLayout() {
    const L = await api("/api/a/layout?plate=" + mapPlate);
    const [nr, nc] = L.dims;
    const kinds = [...new Set(L.wells.map((w) => w.content))];
    const M = Array.from({ length: nr }, () => Array(nc).fill(null));
    const role = (w) => AS.roles[w.content] || "other";
    L.wells.forEach((w) => { M[w.r][w.c] = role(w) === "empty" ? null : 1; });
    const info = L.wells.map((w) => ({ ...w, role: role(w), content: role(w) === "pos" ? "pos" : role(w) === "neg" ? "neg" : w.content }));
    const cat = { sample: C.css("--line2"), neg: C.css("--neg"), pos: C.css("--pos"), other: "#eda100", empty: "transparent" };
    C.plate($("#amLayout"), M, info, { cat, unit: "" });
    $("#amLayout .legend").remove();
    $("#amLegend").innerHTML = `<span><i style="background:${cat.neg}"></i>negative control</span><span><i style="background:${cat.pos}"></i>positive control</span><span><i style="background:${cat.other}"></i>other (not scored)</span><span><i style="background:${cat.sample}"></i>sample</span><span class="faint">${kinds.length} well type(s) on this plate · hover a well for its reagent</span>`;
  }
  $("#amMark").onclick = async () => {
    try {
      AS = await api("/api/a/map", { mark: [{ wells: $("#amWells").value, plates: $("#amPlates").value, content: $("#amContent").value }] });
      toast("Marked."); renderMap();
    } catch (e) { toast(e.message, true); }
  };
  $("#amSave").onclick = async () => {
    const roles = {};
    $$("[data-role]").forEach((s) => (roles[s.dataset.role] = s.value));
    const dir = ($("input[name=adir]:checked") || {}).value;
    try { AS = await api("/api/a/map", { roles, direction: dir, channel: +$("#amChannel").value || 1 }); AQ = AR = null; toast("Saved."); go("quality"); runQC(); }
    catch (e) { toast(e.message, true); }
  };

  /* ---------------- plate QC */
  $("#aqBtn").onclick = () => runQC();
  async function runQC() {
    if (!AS || !AS.screen) return toast("Load a screen first.", true);
    try {
      const j = await api("/api/a/qc", {});
      watchJob(j.id, async () => { AQ = await api("/api/a/qc" + (qcPlate ? "?plate=" + qcPlate : "")); AS.qc_overall = AQ.overall; readout(); renderQC(); });
    } catch (e) { toast(e.message, true); }
  }
  async function renderQC() {
    if (!AQ && AS && AS.has_qc) { try { AQ = await api("/api/a/qc"); } catch (e) {} }
    $("#aqEmpty").classList.toggle("hidden", !!AQ);
    if (!AQ) { $("#aqCards").innerHTML = ""; $("#aqOverall").innerHTML = ""; return; }
    const n = (s) => AQ.checks.filter((c) => c.status === s).length;
    $("#aqOverall").innerHTML = `<div class="card overall">${statusTag(AQ.overall)}<span class="big">${AQ.overall === "pass" ? "The plates pass every check" : AQ.overall === "warn" ? "Usable, with points to look at" : AQ.overall === "fail" ? "Problems that affect the hit calls" : "Not enough information"}</span><span class="spacer"></span><span class="muted small">${n("pass")} pass · ${n("warn")} warn · ${n("fail")} fail · ${n("na")} not applicable</span></div>`;
    const box = $("#aqCards");
    box.innerHTML = "";
    const rows = AQ.dims[0], cols = AQ.dims[1];
    for (const c of AQ.checks) {
      const el = document.createElement("div");
      el.className = "card qc" + (c.id === "plates" ? " wide" : "");
      el.innerHTML = `<h3>${statusTag(c.status)}${esc(c.name)}</h3><div class="why">${esc(c.why)}</div><div class="rule">${esc(c.rule || "")}</div><div class="plot mt"></div>`;
      box.appendChild(el);
      const P = el.querySelector(".plot");
      if (c.id === "plates") {
        const reps = [...new Set(c.per_plate.map((p) => p.replicate))];
        const plates = [...new Set(c.per_plate.map((p) => p.plate))];
        const mx = Math.max(...c.per_plate.map((p) => p.value));
        const items = plates.map((pl) => ({ label: "plate " + pl, value: Math.max(...c.per_plate.filter((p) => p.plate === pl).map((p) => p.value)), status: "na", tip: c.per_plate.filter((p) => p.plate === pl).map((p) => `replicate ${p.replicate}: ${fmt(p.value)}`).join(" · ") }));
        P.innerHTML = `<div class="plot" id="aqPlates"></div>`;
        stripChart($("#aqPlates"), c.per_plate, reps);
      } else if (c.id === "controls") {
        if (c.per_plate && c.per_plate.length) C.bars(P, c.per_plate.map((p) => ({ label: `plate ${p.plate} r${p.replicate}`, value: p.zprime, status: p.status, tip: `Z′ ${fmt(p.zprime)} · SSMD ${fmt(p.ssmd)}` })), { min: Math.min(0, ...c.per_plate.map((p) => p.zprime)), max: 1, lines: [{ v: 0.5, label: "0.5" }] });
        else if (c.pooled) P.innerHTML = `<div class="kv"><span>Z′ (pooled)</span><span>${fmt(c.pooled.zprime)}</span><span>SSMD (pooled)</span><span>${fmt(c.pooled.ssmd)}</span><span>Positive controls</span><span>${c.pooled.n_pos} wells, ${(c.pooled.pos_mean * 100).toFixed(0)}% of plate median</span><span>Negative controls</span><span>${c.pooled.n_neg} wells, ${(c.pooled.neg_mean * 100).toFixed(0)}% of plate median</span></div>`;
      } else if (c.id === "edge") {
        const v = c.per_plate.map((x, i) => ({ label: "", value: x * 100 }));
        P.innerHTML = `<div id="aqRow"></div><div id="aqCol"></div>`;
        C.profile($("#aqRow"), AQ.row_profile, AQ.row_profile.map((_, i) => String.fromCharCode(65 + i)), { xl: "plate row", what: "row" });
        C.profile($("#aqCol"), AQ.col_profile, AQ.col_profile.map((_, i) => String(i + 1)), { xl: "plate column", what: "column" });
      } else if (c.id === "rowcol") {
        const plates = AQ.plates, rep = AQ.replicates[0];
        P.innerHTML = `<div class="row small"><label>Plate <select id="aqPlate">${plates.map((p) => `<option ${p === AQ.plate ? "selected" : ""}>${p}</option>`).join("")}</select></label><span class="muted">raw values, replicate ${AQ.replicates.map((r) => r).join(" and ")}</span></div><div id="aqMaps"></div>`;
        const maps = $("#aqMaps");
        Object.entries(AQ.maps).forEach(([k, M]) => { const d = document.createElement("div"); d.className = "plot mt"; maps.appendChild(d); C.plate(d, M, [], { mode: "seq", unit: `raw signal, replicate ${k.split(":")[0]}` }); });
        $("#aqPlate").onchange = async () => { qcPlate = +$("#aqPlate").value; AQ = await api("/api/a/qc?plate=" + qcPlate); renderQC(); };
      } else if (c.id === "replicates" && AQ.rep_scatter) {
        const s = AQ.rep_scatter;
        C.scatter(P, s.x.map((x, i) => ({ x, y: s.y[i], k: "ns", label: s.labels[i] })), { colors: { ns: C.css("--accent") }, xl: "replicate 1 (÷ plate median)", yl: "replicate 2", height: 300,
          tip: (p) => `<b>${esc(p.label)}</b><br>replicate 1 ${fmt(p.x)} · replicate 2 ${fmt(p.y)}` });
      }
    }
  }
  function stripChart(el, per, reps) {
    const plates = [...new Set(per.map((p) => p.plate))];
    const pts = per.map((p) => ({ x: p.plate + (p.replicate - 1.5) * 0.25, y: p.value, k: p.replicate === reps[0] ? "neg" : "pos", label: `plate ${p.plate}, replicate ${p.replicate}` }));
    C.scatter(el, pts, { colors: { neg: C.css("--accent"), pos: C.css("--accent2"), ns: C.css("--ns") }, xl: "plate", yl: "median raw signal of sample wells", height: 240, yr: [0, Math.max(...per.map((p) => p.value))],
      tip: (p) => `<b>${esc(p.label)}</b><br>median ${fmt(p.y)}` });
    el.insertAdjacentHTML("beforeend", `<div class="legend">${reps.map((r, i) => `<span><i style="background:${i ? C.css("--accent2") : C.css("--accent")}"></i>replicate ${r}</span>`).join("")}</div>`);
  }

  /* ---------------- analyse */
  function renderAnalyse() {
    if (!AS) return;
    $("#aNorm").innerHTML = Object.entries(AS.norm_methods).map(([k, v]) => `<label><input type="radio" name="anorm" value="${k}" ${k === (AS.default_norm || "robust_z") ? "checked" : ""}> <b>${esc(v.split(":")[0])}</b><span class="desc">${esc(v.split(":").slice(1).join(":").trim())}</span></label>`).join("");
    $("#aDirNote").textContent = AS.screen ? `Hits are wells with ${AS.direction === "lower" ? "lower" : AS.direction === "higher" ? "higher" : "lower or higher"} signal than the plate's samples (set on the Plate map page). Scores are oriented so that a hit scores high.` : "Load a screen first.";
  }
  $("#aRunBtn").onclick = async () => {
    const o = { method: ($("input[name=anorm]:checked") || {}).value || AS.default_norm || "robust_z", log2: $("#aLog").checked, use_flagged: $("#aFlagged").checked,
      score: $("#aScore").value, summary: $("#aSumm").value, threshold: +$("#aThr").value, replicate_min: +$("#aRepMin").value, gene_fraction: +$("#aGeneFr").value };
    try { const j = await api("/api/a/analyse", o); watchJob(j.id, async () => { AR = null; resPlate = null; AS = await api("/api/a/state"); toast("Analysis finished."); go("results"); }); }
    catch (e) { toast(e.message, true); }
  };

  /* ---------------- results */
  async function renderResults() {
    const has = AS && AS.has_result;
    $("#arEmpty").classList.toggle("hidden", !!has);
    $("#arBody").classList.toggle("hidden", !has);
    if (!has) return;
    if (!AR || (resPlate && AR.plate !== resPlate)) AR = await api("/api/a/result" + (resPlate ? "?plate=" + resPlate : ""));
    const st = AR.settings;
    $("#arTiles").innerHTML = [["Hit wells", AR.n_hits.toLocaleString(), `of ${AR.n_samples.toLocaleString()} sample wells (${(100 * AR.n_hits / Math.max(1, AR.n_samples)).toFixed(1)}%)`],
      ["Hit genes", AR.n_gene_hits.toLocaleString(), AR.n_genes ? `of ${AR.n_genes.toLocaleString()} genes` : "no gene annotation"],
      ["Normalisation", esc(AS.norm_methods[st.method].split(":")[0]), `${st.log2 ? "log₂, " : ""}${st.score === "ssmd" ? "robust SSMD" : "z-score"}, ${st.direction} signal`],
      ["Hit rule", `≥ ${st.threshold}`, `${st.summary} of replicates, every replicate ≥ ${st.replicate_min}`]].map(([k, v, s]) => `<div class="card tile"><div class="k">${k}</div><div class="v">${v}</div><div class="s">${esc(s)}</div></div>`).join("");
    const cmp = AR.comparison;
    $("#arCompare").classList.toggle("hidden", !cmp);
    if (cmp) $("#arCompare").innerHTML = `<h3>Compared with the published 2003 analysis</h3><p class="small">${cmp.wells.toLocaleString()} wells scored in <code>Analysis2003.txt</code>. Correlation of scores r = <b>${cmp.r.toFixed(4)}</b>${cmp.exact_2dp > 0.5 ? `; ${(cmp.exact_2dp * 100).toFixed(2)}% identical to two decimals` : ""}. Published hits (score ≥ 2): <b>${cmp.published_hits}</b>; this run: <b>${cmp.our_hits}</b>; in both: <b>${cmp.both}</b> (recall ${(cmp.recall * 100).toFixed(1)}%, precision ${(cmp.precision * 100).toFixed(1)}%).</p><p class="small muted">The 2003 analysis used a classic plate z-score (plate mean and SD over all 384 wells), averaged the two replicates, did not require replicate agreement and kept the wells flagged in the screen log. Choose “Classic plate z-score”, keep flagged wells and set the replicate threshold very low to reproduce it exactly.</p>`;
    $("#arPlate").innerHTML = AR.plates.map((p) => `<option ${p === AR.plate ? "selected" : ""}>${p}</option>`).join("");
    $("#arPlate").onchange = async () => { resPlate = +$("#arPlate").value; AR = await api("/api/a/result?plate=" + resPlate); renderResults(); };
    C.plate($("#arMap"), AR.map, AR.map_info, { mode: "div", lim: 4, unit: st.score === "ssmd" ? "SSMD" : "score (z)" });
    if (AR.scatter) {
      const s = AR.scatter;
      C.scatter($("#arScatter"), s.x.map((x, i) => ({ x, y: s.y[i], k: s.hit[i] ? "pos" : "ns", label: s.labels[i] })), { colors: { pos: C.css("--pos"), ns: C.css("--ns") }, xl: "score, replicate 1", yl: "score, replicate 2", height: 320, z: (p) => (p.k === "ns" ? 0 : 1),
        tip: (p) => `<b>${esc(p.label)}</b><br>replicate 1 ${fmt(p.x)} · replicate 2 ${fmt(p.y)}${p.k === "pos" ? "<br><b>hit</b>" : ""}` });
      $("#arScatter").insertAdjacentHTML("beforeend", `<div class="legend"><span><i style="background:${C.css("--pos")}"></i>hit</span><span><i style="background:${C.css("--ns")}"></i>not a hit</span></div>`);
    } else $("#arScatter").innerHTML = `<div class="empty">One replicate only.</div>`;
    C.histBars($("#arHist"), AR.hist.edges, AR.hist.counts, { xl: "combined score", log: true, lines: [{ v: st.threshold, label: "hit ≥ " + st.threshold, dir: "hi" }] });
    C.histBars($("#arPerPlate"), [AR.plates[0] - 0.5, ...AR.plates.map((p) => p + 0.5)], AR.hits_per_plate, { xl: "plate", unit: "hit wells" });
    drawHitTable();
    const G = AR.genes, gc = AR.genes_columns;
    $("#arGenes").innerHTML = G.length ? `<tr>${gc.map((c) => `<th class="${c === "gene" ? "" : "num"}">${esc(c.replace("_", " "))}</th>`).join("")}</tr>` + G.slice(0, 500).map((r) => `<tr>${r.map((v, i) => `<td class="${i ? "num" : ""}" ${gc[i] === "hit" && v ? 'style="color:var(--pos);font-weight:600"' : ""}>${typeof v === "number" ? fmt(v) : v === true ? "hit" : v === false ? "" : esc(v ?? "")}</td>`).join("")}</tr>`).join("") : `<div class="empty">No gene annotation in this screen.</div>`;
  }
  function drawHitTable() {
    const q = $("#arSearch").value.trim().toUpperCase(), cols = AR.hits_columns;
    const rows = AR.hits.filter((r) => !q || String(r[3]).toUpperCase().includes(q) || String(r[4]).toUpperCase().includes(q));
    $("#arHits").innerHTML = `<tr>${cols.map((c) => `<th class="${/score|raw|plate/.test(c) ? "num" : ""}">${esc(c.replace("_", " "))}</th>`).join("")}</tr>` +
      rows.slice(0, 500).map((r) => `<tr class="click" data-pl="${r[0]}">${r.map((v, i) => `<td class="${/score|raw|plate/.test(cols[i]) ? "num" : ""}">${typeof v === "number" ? fmt(v) : esc(v ?? "")}</td>`).join("")}</tr>`).join("") +
      (rows.length > 500 ? `<tr><td colspan="${cols.length}" class="faint">… ${rows.length - 500} more in the export</td></tr>` : "");
    $$("#arHits tr.click").forEach((tr) => (tr.onclick = async () => { resPlate = +tr.dataset.pl; AR = await api("/api/a/result?plate=" + resPlate); renderResults(); window.scrollTo({ top: $("#arMap").getBoundingClientRect().top + scrollY - 120, behavior: "smooth" }); }));
  }
  $("#arSearch").oninput = () => AR && drawHitTable();
  $("#aExpBtn").onclick = (e) => { e.stopPropagation(); $("#aExpMenu").classList.toggle("on"); };
  document.addEventListener("click", (e) => { if (!e.target.closest(".menu")) $("#aExpMenu").classList.remove("on"); });
  $$("[data-aexp]").forEach((a) => (a.onclick = () => {
    const l = document.createElement("a"); l.href = `/api/a/export?what=${a.dataset.aexp}&fmt=${a.dataset.fmt}`; l.download = ""; document.body.appendChild(l); l.click(); l.remove();
    $("#aExpMenu").classList.remove("on");
  }));

  function render(p) {
    if (!AS) return refresh().then(() => render(p));
    if (p === "samples") renderMap();
    if (p === "quality") renderQC();
    if (p === "analyse") renderAnalyse();
    if (p === "results") renderResults();
  }
  return { refresh, readout, render, uploadFiles };
})();
