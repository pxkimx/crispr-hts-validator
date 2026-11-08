/* CRISPR HTS Validator — page logic. */
const $ = (s, r = document) => r.querySelector(s);
const $$ = (s, r = document) => [...r.querySelectorAll(s)];
const { fmt, esc } = C;
let S = null;            // server state summary
let scr = {};            // screen design numbers (cells per guide, MOI, selection)
let QC = null;           // last quality-check result
let R = null;            // current result
let design = null;       // editable copy of the design
let page = "data";
let MODE = (() => { try { return localStorage.getItem("ch-mode") === "arrayed" ? "arrayed" : "pooled"; } catch (e) { return "pooled"; } })();
const SECTION = { pooled: { data: "p-data", samples: "p-samples", quality: "p-quality", analyse: "p-analyse", results: "p-results", methods: "p-methods" },
                  arrayed: { data: "p-data", samples: "p-amap", quality: "p-aquality", analyse: "p-aanalyse", results: "p-aresults", methods: "p-methods" } };
const TAB_LABELS = { pooled: { samples: "Samples", results: "Results" }, arrayed: { samples: "Plate map", results: "Hits" } };
function setMode(m) {
  MODE = m;
  try { localStorage.setItem("ch-mode", m); } catch (e) {}
  $$("#modePill button").forEach((b) => b.classList.toggle("on", b.dataset.m === m));
  $("#pooledData").classList.toggle("hidden", m !== "pooled");
  $("#arrayedData").classList.toggle("hidden", m !== "arrayed");
  $$("#tabs button").forEach((b) => { const l = TAB_LABELS[m][b.dataset.p]; if (l) b.querySelector(".lbl").textContent = l; });
  $("#heroKicker").textContent = m === "arrayed" ? "Arrayed screens · plate reads to hit calls" : "Pooled CRISPR screens · read counts to hit lists";
  $("#heroLede").textContent = m === "arrayed"
    ? "Load plate-reader data from an arrayed CRISPR, RNAi or compound screen. The validator checks plates, controls, edge and row/column effects and replicates, normalises each plate (including B-score), and calls hits that both replicates support. Everything stays on this Mac."
    : "Load a count table or raw FASTQ files. The validator checks depth, evenness, replicates and essential-gene drop-out, normalises, and scores guides and genes the way MAGeCK does — with each result traceable to the guides behind it. Everything stays on this Mac.";
  if (m === "arrayed" && window.A) window.A.refresh(); else if (S) renderReadout();
  if (page !== "data") go(page);
}
$$("#modePill button").forEach((b) => (b.onclick = () => setMode(b.dataset.m)));

/* ---------------- quit when the window closes (server holds this stream) */
try { new EventSource("/api/window"); } catch (e) {}

/* ---------------- theme */
window.setTheme = function (t) {
  document.documentElement.setAttribute("data-theme", t);
  try { localStorage.setItem("ch-theme", t); } catch (e) {}
  $$("#themePill button").forEach((b) => b.classList.toggle("on", b.dataset.t === t));
  redraw();
};
$$("#themePill button").forEach((b) => (b.onclick = () => setTheme(b.dataset.t)));
$$("#themePill button").forEach((b) => b.classList.toggle("on", b.dataset.t === document.documentElement.getAttribute("data-theme")));

/* ---------------- helpers */
async function api(path, body, method) {
  const opt = body === undefined ? { method: method || "GET" } : { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) };
  const r = await fetch(path, opt);
  const j = await r.json().catch(() => ({ error: "The app did not answer." }));
  if (!r.ok || j.error) throw new Error(j.error || r.statusText);
  return j;
}
function toast(msg, err) {
  const t = $("#toast");
  t.textContent = msg;
  t.className = "toast on" + (err ? " err" : "");
  clearTimeout(t._h);
  t._h = setTimeout(() => (t.className = "toast"), err ? 6000 : 2800);
}
const statusTag = (s) => `<span class="status ${s}">${{ pass: "pass", warn: "warn", fail: "fail", na: "n/a" }[s] || s}</span>`;

/* ---------------- navigation */
function go(p) {
  page = p;
  $$("#tabs button").forEach((b) => b.classList.toggle("on", b.dataset.p === p));
  const sid = SECTION[MODE][p];
  $$("section.page").forEach((s) => s.classList.toggle("on", s.id === sid));
  window.scrollTo({ top: 0 });
  if (MODE === "arrayed" && p !== "data" && p !== "methods") { if (window.A) window.A.render(p); return; }
  if (p === "samples") renderSamples();
  if (p === "quality") renderQC();
  if (p === "analyse") renderAnalyse();
  if (p === "results") loadResults();
  if (p === "methods") renderMethods();
}
$$("#tabs button").forEach((b) => (b.onclick = () => go(b.dataset.p)));
document.addEventListener("click", (e) => {
  const g = e.target.closest("[data-go]");
  if (g) go(g.dataset.go);
  if (!e.target.closest(".menu")) $("#expMenu").classList.remove("on");
});

/* ---------------- state */
async function refresh() {
  S = await api("/api/state");
  $("#ver").textContent = "v" + S.version;
  renderReadout();
  design = JSON.parse(JSON.stringify(S.design));
  renderData();
  if (S.jobs.length) watchJob(S.jobs[0].id);
}

function renderReadout() {
  if (MODE === "arrayed") { if (window.A) window.A.readout(); return; }
  const g = S.table ? S.table.guides.toLocaleString() : (S.library ? S.library.guides.toLocaleString() : "—");
  const n = S.table ? S.table.samples.length : (S.fastq.length || "—");
  $("#readout").innerHTML = "<span>" + g + "</span> guides · <span>" + n + "</span> samples";
  $("#qcDot").className = "dot " + (S.qc_overall || "");
}

/* ---------------- DATA page */
function renderData() {
  const ex = $("#examples");
  ex.innerHTML = Object.entries(S.examples).map(([k, v]) => `<div class="card"><div class="kicker">public data</div><h3>${esc(v.title)}</h3><p class="muted small">${esc(v.about)}</p><button class="btn sm" data-ex="${k}" ${S.examples_available ? "" : "disabled"}>Load</button></div>`).join("");
  $$("[data-ex]").forEach((b) => (b.onclick = () => loadExample(b.dataset.ex, b)));
  // FASTQ panel
  const fq = S.fastq.length || S.library;
  $("#fastqPanel").classList.toggle("hidden", !fq);
  if (fq) {
    $("#libInfo").innerHTML = S.library ? `Library: <b>${S.library.guides.toLocaleString()}</b> guides, ${S.library.genes.toLocaleString()} genes${S.library.notes.length ? " · " + esc(S.library.notes.join(" ")) : ""}` : `<span class="status warn">no library yet</span>`;
    $("#fastqTable").innerHTML = `<tr><th>File</th><th class="num">Size</th><th>Sample name</th><th></th></tr>` + S.fastq.map((f) =>
      `<tr><td class="mono">${esc(f.name)}</td><td class="num">${(f.size / 1e6).toFixed(1)} MB</td><td><input type="text" data-fq="${esc(f.name)}" value="${esc(f.sample)}"></td><td><button class="btn sm ghost" data-rmfq="${esc(f.name)}">Remove</button></td></tr>`).join("");
    $$("[data-rmfq]").forEach((b) => (b.onclick = async () => { await api("/api/fastq/remove", { name: b.dataset.rmfq }); refresh(); }));
  }
  if (S.count_info && S.count_info.length) {
    $("#countInfo").innerHTML = `<h3 class="mt">Counting summary</h3><table class="tbl"><tr><th>Sample</th><th>File</th><th class="num">Reads</th><th class="num">Mapped</th><th class="num">Mapped %</th><th class="num">Zero-count guides</th><th class="num">Gini</th><th>Guide found</th></tr>` +
      S.count_info.map((c) => `<tr><td>${esc(c.sample)}</td><td class="mono">${esc(c.file)}</td><td class="num">${c.reads.toLocaleString()}</td><td class="num">${c.mapped.toLocaleString()}</td><td class="num">${(c.fraction * 100).toFixed(2)}%</td><td class="num">${c.zero.toLocaleString()}</td><td class="num">${c.gini.toFixed(4)}</td><td class="small">${c.mode === "scaffold" ? "by scaffold" : "from position " + c.trims.join(", ")}${c.reverse ? ", reverse strand" : ""}${c.mismatch_reads ? `, ${c.mismatch_reads.toLocaleString()} with 1 mismatch` : ""}</td></tr>`).join("") + "</table>";
  } else $("#countInfo").innerHTML = "";
  // loaded
  const t = S.table;
  $("#loaded").classList.toggle("hidden", !t);
  if (t) {
    const reads = t.totals.reduce((a, b) => a + b, 0);
    $("#loadedTiles").innerHTML = [
      ["Format", t.format, t.source], ["Guides", t.guides.toLocaleString(), ""], ["Genes", t.genes.toLocaleString(), ""],
      ["Samples", t.samples.length, (reads >= 1e6 ? (reads / 1e6).toFixed(1) + " M" : reads.toLocaleString()) + " counted reads"]].map(([k, v, s]) => `<div class="card tile"><div class="k">${k}</div><div class="v">${esc(v)}</div><div class="s">${esc(s)}</div></div>`).join("");
    $("#loadedName").textContent = S.name;
    $("#loadedNotes").innerHTML = t.notes.length ? `<div class="notes">${t.notes.map(esc).join("<br>")}</div>` : "";
    $("#preview").innerHTML = `<tr><th>guide</th><th>gene</th>${t.samples.map((s) => `<th class="num">${esc(s)}</th>`).join("")}</tr>` +
      t.preview.map((r) => `<tr>${r.map((x, i) => `<td class="${i > 1 ? "num" : ""}">${i > 1 ? fmt(x) : esc(x)}</td>`).join("")}</tr>`).join("");
  }
}
async function loadExample(k, btn) {
  btn.disabled = true;
  try { await api("/api/example", { key: k }); QC = null; R = null; await refresh(); toast("Example loaded."); if (S.table) window.scrollTo({ top: $("#loaded").offsetTop - 70, behavior: "smooth" }); else window.scrollTo({ top: $("#fastqPanel").offsetTop - 70, behavior: "smooth" }); }
  catch (e) { toast(e.message, true); }
  btn.disabled = false;
}

/* uploads: browse + drag and drop */
let browseKind = "counts";
$$("[data-browse]").forEach((b) => (b.onclick = () => { browseKind = b.dataset.browse; $("#fileIn").multiple = browseKind === "fastq" || browseKind === "arrayed"; $("#fileIn").click(); }));
$("#fileIn").onchange = (e) => { uploadFiles([...e.target.files], browseKind); e.target.value = ""; };
$$(".drop").forEach((d) => {
  d.ondragover = (e) => { e.preventDefault(); d.classList.add("over"); };
  d.ondragleave = () => d.classList.remove("over");
  d.ondrop = (e) => { e.preventDefault(); d.classList.remove("over"); uploadFiles([...e.dataTransfer.files], d.dataset.kind); };
});
document.addEventListener("dragover", (e) => e.preventDefault());
document.addEventListener("drop", (e) => e.preventDefault());
function guessKind(f, hint) {
  const n = f.name.toLowerCase();
  if (/\.(fastq|fq)(\.gz)?$/.test(n)) return "fastq";
  if (hint === "screeninfo") return "screeninfo";
  if (hint === "positions") return "positions";
  if (hint === "library" || /librar|sgrna_lib|guides?\./.test(n)) return "library";
  if (hint === "fastq") return "library";      // non-FASTQ file dropped with FASTQs = the library
  return "counts";
}
async function uploadFiles(files, hint) {
  if (!files.length) return;
  if (hint === "arrayed") return window.A.uploadFiles(files);
  // library first, then FASTQs
  const items = files.map((f) => ({ f, k: guessKind(f, hint) })).sort((a, b) => (a.k === "library" ? -1 : 0) - (b.k === "library" ? -1 : 0));
  for (const { f, k } of items) {
    showJob({ message: `Reading ${f.name}…`, progress: 0.05 });
    try {
      const r = await uploadOne(f, k);
      if (k === "counts") { QC = null; R = null; toast(`${f.name}: ${r.table.guides.toLocaleString()} guides, ${r.table.samples.length} samples (${r.table.format}).`); }
      else if (k === "screeninfo") toast(`Imported ${r.design.comparisons.length} comparison(s) from ${f.name}.`);
      else toast(`${f.name} added.`);
    } catch (e) { toast(`${f.name}: ${e.message}`, true); }
  }
  hideJob();
  await refresh();
  if (page === "samples") renderSamples();
}
function uploadOne(f, kind) {
  return new Promise((res, rej) => {
    const x = new XMLHttpRequest();
    x.open("POST", `/api/upload?kind=${kind}&name=${encodeURIComponent(f.name)}`);
    x.upload.onprogress = (e) => e.lengthComputable && showJob({ message: `Copying ${f.name}… ${Math.round(100 * e.loaded / e.total)}%`, progress: e.loaded / e.total });
    x.onload = () => { let j = {}; try { j = JSON.parse(x.responseText); } catch (e) {} x.status === 200 && !j.error ? res(j) : rej(new Error(j.error || "Upload failed.")); };
    x.onerror = () => rej(new Error("Upload failed."));
    x.send(f);
  });
}

/* counting */
$("#cMode").onchange = () => { $(".cTrimL").classList.toggle("hidden", $("#cMode").value !== "offset"); $(".cFlankL").classList.toggle("hidden", $("#cMode").value !== "scaffold"); };
$("#countBtn").onclick = async () => {
  const samples = {};
  $$("[data-fq]").forEach((i) => (samples[i.dataset.fq] = i.value));
  try {
    const j = await api("/api/count", { samples, mode: $("#cMode").value, trim: Math.max(0, parseInt($("#cTrim").value || "1") - 1), mismatch: $("#cMismatch").checked, flank: $("#cFlank").value });
    watchJob(j.id, async () => { await refresh(); toast("Counting finished."); });
  } catch (e) { toast(e.message, true); }
};

/* ---------------- jobs */
function showJob(j) {
  $("#jobbar").classList.add("on");
  $("#jobfill").style.width = Math.round((j.progress || 0) * 100) + "%";
  $("#jobmsg").textContent = j.message || "";
  $("#jobCancel").style.display = j.id ? "" : "none";
  $("#jobCancel").onclick = () => j.id && api("/api/cancel", { id: j.id });
}
function hideJob() { $("#jobbar").classList.remove("on"); }
function watchJob(id, done) {
  const tick = async () => {
    let j;
    try { j = await api("/api/job?id=" + id); } catch (e) { hideJob(); return; }
    showJob(j);
    if (j.state === "running") return setTimeout(tick, 500);
    hideJob();
    if (j.state === "error") toast(j.error, true);
    else if (done) done(j);
  };
  tick();
}

/* ---------------- SAMPLES page */
function renderSamples() {
  const has = S && S.table;
  $("#samplesEmpty").classList.toggle("hidden", !!has);
  $("#samplesBody").classList.toggle("hidden", !has);
  if (!has) return;
  $("#sampleTable").innerHTML = `<tr><th>Sample</th><th>Condition</th><th>Reference</th><th class="num">Rep.</th><th class="num">Reads</th></tr>` +
    design.samples.map((r, i) => `<tr><td class="mono">${esc(r.sample)}</td><td><input type="text" data-cond="${i}" value="${esc(r.condition)}"></td><td><label class="chk"><input type="checkbox" data-ref="${i}" ${r.reference ? "checked" : ""}> plasmid / T0</label></td><td class="num">${r.replicate}</td><td class="num">${fmt(S.table.totals[S.table.samples.indexOf(r.sample)])}</td></tr>`).join("");
  $$("[data-cond]").forEach((i) => (i.oninput = () => { design.samples[+i.dataset.cond].condition = i.value; renumber(); }));
  $$("[data-ref]").forEach((i) => (i.onchange = () => (design.samples[+i.dataset.ref].reference = i.checked)));
  renderComps();
  renderScreen();
}
function renderScreen() {
  const F = (S && S.screen_fields) || {};
  scr = { ...(S.screen || {}), ...scr };
  $("#screenForm").innerHTML = Object.entries(F).map(([k, f]) =>
    `<label><b>${esc(f.label)}</b><input type="text" inputmode="${f.type === "text" ? "text" : "decimal"}" data-sc="${k}" value="${esc(scr[k] == null ? "" : scr[k])}"><span class="hint">${esc(f.hint || "")}</span></label>`).join("");
  $$("[data-sc]").forEach((i) => (i.oninput = () => {
    const v = i.value.trim();
    if (v === "") delete scr[i.dataset.sc]; else scr[i.dataset.sc] = v;
    screenCalc();
  }));
  screenCalc();
}
function screenCalc() {
  const n = (k) => { const v = parseFloat(String(scr[k] == null ? "" : scr[k]).replace(/[, ]/g, "")); return v > 0 ? v : null; };
  const lib = (S.table && S.table.guides) || 0;
  const out = [];
  if (lib) {
    for (const [k, what] of [["cells_transduced", "after selection"], ["cells_bottleneck", "at the tightest passage"], ["cells_harvested", "at harvest"]]) {
      const c = n(k);
      if (c) { const x = c / lib; out.push(`<b class="${x >= 500 ? "pass" : x >= 200 ? "warn" : "fail"}">${fmt(x)}×</b> ${what}`); }
    }
  }
  const m = n("moi");
  if (m) { const one = 100 * (m * Math.exp(-m)) / (1 - Math.exp(-m)); out.push(`<b class="${m <= 0.4 && m >= 0.05 ? "pass" : m > 0.6 ? "fail" : "warn"}">MOI ${m}</b> — about ${one.toFixed(0)}% of transduced cells carry one guide`); }
  $("#screenCalc").innerHTML = out.length ? out.join('<span class="faint"> · </span>')
    : `<span class="muted">Fill these in and the quality check compares them with 500 cells per guide and MOI 0.3.</span>`;
}
function renumber() {
  const seen = {};
  design.samples.forEach((r) => (r.replicate = seen[r.condition] = (seen[r.condition] || 0) + 1));
  $$("#sampleTable tr").slice(1).forEach((tr, i) => (tr.children[3].textContent = design.samples[i].replicate));
}
function renderComps() {
  const samples = S.table.samples;
  $("#compList").innerHTML = design.comparisons.length ? design.comparisons.map((c, i) => `<div class="comp" data-ci="${i}">
    <div class="row"><input type="text" data-cname="${i}" value="${esc(c.name)}" style="flex:1"><button class="btn sm ghost" data-cdel="${i}">Remove</button></div>
    <div class="grp"><b>treatment</b>${samples.map((s) => `<span class="chip ${c.treatment.includes(s) ? "t" : ""}" data-ct="${i}" data-s="${esc(s)}">${esc(s)}</span>`).join("")}</div>
    <div class="grp"><b>control</b>${samples.map((s) => `<span class="chip ${c.control.includes(s) ? "c" : ""}" data-cc="${i}" data-s="${esc(s)}">${esc(s)}</span>`).join("")}</div></div>`).join("")
    : `<div class="empty">No comparisons yet. Mark a reference sample and press + Add.</div>`;
  $$("[data-cname]").forEach((i) => (i.oninput = () => (design.comparisons[+i.dataset.cname].name = i.value)));
  $$("[data-cdel]").forEach((b) => (b.onclick = () => { design.comparisons.splice(+b.dataset.cdel, 1); renderComps(); }));
  const tog = (arr, s) => (arr.includes(s) ? arr.splice(arr.indexOf(s), 1) : arr.push(s));
  $$("[data-ct]").forEach((c) => (c.onclick = () => { const k = design.comparisons[+c.dataset.ct]; tog(k.treatment, c.dataset.s); k.control = k.control.filter((x) => x !== c.dataset.s); renderComps(); }));
  $$("[data-cc]").forEach((c) => (c.onclick = () => { const k = design.comparisons[+c.dataset.cc]; tog(k.control, c.dataset.s); k.treatment = k.treatment.filter((x) => x !== c.dataset.s); renderComps(); }));
}
$("#addComp").onclick = () => {
  const refs = design.samples.filter((r) => r.reference).map((r) => r.sample);
  const used = new Set(design.comparisons.map((c) => c.name));
  const conds = [...new Set(design.samples.filter((r) => !r.reference).map((r) => r.condition))];
  const c = conds.find((x) => !used.has(`${x} vs reference`)) || "new";
  design.comparisons.push({ name: `${c} vs reference`, treatment: design.samples.filter((r) => r.condition === c && !r.reference).map((r) => r.sample), control: refs });
  renderComps();
};
$("#saveDesign").onclick = async () => {
  try { const d = await api("/api/design", { ...design, screen: scr }); design = d; S.design = d; QC = null; $("#qcDot").className = "dot"; toast("Saved."); go("quality"); runQC(); }
  catch (e) { toast(e.message, true); }
};

/* ---------------- QUALITY page */
const TH_LABELS = {
  coverage_pass: "Reads per guide — pass at", coverage_fail: "Reads per guide — fail below", mapped_pass: "Mapped fraction — pass at", mapped_fail: "Mapped fraction — fail below",
  zero_ref_pass: "Zero-count (reference) — pass up to", zero_ref_fail: "Zero-count (reference) — fail above", zero_late_pass: "Zero-count (later) — pass up to", zero_late_fail: "Zero-count (later) — fail above",
  gini_ref_pass: "Gini (reference) — pass up to", gini_ref_fail: "Gini (reference) — fail above", gini_late_pass: "Gini (later) — pass up to", gini_late_fail: "Gini (later) — fail above",
  skew_pass: "90/10 skew — pass up to", skew_fail: "90/10 skew — fail above", rep_r_pass: "Replicate r — pass at", rep_r_fail: "Replicate r — fail below", auc_pass: "Essential AUC — pass at", auc_fail: "Essential AUC — fail below",
  cells_pass: "Cells per guide — pass at", cells_fail: "Cells per guide — fail below", moi_pass_hi: "MOI — pass up to", moi_fail_hi: "MOI — fail above",
  guides_gene_pass: "Guides per gene — pass at", guides_gene_fail: "Guides per gene — fail below", ntc_pass: "Control guides — pass at", ntc_fail: "Control guides — fail below",
};
$("#thBtn").onclick = () => { $("#thPanel").classList.toggle("hidden"); renderThresholds(); };
function renderThresholds() {
  const th = (QC && QC.thresholds) || {};
  const def = (QC && QC.threshold_defaults) || {};
  $("#thGrid").innerHTML = Object.entries(TH_LABELS).map(([k, l]) => {
    const v = th[k], d = def[k];
    const changed = v != null && d != null && +v !== +d;
    return `<label>${l}<input type="text" class="${changed ? "set" : ""}" data-th="${k}" value="${v == null ? "" : v}" placeholder="${d == null ? "default" : d}"></label>`;
  }).join("");
  const n = Object.keys(TH_LABELS).filter((k) => th[k] != null && def[k] != null && +th[k] !== +def[k]).length;
  $("#thSaved").textContent = n ? `${n} changed from the default, remembered` : "all at the default";
  $$("[data-th]").forEach((i) => (i.oninput = () => i.classList.toggle("set", i.value.trim() !== "" && +i.value !== +def[i.dataset.th])));
}
$("#thReset").onclick = async () => {
  try {
    const r = await api("/api/thresholds/reset", { which: "pooled" });
    if (QC) { QC.thresholds = r.thresholds; }
    renderThresholds();
    toast("Thresholds back to the defaults \u2014 run the checks to apply them.");
  } catch (e) { toast(e.message, true); }
};
$("#qcBtn").onclick = () => runQC();
async function runQC() {
  if (!S.table) return toast("Load data first.", true);
  const thresholds = {};
  $$("[data-th]").forEach((i) => { if (i.value.trim() !== "" && isFinite(+i.value)) thresholds[i.dataset.th] = +i.value; });
  try {
    const j = await api("/api/qc", { thresholds });
    watchJob(j.id, async () => { QC = await api("/api/qc"); S.qc_overall = QC.overall; $("#qcDot").className = "dot " + QC.overall; renderQC(); });
  } catch (e) { toast(e.message, true); }
}
async function renderQC() {
  if (!QC && S && S.has_qc) { try { QC = await api("/api/qc"); } catch (e) {} }
  $("#qcEmpty").classList.toggle("hidden", !!QC);
  if (!QC) { $("#qcCards").innerHTML = ""; $("#qcOverall").innerHTML = ""; return; }
  const n = (s) => QC.checks.filter((c) => c.status === s).length;
  $("#qcOverall").innerHTML = `<div class="card overall">${statusTag(QC.overall)}<span class="big">${QC.overall === "pass" ? "The screen passes every check" : QC.overall === "warn" ? "Usable, with points to look at" : QC.overall === "fail" ? "Problems that affect the results" : "Not enough information"}</span><span class="spacer"></span><span class="muted small">${n("pass")} pass · ${n("warn")} warn · ${n("fail")} fail · ${n("na")} not applicable</span></div>`;
  const refs = (S.design.samples || []).filter((r) => r.reference).map((r) => r.sample);
  const condOf = Object.fromEntries((S.design.samples || []).map((r) => [r.sample, r.condition]));
  const box = $("#qcCards");
  box.innerHTML = "";
  for (const c of QC.checks) {
    const wide = c.id === "essentials" || c.id === "pca";
    const el = document.createElement("div");
    el.className = "card qc" + (wide ? " wide" : "");
    el.innerHTML = `<h3>${statusTag(c.status)}${esc(c.name)}</h3><div class="why">${esc(c.why)}</div><div class="rule">${esc(c.rule || "")}</div><div class="plot mt" id="qc-${c.id}"></div>`;
    box.appendChild(el);
    const P = el.querySelector(".plot");
    const th = QC.thresholds;
    if (c.id === "reads" && c.per_sample.length) C.bars(P, c.per_sample.map((p) => ({ label: p.sample, value: p.coverage, status: p.status, tip: `${fmt(p.value)} reads · ${fmt(p.coverage)} per guide` })), { lines: [{ v: th.coverage_pass, label: "pass" }, { v: th.coverage_fail, label: "fail" }], vf: (v) => fmt(v) + "×" });
    else if (c.id === "mapped" && c.per_sample.length) C.bars(P, c.per_sample.map((p) => ({ label: p.sample, value: p.value * 100, status: p.status })), { max: 100, lines: [{ v: th.mapped_pass * 100, label: "pass" }], vf: (v) => fmt(v) + "%" });
    else if (c.id === "zeros") C.bars(P, c.per_sample.map((p) => ({ label: p.sample + (p.reference ? " ●" : ""), value: p.value * 100, status: p.status, tip: `${p.count.toLocaleString()} guides with no reads (${(p.value * 100).toFixed(2)}%)${p.reference ? " — reference sample" : ""}` })), { vf: (v) => fmt(v) + "%" });
    else if (c.id === "gini") C.bars(P, c.per_sample.map((p) => ({ label: p.sample + (p.reference ? " ●" : ""), value: p.value, status: p.status, tip: `MAGeCK Gini (log counts) ${p.value.toFixed(4)} · Gini on raw counts ${p.raw.toFixed(3)}` })), { min: 0, vf: (v) => v.toFixed(3) });
    else if (c.id === "coverage") { if (QC.distributions) C.boxes(P, QC.distributions.samples, QC.distributions.quantiles, refs); }
    else if (c.id === "replicates" && c.pairs.length) P.innerHTML = `<table class="tbl"><tr><th>Condition</th><th>Replicates</th><th class="num">r (guide LFC)</th><th class="num">ρ (Spearman)</th><th class="num">r (log counts)</th><th></th></tr>${c.pairs.map((p) => `<tr><td>${esc(p.condition)}</td><td>${esc(p.a)} · ${esc(p.b)}</td><td class="num">${fmt(p.value)}</td><td class="num">${fmt(p.spearman)}</td><td class="num">${fmt(p.r_counts)}</td><td>${statusTag(p.status)}</td></tr>`).join("")}</table>`;
    else if (c.id === "pca" && c.pca) C.pca(P, c.pca, condOf, refs);
    else if (c.id === "essentials" && c.per_condition && c.per_condition.length) renderEss(P, c);
    else if (c.id === "representation") renderRep(P, c);
    else if (c.id === "guides_gene" && c.hist && c.hist.length) C.bars(P, c.hist.map((h) => ({ label: String(h.n), value: h.genes, status: h.n >= th.guides_gene_pass ? "pass" : h.n >= th.guides_gene_fail ? "warn" : "fail", tip: `${h.genes.toLocaleString()} genes with ${h.n} guide${h.n === 1 ? "" : "s"}` })), { vf: (v) => fmt(v) });
    else if (c.id === "ntc" && c.per_condition && c.per_condition.length) P.innerHTML = `<table class="tbl"><tr><th>Condition</th><th class="num">Median LFC, controls</th><th class="num">Spread (robust SD)</th><th class="num">Median LFC, targeting</th><th class="num">Smallest effect above noise</th><th></th></tr>${c.per_condition.map((p) => `<tr><td>${esc(p.condition)}</td><td class="num">${fmt(p.median)}</td><td class="num">${fmt(p.sd)}</td><td class="num">${fmt(p.median_targeting)}</td><td class="num">±${fmt(p.floor)}</td><td>${statusTag(p.status)}</td></tr>`).join("")}</table>`;
  }
}
function renderRep(P, c) {
  if (!c.rows.length && !c.moi) { P.innerHTML = ""; return; }
  const cells = c.rows.map((r) => `<div class="stat"><b class="${r.status}">${fmt(r.per_guide)}×</b><span>${esc(r.step)}<br>${fmt(r.cells)} cells</span></div>`).join("");
  const moi = c.moi ? `<div class="stat"><b class="${c.moi.status}">${c.moi.moi}</b><span>MOI<br>about ${c.moi.one_guide.toFixed(0)}% carry one guide</span></div>` : "";
  P.innerHTML = `<div class="repcard">${cells}${moi}</div><div class="small muted mt">Cells per guide = cells ÷ ${fmt(c.library_guides)} guides in the library.</div>`;
}
function renderEss(P, c) {
  P.innerHTML = `<div class="row small muted">${c.n_ceg} core-essential (CEGv2) and ${c.n_neg} non-essential (NEGv1) genes are in this library. Gene LFC = mean of its guides, replicates averaged, against the reference.</div>
    <table class="tbl mt"><tr><th>Condition</th><th class="num">AUC</th><th class="num">Recall at 5% FDR</th><th class="num">Median LFC essential</th><th class="num">Median LFC non-essential</th><th></th></tr>
    ${c.per_condition.map((p, i) => `<tr class="click ${i === 0 ? "sel" : ""}" data-pc="${i}"><td>${esc(p.condition)}</td><td class="num">${p.value.toFixed(3)}</td><td class="num">${(p.recall_at_fdr5 * 100).toFixed(1)}%</td><td class="num">${fmt(p.median_ess)}</td><td class="num">${fmt(p.median_non)}</td><td>${statusTag(p.status)}</td></tr>`).join("")}</table>
    <div class="grid g3 mt"><div><h3 class="small">ROC</h3><div id="essRoc"></div></div><div><h3 class="small">Precision–recall</h3><div id="essPr"></div></div><div><h3 class="small">Distribution <span class="faint" id="essWhich"></span></h3><div id="essHist"></div></div></div>`;
  const pal = [C.css("--accent"), C.css("--ess"), C.css("--accent2"), "#e87ba4", "#008300", "#6250d6"];
  C.lines($("#essRoc"), c.per_condition.map((p, i) => ({ name: p.condition, x: p.curves.fpr, y: p.curves.tpr, color: pal[i % pal.length] })), { xl: "FPR", yl: "TPR", diag: true });
  C.lines($("#essPr"), c.per_condition.map((p, i) => ({ name: p.condition, x: p.curves.recall, y: p.curves.precision, color: pal[i % pal.length] })), { xl: "recall", yl: "precision" });
  const legend = `<div class="legend">${c.per_condition.map((p, i) => `<span><i style="background:${pal[i % pal.length]}"></i>${esc(p.condition)}</span>`).join("")}</div>`;
  $("#essRoc").insertAdjacentHTML("beforeend", legend);
  const show = (i) => { $("#essWhich").textContent = "· " + c.per_condition[i].condition; C.hist2($("#essHist"), c.per_condition[i].hist, ["core-essential", "non-essential"]); $$("[data-pc]").forEach((r) => r.classList.toggle("sel", +r.dataset.pc === i)); };
  $$("[data-pc]").forEach((r) => (r.onclick = () => show(+r.dataset.pc)));
  show(c.per_condition.length - 1 >= 0 ? c.per_condition.map((p) => p.value).indexOf(Math.max(...c.per_condition.map((p) => p.value))) : 0);
}

/* ---------------- ANALYSE page */
function renderAnalyse() {
  if (!S) return;
  $("#normOpts").innerHTML = Object.entries(S.norm_methods).map(([k, v]) => `<label><input type="radio" name="norm" value="${k}" ${k === "median" ? "checked" : ""}> <b>${esc(v.split(":")[0])}</b><span class="desc">${esc((t => t.charAt(0).toUpperCase() + t.slice(1))(v.split(":").slice(1).join(":").trim()))}</span></label>`).join("");
  if (!$("#ctlPattern").value) $("#ctlPattern").value = S.control_pattern;
  const comps = (S.design && S.design.comparisons) || [];
  $("#aComps").innerHTML = comps.length ? comps.map((c) => `<label><input type="checkbox" data-ac="${esc(c.name)}" checked> ${esc(c.name)} <span class="faint small">(${c.treatment.length} vs ${c.control.length})</span></label>`).join("") : `<div class="empty">Set up comparisons on the Samples page first.</div>`;
  ctlPreview();
}
async function ctlPreview() {
  if (!S || !S.table) return;
  const how = ($("input[name=ctl]:checked") || {}).value || "none";
  if (how === "none") { $("#ctlInfo").textContent = "No control guides used."; return; }
  try {
    const j = await api("/api/controls", { controls: how, control_pattern: $("#ctlPattern").value, control_list: $("#ctlList").value });
    $("#ctlInfo").textContent = j.guides ? `${j.guides.toLocaleString()} control guides in ${j.n_genes} gene name(s): ${j.genes.join(", ")}${j.n_genes > 30 ? "…" : ""}` : "No guides match.";
  } catch (e) { $("#ctlInfo").textContent = e.message; }
}
$$("input[name=ctl]").forEach((r) => (r.onchange = ctlPreview));
$("#ctlPattern").oninput = () => { clearTimeout(window._cp); window._cp = setTimeout(ctlPreview, 300); };
$("#ctlList").oninput = () => { clearTimeout(window._cp); window._cp = setTimeout(ctlPreview, 300); };
$("#runBtn").onclick = async () => {
  const opts = {
    norm: ($("input[name=norm]:checked") || {}).value || "median", controls: ($("input[name=ctl]:checked") || {}).value,
    control_pattern: $("#ctlPattern").value, control_list: $("#ctlList").value, control_null: $("#ctlNull").checked,
    remove_zero: $("#aZero").value, gene_lfc: $("#aLfc").value, gene_test_fdr_threshold: +$("#aAlpha").value || 0.25,
    comparisons: $$("[data-ac]:checked").map((i) => i.dataset.ac),
  };
  if (!opts.comparisons.length) return toast("Choose at least one comparison.", true);
  try {
    const j = await api("/api/analyse", opts);
    watchJob(j.id, async () => { await refresh(); toast("Analysis finished."); go("results"); });
  } catch (e) { toast(e.message, true); }
};

/* ---------------- RESULTS page */
let rankDir = "neg", sortCol = "neg_rank", sortAsc = true, gpage = 0, selGene = null;
async function loadResults(name) {
  const list = (S && S.results) || [];
  $("#resEmpty").classList.toggle("hidden", list.length > 0);
  $("#resBody").classList.toggle("hidden", !list.length);
  if (!list.length) return;
  const sel = $("#rComp");
  const cur = name || sel.value || (R && R.name) || list[0].name;
  sel.innerHTML = list.map((r) => `<option ${r.name === cur ? "selected" : ""}>${esc(r.name)}</option>`).join("");
  if (!R || R.name !== cur) { R = await api("/api/result?name=" + encodeURIComponent(cur)); R.idx = Object.fromEntries(R.columns.map((c, i) => [c, i])); selGene = null; gpage = 0; }
  renderResults();
}
$("#rComp").onchange = () => loadResults($("#rComp").value);
$$("#rankDir button").forEach((b) => (b.onclick = () => { rankDir = b.dataset.d; $$("#rankDir button").forEach((x) => x.classList.toggle("on", x === b)); drawRank(); }));
$("#fdrCut").onchange = () => renderResults();
function cut() { const v = +$("#fdrCut").value; return v > 0 && v < 1 ? v : 0.05; }
function renderResults() {
  const I = R.idx, rows = R.rows, f = cut();
  const nneg = rows.filter((r) => r[I.neg_fdr] <= f).length, npos = rows.filter((r) => r[I.pos_fdr] <= f).length;
  const sf = Object.entries(R.size_factors).map(([k, v]) => `${k} ${v.toFixed(3)}`).join(" · ");
  $("#resTiles").innerHTML = [
    ["Depleted genes", nneg, `FDR ≤ ${f} (negative selection)`], ["Enriched genes", npos, `FDR ≤ ${f} (positive selection)`],
    ["Genes tested", rows.length.toLocaleString(), `${R.comparison.treatment.length} treatment vs ${R.comparison.control.length} control`],
    ["α cut-off", `${R.alpha_neg.toFixed(3)} / ${R.alpha_pos.toFixed(3)}`, "depleted / enriched (fraction of guides used by RRA)"]].map(([k, v, s]) => `<div class="card tile"><div class="k">${k}</div><div class="v">${v}</div><div class="s">${esc(s)}</div></div>`).join("");
  $("#resNotes").innerHTML = `<b>${esc(R.name)}</b>: treatment ${R.comparison.treatment.map(esc).join(", ")} vs control ${R.comparison.control.map(esc).join(", ")} · normalisation ${esc(R.settings.norm)}${R.settings.controls ? ` (${R.settings.controls} control guides)` : ""} · size factors ${esc(sf)} · variance model k=${R.model[0].toFixed(3)}, b=${R.model[1].toFixed(3)} · run ${esc(R.when)}${R.norm_warnings.length ? "<br>" + R.norm_warnings.map(esc).join("<br>") : ""}`;
  drawVolcano(); drawRank(); drawTable();
}
const colors = () => ({ neg: C.css("--neg"), pos: C.css("--pos"), ns: C.css("--ns") });
function drawVolcano() {
  const I = R.idx, f = cut();
  const pts = R.rows.map((r) => {
    const dep = r[I.lfc] < 0, p = dep ? r[I.neg_p] : r[I.pos_p], q = dep ? r[I.neg_fdr] : r[I.pos_fdr];
    return { x: r[I.lfc], y: -Math.log10(Math.max(p, 1e-300)), k: q <= f ? (dep ? "neg" : "pos") : "ns", label: r[0], r };
  });
  const top = [...pts].filter((p) => p.k !== "ns").sort((a, b) => b.y - a.y || Math.abs(b.x) - Math.abs(a.x)).slice(0, 12);
  C.scatter($("#volcano"), pts, {
    colors: colors(), xl: "gene log₂ fold change (median of guides)", yl: "−log₁₀ p (RRA)", vline: 0, labels: top, selected: selGene, z: (p) => (p.k === "ns" ? 0 : 1),
    tip: (p) => `<b>${esc(p.label)}</b><br>LFC ${fmt(p.x)}<br>depleted: p ${fmt(p.r[I.neg_p])}, FDR ${fmt(p.r[I.neg_fdr])}, rank ${p.r[I.neg_rank]}<br>enriched: p ${fmt(p.r[I.pos_p])}, FDR ${fmt(p.r[I.pos_fdr])}, rank ${p.r[I.pos_rank]}<br>${p.r[I.guides]} guides`,
    click: (p) => showGene(p.label),
  });
  $("#volcano").insertAdjacentHTML("beforeend", `<div class="legend"><span><i style="background:${colors().neg}"></i>depleted, FDR ≤ ${f}</span><span><i style="background:${colors().pos}"></i>enriched, FDR ≤ ${f}</span><span><i style="background:${colors().ns}"></i>not significant</span></div>`);
}
function drawRank() {
  const I = R.idx, f = cut(), d = rankDir;
  const pts = R.rows.map((r) => ({ x: r[I[d + "_rank"]], y: -Math.log10(Math.max(r[I[d + "_score"]], 1e-300)), k: r[I[d + "_fdr"]] <= f ? d : "ns", label: r[0], r }));
  const top = pts.filter((p) => p.x <= 10);
  C.scatter($("#rankplot"), pts, {
    colors: colors(), xl: `gene rank (${d === "neg" ? "depleted" : "enriched"})`, yl: "−log₁₀ RRA score", labels: top, selected: selGene, z: (p) => (p.k === "ns" ? 0 : 1),
    tip: (p) => `<b>${esc(p.label)}</b><br>rank ${p.x} · RRA score ${fmt(p.r[I[d + "_score"]])}<br>p ${fmt(p.r[I[d + "_p"]])} · FDR ${fmt(p.r[I[d + "_fdr"]])}<br>${p.r[I[d + "_good"]]} of ${p.r[I.guides]} guides under α`,
    click: (p) => showGene(p.label),
  });
}
const COLS = [["gene", "Gene"], ["guides", "Guides"], ["lfc", "LFC"], ["neg_p", "p depleted"], ["neg_fdr", "FDR depl."], ["neg_rank", "Rank depl."], ["pos_p", "p enriched"], ["pos_fdr", "FDR enr."], ["pos_rank", "Rank enr."]];
function drawTable() {
  const I = R.idx, q = $("#gSearch").value.trim().toUpperCase(), f = cut();
  let rows = R.rows;
  if (q) rows = rows.filter((r) => String(r[0]).toUpperCase().includes(q));
  const ci = I[sortCol];
  rows = [...rows].sort((a, b) => (a[ci] < b[ci] ? -1 : a[ci] > b[ci] ? 1 : 0) * (sortAsc ? 1 : -1));
  const per = 100, pages = Math.max(1, Math.ceil(rows.length / per));
  gpage = Math.min(gpage, pages - 1);
  const show = rows.slice(gpage * per, gpage * per + per);
  $("#geneTable").innerHTML = `<tr>${COLS.map(([k, l]) => `<th class="${k === "gene" ? "" : "num"} ${k === sortCol ? "sorted" + (sortAsc ? " asc" : "") : ""}" data-sort="${k}">${l}</th>`).join("")}</tr>` +
    show.map((r) => `<tr class="click ${r[0] === selGene ? "sel" : ""}" data-g="${esc(r[0])}">${COLS.map(([k]) => {
      const v = r[I[k]];
      if (k === "gene") return `<td><b>${esc(v)}</b></td>`;
      const hit = (k === "neg_fdr" || k === "pos_fdr") && v <= f;
      return `<td class="num" ${hit ? `style="color:${k === "neg_fdr" ? "var(--neg)" : "var(--pos)"};font-weight:600"` : ""}>${typeof v === "number" ? fmt(v) : esc(v)}</td>`;
    }).join("")}</tr>`).join("");
  $("#gPager").innerHTML = `${rows.length.toLocaleString()} genes · page ${gpage + 1} of ${pages} <span class="spacer"></span><button class="btn sm" id="gPrev" ${gpage ? "" : "disabled"}>‹</button><button class="btn sm" id="gNext" ${gpage < pages - 1 ? "" : "disabled"}>›</button>`;
  $("#gPrev").onclick = () => { gpage--; drawTable(); };
  $("#gNext").onclick = () => { gpage++; drawTable(); };
  $$("#geneTable th").forEach((th) => (th.onclick = () => { const k = th.dataset.sort; if (k === sortCol) sortAsc = !sortAsc; else { sortCol = k; sortAsc = !/fdr|_p|rank|gene/.test(k) ? false : true; } drawTable(); }));
  $$("#geneTable tr.click").forEach((tr) => (tr.onclick = () => showGene(tr.dataset.g)));
}
$("#gSearch").oninput = () => { gpage = 0; drawTable(); };
async function showGene(g) {
  selGene = g;
  drawTable(); drawVolcano(); drawRank();
  try {
    const d = await api(`/api/gene?name=${encodeURIComponent(R.name)}&gene=${encodeURIComponent(g)}`);
    const s = d.summary || {};
    const nc = d.n_control;
    $("#geneView").innerHTML = `<div class="gv-head"><h2>${esc(d.gene)}</h2><span class="muted small">${d.guides.length} guides · median LFC ${fmt(s.lfc)}</span></div>
      <div class="kv"><span>Depleted</span><span>RRA score ${fmt(s.neg_score)} · p ${fmt(s.neg_p)} · FDR ${fmt(s.neg_fdr)} · rank ${s.neg_rank} · ${s.neg_good} guides under α</span>
      <span>Enriched</span><span>RRA score ${fmt(s.pos_score)} · p ${fmt(s.pos_p)} · FDR ${fmt(s.pos_fdr)} · rank ${s.pos_rank} · ${s.pos_good} guides under α</span></div>
      <div class="plot" id="gvBars"></div>
      <h3 class="mt">Normalised counts</h3>
      <div class="scroll"><table class="tbl"><tr><th>Guide</th>${d.samples.map((x, i) => `<th class="num" style="color:${i < nc ? "var(--neg)" : "var(--pos)"}">${esc(x)}</th>`).join("")}<th class="num">LFC</th><th class="num">p low</th><th class="num">p high</th><th class="num">FDR</th></tr>
      ${d.guides.map((r) => `<tr style="${r.valid ? "" : "opacity:.5"}"><td class="mono small">${esc(r.guide)}</td>${r.norm.map((x, i) => `<td class="num" title="raw ${r.raw[i]}">${fmt(x)}</td>`).join("")}<td class="num">${fmt(r.lfc)}</td><td class="num">${fmt(r.p_low)}</td><td class="num">${fmt(r.p_high)}</td><td class="num">${fmt(r.fdr)}</td></tr>`).join("")}</table></div>
      <p class="faint small">Blue columns are control samples, red are treatment. Hover a count for the raw value. Faded guides were left out of the gene test (zero counts in both groups).</p>`;
    C.guideBars($("#gvBars"), d.guides);
  } catch (e) { toast(e.message, true); }
}
/* export */
$("#expBtn").onclick = (e) => { e.stopPropagation(); $("#expMenu").classList.toggle("on"); };
$$("[data-exp]").forEach((a) => (a.onclick = () => {
  const what = a.dataset.exp, name = R ? R.name : "";
  if ((what === "genes" || what === "guides") && !R) return toast("Run the analysis first.", true);
  const url = `/api/export?what=${what}&fmt=${a.dataset.fmt}&name=${encodeURIComponent(name)}`;
  const l = document.createElement("a"); l.href = url; l.download = ""; document.body.appendChild(l); l.click(); l.remove();
  $("#expMenu").classList.remove("on");
}));

/* ---------------- METHODS page */
function renderMethods() {
  const th = (QC && QC.threshold_text) || null;
  const thList = th ? Object.values(th).map((t) => `<li>${esc(t)}</li>`).join("") : "<li>Run the quality checks once to list the thresholds in use.</li>";
  $("#methodsBody").innerHTML = `
  <div class="card"><h3>Read counting</h3><p>As <code>mageck count</code> (0.5.9.5): the guide's start in the read is found from the first 100,000 reads, then every read is matched exactly at those positions, longest guide length first. Guides with repeated sequences in the library keep their first record. Options: one mismatch (assigned only when a single guide is one base away), reverse-complement reads, or locating the guide by the scaffold that follows it.</p></div>
  <div class="card"><h3>Normalisation</h3><ul><li><b>Median ratio</b> (default, as MAGeCK): factor = 1 / median over guides of count ÷ geometric mean of (count + 1) across the samples of the comparison. Falls back to total count if a sample has more than 45% zero-count guides.</li><li><b>Total count</b>: same total reads per sample.</li><li><b>Control guides</b>: median ratio over control guides only.</li></ul></div>
  <div class="card"><h3>Guide test (MAGeCK-style)</h3><p>Control and treatment means are the medians of normalised replicate counts. Variance comes from a mean–variance model, log₂(var − mean + 1) = k·log₂(mean + 1) + b, fitted on the control replicates (all samples when there is one control). The p-value for depletion or enrichment uses a normal distribution truncated at zero, as MAGeCK 0.5.9.5 ships; guide scores are the normal quantiles of the depletion p-value. Guide LFC = log₂(treatment + 1) − log₂(control + 1).</p></div>
  <div class="card"><h3>Gene test (α-RRA)</h3><p>Guides are ranked by score; each gene takes the smallest beta-distribution probability of its top guides that fall under α (α = fraction of guides with p ≤ 0.25, kept between 0.05 and 0.5 for depletion). Gene p-values come from scoring random guide sets of the same size (≥ 100,000 draws, using MAGeCK's own random-number generator and seed, so results are reproducible); FDR by Benjamini–Hochberg. Gene LFC = median of its guides.</p></div>
  <div class="card"><h3>Quality checks — defaults (to confirm)</h3><ul>${thList}</ul><p class="small muted">Essential genes: CEGv2 (684) and non-essential NEGv1 (927) human gene lists from the Hart lab BAGEL repository (MIT licence).</p></div>
  <div class="card"><h3>Thresholds</h3><p>Every pass / warn / fail line is a scientific choice, not a fact. The defaults follow published guidance and were confirmed for this tool on 7 November 2026: <b>300 / 100</b> reads per guide; <b>80% / 65%</b> mapped reads; Gini <b>0.10 / 0.20</b> at baseline and 0.20 / 0.30 later; zero-count guides <b>1% / 5%</b> in the reference and 5% / 15% later; library skew (90/10) <b>10 / 20</b>; replicate Pearson r <b>0.60 / 0.40</b>, with Spearman shown beside it; essential-gene AUC <b>0.90 / 0.80</b>; representation <b>500× / 200×</b> cells per guide; MOI <b>0.4 / 0.6</b>; <b>4 / 2</b> guides per gene; <b>100 / 20</b> non-targeting controls. For arrayed screens: plate median ratio 2 / 5, Z′ 0.5 / 0, |SSMD| 3 / 2, edge effect 10% / 20%, row and column share 10% / 25%, replicate r 0.70 / 0.40; the default normalisation is <b>robust z</b> and a well is a hit when its combined score passes <b>2</b> and every replicate also passes <b>1</b>. Change any of them under <i>Thresholds</i> on the Quality page (pooled and arrayed screens each have their own panel): a box you change is kept for next time on this Mac, and <i>Restore defaults</i> puts them all back.</p></div>
  <div class="card"><h3>Screen design: representation and MOI</h3><p>Cells per guide = the cells you enter for each step ÷ the guides in the library. At least 500 cells per guide is expected at every step (transduction, the tightest passage, harvest): below that, guides disappear because too few cells carried them, which looks exactly like a real dropout afterwards and cannot be fixed by sequencing more deeply. MOI about 0.3 keeps one guide per cell; the app shows the Poisson fraction of transduced cells carrying exactly one guide, so you can see what a higher MOI costs. These are the numbers you enter on the Samples page — the app cannot read them from the counts (to confirm: 500× pass, 200× fail; MOI 0.4 pass, 0.6 fail).</p></div>
  <div class="card"><h3>Guides per gene and non-targeting controls</h3><p>Guides per gene is counted over targeting guides only; genes with fewer than four are listed, because a gene score built on one or two guides cannot tell a real effect from one bad guide. Non-targeting controls are found by name (non-targeting, NTC, safe-harbour, EGFP, LacZ, luciferase, AAVS1, or your own pattern). They should sit at a log₂ fold change of zero: a shift means the whole library moved, which is a normalisation problem rather than biology. Their robust spread (IQR ÷ 1.349) is the noise floor, reported as the smallest gene effect that can be told apart from no effect.</p></div>
  <div class="card"><h3>BAGEL2 Bayes factors</h3><p>Ported from BAGEL.py (Kim &amp; Hart 2021, MIT): counts + 5 scaled to 10 million reads, log₂ fold change against the summed controls; in 10-fold cross-validation over genes, kernel densities of guide fold changes of reference essential (CEGv2) and non-essential (NEGv1) genes give a log-likelihood ratio, fitted by a straight line that scores the held-out genes; a gene's Bayes factor is the sum over its guides and replicates. BF &gt; 5 is shown as essential: the value BAGEL2's tutorial uses, confirmed here because on HAP1 TKOv3 T18 it calls 1,796 genes at 94.0% recall with no false positives, and still holds 98% precision on the too-early T3 timepoint, where BF &gt; 1 would admit 55 non-essential genes. You can type any cut-off on the results page to see how the count moves; exports use 5. Mouse screens use the Hart lab mouse CEG/NEG lists.</p></div>
  <div class="card"><h3>drugZ</h3><p>Ported from drugz.py (Colic et al. 2019, MIT): treated vs untreated fold changes per replicate pair (or of the means, unpaired), guide z-scores from a standard deviation estimated in a sliding window of 1,000 guides sorted by control reads and never allowed to decrease, gene normZ = z-scored sum of guide z over √observations; p-values and FDR for synthetic sick/lethal (negative) and suppressor (positive) genes.</p></div>
  <div class="card"><h3>Copy-number correction</h3><p>After CRISPRcleanR (Iorio et al. 2018), without copy-number data: guide fold changes ordered along each chromosome are segmented (binary segmentation, permutation test, α = 0.01, segments of ≥ 10 guides); segments covering ≥ 3 genes have their mean removed (for MAGeCK, the treatment counts are rescaled accordingly). On HAP1 TKOv3 it removes planted amplicons of 6–81 genes to within ±12%, but also shifts clusters of neighbouring essential genes (essential-gene AUC 0.998 → 0.994), so it is off by default (confirmed) and meant for cancer lines with copy-number changes. Segments must cover at least 3 genes and 10 guides, CRISPRcleanR's own defaults, confirmed so results stay comparable with that tool.</p></div>
  <div class="card"><h3>Pathways and gene sets</h3><p>Offline gene sets: Gene Ontology (biological process, molecular function, cellular component; human and mouse; annotations propagated to parent terms) and Reactome pathways (human), both CC BY 4.0, sets of 10–500 genes; or your own .gmt. Over-representation: one-sided hypergeometric test of the hit genes against all genes in the screen. Ranked: Mann–Whitney test of the set's gene scores against all other genes, no cut-off needed. Benjamini–Hochberg FDR.</p></div>
  <div class="card"><h3>Arrayed screens — plate normalisation</h3><ul><li><b>Robust z-score</b> (default, to confirm): (value − plate median of sample wells) ÷ (1.4826 × MAD).</li><li><b>B-score</b> (Brideau et al. 2003): a two-way median polish removes row and column effects within each plate; residual ÷ (1.4826 × MAD).</li><li><b>Plate median</b>, <b>percent of control</b> (÷ mean of negative controls), <b>NPI</b> (0% = negative, 100% = positive controls), and the <b>classic plate z-score</b> (plate mean and SD of all wells, as in the 2003 KcViab analysis).</li></ul><p>Non-z methods are converted to robust z-scores over all sample wells. For “lower signal” screens the sign is flipped so a hit always scores high.</p></div>
  <div class="card"><h3>Arrayed screens — checks and hits</h3><p>Plate-to-plate medians; Z′ and SSMD of the controls per plate when every plate has at least two wells of each (otherwise pooled across plates, and the page says so); edge wells against interior wells; the share of within-plate variance explained by rows and columns; replicate correlation. A well is a hit when the combined score (mean of replicates by default) reaches the threshold <i>and</i> every replicate reaches the replicate threshold (defaults 2 and 1, to confirm). A gene is a hit when at least half of its wells are hits. SSMD scores use Zhang's robust form for screens without replicates.</p></div>
  <div class="card"><h3>Validation</h3><p>Automated tests compare this tool with the official programs on public data: read counts on the MAGeCK demo FASTQ files against <code>mageck count</code>, gene rankings on the MAGeCK demos and the HAP1 TKOv3 screen against <code>mageck test</code> (identical counts; statistics equal to MAGeCK's printed precision; identical hit lists). The KcViab arrayed screen is compared with its published 2003 scores: the classic plate z-score reproduces 99.98% of them to two decimals. BAGEL2 Bayes factors and drugZ scores on HAP1 TKOv3 are tested against the official Hart lab programs: identical to their printed precision for every gene, with the same essential, synthetic-lethal and suppressor genes called. See <code>tests/</code> and <code>test-data/reference/</code> in the tool folder.</p><p class="small muted">Methods follow: Li W et al. MAGeCK, Genome Biology 15:554 (2014); Hart T et al. G3 7:2719 (2017); Kim E &amp; Hart T, Genome Medicine 13:2 (2021); Colic M et al., Genome Medicine 11:52 (2019); Iorio F et al., BMC Genomics 19:604 (2018); Brideau C et al., J Biomol Screen 8:634 (2003); Boutros M et al., Genome Biology 7:R66 (2006).</p></div>`;
}

/* ---------------- hero read animation */
(function heroRead() {
  const el = $("#heroRead"), stat = $("#hrStat");
  const NB = "ACGT", SC = "GTTTTAGAGCTA", r = (n) => Array.from({ length: n }, () => NB[Math.floor(Math.random() * 4)]).join("");
  let lines = [];
  function line() {
    const lead = r(Math.random() < 0.5 ? 0 : 3), g = r(20), ok = Math.random() < 0.85;
    const gm = ok ? g : g.slice(0, 7) + "N" + g.slice(8);
    return `<span class="o">${lead}</span><span class="${ok ? "g" : "x"}">${gm}</span><span class="s">${SC}</span>`;
  }
  let n = 0, hit = 0;
  setInterval(() => {
    lines.push(line()); if (lines.length > 5) lines.shift();
    n++; if (!/class="x"/.test(lines[lines.length - 1])) hit++;
    el.innerHTML = lines.join("\n");
    stat.textContent = `${hit}/${n} reads matched · ${(100 * hit / n).toFixed(0)}%`;
  }, 900);
})();

function redraw() {
  if (MODE === "arrayed") { if (window.A && page !== "data") window.A.render(page); return; }
  if (page === "quality" && QC) renderQC(); if (page === "results" && R) renderResults();
}
let rz; window.addEventListener("resize", () => { clearTimeout(rz); rz = setTimeout(() => { if (page === "results" && R) { drawVolcano(); drawRank(); } }, 200); });

refresh().then(() => setMode(MODE)).catch((e) => toast(e.message, true));
