/* Small chart kit (no libraries): SVG for few marks, canvas for thousands of genes. Colours come from CSS
   variables so charts follow the Light | Dark switch; every chart has a hover tooltip. */
const C = (() => {
  const css = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
  const tipEl = () => document.getElementById("tip");
  function tip(html, e) {
    const t = tipEl();
    if (!html) { t.classList.remove("on"); return; }
    t.innerHTML = html;
    t.classList.add("on");
    const w = t.offsetWidth, h = t.offsetHeight;
    let x = e.clientX + 14, y = e.clientY + 14;
    if (x + w > innerWidth - 8) x = e.clientX - w - 14;
    if (y + h > innerHeight - 8) y = e.clientY - h - 14;
    t.style.left = x + "px"; t.style.top = y + "px";
  }
  function ticks(lo, hi, n = 5) {
    if (!isFinite(lo) || !isFinite(hi)) return [0];
    if (hi === lo) { hi = lo + 1; }
    const span = hi - lo, step0 = Math.pow(10, Math.floor(Math.log10(span / n)));
    const err = (span / n) / step0;
    const step = step0 * (err >= 7.5 ? 10 : err >= 3.5 ? 5 : err >= 1.5 ? 2 : 1);
    const out = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(+v.toFixed(10));
    return out;
  }
  const fmt = (v) => {
    if (v === null || v === undefined || !isFinite(v)) return "–";
    const a = Math.abs(v);
    if (a >= 1e6 && a < 1e12) return (v / 1e6).toFixed(a >= 1e8 ? 0 : 1) + " M";
    if (a !== 0 && (a < 1e-3 || a >= 1e12)) return v.toExponential(1).replace("e+", "e");
    if (a >= 1000) return v.toLocaleString("en-US", { maximumFractionDigits: 0 });
    return +v.toPrecision(3) + "";
  };
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const svgEl = (w, h) => `<svg viewBox="0 0 ${w} ${h}" preserveAspectRatio="xMidYMid meet">`;

  function frame(o) {
    // returns svg pieces for axes/grid given scales
    const { W, H, m, xt, yt, sx, sy, xl, yl, xf = fmt, yf = fmt } = o;
    let s = "";
    for (const v of yt) s += `<line class="gl" x1="${m.l}" x2="${W - m.r}" y1="${sy(v)}" y2="${sy(v)}"/><text x="${m.l - 6}" y="${sy(v) + 3.5}" text-anchor="end">${yf(v)}</text>`;
    for (const v of xt) s += `<text x="${sx(v)}" y="${H - m.b + 15}" text-anchor="middle">${xf(v)}</text>`;
    s += `<line class="ax" x1="${m.l}" x2="${W - m.r}" y1="${H - m.b}" y2="${H - m.b}"/>`;
    if (xl) s += `<text class="lab" x="${(m.l + W - m.r) / 2}" y="${H - 4}" text-anchor="middle">${xl}</text>`;
    if (yl) s += `<text class="lab" transform="translate(12 ${(m.t + H - m.b) / 2}) rotate(-90)" text-anchor="middle">${yl}</text>`;
    return s;
  }

  /* Horizontal bars, one per sample, coloured by status, with threshold marks. */
  function bars(el, items, o = {}) {
    const W = Math.max(360, Math.min(900, el.clientWidth || 520)), rowH = 22, m = { l: 120, r: 50, t: 8, b: 30 }, H = m.t + m.b + rowH * items.length;
    const vals = items.map((d) => d.value).filter((v) => isFinite(v));
    let lo = o.min ?? 0, hi = o.max ?? Math.max(...vals, ...(o.lines || []).map((l) => l.v)) * 1.08;
    const sx = (v) => m.l + (W - m.l - m.r) * (Math.max(lo, Math.min(hi, v)) - lo) / (hi - lo || 1);
    let s = svgEl(W, H);
    const xt = ticks(lo, hi, 5);
    for (const v of xt) s += `<line class="gl" x1="${sx(v)}" x2="${sx(v)}" y1="${m.t}" y2="${H - m.b}"/><text x="${sx(v)}" y="${H - m.b + 14}" text-anchor="middle">${(o.xf || fmt)(v)}</text>`;
    items.forEach((d, i) => {
      const y = m.t + i * rowH + 4, col = d.status === "pass" ? css("--good") : d.status === "warn" ? css("--warning") : d.status === "fail" ? css("--critical") : css("--ns");
      const w = Math.max(2, sx(d.value) - m.l);
      s += `<text x="${m.l - 8}" y="${y + 11}" text-anchor="end" class="lab">${esc(d.label.length > 18 ? d.label.slice(0, 17) + "…" : d.label)}</text>`;
      s += `<rect x="${m.l}" y="${y}" width="${w}" height="${rowH - 9}" rx="4" fill="${col}" opacity=".85" data-i="${i}"/>`;
      s += `<text x="${m.l + w + 5}" y="${y + 10.5}">${(o.vf || fmt)(d.value)}</text>`;
    });
    for (const l of o.lines || []) s += `<line x1="${sx(l.v)}" x2="${sx(l.v)}" y1="${m.t - 2}" y2="${H - m.b}" stroke="${css("--ink2")}" stroke-dasharray="3 3" opacity=".6"/><text x="${sx(l.v)}" y="${H - 3}" text-anchor="middle" class="lab">${esc(l.label)}</text>`;
    s += `</svg>`;
    el.innerHTML = s;
    el.querySelectorAll("rect[data-i]").forEach((r) => {
      const d = items[+r.dataset.i];
      r.onmousemove = (e) => tip(`<b>${esc(d.label)}</b><br>${esc(d.tip || (o.vf || fmt)(d.value))}`, e);
      r.onmouseleave = () => tip(null);
    });
  }

  /* Box plots of log2 counts per sample (quantiles 1,10,25,50,75,90,99). */
  function boxes(el, samples, q, refs) {
    const W = Math.max(360, Math.min(900, el.clientWidth || 520)), m = { l: 44, r: 10, t: 10, b: 62 }, H = 230;
    const all = q.flat(), lo = Math.floor(Math.min(...all)), hi = Math.ceil(Math.max(...all));
    const sy = (v) => H - m.b - (H - m.t - m.b) * (v - lo) / (hi - lo || 1);
    const bw = (W - m.l - m.r) / samples.length, sx = (i) => m.l + bw * (i + 0.5);
    let s = svgEl(W, H) + frame({ W, H, m, xt: [], yt: ticks(lo, hi, 5), sx, sy, yl: "log₂(count + 1)" });
    samples.forEach((n, i) => {
      const [p1, p10, p25, p50, p75, p90, p99] = q[i], x = sx(i), w = Math.min(26, bw * 0.55);
      const col = refs.includes(n) ? css("--accent2") : css("--accent");
      s += `<g data-i="${i}"><line x1="${x}" x2="${x}" y1="${sy(p10)}" y2="${sy(p90)}" stroke="${col}" stroke-width="1.5"/>`;
      s += `<rect x="${x - w / 2}" y="${sy(p75)}" width="${w}" height="${Math.max(1, sy(p25) - sy(p75))}" rx="4" fill="${col}" fill-opacity=".35" stroke="${col}"/>`;
      s += `<line x1="${x - w / 2}" x2="${x + w / 2}" y1="${sy(p50)}" y2="${sy(p50)}" stroke="${css("--ink")}" stroke-width="2"/>`;
      s += `<rect x="${x - bw / 2}" y="${m.t}" width="${bw}" height="${H - m.t - m.b}" fill="transparent"/></g>`;
      s += `<text transform="translate(${x + 3} ${H - m.b + 10}) rotate(40)" class="lab">${esc(n.length > 14 ? n.slice(0, 13) + "…" : n)}</text>`;
    });
    el.innerHTML = s + "</svg>";
    el.querySelectorAll("g[data-i]").forEach((g) => {
      const i = +g.dataset.i, v = q[i];
      g.onmousemove = (e) => tip(`<b>${esc(samples[i])}</b>${refs.includes(samples[i]) ? " (reference)" : ""}<br>median ${fmt(2 ** v[3] - 1)} reads<br>middle 50%: ${fmt(2 ** v[2] - 1)}–${fmt(2 ** v[4] - 1)}`, e);
      g.onmouseleave = () => tip(null);
    });
  }

  /* Line chart (ROC / PR). series: [{name, x, y, color}] */
  function lines(el, series, o = {}) {
    const W = Math.max(260, Math.min(520, el.clientWidth || 300)), H = 250, m = { l: 42, r: 10, t: 10, b: 38 };
    const sx = (v) => m.l + (W - m.l - m.r) * v, sy = (v) => H - m.b - (H - m.t - m.b) * v;
    let s = svgEl(W, H) + frame({ W, H, m, xt: [0, 0.25, 0.5, 0.75, 1], yt: [0, 0.25, 0.5, 0.75, 1], sx, sy, xl: o.xl, yl: o.yl });
    if (o.diag) s += `<line x1="${sx(0)}" y1="${sy(0)}" x2="${sx(1)}" y2="${sy(1)}" stroke="${css("--axis")}" stroke-dasharray="3 3"/>`;
    series.forEach((se, k) => {
      const d = se.x.map((x, i) => `${i ? "L" : "M"}${sx(x).toFixed(1)},${sy(se.y[i]).toFixed(1)}`).join("");
      s += `<path d="${d}" fill="none" stroke="${se.color}" stroke-width="2" stroke-linejoin="round"/>`;
    });
    s += `<rect x="${m.l}" y="${m.t}" width="${W - m.l - m.r}" height="${H - m.t - m.b}" fill="transparent" class="hov"/>`;
    el.innerHTML = s + "</svg>";
    const r = el.querySelector(".hov");
    r.onmousemove = (e) => {
      const bb = el.querySelector("svg").getBoundingClientRect(), xv = ((e.clientX - bb.left) * W / bb.width - m.l) / (W - m.l - m.r);
      const rows = series.map((se) => {
        let best = 0;
        se.x.forEach((x, i) => { if (Math.abs(x - xv) < Math.abs(se.x[best] - xv)) best = i; });
        return `<i style="display:inline-block;width:8px;height:8px;border-radius:50%;background:${se.color}"></i> ${esc(se.name)}: ${o.xl} ${fmt(se.x[best])}, ${o.yl} ${fmt(se.y[best])}`;
      });
      tip(rows.join("<br>"), e);
    };
    r.onmouseleave = () => tip(null);
  }

  /* Overlaid density histograms (essential vs non-essential gene LFC). */
  function hist2(el, h, names) {
    const W = Math.max(260, Math.min(520, el.clientWidth || 300)), H = 250, m = { l: 36, r: 10, t: 10, b: 38 };
    const e = h.edges, n = e.length - 1, mx = Math.max(...h.ess, ...h.non);
    const sx = (v) => m.l + (W - m.l - m.r) * (v - e[0]) / (e[n] - e[0]), sy = (v) => H - m.b - (H - m.t - m.b) * v / mx;
    let s = svgEl(W, H) + frame({ W, H, m, xt: ticks(e[0], e[n], 4), yt: [], sx, sy, xl: "gene log₂ fold change" });
    const path = (arr) => "M" + sx(e[0]) + "," + sy(0) + arr.map((v, i) => `L${sx(e[i])},${sy(v)}L${sx(e[i + 1])},${sy(v)}`).join("") + `L${sx(e[n])},${sy(0)}Z`;
    s += `<path d="${path(h.non)}" fill="${css("--non")}" fill-opacity=".28" stroke="${css("--non")}" stroke-width="1.5"/>`;
    s += `<path d="${path(h.ess)}" fill="${css("--ess")}" fill-opacity=".28" stroke="${css("--ess")}" stroke-width="1.5"/>`;
    el.innerHTML = s + "</svg>" + `<div class="legend"><span><i style="background:${css("--ess")}"></i>${names[0]}</span><span><i style="background:${css("--non")}"></i>${names[1]}</span></div>`;
  }

  /* Labelled scatter in SVG (PCA). */
  function pca(el, p, condOf, refs) {
    const W = Math.max(360, Math.min(1000, el.clientWidth || 520)), H = 280, m = { l: 46, r: 16, t: 12, b: 38 };
    const xs = p.pc1, ys = p.pc2;
    const pad = (a) => { const lo = Math.min(...a), hi = Math.max(...a), d = (hi - lo) * 0.15 || 1; return [lo - d, hi + d]; };
    const [x0, x1] = pad(xs), [y0, y1] = pad(ys);
    const sx = (v) => m.l + (W - m.l - m.r) * (v - x0) / (x1 - x0), sy = (v) => H - m.b - (H - m.t - m.b) * (v - y0) / (y1 - y0);
    const conds = [...new Set(p.samples.map((s) => condOf[s] || s))];
    const pal = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"];
    const pald = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300", "#9085e9", "#e66767"];
    const P = document.documentElement.dataset.theme === "dark" ? pald : pal;
    const col = (s) => conds.length > 8 ? css("--accent") : P[conds.indexOf(condOf[s] || s) % 8];
    let s = svgEl(W, H) + frame({ W, H, m, xt: ticks(x0, x1, 5), yt: ticks(y0, y1, 4), sx, sy,
      xl: `PC1 (${Math.round(p.explained[0] * 100)}%)`, yl: `PC2 (${Math.round((p.explained[1] || 0) * 100)}%)` });
    p.samples.forEach((n, i) => {
      s += `<circle cx="${sx(xs[i])}" cy="${sy(ys[i])}" r="${refs.includes(n) ? 7 : 6}" fill="${col(n)}" stroke="${css("--surface")}" stroke-width="2" data-i="${i}" ${refs.includes(n) ? 'stroke-dasharray="2 1.5"' : ""}/>`;
    });
    // one label per condition, at the centre of its replicates (sample names are in the tooltip)
    conds.forEach((c) => {
      const ii = p.samples.map((n, i) => ((condOf[n] || n) === c ? i : -1)).filter((i) => i >= 0);
      const cx = ii.reduce((a, i) => a + sx(xs[i]), 0) / ii.length, cy = Math.max(m.t + 22, Math.min(...ii.map((i) => sy(ys[i]))));
      s += `<text x="${cx}" y="${cy - 11}" text-anchor="middle" class="ttl">${esc(c)}</text>`;
    });
    el.innerHTML = s + "</svg>" + (conds.length <= 8 ? `<div class="legend">${conds.map((c) => `<span><i style="background:${P[conds.indexOf(c)]}"></i>${esc(c)}</span>`).join("")}</div>` : "");
    el.querySelectorAll("circle").forEach((c) => {
      const i = +c.dataset.i;
      c.onmousemove = (e) => tip(`<b>${esc(p.samples[i])}</b><br>condition ${esc(condOf[p.samples[i]] || p.samples[i])}<br>PC1 ${fmt(xs[i])}, PC2 ${fmt(ys[i])}`, e);
      c.onmouseleave = () => tip(null);
    });
  }

  /* Canvas scatter for thousands of points, with hover + click. pts: {x, y, k (colour key), label, i} */
  function scatter(el, pts, o) {
    el.innerHTML = "";
    const cv = document.createElement("canvas");
    el.appendChild(cv);
    const W = el.clientWidth || 600, H = o.height || 340, dpr = Math.min(2, devicePixelRatio || 1);
    cv.width = W * dpr; cv.height = H * dpr; cv.style.height = H + "px";
    const ctx = cv.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    const m = { l: 52, r: 14, t: 12, b: 40 };
    const xs = pts.map((p) => p.x).filter(isFinite), ys = pts.map((p) => p.y).filter(isFinite);
    let [x0, x1] = o.xr || [Math.min(...xs), Math.max(...xs)], [y0, y1] = o.yr || [Math.min(0, ...ys), Math.max(...ys)];
    if (x0 === x1) { x0 -= 1; x1 += 1; }
    if (y0 === y1) { y1 += 1; }
    const px = (x1 - x0) * 0.04, py = (y1 - y0) * 0.05;
    x0 -= px; x1 += px; y1 += py;
    const sx = (v) => m.l + (W - m.l - m.r) * (v - x0) / (x1 - x0), sy = (v) => H - m.b - (H - m.t - m.b) * (v - y0) / (y1 - y0);
    const txt = css("--muted"), grid = css("--grid"), axis = css("--axis");
    ctx.font = "11px -apple-system, sans-serif";
    ctx.fillStyle = txt; ctx.strokeStyle = grid; ctx.lineWidth = 1;
    for (const v of ticks(y0, y1, 5)) { ctx.beginPath(); ctx.moveTo(m.l, sy(v) + .5); ctx.lineTo(W - m.r, sy(v) + .5); ctx.stroke(); ctx.textAlign = "right"; ctx.fillText(fmt(v), m.l - 6, sy(v) + 4); }
    ctx.textAlign = "center";
    for (const v of ticks(x0, x1, 7)) ctx.fillText(fmt(v), sx(v), H - m.b + 15);
    ctx.strokeStyle = axis; ctx.beginPath(); ctx.moveTo(m.l, H - m.b + .5); ctx.lineTo(W - m.r, H - m.b + .5); ctx.stroke();
    if (o.vline !== undefined) { ctx.setLineDash([3, 3]); ctx.beginPath(); ctx.moveTo(sx(o.vline), m.t); ctx.lineTo(sx(o.vline), H - m.b); ctx.stroke(); ctx.setLineDash([]); }
    if (o.hline !== undefined && o.hline > y0 && o.hline < y1) { ctx.setLineDash([3, 3]); ctx.beginPath(); ctx.moveTo(m.l, sy(o.hline)); ctx.lineTo(W - m.r, sy(o.hline)); ctx.stroke(); ctx.setLineDash([]); }
    ctx.fillStyle = css("--ink2");
    ctx.fillText(o.xl || "", (m.l + W - m.r) / 2, H - 6);
    ctx.save(); ctx.translate(13, (m.t + H - m.b) / 2); ctx.rotate(-Math.PI / 2); ctx.fillText(o.yl || "", 0, 0); ctx.restore();
    const colors = o.colors;
    const order = [...pts].sort((a, b) => (o.z ? o.z(a) - o.z(b) : 0));
    const ring = css("--surface");
    for (const p of order) {
      if (!isFinite(p.x) || !isFinite(p.y)) continue;
      if (o.yr && (p.y < o.yr[0] || p.y > o.yr[1] + py)) continue;      // outside the chosen range (named in the axis title)
      const r = p.k === "ns" ? 2.2 : 3.6;
      ctx.beginPath(); ctx.arc(sx(p.x), sy(p.y), r, 0, 6.283);
      ctx.fillStyle = colors[p.k] || colors.ns;
      ctx.globalAlpha = p.k === "ns" ? 0.55 : 0.92;
      ctx.fill();
      if (p.k !== "ns") { ctx.globalAlpha = 1; ctx.lineWidth = 1; ctx.strokeStyle = ring; ctx.stroke(); }
    }
    ctx.globalAlpha = 1;
    // labels for top points
    ctx.fillStyle = css("--ink"); ctx.font = "600 11px -apple-system, sans-serif"; ctx.textAlign = "left";
    // labels: try right, left, above, below; skip a label that would overlap another or leave the plot
    const placed = [];
    const overlaps = (b) => placed.some((q) => b.x < q.x + q.w && q.x < b.x + b.w && b.y < q.y + q.h && q.y < b.y + b.h);
    for (const p of o.labels || []) {
      if (!isFinite(p.x) || !isFinite(p.y)) continue;
      const w = ctx.measureText(p.label).width + 4, h = 13, X = sx(p.x), Y = sy(p.y);
      const cands = [[X + 6, Y - h / 2], [X - w - 6, Y - h / 2], [X - w / 2, Y - h - 6], [X - w / 2, Y + 6]];
      for (const [bx, by] of cands) {
        const b = { x: bx, y: by, w, h };
        if (bx < m.l || bx + w > W - m.r || by < 0 || by + h > H - m.b || overlaps(b)) continue;
        placed.push(b);
        ctx.fillText(p.label, bx + 2, by + 10);
        break;
      }
    }
    const sel = o.selected && pts.find((p) => p.label === o.selected);
    if (sel) { ctx.strokeStyle = css("--ink"); ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(sx(sel.x), sy(sel.y), 7, 0, 6.283); ctx.stroke(); }
    const near = (e) => {
      const b = cv.getBoundingClientRect(), x = e.clientX - b.left, y = e.clientY - b.top;
      let best = null, bd = 64;
      for (const p of pts) { const d = (sx(p.x) - x) ** 2 + (sy(p.y) - y) ** 2; if (d < bd) { bd = d; best = p; } }
      return best;
    };
    cv.onmousemove = (e) => { const p = near(e); cv.style.cursor = p ? "pointer" : "default"; tip(p ? o.tip(p) : null, e); };
    cv.onmouseleave = () => tip(null);
    cv.onclick = (e) => { const p = near(e); if (p && o.click) o.click(p); };
  }

  /* Per-guide bars of log fold change (gene view). */
  function guideBars(el, rows) {
    const W = Math.max(360, Math.min(900, el.clientWidth || 520)), rowH = 24, m = { l: 150, r: 46, t: 8, b: 30 }, H = m.t + m.b + rowH * rows.length;
    const v = rows.map((r) => r.lfc), lim = Math.max(1, ...v.map(Math.abs)) * 1.1;
    const sx = (x) => m.l + (W - m.l - m.r) * (x + lim) / (2 * lim);
    let s = svgEl(W, H);
    for (const t of ticks(-lim, lim, 6)) s += `<line class="gl" x1="${sx(t)}" x2="${sx(t)}" y1="${m.t}" y2="${H - m.b}"/><text x="${sx(t)}" y="${H - m.b + 14}" text-anchor="middle">${fmt(t)}</text>`;
    s += `<line class="ax" x1="${sx(0)}" x2="${sx(0)}" y1="${m.t}" y2="${H - m.b}"/>`;
    rows.forEach((r, i) => {
      const y = m.t + i * rowH + 4, x = sx(Math.min(0, r.lfc)), w = Math.max(2, Math.abs(sx(r.lfc) - sx(0)));
      const col = !r.valid ? css("--ns") : r.lfc < 0 ? css("--neg") : css("--pos");
      s += `<text x="${m.l - 8}" y="${y + 11}" text-anchor="end" class="lab" style="font-family:var(--mono);font-size:10.5px">${esc(r.guide.length > 22 ? r.guide.slice(0, 21) + "…" : r.guide)}</text>`;
      s += `<rect x="${x}" y="${y}" width="${w}" height="${rowH - 10}" rx="4" fill="${col}" data-i="${i}"/>`;
    });
    s += `<text class="lab" x="${(m.l + W - m.r) / 2}" y="${H - 3}" text-anchor="middle">guide log₂ fold change</text>`;
    el.innerHTML = s + "</svg>";
    el.querySelectorAll("rect[data-i]").forEach((b) => {
      const r = rows[+b.dataset.i];
      b.onmousemove = (e) => tip(`<b>${esc(r.guide)}</b><br>LFC ${fmt(r.lfc)} · p(low) ${fmt(r.p_low)} · p(high) ${fmt(r.p_high)}<br>control ${fmt(r.control_mean)} → treatment ${fmt(r.treat_mean)} (normalised)${r.valid ? "" : "<br>left out of the gene test (zero counts)"}`, e);
      b.onmouseleave = () => tip(null);
    });
  }

  /* Plate heatmap. M: rows x cols (null = no value). info: [{r, c, well, content, gene, hit}].
     mode "seq" (raw values, one-hue ramp) or "div" (scores, blue-grey-red around 0, clipped at ±lim). */
  function plate(el, M, info, o = {}) {
    const nr = M.length, nc = M[0].length;
    const W = Math.max(360, Math.min(980, el.clientWidth || 600)), m = { l: 24, r: 8, t: 18, b: 8 };
    const cell = Math.min((W - m.l - m.r) / nc, 30), H = m.t + m.b + cell * nr;
    const vals = M.flat().filter((v) => v !== null && isFinite(v));
    const dk = document.documentElement.dataset.theme === "dark";
    const sorted = [...vals].sort((a, b) => a - b), q = (f) => sorted[Math.min(sorted.length - 1, Math.max(0, Math.round(f * (sorted.length - 1))))];
    const lo = o.lo ?? q(0.02), hi = o.hi ?? q(0.98), lim = o.lim ?? 4;   // robust range: one extreme well should not wash out the plate
    const mix = (a, b, t) => a.map((x, i) => Math.round(x + (b[i] - x) * t));
    const hex = (h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16));
    const BLUE = hex(dk ? "#3987e5" : "#2a78d6"), RED = hex(dk ? "#e66767" : "#e34948"), MID = hex(dk ? "#383835" : "#f0efec");
    const SEQ0 = hex(dk ? "#104281" : "#cde2fb"), SEQ1 = hex(dk ? "#9ec5f4" : "#104281");
    const color = (v) => {
      if (v === null || !isFinite(v)) return "transparent";
      let c;
      if (o.mode === "div") { const t = Math.max(-1, Math.min(1, v / lim)); c = t < 0 ? mix(MID, BLUE, -t) : mix(MID, RED, t); }
      else { const t = hi > lo ? (v - lo) / (hi - lo) : 0.5; c = mix(SEQ0, SEQ1, Math.max(0, Math.min(1, t))); }
      return `rgb(${c.join(",")})`;
    };
    const at = {};
    (info || []).forEach((d) => (at[d.r + ":" + d.c] = d));
    let s = svgEl(W, H);
    for (let c = 0; c < nc; c++) if (nc <= 24 || c % 2 === 0) s += `<text x="${m.l + cell * (c + 0.5)}" y="${m.t - 5}" text-anchor="middle" style="font-size:9.5px">${c + 1}</text>`;
    for (let r = 0; r < nr; r++) s += `<text x="${m.l - 6}" y="${m.t + cell * (r + 0.5) + 3.5}" text-anchor="end" style="font-size:9.5px">${String.fromCharCode(65 + (r % 26))}</text>`;
    for (let r = 0; r < nr; r++) for (let c = 0; c < nc; c++) {
      const v = M[r][c], d = at[r + ":" + c] || {};
      const x = m.l + cell * c, y = m.t + cell * r;
      const fill = o.cat ? (o.cat[d.role] || "transparent") : color(v);
      s += `<rect x="${x + 0.75}" y="${y + 0.75}" width="${cell - 1.5}" height="${cell - 1.5}" rx="${Math.min(4, cell / 5)}" fill="${fill}" stroke="${v === null ? css("--line2") : "none"}" data-k="${r}:${c}"/>`;
      if (d.content === "pos" || d.content === "neg") s += `<path d="M${x + cell / 2},${y + cell * 0.22} L${x + cell * 0.78},${y + cell / 2} L${x + cell / 2},${y + cell * 0.78} L${x + cell * 0.22},${y + cell / 2}Z" fill="${d.content === "pos" ? css("--ink") : "none"}" stroke="${css("--ink")}" stroke-width="1.2" pointer-events="none"/>`;
      if (d.hit) s += `<circle cx="${x + cell / 2}" cy="${y + cell / 2}" r="${cell * 0.36}" fill="none" stroke="${css("--ink")}" stroke-width="2" pointer-events="none"/>`;
    }
    el.innerHTML = s + "</svg>" + (o.mode === "div"
      ? `<div class="legend"><span><i style="background:${color(-lim)}"></i>≤ −${lim}</span><span><i style="background:${color(0)};border:1px solid var(--line2)"></i>0</span><span><i style="background:${color(lim)}"></i>≥ ${lim}</span><span class="faint">${esc(o.unit || "")}</span></div>`
      : `<div class="legend"><span><i style="background:${color(lo)}"></i>${fmt(lo)}</span><span><i style="background:${color(hi)}"></i>${fmt(hi)}</span><span class="faint">${esc(o.unit || "")}</span></div>`);
    el.querySelectorAll("rect[data-k]").forEach((rc) => {
      const [r, c] = rc.dataset.k.split(":").map(Number), d = at[r + ":" + c] || {};
      rc.onmousemove = (e) => tip(`<b>${esc(d.well || String.fromCharCode(65 + r) + String(c + 1).padStart(2, "0"))}</b> · ${esc(d.content || "")}${d.gene ? " · " + esc(d.gene) : ""}${d.reagent ? " · " + esc(d.reagent) : ""}<br>${esc(o.unit || "value")}: ${fmt(M[r][c])}${d.hit ? "<br><b>hit</b>" : ""}`, e);
      rc.onmouseleave = () => tip(null);
    });
  }

  /* Simple line profile (row or column medians). */
  function profile(el, vals, labels, o = {}) {
    const W = Math.max(300, Math.min(700, el.clientWidth || 400)), H = 170, m = { l: 44, r: 10, t: 10, b: 30 };
    const v = vals.map((x) => (x === null ? NaN : x)), ok = v.filter(isFinite);
    const lo = Math.min(...ok, 1) - 0.02, hi = Math.max(...ok, 1) + 0.02;
    const sx = (i) => m.l + (W - m.l - m.r) * (v.length > 1 ? i / (v.length - 1) : 0.5), sy = (y) => H - m.b - (H - m.t - m.b) * (y - lo) / (hi - lo);
    let s = svgEl(W, H) + frame({ W, H, m, xt: [], yt: ticks(lo, hi, 4), sx, sy, xl: o.xl });
    labels.forEach((l, i) => { if (labels.length <= 24 || i % 2 === 0) s += `<text x="${sx(i)}" y="${H - m.b + 13}" text-anchor="middle" style="font-size:9.5px">${l}</text>`; });
    s += `<line x1="${m.l}" x2="${W - m.r}" y1="${sy(1)}" y2="${sy(1)}" stroke="${css("--axis")}" stroke-dasharray="3 3"/>`;
    s += `<path d="${v.map((y, i) => (isFinite(y) ? `${i ? "L" : "M"}${sx(i)},${sy(y)}` : "")).join("")}" fill="none" stroke="${css("--accent")}" stroke-width="2"/>`;
    v.forEach((y, i) => { if (isFinite(y)) s += `<circle cx="${sx(i)}" cy="${sy(y)}" r="3.5" fill="${css("--accent")}" stroke="${css("--surface")}" stroke-width="1.5" data-i="${i}"/>`; });
    el.innerHTML = s + "</svg>";
    el.querySelectorAll("circle").forEach((c) => {
      const i = +c.dataset.i;
      c.onmousemove = (e) => tip(`${o.what || ""} <b>${labels[i]}</b>: ${(v[i] * 100).toFixed(1)}% of plate median`, e);
      c.onmouseleave = () => tip(null);
    });
  }

  /* Histogram with threshold lines. */
  function histBars(el, edges, counts, o = {}) {
    const W = Math.max(300, Math.min(700, el.clientWidth || 400)), H = 210, m = { l: 48, r: 10, t: 10, b: 34 };
    const x0 = edges[0], x1 = edges[edges.length - 1], mx = Math.max(...counts, 1);
    const sx = (x) => m.l + (W - m.l - m.r) * (x - x0) / (x1 - x0), sy = (y) => H - m.b - (H - m.t - m.b) * (o.log ? Math.log10(y + 1) / Math.log10(mx + 1) : y / mx);
    let s = svgEl(W, H) + frame({ W, H, m, xt: ticks(x0, x1, 6), yt: [], sx, sy, xl: o.xl });
    counts.forEach((c, i) => {
      const hit = o.lines && o.lines.some((l) => (l.dir === "lo" ? edges[i + 1] <= l.v : edges[i] >= l.v));
      s += `<rect x="${sx(edges[i]) + 0.5}" y="${sy(c)}" width="${Math.max(1, sx(edges[i + 1]) - sx(edges[i]) - 1)}" height="${H - m.b - sy(c)}" rx="1.5" fill="${hit ? css("--pos") : css("--accent")}" opacity=".8" data-i="${i}"/>`;
    });
    (o.lines || []).forEach((l) => (s += `<line x1="${sx(l.v)}" x2="${sx(l.v)}" y1="${m.t}" y2="${H - m.b}" stroke="${css("--ink2")}" stroke-dasharray="3 3"/><text x="${sx(l.v) + 4}" y="${m.t + 10}" class="lab">${esc(l.label)}</text>`));
    el.innerHTML = s + "</svg>";
    el.querySelectorAll("rect[data-i]").forEach((r) => {
      const i = +r.dataset.i;
      r.onmousemove = (e) => tip(o.unit ? `${esc(o.xl)} ${fmt((edges[i] + edges[i + 1]) / 2)}: <b>${counts[i].toLocaleString()}</b> ${esc(o.unit)}` : `${fmt(edges[i])} to ${fmt(edges[i + 1])}: <b>${counts[i].toLocaleString()}</b> wells`, e);
      r.onmouseleave = () => tip(null);
    });
  }

  return { css, tip, ticks, fmt, esc, bars, boxes, lines, hist2, pca, scatter, guideBars, plate, profile, histBars };
})();
