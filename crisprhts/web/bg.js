/* Background: sequencing reads drift across the window the way they come off the sequencer. In each read the
   20-nt guide (teal) sits in front of the Cas9 scaffold (blue). Now and then a read is matched to its guide and
   adds one count to the guide-count histogram along the bottom, which slowly fills into the skewed, roughly
   log-normal shape real libraries have. Kept faint, under a veil, so it never competes with the text. */
(function () {
  const cv = document.getElementById("bg");
  if (!cv) return;
  const ctx = cv.getContext("2d");
  const reduce = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  const SCAF = "GTTTTAGAGCTAGAAATAGCAAGTT";
  const NB = "ACGT";
  let W = 0, H = 0, dpr = 1, last = 0, reads = [], bins = [], flyers = [];
  const rnd = (n) => Math.floor(Math.random() * n);
  const seq = (n) => Array.from({ length: n }, () => NB[rnd(4)]).join("");
  const dark = () => document.documentElement.getAttribute("data-theme") === "dark";
  // the Contrast style hides the canvas: stop animating rather than draw behind display:none
  const off = () => document.documentElement.getAttribute("data-style") === "contrast";
  // the guide, the scaffold and the unmatched bases take their colour from whichever style is on
  let cG = "63,208,192", cS = "57,135,229", cO = "160,180,176";
  function readColours() {
    const cs = getComputedStyle(document.documentElement);
    const rgb = (name, fallback) => {
      const v = (cs.getPropertyValue(name) || "").trim();
      const m = /^#([0-9a-f]{3}|[0-9a-f]{6})$/i.exec(v);
      if (!m) return fallback;
      const h = m[1].length === 3 ? m[1].replace(/./g, (c) => c + c) : m[1];
      return `${parseInt(h.slice(0, 2), 16)},${parseInt(h.slice(2, 4), 16)},${parseInt(h.slice(4, 6), 16)}`;
    };
    cG = rgb("--accent", cG); cS = rgb("--accent2", cS); cO = rgb("--faint", cO);
  }
  const NBIN = 64;

  function newRead(anywhere) {
    const lead = rnd(3) === 0 ? seq(2 + rnd(6)) : "";
    const z = 0.35 + Math.random() * 0.65;          // depth: far reads are smaller, slower and fainter
    const text = lead + seq(20) + SCAF.slice(0, 10 + rnd(14));
    const fs = 9 + z * 6;
    return { text, lead: lead.length, z, fs, y: 40 + Math.random() * (H - 140), x: anywhere ? Math.random() * W : -text.length * fs * 0.62,
             v: 8 + z * 22, mapped: Math.random() < 0.82, hit: false };
  }
  function resize() {
    dpr = Math.min(2, window.devicePixelRatio || 1);
    W = window.innerWidth; H = window.innerHeight;
    cv.width = W * dpr; cv.height = H * dpr;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    const n = Math.round(Math.min(34, W * H / 42000));
    reads = Array.from({ length: n }, () => newRead(true));
    bins = Array.from({ length: NBIN }, (_, i) => 0);
    seedBins(400);
  }
  function seedBins(n) { for (let k = 0; k < n; k++) addCount(); }
  function addCount() {
    // log-normal guide abundance -> bin index
    const u1 = Math.random() || 1e-9, u2 = Math.random();
    const z = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
    const b = Math.max(0, Math.min(NBIN - 1, Math.round(NBIN * 0.42 + z * NBIN * 0.14)));
    bins[b] += 1;
    if (Math.max(...bins) > 260) bins = bins.map((x) => x * 0.6);
  }

  function draw(now, dt) {
    const dk = dark();
    ctx.clearRect(0, 0, W, H);
    for (const r of reads) {
      r.x += r.v * dt;
      const cw = r.fs * 0.62;
      ctx.font = `${r.fs}px ui-monospace, SFMono-Regular, Menlo, monospace`;
      const a = (dk ? 0.38 : 0.30) * r.z;
      for (let i = 0; i < r.text.length; i++) {
        const x = r.x + i * cw;
        if (x < -20 || x > W + 20) continue;
        const inG = i >= r.lead && i < r.lead + 20, inS = i >= r.lead + 20;
        const c = inG ? cG : inS ? cS : cO;
        ctx.fillStyle = `rgba(${c},${inG ? a * 1.25 : a * 0.8})`;
        ctx.fillText(r.text[i], x, r.y);
      }
      // a read is "matched" when its guide crosses the middle third of the window
      const gx = r.x + (r.lead + 10) * cw;
      if (!r.hit && r.mapped && gx > W * 0.45) {
        r.hit = true;
        flyers.push({ x: gx, y: r.y - r.fs * 0.35, tx: 0, ty: 0, t: 0, z: r.z });
      }
      if (r.x > W + 40) Object.assign(r, newRead(false));
    }
    // flyers: one count travels down to the histogram
    const bw = Math.min(14, W / NBIN), x0 = (W - bw * NBIN) / 2, base = H - 18;
    for (const f of flyers) {
      if (!f.b) { addCount(); f.b = 1; f.tx = x0 + bw * (NBIN * 0.42 + (Math.random() - 0.5) * NBIN * 0.4); f.ty = base - 10; }
      f.t += dt * 0.9;
      const e = Math.min(1, f.t), ix = f.x + (f.tx - f.x) * e, iy = f.y + (f.ty - f.y) * e * e;
      ctx.fillStyle = `rgba(${cG},${0.55 * (1 - e) * f.z})`;
      ctx.beginPath(); ctx.arc(ix, iy, 2.2, 0, 6.283); ctx.fill();
    }
    flyers = flyers.filter((f) => f.t < 1);
    // guide-count histogram
    const mx = Math.max(1, ...bins), hmax = Math.min(110, H * 0.14);
    for (let i = 0; i < NBIN; i++) {
      const h = (bins[i] / mx) * hmax;
      ctx.fillStyle = `rgba(${cG},${dk ? 0.16 : 0.12})`;
      ctx.fillRect(x0 + i * bw + 1, base - h, bw - 2, h);
    }
    ctx.fillStyle = `rgba(${cO},${dk ? 0.35 : 0.3})`;
    ctx.fillRect(x0, base, bw * NBIN, 1);
    ctx.font = "10px ui-monospace, Menlo, monospace";
    ctx.fillText("log₂ guide count", x0, base + 13);
  }

  function frame(now) {
    requestAnimationFrame(frame);
    if (document.hidden || off() || now - last < 33) return;
    const dt = Math.min(0.1, (now - (last || now)) / 1000);
    last = now;
    draw(now, dt);
  }
  window.addEventListener("resize", resize);
  window.addEventListener("ch-style", readColours);     // the style or theme changed
  readColours();
  resize();
  if (reduce) draw(0, 0); else requestAnimationFrame(frame);
})();
