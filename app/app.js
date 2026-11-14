// Fire Activity Explorer: harmonized MODIS + VIIRS hotspots -> region/date analytics + downloads.
// Input: data/hotspots.csv.gz + data/meta.json (written by ../hotspots.py), or any CSV with the same columns.
'use strict';

const $ = id => document.getElementById(id);
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const SENS = ['MODIS', 'VIIRS'], SKEY = ['--modis', '--viirs'];
const DAY = 864e5;
const toDay = s => Date.UTC(+s.slice(0, 4), +s.slice(5, 7) - 1, +s.slice(8, 10)) / DAY;
const dayStr = d => new Date(d * DAY).toISOString().slice(0, 10);
const fmt = n => n == null || !isFinite(n) ? '–' : Math.abs(n) >= 100 ? Math.round(n).toLocaleString() : Math.abs(n) >= 10 ? n.toFixed(1) : n.toFixed(2);
const MON = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
const lonS = v => `${Math.abs(v).toFixed(2)}°${v < 0 ? 'W' : 'E'}`, latS = v => `${Math.abs(v).toFixed(2)}°${v < 0 ? 'S' : 'N'}`;
const boxS = () => `${latS(st.s)}–${latS(st.n)}, ${lonS(st.w)}–${lonS(st.e)}`;
const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

// Approximate boxes for common areas (lon/lat, degrees).
const PRESETS = [
  ['North Coast & Klamath', [-124.5, 39.5, -122.0, 42.0]], ['Sacramento Valley', [-122.6, 38.3, -121.2, 40.6]],
  ['Bay Area', [-123.1, 36.9, -121.5, 38.6]], ['Sierra Nevada', [-121.3, 35.5, -118.0, 40.2]],
  ['Central Coast', [-122.0, 34.4, -119.5, 37.0]], ['San Joaquin Valley', [-121.2, 34.9, -118.6, 37.9]],
  ['Southern California', [-119.5, 32.5, -114.5, 35.0]],
];

let H = null;          // hotspot columns
let META = null;       // dataset metadata
let SAMPLED = [[], []];// sorted sampled days per sensor
const st = { w: 0, s: 0, e: 0, n: 0, from: 0, to: 0, sens: [true, true], units: 'harm', conf: 0 };
let A = null;          // last analytics

// ---------------------------------------------------------------- loading
async function gunzipText(blobOrResp) {
  const body = blobOrResp.body || blobOrResp.stream();
  return await new Response(body.pipeThrough(new DecompressionStream('gzip'))).text();
}
async function readFile(file) {
  const head = new Uint8Array(await file.slice(0, 2).arrayBuffer());
  return head[0] === 0x1f && head[1] === 0x8b ? gunzipText(file) : file.text();
}

function parseCSV(text) {
  const lines = text.split(/\r?\n/).filter(l => l && !l.startsWith('#'));
  const cols = lines[0].split(',').map(c => c.trim().toLowerCase());
  const ix = n => cols.indexOf(n);
  const need = ['date', 'lat', 'lon', 'sensor'].filter(c => ix(c) < 0);
  if (need.length) throw new Error(`missing column(s): ${need.join(', ')}. Expected: ${'date,time,lat,lon,sensor,satellite,daynight,confidence,frp,pixels,modis_equiv'}`);
  const n = lines.length - 1;
  const h = { n, day: new Int32Array(n), time: new Uint16Array(n), lat: new Float32Array(n), lon: new Float32Array(n),
    sensor: new Uint8Array(n), sat: [], night: new Uint8Array(n), conf: new Uint8Array(n), frp: new Float32Array(n),
    pix: new Uint16Array(n), weq: new Float32Array(n), hasWeq: ix('modis_equiv') >= 0 };
  const c = { date: ix('date'), time: ix('time'), lat: ix('lat'), lon: ix('lon'), sensor: ix('sensor'), sat: ix('satellite'),
    dn: ix('daynight'), conf: ix('confidence'), frp: ix('frp'), pix: ix('pixels'), weq: ix('modis_equiv') };
  const CONF = { low: 0, l: 0, nominal: 1, n: 1, high: 2, h: 2 };
  const dayCache = new Map();
  let k = 0;
  for (let i = 1; i <= n; i++) {
    const f = lines[i].split(',');
    const ds = f[c.date]; let d = dayCache.get(ds);
    if (d === undefined) { d = toDay(ds); dayCache.set(ds, d); }
    const sen = /^v/i.test(f[c.sensor]) ? 1 : 0;
    h.day[k] = d; h.time[k] = c.time >= 0 ? +f[c.time] || 0 : 0;
    h.lat[k] = +f[c.lat]; h.lon[k] = +f[c.lon]; h.sensor[k] = sen;
    h.sat[k] = c.sat >= 0 ? f[c.sat] : '';
    h.night[k] = c.dn >= 0 && f[c.dn] === 'N' ? 1 : 0;
    const cv = c.conf >= 0 ? f[c.conf] : 'nominal';
    h.conf[k] = cv in CONF ? CONF[cv] : (+cv < 30 ? 0 : +cv < 80 ? 1 : 2);
    h.frp[k] = c.frp >= 0 ? +f[c.frp] || 0 : 0;
    h.pix[k] = c.pix >= 0 ? +f[c.pix] || 1 : 1;
    h.weq[k] = c.weq >= 0 ? +f[c.weq] : 1;
    if (isFinite(d) && isFinite(h.lat[k]) && isFinite(h.lon[k])) k++;
  }
  for (const key of ['day', 'time', 'lat', 'lon', 'sensor', 'night', 'conf', 'frp', 'pix', 'weq']) h[key] = h[key].subarray(0, k);
  h.n = k;
  return h;
}

function useData(h, meta, label) {
  H = h;
  let w = 180, s = 90, e = -180, n = -90;
  for (let i = 0; i < h.n; i++) { w = Math.min(w, h.lon[i]); e = Math.max(e, h.lon[i]); s = Math.min(s, h.lat[i]); n = Math.max(n, h.lat[i]); }
  const days = [new Set(), new Set()];
  for (let i = 0; i < h.n; i++) days[h.sensor[i]].add(h.day[i]);
  META = meta || { name: label, bbox: [w, s, e, n], sampled_days: null, harmonization: null };
  SAMPLED = SENS.map((k, j) => META.sampled_days && META.sampled_days[k] ? META.sampled_days[k].map(toDay) : [...days[j]]).map(a => Int32Array.from(a).sort());
  $('region').textContent = META.name;
  const bb = META.bbox;
  PRESETS.splice(0, 0, ['Whole dataset', bb]);
  $('preset').innerHTML = PRESETS.map((p, i) => `<option value="${i}">${esc(p[0])}</option>`).join('') + '<option value="custom" hidden>Custom</option>';
  const dmin = Math.min(...SAMPLED.filter(a => a.length).map(a => a[0])), dmax = Math.max(...SAMPLED.filter(a => a.length).map(a => a[a.length - 1]));
  $('from').min = $('to').min = dayStr(dmin); $('from').max = $('to').max = dayStr(dmax);
  const y0 = +dayStr(dmin).slice(0, 4), y1 = +dayStr(dmax).slice(0, 4);
  const q = [['all', 'All dates'], ['last12', 'Last 12 months'], ['modis', 'MODIS-only era'], ['overlap', 'Both sensors (VIIRS era)']];
  for (let y = y1; y >= y0; y--) q.push(['y' + y, String(y)]);
  $('quick').innerHTML = '<option value="">Quick range…</option>' + q.map(([v, t]) => `<option value="${v}">${t}</option>`).join('');
  $('quick').dataset.min = dmin; $('quick').dataset.max = dmax;
  const hm = META.harmonization;
  $('banner').innerHTML = `${META.sampling ? `Sampled data: ${esc(META.sampling)}, so counts cover those days only and rates are <b>per sampled day</b>. ` : 'Sampled days are taken to be the days present in the file. '}` +
    (hm ? `VIIRS is in MODIS units: pixels merged to ${hm.merge_km} km per overpass, then calibrated (${esc(hm.model)}, typical monthly error ×/÷ ${hm.typical_error}). ${esc(hm.note)}` :
      h.hasWeq ? '' : '<span class="err">This file has no modis_equiv column, so VIIRS can only be shown as raw detections.</span>');
  $('foot').innerHTML = `${esc(META.source || label)} · ${h.n.toLocaleString()} hotspots` + (META.generated ? ` · built ${esc(META.generated)}` : '') +
    ` · Columns: date, time, lat, lon, sensor, satellite, daynight, confidence, frp, pixels, modis_equiv.`;
  if (!h.hasWeq) setUnits('raw');
  initMap();
  if (!readHash()) { setBox(bb); st.from = dmin; st.to = dmax; $('preset').value = '0'; }
  syncInputs(); update();
  $('loader').hidden = true;
}

async function boot() {
  try {
    const [r, m] = await Promise.all([fetch('data/hotspots.csv.gz'), fetch('data/meta.json')]);
    if (!r.ok) throw new Error(r.status);
    const meta = m.ok ? await m.json() : null;
    // Some static hosts serve .gz with Content-Encoding: gzip (already decoded); detect by peeking.
    const buf = new Uint8Array(await r.arrayBuffer());
    const text = buf[0] === 0x1f && buf[1] === 0x8b ? await gunzipText(new Blob([buf])) : new TextDecoder().decode(buf);
    useData(parseCSV(text), meta, 'data/hotspots.csv.gz');
  } catch (err) {
    $('ltitle').textContent = 'Load a hotspot file';
    $('lmsg').innerHTML = location.protocol === 'file:' ?
      'Browsers block opened-from-disk pages from reading files next to them. Pick <code>app/data/hotspots.csv.gz</code> below, or serve the folder (for example <code>python3 -m http.server</code> in the repo, then open <code>/app/</code>).' :
      `Couldn't load data/hotspots.csv.gz (${esc(err.message)}). Pick a harmonized hotspot CSV (or .csv.gz) instead.`;
    $('lpick').hidden = false;
  }
}
async function onFile(ev) {
  const f = ev.target.files[0]; if (!f) return;
  try {
    $('loader').hidden = false; $('ltitle').textContent = 'Reading ' + f.name + '…'; $('lpick').hidden = true; $('lmsg').textContent = '';
    const h = parseCSV(await readFile(f));
    if (!h.n) throw new Error('no rows');
    PRESETS.splice(0, PRESETS.length - 7 > 0 ? 1 : 0);
    useData(h, f.name === 'hotspots.csv.gz' && META && META.sampled_days ? META : null, f.name);
  } catch (err) {
    $('ltitle').textContent = 'Could not read ' + f.name; $('lmsg').innerHTML = `<span class="err">${esc(err.message)}</span>`; $('lpick').hidden = false;
  }
}
$('file').addEventListener('change', onFile); $('file2').addEventListener('change', onFile);

// ---------------------------------------------------------------- selection state
function setBox([w, s, e, n]) { Object.assign(st, { w: Math.min(w, e), e: Math.max(w, e), s: Math.min(s, n), n: Math.max(s, n) }); }
function syncInputs() {
  for (const k of ['w', 's', 'e', 'n']) $(k).value = (+st[k]).toFixed(2);
  $('from').value = dayStr(st.from); $('to').value = dayStr(st.to);
  $('sM').checked = st.sens[0]; $('sV').checked = st.sens[1]; $('conf').value = st.conf;
  [...$('units').children].forEach(b => b.setAttribute('aria-pressed', b.dataset.u === st.units));
}
function writeHash() {
  const p = new URLSearchParams({ bbox: [st.w, st.s, st.e, st.n].map(v => (+v).toFixed(3)).join(','), from: dayStr(st.from), to: dayStr(st.to),
    sensors: SENS.filter((_, j) => st.sens[j]).join(','), units: st.units, conf: st.conf });
  history.replaceState(null, '', '#' + p);
}
function readHash() {
  if (!location.hash) return false;
  const p = new URLSearchParams(location.hash.slice(1));
  const bb = (p.get('bbox') || '').split(',').map(Number);
  if (bb.length !== 4 || bb.some(v => !isFinite(v)) || !p.get('from') || !p.get('to')) return false;
  setBox(bb); st.from = toDay(p.get('from')); st.to = toDay(p.get('to'));
  const s = (p.get('sensors') || 'MODIS,VIIRS').split(','); st.sens = SENS.map(k => s.includes(k));
  st.units = p.get('units') === 'raw' || !H.hasWeq ? 'raw' : 'harm'; st.conf = +p.get('conf') || 0;
  $('preset').value = 'custom';
  return true;
}
function setUnits(u) { st.units = u; [...$('units').children].forEach(b => b.setAttribute('aria-pressed', b.dataset.u === u)); }

let timer;
const schedule = () => { clearTimeout(timer); timer = setTimeout(update, 120); };
$('preset').addEventListener('change', () => { const p = PRESETS[+$('preset').value]; if (p) { setBox(p[1]); syncInputs(); fitRegion(); schedule(); } });
for (const k of ['w', 's', 'e', 'n']) $(k).addEventListener('change', () => {
  const v = ['w', 's', 'e', 'n'].map(x => +$(x).value); if (v.some(x => !isFinite(x))) return;
  setBox(v); $('preset').value = 'custom'; syncInputs(); schedule();
});
for (const k of ['from', 'to']) $(k).addEventListener('change', () => {
  const f = $('from').value, t = $('to').value; if (!f || !t) return;
  st.from = Math.min(toDay(f), toDay(t)); st.to = Math.max(toDay(f), toDay(t)); $('quick').value = ''; syncInputs(); schedule();
});
$('quick').addEventListener('change', () => {
  const v = $('quick').value, dmin = +$('quick').dataset.min, dmax = +$('quick').dataset.max;
  const vStart = SAMPLED[1].length ? SAMPLED[1][0] : dmin;
  if (v === 'all') { st.from = dmin; st.to = dmax; }
  else if (v === 'last12') { st.to = dmax; st.from = dmax - 365; }
  else if (v === 'modis') { st.from = dmin; st.to = vStart - 1; }
  else if (v === 'overlap') { st.from = vStart; st.to = dmax; }
  else if (v[0] === 'y') { st.from = Math.max(dmin, toDay(v.slice(1) + '-01-01')); st.to = Math.min(dmax, toDay(v.slice(1) + '-12-31')); }
  syncInputs(); schedule();
});
$('sM').addEventListener('change', () => { st.sens[0] = $('sM').checked; schedule(); });
$('sV').addEventListener('change', () => { st.sens[1] = $('sV').checked; schedule(); });
$('conf').addEventListener('change', () => { st.conf = +$('conf').value; schedule(); });
$('units').addEventListener('click', ev => { const u = ev.target.dataset.u; if (!u || (u === 'harm' && !H.hasWeq)) return; setUnits(u); schedule(); });

// ---------------------------------------------------------------- analytics
const lowerBound = (a, x) => { let lo = 0, hi = a.length; while (lo < hi) { const m = (lo + hi) >> 1; if (a[m] < x) lo = m + 1; else hi = m; } return lo; };
const sampledIn = (j, from, to) => lowerBound(SAMPLED[j], to + 1) - lowerBound(SAMPLED[j], from);
const monthKey = d => dayStr(d).slice(0, 7);

function analyse() {
  const harm = st.units === 'harm';
  const idx = [];
  for (let i = 0; i < H.n; i++) {
    const d = H.day[i];
    if (d < st.from || d > st.to || !st.sens[H.sensor[i]] || H.conf[i] < st.conf) continue;
    const la = H.lat[i], lo = H.lon[i];
    if (la < st.s || la > st.n || lo < st.w || lo > st.e) continue;
    idx.push(i);
  }
  const wt = i => H.sensor[i] === 0 ? 1 : harm ? H.weq[i] : H.pix[i];
  const per = SENS.map((name, j) => ({ name, on: st.sens[j], hotspots: 0, count: 0, frp: 0, night: 0, high: 0, days: sampledIn(j, st.from, st.to), frps: [] }));
  const months = new Map();      // 'YYYY-MM' -> [countM, countV]
  for (const i of idx) {
    const j = H.sensor[i], p = per[j], w = wt(i);
    p.hotspots++; p.count += w; p.frp += H.frp[i]; p.frps.push(H.frp[i]);
    if (H.night[i]) p.night += w; if (H.conf[i] === 2) p.high += w;
    const mk = monthKey(H.day[i]); let m = months.get(mk); if (!m) months.set(mk, m = [0, 0]); m[j] += w;
  }
  for (const p of per) {
    p.rate = p.days ? p.count / p.days : null;
    p.frps.sort((a, b) => a - b); p.frpMedian = p.frps.length ? p.frps[p.frps.length >> 1] : null;
    p.frpP90 = p.frps.length ? p.frps[Math.floor(p.frps.length * .9)] : null;
  }
  // monthly series over the selected range, with sampled days per month (clipped to the range)
  const series = [];
  const sDays = SAMPLED.map(a => { const m = new Map(); for (const d of a) if (d >= st.from && d <= st.to) { const k = monthKey(d); m.set(k, (m.get(k) || 0) + 1); } return m; });
  let [y, mo] = dayStr(st.from).slice(0, 7).split('-').map(Number); const end = monthKey(st.to);
  for (;;) {
    const k = `${y}-${String(mo).padStart(2, '0')}`, c = months.get(k) || [0, 0];
    series.push({ month: k, t: y + (mo - .5) / 12, count: c, days: [sDays[0].get(k) || 0, sDays[1].get(k) || 0],
      rate: [0, 1].map(j => per[j].on && sDays[j].get(k) ? c[j] / sDays[j].get(k) : null) });
    if (k >= end) break; mo++; if (mo > 12) { mo = 1; y++; }
  }
  // seasonal cycle
  const season = MON.map((name, m) => {
    const r = [0, 1].map(j => { let c = 0, d = 0; for (const s of series) if (+s.month.slice(5) === m + 1) { c += s.count[j]; d += s.days[j]; } return per[j].on && d ? c / d : null; });
    return { name, rate: r };
  });
  // FRP distribution (log bins), share of hotspots
  const edges = [0, 5, 10, 20, 50, 100, 200, 500, 1000, Infinity];
  const frp = edges.slice(0, -1).map((lo, b) => ({ name: edges[b + 1] === Infinity ? `${lo}+` : `${lo}–${edges[b + 1]}`,
    share: per.map(p => p.on && p.frps.length ? (lowerBound(p.frps, edges[b + 1]) - lowerBound(p.frps, lo)) / p.frps.length : null) }));
  // top 0.1° cells and map grid
  const span = Math.max(st.e - st.w, st.n - st.s);
  const mc = span > 6 ? 0.1 : span > 2 ? 0.05 : span > 0.6 ? 0.02 : 0.01;
  const grid = new Map(), top = new Map();
  for (const i of idx) {
    const j = H.sensor[i], w = wt(i);
    const gk = Math.floor(H.lon[i] / mc) + ',' + Math.floor(H.lat[i] / mc); let g = grid.get(gk); if (!g) grid.set(gk, g = [0, 0]); g[j] += w;
    const tk = Math.floor(H.lon[i] * 10) + ',' + Math.floor(H.lat[i] * 10);
    let t = top.get(tk); if (!t) top.set(tk, t = { c: [0, 0], frp: 0, first: H.day[i], last: H.day[i] });
    t.c[j] += w; t.frp = Math.max(t.frp, H.frp[i]); t.first = Math.min(t.first, H.day[i]); t.last = Math.max(t.last, H.day[i]);
  }
  const rateOf = c => { const r = [0, 1].filter(j => per[j].on && per[j].days).map(j => c[j] / per[j].days); return r.length ? r : [0]; };
  const cells = [...grid].map(([k, c]) => { const [x, yy] = k.split(',').map(Number); const r = rateOf(c); return [x * mc, yy * mc, r.reduce((a, b) => a + b) / r.length]; }).filter(c => c[2] > 0);
  const tops = [...top].map(([k, t]) => { const [x, yy] = k.split(',').map(Number); return { lon: x / 10 + .05, lat: yy / 10 + .05, ...t, score: Math.max(...rateOf(t.c)) }; })
    .sort((a, b) => b.score - a.score).slice(0, 10);
  const peak = [0, 1].map(j => series.reduce((b, s) => s.rate[j] != null && (b == null || s.rate[j] > b.rate[j]) ? s : b, null));
  return { idx, per, series, season, frp, cells, mc, tops, peak, harm, unit: harm ? 'MODIS-equivalent hotspots' : 'raw detections' };
}

// ---------------------------------------------------------------- charts (SVG with explicit colours, so exports are self-contained)
const NS = 'http://www.w3.org/2000/svg';
const el = (n, a, txt) => { const e = document.createElementNS(NS, n); for (const k in a) e.setAttribute(k, a[k]); if (txt != null) e.textContent = txt; return e; };
const pal = () => ({ bg: css('--surface'), grid: css('--grid'), text: css('--text'), text2: css('--text-2'), text3: css('--text-3'), band: css('--band'), s: [css('--modis'), css('--viirs')], heat: css('--heat') });
const LIGHT = { bg: '#fcfcfb', grid: '#e9e7e2', text: '#0b0b0b', text2: '#52514e', text3: '#85837c', band: 'rgba(11,11,11,.07)', s: ['#2a78d6', '#eb6834'], heat: '#b3261e' };
const niceStep = max => [1, 2, 5].map(k => k * 10 ** Math.floor(Math.log10(max / 4))).concat(10 ** Math.ceil(Math.log10(max / 4))).find(s => max / s <= 5);
function axisY(svg, P, m, W, H, max, y) {
  const step = niceStep(max || 1), top = Math.ceil((max || 1) / step) * step;
  for (let v = 0; v <= top + 1e-9; v += step) {
    svg.append(el('line', { x1: m.l, x2: W - m.r, y1: y(v, top), y2: y(v, top), stroke: P.grid, 'stroke-width': 1 }));
    svg.append(el('text', { x: m.l - 7, y: y(v, top) + 4, 'text-anchor': 'end', fill: P.text3, 'font-size': 11 }, (+v.toFixed(6)).toLocaleString(undefined, { maximumFractionDigits: 3 })));
  }
  return top;
}
function svgBase(W, H, P, title) {
  const svg = el('svg', { xmlns: NS, viewBox: `0 0 ${W} ${H}`, width: W, height: H, 'font-family': 'system-ui, -apple-system, Segoe UI, sans-serif', role: 'img', 'aria-label': title });
  svg.append(el('title', {}, title));
  return svg;
}

function monthlyChart(a, W, P, interactive) {
  const Hh = 300, m = { l: 48, r: 14, t: 12, b: 26 }, S = a.series;
  const svg = svgBase(W, Hh, P, `Fire activity per sampled day by month (${a.unit})`);
  const t0 = S[0].t - 1 / 24, t1 = S[S.length - 1].t + 1 / 24;
  const x = t => m.l + (t - t0) / Math.max(t1 - t0, 1e-9) * (W - m.l - m.r);
  const max = Math.max(0, ...S.flatMap(s => s.rate.filter(v => v != null)));
  const y = (v, top) => Hh - m.b - v / top * (Hh - m.t - m.b);
  const top = axisY(svg, P, m, W, Hh, max, y);
  const yrs = t1 - t0, every = yrs > 16 ? 4 : yrs > 8 ? 2 : 1;
  if (yrs > 2) { for (let yr = Math.ceil(t0); yr <= t1; yr += every) svg.append(el('text', { x: x(yr), y: Hh - 8, 'text-anchor': 'middle', fill: P.text3, 'font-size': 11 }, yr)); }
  else S.forEach((s, i) => { if (S.length <= 12 || i % 3 === 0) svg.append(el('text', { x: x(s.t), y: Hh - 8, 'text-anchor': 'middle', fill: P.text3, 'font-size': 11 }, MON[+s.month.slice(5) - 1] + (s.month.endsWith('01') || i === 0 ? ' ' + s.month.slice(0, 4) : ''))); });
  for (const j of [0, 1]) {
    if (!a.per[j].on) continue;
    let d = '', pen = false, pts = [];
    S.forEach(s => { const v = s.rate[j]; if (v == null) { pen = false; return; } d += (pen ? 'L' : 'M') + x(s.t).toFixed(1) + ' ' + y(v, top).toFixed(1); pen = true; pts.push(s); });
    svg.append(el('path', { d, fill: 'none', stroke: P.s[j], 'stroke-width': 2, 'stroke-linejoin': 'round', 'stroke-linecap': 'round' }));
    if (S.length <= 36) pts.forEach(s => svg.append(el('circle', { cx: x(s.t), cy: y(s.rate[j], top), r: 3, fill: P.s[j], stroke: P.bg, 'stroke-width': 1.5 })));
  }
  if (!S.some(s => s.rate.some(v => v != null))) svg.append(el('text', { x: W / 2, y: Hh / 2, 'text-anchor': 'middle', fill: P.text3, 'font-size': 13 }, 'No sampled days in this range'));
  if (interactive) hoverLine(svg, a, x, y, top, m, Hh, W);
  return svg;
}
function hoverLine(svg, a, x, y, top, m, Hh, W) {
  const g = el('g'); svg.append(g);
  const hit = el('rect', { x: m.l, y: m.t, width: W - m.l - m.r, height: Hh - m.t - m.b, fill: 'transparent' }); svg.append(hit);
  const tip = $('monthly').querySelector('.tip');
  hit.addEventListener('mousemove', ev => {
    const r = svg.getBoundingClientRect(), px = (ev.clientX - r.left) * W / r.width;
    const s = a.series.reduce((b, q) => Math.abs(x(q.t) - px) < Math.abs(x(b.t) - px) ? q : b);
    g.innerHTML = ''; g.append(el('line', { x1: x(s.t), x2: x(s.t), y1: m.t, y2: Hh - m.b, stroke: css('--text-3'), 'stroke-width': 1 }));
    [0, 1].forEach(j => s.rate[j] != null && g.append(el('circle', { cx: x(s.t), cy: y(s.rate[j], top), r: 4.5, fill: css(SKEY[j]), stroke: css('--surface'), 'stroke-width': 2 })));
    tip.innerHTML = `<b>${MON[+s.month.slice(5) - 1]} ${s.month.slice(0, 4)}</b>` + [0, 1].filter(j => a.per[j].on).map(j =>
      `<i class="sw" style="background:var(${SKEY[j]})"></i> ${SENS[j]} ${s.rate[j] == null ? '<span class="m">no sampled days</span>' : `${fmt(s.rate[j])} / day <span class="m">(${fmt(s.count[j])} in ${s.days[j]} days)</span>`}`).join('<br>');
    tip.style.display = 'block';
    const sx = x(s.t) * r.width / W; tip.style.left = (sx + 14 + tip.offsetWidth > r.width ? sx - tip.offsetWidth - 14 : sx + 14) + 'px'; tip.style.top = '8px';
  });
  hit.addEventListener('mouseleave', () => { g.innerHTML = ''; tip.style.display = 'none'; });
}

function barChart(id, rows, key, W, P, interactive, opts) {
  const Hh = 250, m = { l: 44, r: 10, t: 10, b: 30 }, on = [0, 1].filter(j => A.per[j].on);
  const svg = svgBase(W, Hh, P, opts.title);
  const max = Math.max(0, ...rows.flatMap(r => r[key].filter(v => v != null)));
  const y = (v, top) => Hh - m.b - v / top * (Hh - m.t - m.b);
  const top = axisY(svg, P, m, W, Hh, opts.pct ? Math.max(max, .01) : max, y);
  if (opts.pct) svg.querySelectorAll('text').forEach(t => { if (t.getAttribute('text-anchor') === 'end') t.textContent = Math.round(parseFloat(t.textContent) * 100) + '%'; });
  const bw = (W - m.l - m.r) / rows.length, gap = 2, inner = Math.min(18, (bw - 8 - gap * (on.length - 1)) / Math.max(on.length, 1));
  const tip = interactive ? $(id).querySelector('.tip') : null;
  rows.forEach((r, i) => {
    const cx = m.l + bw * (i + .5), x0 = cx - (inner * on.length + gap * (on.length - 1)) / 2;
    svg.append(el('text', { x: cx, y: Hh - 10, 'text-anchor': 'middle', fill: P.text3, 'font-size': 10.5 }, r.name));
    on.forEach((j, k) => {
      const v = r[key][j]; if (v == null) return;
      const h = Math.max(0, Hh - m.b - y(v, top)), bx = x0 + k * (inner + gap);
      svg.append(el('path', { d: barPath(bx, Hh - m.b, inner, h, Math.min(3, inner / 2, h)), fill: P.s[j] }));
    });
    if (interactive) {
      const hit = el('rect', { x: m.l + bw * i, y: m.t, width: bw, height: Hh - m.t - m.b, fill: 'transparent' }); svg.append(hit);
      hit.addEventListener('mousemove', ev => {
        const rr = svg.getBoundingClientRect();
        tip.innerHTML = `<b>${opts.label(r)}</b>` + on.map(j => `<i class="sw" style="background:var(${SKEY[j]})"></i> ${SENS[j]} ${r[key][j] == null ? '–' : opts.value(r[key][j])}`).join('<br>');
        tip.style.display = 'block'; const sx = (ev.clientX - rr.left);
        tip.style.left = (sx + 14 + tip.offsetWidth > rr.width ? sx - tip.offsetWidth - 14 : sx + 14) + 'px'; tip.style.top = '8px';
        hit.setAttribute('fill', css('--band'));
      });
      hit.addEventListener('mouseleave', () => { tip.style.display = 'none'; hit.setAttribute('fill', 'transparent'); });
    }
  });
  return svg;
}
const barPath = (x, base, w, h, r) => h <= 0 ? '' : `M${x},${base}V${base - h + r}Q${x},${base - h} ${x + r},${base - h}H${x + w - r}Q${x + w},${base - h} ${x + w},${base - h + r}V${base}Z`;

function mapSVG(a, W, P) {   // static density map for the report (no basemap)
  const asp = Math.cos((st.s + st.n) / 2 * Math.PI / 180), dw = (st.e - st.w) * asp, dh = st.n - st.s;
  const m = { l: 44, r: 10, t: 10, b: 26 }, iw = W - m.l - m.r, ih = Math.min(520, iw * dh / dw), Hh = ih + m.t + m.b;
  const svg = svgBase(W, Hh, P, 'Hotspot density in the selected region');
  const sx = iw / dw, x = lon => m.l + (lon - st.w) * asp * sx, y = lat => m.t + (st.n - lat) * sx;
  svg.append(el('rect', { x: m.l, y: m.t, width: iw, height: ih, fill: P.bg, stroke: P.text3, 'stroke-width': 1 }));
  const clip = el('clipPath', { id: 'mapclip' }); clip.append(el('rect', { x: m.l, y: m.t, width: iw, height: ih })); svg.append(clip);
  const cg = el('g', { 'clip-path': 'url(#mapclip)' }); svg.append(cg);
  const max = Math.max(1e-9, ...a.cells.map(c => c[2]));
  for (const [lo, la, v] of a.cells) {
    if (lo + a.mc < st.w || lo > st.e || la + a.mc < st.s || la > st.n) continue;
    cg.append(el('rect', { x: x(lo), y: y(la + a.mc), width: Math.max(1, a.mc * asp * sx), height: Math.max(1, a.mc * sx), fill: P.heat, 'fill-opacity': (.12 + .8 * Math.log1p(v) / Math.log1p(max)).toFixed(2) }));
  }
  const tick = (lo, hi) => { const s = niceStep((hi - lo) * 1.2) ; const out = []; for (let v = Math.ceil(lo / s) * s; v <= hi; v += s) out.push(+v.toFixed(3)); return out; };
  tick(st.w, st.e).forEach(v => svg.append(el('text', { x: x(v), y: Hh - 8, 'text-anchor': 'middle', fill: P.text3, 'font-size': 11 }, lonS(v).replace('.00', ''))));
  tick(st.s, st.n).forEach(v => svg.append(el('text', { x: m.l - 6, y: y(v) + 4, 'text-anchor': 'end', fill: P.text3, 'font-size': 11 }, latS(v).replace('.00', ''))));
  return svg;
}

// ---------------------------------------------------------------- map
let map, cellLayer, boxLayer, drawing = false;
function initMap() {
  if (map) { fitRegion(); return; }
  map = L.map('map', { renderer: L.canvas(), zoomSnap: .25 });
  // don't hijack page scrolling: wheel-zoom only after the map is clicked, off again when the pointer leaves
  map.scrollWheelZoom.disable(); map.on('click', () => map.scrollWheelZoom.enable()); map.on('mouseout', () => map.scrollWheelZoom.disable());
  L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Base/MapServer/tile/{z}/{y}/{x}',
    { attribution: 'Basemap © Esri · Fire data: NASA FIRMS', maxZoom: 13, className: 'base' }).addTo(map);
  L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Reference/MapServer/tile/{z}/{y}/{x}', { maxZoom: 13, className: 'base' }).addTo(map);
  cellLayer = L.layerGroup().addTo(map); boxLayer = L.layerGroup().addTo(map);
  const bb = META.bbox; map.fitBounds([[bb[1], bb[0]], [bb[3], bb[2]]]);
  let start = null, rect = null;
  map.on('mousedown', ev => { if (!drawing) return; start = ev.latlng; rect = L.rectangle([start, start], { color: css('--text'), weight: 1.5, dashArray: '4 3', fill: false }).addTo(boxLayer); });
  map.on('mousemove', ev => { if (drawing && start) rect.setBounds([start, ev.latlng]); });
  map.on('mouseup', ev => {
    if (!drawing || !start) return;
    const b = L.latLngBounds(start, ev.latlng); start = null; setDrawing(false);
    if (Math.abs(b.getEast() - b.getWest()) < .01 || Math.abs(b.getNorth() - b.getSouth()) < .01) { drawBox(); return; }
    setBox([b.getWest(), b.getSouth(), b.getEast(), b.getNorth()]); $('preset').value = 'custom'; syncInputs(); update();
  });
}
function setDrawing(on) {
  drawing = on; $('draw').setAttribute('aria-pressed', on); $('map').classList.toggle('drawing', on);
  if (on) { map.dragging.disable(); $('draw').textContent = 'Drag a box on the map…'; } else { map.dragging.enable(); $('draw').textContent = 'Draw on map'; }
}
$('draw').addEventListener('click', () => setDrawing(!drawing));
addEventListener('keydown', ev => { if (ev.key === 'Escape' && drawing) setDrawing(false); });
function drawBox() {
  boxLayer.clearLayers();
  L.rectangle([[st.s, st.w], [st.n, st.e]], { color: css('--text'), weight: 1.5, dashArray: '5 4', fill: false, interactive: false }).addTo(boxLayer);
}
function fitRegion() { if (map) map.fitBounds([[st.s, st.w], [st.n, st.e]], { padding: [16, 16] }); }
function drawCells(a) {
  cellLayer.clearLayers();
  const max = Math.max(1e-9, ...a.cells.map(c => c[2])), col = css('--heat');
  for (const [lo, la, v] of a.cells)
    L.rectangle([[la, lo], [la + a.mc, lo + a.mc]], { stroke: false, fillColor: col, interactive: false, fillOpacity: .12 + .8 * Math.log1p(v) / Math.log1p(max) }).addTo(cellLayer);
  drawBox();
}

// ---------------------------------------------------------------- render
function tile(k, v, n, wide) { return `<div class="tile${wide ? ' wide' : ''}"><div class="k">${k}</div><div class="v">${v}</div><div class="n">${n}</div></div>`; }
function tilesHTML(a) {
  const on = a.per.filter(p => p.on), out = [];
  if (!on.length) return '<div class="tile wide"><div class="k">No sensor selected</div></div>';
  for (const p of on) {
    const j = SENS.indexOf(p.name);
    out.push(tile(`<i class="sw" style="background:var(${SKEY[j]})"></i> ${p.name} ${a.harm && j ? '(MODIS units)' : ''}`,
      p.days ? fmt(p.rate) + '<span style="font-size:13px;font-weight:500;color:var(--text-2)"> / day</span>' : '–',
      p.days ? `${fmt(p.count)} ${a.harm ? 'MODIS-eq.' : 'detections'} in ${p.days} sampled days · ${p.hotspots.toLocaleString()} hotspots` : 'no sampled days in range'));
  }
  const both = on.length === 2 && a.per[0].days && a.per[1].days && a.per[0].rate;
  if (both) {
    const ov = overlapRatio(a);
    out.push(tile(a.harm ? 'VIIRS ÷ MODIS (same days)' : 'VIIRS per MODIS detection', ov ? ov.toFixed(2) + '×' : '–',
      a.harm ? 'near 1× when harmonization holds for this region' : 'raw ratio on days both sensors were sampled', true));
  }
  for (const p of on) {
    const j = SENS.indexOf(p.name), pk = a.peak[j];
    out.push(tile(`${p.name} peak month`, pk ? MON[+pk.month.slice(5) - 1] + ' ' + pk.month.slice(0, 4) : '–', pk ? fmt(pk.rate[j]) + ' per sampled day' : ''));
  }
  for (const p of on) out.push(tile(`${p.name} fire intensity`, p.frpMedian == null ? '–' : fmt(p.frpMedian) + '<span style="font-size:13px;font-weight:500;color:var(--text-2)"> MW</span>',
    p.frpMedian == null ? '' : `median FRP · 90th pct ${fmt(p.frpP90)} MW · ${Math.round(p.count ? p.night / p.count * 100 : 0)}% at night`));
  return out.join('');
}
function overlapRatio(a) {   // compare sensors only over months both were sampled
  let m = 0, v = 0;
  for (const s of a.series) if (s.days[0] && s.days[1]) { m += s.count[0] / s.days[0]; v += s.count[1] / s.days[1]; }
  return m ? v / m : null;
}
function topHTML(a) {
  const on = [0, 1].filter(j => a.per[j].on);
  if (!a.tops.length) return '<tbody><tr><td class="empty">No hotspots in this selection</td></tr></tbody>';
  return `<thead><tr><th>Place (cell centre)</th>${on.map(j => `<th>${SENS[j]} / sampled day</th>`).join('')}<th>Max FRP (MW)</th><th>First seen</th><th>Last seen</th></tr></thead><tbody>` +
    a.tops.map(t => `<tr><td>${latS(t.lat)}, ${lonS(t.lon)}</td>${on.map(j => `<td>${a.per[j].days ? fmt(t.c[j] / a.per[j].days) : '–'}</td>`).join('')}<td>${fmt(t.frp)}</td><td>${dayStr(t.first)}</td><td>${dayStr(t.last)}</td></tr>`).join('') + '</tbody>';
}
const chartW = id => Math.max(320, Math.round($(id).clientWidth || 600));
function renderCharts() {
  const P = pal();
  for (const id of ['monthly', 'season', 'frp']) { $(id).innerHTML = '<div class="tip"></div>'; }
  $('monthly').prepend(monthlyChart(A, chartW('monthly'), P, true));
  $('season').prepend(barChart('season', A.season, 'rate', chartW('season'), P, true, { title: 'Seasonal cycle', label: r => r.name, value: v => fmt(v) + ' / day' }));
  $('frp').prepend(barChart('frp', A.frp, 'share', chartW('frp'), P, true, { title: 'Fire radiative power distribution', pct: true, label: r => r.name + ' MW', value: v => Math.round(v * 100) + '% of hotspots' }));
}
function update() {
  if (!H) return;
  A = analyse(); writeHash();
  const nd = [0, 1].map(j => sampledIn(j, st.from, st.to));
  $('daysnote').textContent = `${nd[0]} MODIS and ${nd[1]} VIIRS sampled days in range · ${A.idx.length.toLocaleString()} hotspots selected`;
  $('tiles').innerHTML = tilesHTML(A);
  $('leg1').innerHTML = [0, 1].filter(j => A.per[j].on).map(j => `<span><i class="sw" style="background:var(${SKEY[j]})"></i>${SENS[j]}${A.harm && j ? ' (MODIS units)' : ''}</span>`).join('');
  renderCharts(); drawCells(A);
  $('mapnote').textContent = `Cells of ${A.mc}°; shade = ${A.unit} per sampled day${[0, 1].filter(j => A.per[j].on).length > 1 ? ', averaged over the two sensors' : ''}. Dashed box = selected region.`;
  $('top').innerHTML = topHTML(A);
}
addEventListener('hashchange', () => { if (H && readHash()) { syncInputs(); fitRegion(); update(); } });
let rz; addEventListener('resize', () => { clearTimeout(rz); rz = setTimeout(() => A && renderCharts(), 150); });

// ---------------------------------------------------------------- downloads
function save(name, data, type) {
  const url = URL.createObjectURL(data instanceof Blob ? data : new Blob([data], { type }));
  const a = document.createElement('a'); a.href = url; a.download = name; document.body.append(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 2000);
}
const stem = () => `fires_${(META.name || 'region').toLowerCase().replace(/\W+/g, '-')}_${dayStr(st.from)}_${dayStr(st.to)}`;
const svgText = svg => '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(svg);
function chartFor(name, W, P) {
  return name === 'monthly' ? monthlyChart(A, W, P, false) : name === 'season' ? barChart('season', A.season, 'rate', W, P, false, { title: 'Seasonal cycle' }) :
    barChart('frp', A.frp, 'share', W, P, false, { title: 'Fire radiative power distribution', pct: true });
}
function withTitle(svg, P, title, sub) {  // standalone chart export: add a heading and legend
  const W = +svg.getAttribute('width'), H0 = +svg.getAttribute('height'), pad = 46;
  const out = svgBase(W, H0 + pad + 18, P, title);
  out.append(el('rect', { width: W, height: H0 + pad + 18, fill: P.bg }));
  out.append(el('text', { x: 12, y: 20, fill: P.text, 'font-size': 14, 'font-weight': 650 }, title));
  out.append(el('text', { x: 12, y: 37, fill: P.text2, 'font-size': 11.5 }, sub));
  const g = el('g', { transform: `translate(0 ${pad})` }); [...svg.childNodes].forEach(n => n.nodeName !== 'title' && g.append(n)); out.append(g);
  let lx = 12; [0, 1].filter(j => A.per[j].on).forEach(j => {
    out.append(el('rect', { x: lx, y: H0 + pad + 6, width: 14, height: 3, rx: 1.5, fill: P.s[j] }));
    const label = SENS[j] + (A.harm && j ? ' (MODIS units)' : ''); out.append(el('text', { x: lx + 19, y: H0 + pad + 11, fill: P.text2, 'font-size': 11.5 }, label)); lx += 30 + label.length * 6.5;
  });
  return out;
}
const TITLES = { monthly: 'Fire activity per sampled day, by month', season: 'Seasonal cycle (per sampled day)', frp: 'Fire radiative power (share of hotspots)' };
const subtitle = () => `${META.name} · ${boxS()} · ${dayStr(st.from)} to ${dayStr(st.to)} · ${A.unit}`;
document.querySelectorAll('.dl').forEach(box => box.addEventListener('click', async ev => {
  const f = ev.target.dataset.f; if (!f) return;
  const name = box.dataset.chart, svg = withTitle(chartFor(name, 900, LIGHT), LIGHT, TITLES[name], subtitle());
  if (f === 'svg') return save(`${stem()}_${name}.svg`, svgText(svg), 'image/svg+xml');
  const img = new Image(), url = URL.createObjectURL(new Blob([svgText(svg)], { type: 'image/svg+xml' }));
  await new Promise((ok, no) => { img.onload = ok; img.onerror = no; img.src = url; });
  const c = document.createElement('canvas'), s = 2; c.width = svg.getAttribute('width') * s; c.height = svg.getAttribute('height') * s;
  const g = c.getContext('2d'); g.scale(s, s); g.drawImage(img, 0, 0); URL.revokeObjectURL(url);
  c.toBlob(b => save(`${stem()}_${name}.png`, b), 'image/png');
}));

function hotspotRows() {
  const harm = A.harm;
  return A.idx.map(i => [dayStr(H.day[i]), String(H.time[i]).padStart(4, '0'), H.lat[i].toFixed(4), H.lon[i].toFixed(4), SENS[H.sensor[i]], H.sat[i],
    H.night[i] ? 'N' : 'D', ['low', 'nominal', 'high'][H.conf[i]], H.frp[i].toFixed(1), H.pix[i], H.sensor[i] ? +H.weq[i].toFixed(4) : 1]);
}
function summaryJSON() {
  return {
    dataset: { name: META.name, source: META.source, sampling: META.sampling, harmonization: META.harmonization, generated: META.generated },
    selection: { bbox_wsen: [st.w, st.s, st.e, st.n], from: dayStr(st.from), to: dayStr(st.to), sensors: SENS.filter((_, j) => st.sens[j]),
      units: A.unit, min_confidence: ['all', 'nominal', 'high'][st.conf], link: location.href },
    summary: A.per.filter(p => p.on).map(p => ({ sensor: p.name, hotspots: p.hotspots, count: +p.count.toFixed(2), sampled_days: p.days, per_sampled_day: p.rate && +p.rate.toFixed(4),
      frp_total_mw: +p.frp.toFixed(1), frp_median_mw: p.frpMedian, frp_p90_mw: p.frpP90, night_share: p.count ? +(p.night / p.count).toFixed(3) : null, high_conf_share: p.count ? +(p.high / p.count).toFixed(3) : null })),
    viirs_over_modis_same_months: overlapRatio(A),
    monthly: A.series.map(s => ({ month: s.month, modis: s.rate[0], viirs: s.rate[1], modis_count: s.count[0], viirs_count: s.count[1], modis_days: s.days[0], viirs_days: s.days[1] })),
    seasonal: A.season.map(s => ({ month: s.name, modis: s.rate[0], viirs: s.rate[1] })),
    frp_distribution: A.frp.map(b => ({ bin_mw: b.name, modis_share: b.share[0], viirs_share: b.share[1] })),
    top_places: A.tops.map(t => ({ lat: +t.lat.toFixed(2), lon: +t.lon.toFixed(2), modis: t.c[0], viirs: t.c[1], max_frp_mw: t.frp, first: dayStr(t.first), last: dayStr(t.last) })),
  };
}
function reportHTML() {
  const P = LIGHT, W = 860, sj = summaryJSON(), on = A.per.filter(p => p.on);
  const sv = (name, t) => `<h2>${t}</h2>` + withTitle(chartFor(name, W, P), P, TITLES[name], subtitle()).outerHTML;
  const rows = on.map(p => `<tr><td>${p.name}${A.harm && p.name === 'VIIRS' ? ' (MODIS units)' : ''}</td><td>${p.days}</td><td>${p.hotspots.toLocaleString()}</td><td>${fmt(p.count)}</td><td>${fmt(p.rate)}</td><td>${fmt(p.frpMedian)}</td><td>${fmt(p.frpP90)}</td><td>${p.count ? Math.round(p.night / p.count * 100) : 0}%</td></tr>`).join('');
  const hm = META.harmonization;
  return `<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Fire activity report — ${esc(META.name)}</title>
<style>body{margin:0;background:#f4f3f0;color:#0b0b0b;font:14px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif}main{max-width:900px;margin:0 auto;padding:28px 18px 48px;background:#fcfcfb}
h1{font-size:22px;margin:0 0 4px}h2{font-size:15px;margin:26px 0 8px}p,li{color:#52514e}table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums;font-size:13px}
th,td{padding:5px 8px;text-align:right;border-bottom:1px solid #e9e7e2}th:first-child,td:first-child{text-align:left}thead th{color:#52514e;border-bottom:1px solid #e2e0da}svg{max-width:100%;height:auto;display:block}.k{color:#85837c;font-size:12.5px}
@media print{body{background:#fff}main{padding:0}h2{break-after:avoid}svg{break-inside:avoid}}</style></head><body><main>
<h1>Fire activity — ${esc(META.name)}</h1>
<p>${boxS()} · ${dayStr(st.from)} to ${dayStr(st.to)} · sensors: ${on.map(p => p.name).join(' + ') || 'none'} · units: ${A.unit} · confidence: ${['all', 'nominal and high', 'high only'][st.conf]}</p>
<p class="k">Generated ${new Date().toISOString().slice(0, 16).replace('T', ' ')} UTC from ${esc(META.source || 'a user-supplied file')}. Recreate this view: <a href="${esc(location.href)}">${esc(location.href)}</a></p>
<h2>Summary</h2><table><thead><tr><th>Sensor</th><th>Sampled days</th><th>Hotspots</th><th>${A.harm ? 'MODIS-eq. count' : 'Detections'}</th><th>Per sampled day</th><th>Median FRP (MW)</th><th>90th pct FRP</th><th>Night</th></tr></thead><tbody>${rows}</tbody></table>
${sj.viirs_over_modis_same_months ? `<p>VIIRS ÷ MODIS over months both were sampled: <b>${sj.viirs_over_modis_same_months.toFixed(2)}×</b>${A.harm ? ' (near 1× means the harmonization holds for this region).' : '.'}</p>` : ''}
${sv('monthly', 'Activity by month')}${sv('season', 'Seasonal cycle')}${sv('frp', 'Fire intensity')}
<h2>Hotspot density</h2>${mapSVG(A, W, P).outerHTML}<p class="k">Cells of ${A.mc}°; darker = more ${A.unit} per sampled day. No basemap.</p>
<h2>Most active places</h2><table>${topHTML(A).replace(/var\(--[a-z]+\)/g, '')}</table>
<h2>Method and caveats</h2><ul>
<li>Source: NASA FIRMS active-fire detections${META.sampling ? `, sampled ${esc(META.sampling)}. Counts cover the sampled days only; rates are per sampled day` : ''}.</li>
${hm ? `<li>VIIRS (375 m) is put in MODIS (1 km) units: pixels in the same ${hm.merge_km} km cell within one overpass are merged into one hotspot, then calibrated to MODIS (${esc(hm.model)}; typical monthly error ×/÷ ${hm.typical_error}). ${esc(hm.note)}</li>` : ''}
<li>MODIS counts combine Terra and Aqua (Aqua from mid-2002); VIIRS is Suomi NPP only (from January 2012).</li>
<li>Detections are thermal anomalies, mostly fires; they also include some industrial heat sources and can miss fires under cloud or smoke.</li></ul>
</main></body></html>`;
}
document.querySelector('.downloads').addEventListener('click', ev => {
  const b = ev.target.closest('button'); if (!b || !A) return;
  const k = b.dataset.dl, s = stem();
  if (k === 'csv') save(s + '_hotspots.csv', 'date,time,lat,lon,sensor,satellite,daynight,confidence,frp,pixels,modis_equiv\n' + hotspotRows().map(r => r.join(',')).join('\n') + '\n', 'text/csv');
  else if (k === 'geojson') save(s + '_hotspots.geojson', JSON.stringify({ type: 'FeatureCollection', features: hotspotRows().map(r => ({ type: 'Feature', geometry: { type: 'Point', coordinates: [+r[3], +r[2]] },
    properties: { date: r[0], time: r[1], sensor: r[4], satellite: r[5], daynight: r[6], confidence: r[7], frp: +r[8], pixels: r[9], modis_equiv: r[10] } })) }), 'application/geo+json');
  else if (k === 'monthly') save(s + '_monthly.csv', 'month,sensor,count,units,sampled_days,per_sampled_day\n' + A.series.flatMap(m => [0, 1].filter(j => A.per[j].on).map(j =>
    [m.month, SENS[j], +m.count[j].toFixed(3), A.harm ? 'modis_equiv' : 'detections', m.days[j], m.rate[j] == null ? '' : +m.rate[j].toFixed(4)].join(','))).join('\n') + '\n', 'text/csv');
  else if (k === 'json') save(s + '_analytics.json', JSON.stringify(summaryJSON(), null, 1), 'application/json');
  else if (k === 'report') save(s + '_report.html', reportHTML(), 'text/html');
});

boot();
