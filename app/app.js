// Fire Activity Explorer: harmonized MODIS + VIIRS hotspots -> region/date analytics, recorded fires, downloads.
// Inputs: data/hotspots.csv.gz + data/meta.json (../hotspots.py), data/fires.json (../known_fires.py, optional),
// or any CSV with the same columns via "Load hotspot file…".
'use strict';

const $ = id => document.getElementById(id);
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const SENS = ['MODIS', 'VIIRS'], SKEY = ['--modis', '--viirs'];
const DAY = 864e5, MON = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
const toDay = s => Date.UTC(+s.slice(0, 4), +s.slice(5, 7) - 1, +s.slice(8, 10)) / DAY;
const dayStr = d => new Date(d * DAY).toISOString().slice(0, 10);
const monOf = d => { const s = dayStr(d); return (+s.slice(0, 4) - 2000) * 12 + +s.slice(5, 7) - 1; };
const monKey = m => `${2000 + Math.floor(m / 12)}-${String(m % 12 + 1).padStart(2, '0')}`;
const monLabel = m => `${MON[m % 12]} ${2000 + Math.floor(m / 12)}`;
const fmt = n => n == null || !isFinite(n) ? '–' : Math.abs(n) >= 100 ? Math.round(n).toLocaleString() : Math.abs(n) >= 10 ? n.toFixed(1) : n.toFixed(2);
const acresS = a => a >= 1e4 ? Math.round(a / 1e3) + 'k' : a >= 1e3 ? (a / 1e3).toFixed(1) + 'k' : String(a);
const lonS = v => `${Math.abs(v).toFixed(2)}°${v < 0 ? 'W' : 'E'}`, latS = v => `${Math.abs(v).toFixed(2)}°${v < 0 ? 'S' : 'N'}`;
const boxS = () => `${latS(st.s)}–${latS(st.n)}, ${lonS(st.w)}–${lonS(st.e)}`;
const causeS = c => /^missing|undetermined/i.test(c || '') ? 'Unknown' : String(c).replace(/\/.*$/, '');
const esc = s => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

const PRESETS = [   // approximate boxes (W, S, E, N)
  ['North Coast & Klamath', [-124.5, 39.5, -122.0, 42.0]], ['Sacramento Valley', [-122.6, 38.3, -121.2, 40.6]],
  ['Bay Area', [-123.1, 36.9, -121.5, 38.6]], ['Sierra Nevada', [-121.3, 35.5, -118.0, 40.2]],
  ['Central Coast', [-122.0, 34.4, -119.5, 37.0]], ['San Joaquin Valley', [-121.2, 34.9, -118.6, 37.9]],
  ['Southern California', [-119.5, 32.5, -114.5, 35.0]],
];
const N_PRESETS = PRESETS.length;

let H = null, META = null, FIRES = [], FIRESRC = null;
let SAMPLED = [[], []], SMON = [new Map(), new Map()];   // sampled days per sensor; per month
const st = { w: 0, s: 0, e: 0, n: 0, from: 0, to: 0, sens: [true, true], units: 'harm', conf: 0, view: 'one' };
let A = null;

// ---------------------------------------------------------------- loading
const gunzipText = src => new Response((src.body || src.stream()).pipeThrough(new DecompressionStream('gzip'))).text();
async function readFile(file) {
  const h = new Uint8Array(await file.slice(0, 2).arrayBuffer());
  return h[0] === 0x1f && h[1] === 0x8b ? gunzipText(file) : file.text();
}
function parseCSV(text) {
  const lines = text.split(/\r?\n/).filter(l => l && !l.startsWith('#'));
  const cols = lines[0].split(',').map(c => c.trim().toLowerCase()), ix = n => cols.indexOf(n);
  const need = ['date', 'lat', 'lon', 'sensor'].filter(c => ix(c) < 0);
  if (need.length) throw new Error(`missing column(s): ${need.join(', ')}. Expected date,time,lat,lon,sensor,satellite,daynight,confidence,frp,pixels,modis_equiv`);
  const n = lines.length - 1;
  const h = { n, day: new Int32Array(n), mon: new Int16Array(n), time: new Uint16Array(n), lat: new Float32Array(n), lon: new Float32Array(n),
    sensor: new Uint8Array(n), sat: [], night: new Uint8Array(n), conf: new Uint8Array(n), frp: new Float32Array(n),
    pix: new Uint16Array(n), weq: new Float32Array(n), hasWeq: ix('modis_equiv') >= 0 };
  const c = { date: ix('date'), time: ix('time'), lat: ix('lat'), lon: ix('lon'), sensor: ix('sensor'), sat: ix('satellite'),
    dn: ix('daynight'), conf: ix('confidence'), frp: ix('frp'), pix: ix('pixels'), weq: ix('modis_equiv') };
  const CONF = { low: 0, l: 0, nominal: 1, n: 1, high: 2, h: 2 }, cache = new Map();
  let k = 0;
  for (let i = 1; i <= n; i++) {
    const f = lines[i].split(','), ds = f[c.date];
    let dm = cache.get(ds); if (!dm) { const d = toDay(ds); dm = [d, isFinite(d) ? monOf(d) : 0]; cache.set(ds, dm); }
    h.day[k] = dm[0]; h.mon[k] = dm[1]; h.time[k] = c.time >= 0 ? +f[c.time] || 0 : 0;
    h.lat[k] = +f[c.lat]; h.lon[k] = +f[c.lon]; h.sensor[k] = /^v/i.test(f[c.sensor]) ? 1 : 0;
    h.sat[k] = c.sat >= 0 ? f[c.sat] : ''; h.night[k] = c.dn >= 0 && f[c.dn] === 'N' ? 1 : 0;
    const cv = c.conf >= 0 ? f[c.conf] : 'nominal'; h.conf[k] = cv in CONF ? CONF[cv] : (+cv < 30 ? 0 : +cv < 80 ? 1 : 2);
    h.frp[k] = c.frp >= 0 ? +f[c.frp] || 0 : 0; h.pix[k] = c.pix >= 0 ? +f[c.pix] || 1 : 1; h.weq[k] = c.weq >= 0 ? +f[c.weq] : 1;
    if (isFinite(dm[0]) && isFinite(h.lat[k]) && isFinite(h.lon[k])) k++;
  }
  for (const key of ['day', 'mon', 'time', 'lat', 'lon', 'sensor', 'night', 'conf', 'frp', 'pix', 'weq']) h[key] = h[key].subarray(0, k);
  h.n = k; return h;
}

function useData(h, meta, label) {
  H = h;
  let w = 180, s = 90, e = -180, n = -90; const days = [new Set(), new Set()];
  for (let i = 0; i < h.n; i++) { w = Math.min(w, h.lon[i]); e = Math.max(e, h.lon[i]); s = Math.min(s, h.lat[i]); n = Math.max(n, h.lat[i]); days[h.sensor[i]].add(h.day[i]); }
  META = meta || { name: label, bbox: [w, s, e, n] };
  SAMPLED = SENS.map((k, j) => Int32Array.from(META.sampled_days && META.sampled_days[k] ? META.sampled_days[k].map(toDay) : [...days[j]]).sort());
  SMON = SAMPLED.map(a => { const m = new Map(); for (const d of a) { const k = monOf(d); m.set(k, (m.get(k) || 0) + 1); } return m; });
  $('region').textContent = META.name;
  PRESETS.splice(0, PRESETS.length - N_PRESETS, ['Whole dataset', META.bbox]);
  $('preset').innerHTML = PRESETS.map((p, i) => `<option value="${i}">${esc(p[0])}</option>`).join('') + '<option value="custom" hidden>Custom</option>';
  const ok = SAMPLED.filter(a => a.length), dmin = Math.min(...ok.map(a => a[0])), dmax = Math.max(...ok.map(a => a[a.length - 1]));
  for (const k of ['from', 'to']) { $(k).min = dayStr(dmin); $(k).max = dayStr(dmax); }
  const y0 = +dayStr(dmin).slice(0, 4), y1 = +dayStr(dmax).slice(0, 4), q = [['all', 'All dates'], ['last12', 'Last 12 months']];
  for (let y = y1; y >= y0; y--) q.push(['y' + y, String(y)]);
  $('quick').innerHTML = '<option value="">Quick range…</option>' + q.map(([v, t]) => `<option value="${v}">${t}</option>`).join('');
  Object.assign($('quick').dataset, { min: dmin, max: dmax });
  const hm = META.harmonization;
  $('how').innerHTML = (META.sampling ? `<p><b>Sampled data.</b> ${esc(META.sampling)}, so counts cover those days only and every rate is <i>per sampled day</i>. A date range with no sampled days shows nothing.</p>` :
      '<p>Sampled days are taken to be the days present in the file.</p>') +
    (hm ? `<p><b>Harmonized.</b> MODIS (1 km) is the reference. VIIRS (375 m) pixels in the same ${hm.merge_km} km cell within one overpass are merged into one hotspot, then calibrated to MODIS (${esc(hm.model)}; typical monthly error ×/÷ ${hm.typical_error}). The default line averages MODIS and harmonized VIIRS where both flew. ${esc(hm.note)}</p>` :
      h.hasWeq ? '' : '<p class="err">This file has no modis_equiv column, so VIIRS can only be shown as raw detections.</p>') +
    (FIRESRC ? `<p><b>Recorded fires</b> are from ${esc(FIRESRC.source)}. Their location is ${esc(FIRESRC.location)}. Wikipedia links: ${esc(FIRESRC.wikipedia)}; other fires link to a Wikipedia search.</p>` : '') +
    '<p><a href="../fires.html">Full method and sensor comparison →</a></p>';
  $('foot').innerHTML = `${esc(META.source || label)} · ${h.n.toLocaleString()} hotspots` + (META.generated ? ` · built ${esc(META.generated)}` : '');
  if (!h.hasWeq) st.units = 'raw';
  initMap();
  if (!readHash()) { setBox(META.bbox); st.from = dmin; st.to = dmax; $('preset').value = '0'; }
  syncInputs(); fitRegion(); update();
  $('loader').hidden = true;
}

async function boot() {
  try {
    const [r, m, f] = await Promise.all([fetch('data/hotspots.csv.gz'), fetch('data/meta.json'), fetch('data/fires.json').catch(() => null)]);
    if (!r.ok) throw new Error(r.status);
    if (f && f.ok) setFires(await f.json());
    const meta = m.ok ? await m.json() : null, buf = new Uint8Array(await r.arrayBuffer());
    const text = buf[0] === 0x1f && buf[1] === 0x8b ? await gunzipText(new Blob([buf])) : new TextDecoder().decode(buf);   // some hosts pre-decode .gz
    useData(parseCSV(text), meta, 'data/hotspots.csv.gz');
  } catch (err) {
    $('ltitle').textContent = 'Load a hotspot file';
    $('lmsg').innerHTML = location.protocol === 'file:' ?
      'Pages opened from disk can\'t read the files next to them. Pick <code>app/data/hotspots.csv.gz</code> below, or run <code>python3 -m http.server</code> in the repo and open <code>http://localhost:8000/app/</code>.' :
      `Couldn't load data/hotspots.csv.gz (${esc(err.message)}). Pick a harmonized hotspot CSV (or .csv.gz) instead.`;
    $('lpick').hidden = false;
  }
}
function setFires(j) {
  FIRESRC = j;
  FIRES = (j.fires || []).map(f => { const a = toDay(f.start), b = Math.max(a, toDay(f.end || f.start)); return { ...f, a, b, ma: monOf(a), mb: monOf(b) }; });
}
async function onFile(ev) {
  const f = ev.target.files[0]; if (!f) return;
  try {
    $('loader').hidden = false; $('ltitle').textContent = 'Reading ' + f.name + '…'; $('lpick').hidden = true; $('lmsg').textContent = '';
    const h = parseCSV(await readFile(f)); if (!h.n) throw new Error('no rows');
    useData(h, null, f.name);
  } catch (err) { $('ltitle').textContent = 'Could not read ' + f.name; $('lmsg').innerHTML = `<span class="err">${esc(err.message)}</span>`; $('lpick').hidden = false; }
}
$('file').addEventListener('change', onFile); $('file2').addEventListener('change', onFile);

// ---------------------------------------------------------------- selection state
function setBox([w, s, e, n]) { Object.assign(st, { w: Math.min(w, e), e: Math.max(w, e), s: Math.min(s, n), n: Math.max(s, n) }); }
function syncInputs() {
  for (const k of ['w', 's', 'e', 'n']) $(k).value = (+st[k]).toFixed(2);
  $('from').value = dayStr(st.from); $('to').value = dayStr(st.to);
  $('sM').checked = st.sens[0]; $('sV').checked = st.sens[1]; $('conf').value = st.conf;
  [...$('units').children].forEach(b => b.setAttribute('aria-pressed', b.dataset.u === st.units));
  if (st.units !== 'harm') st.view = 'sensors';
  $('view').children[0].disabled = st.units !== 'harm';
  [...$('view').children].forEach(b => b.setAttribute('aria-pressed', b.dataset.v === st.view));
}
function writeHash() {
  const p = new URLSearchParams({ bbox: [st.w, st.s, st.e, st.n].map(v => (+v).toFixed(3)).join(','), from: dayStr(st.from), to: dayStr(st.to) });
  if (!st.sens[0] || !st.sens[1]) p.set('sensors', SENS.filter((_, j) => st.sens[j]).join(','));
  if (st.units !== 'harm') p.set('units', st.units);
  if (st.conf) p.set('conf', st.conf);
  if (st.view !== 'one') p.set('view', st.view);
  history.replaceState(null, '', '#' + p);
}
function readHash() {
  if (!location.hash) return false;
  const p = new URLSearchParams(location.hash.slice(1)), bb = (p.get('bbox') || '').split(',').map(Number);
  if (bb.length !== 4 || bb.some(v => !isFinite(v)) || !p.get('from') || !p.get('to')) return false;
  setBox(bb); st.from = toDay(p.get('from')); st.to = toDay(p.get('to'));
  const s = (p.get('sensors') || 'MODIS,VIIRS').split(','); st.sens = SENS.map(k => s.includes(k));
  st.units = p.get('units') === 'raw' || !H.hasWeq ? 'raw' : 'harm'; st.conf = +p.get('conf') || 0; st.view = p.get('view') === 'sensors' ? 'sensors' : 'one';
  const pi = PRESETS.findIndex(([, b]) => b.every((v, i) => Math.abs(v - bb[i]) < 1e-3));
  $('preset').value = pi >= 0 ? String(pi) : 'custom';
  return true;
}

let timer; const schedule = () => { clearTimeout(timer); timer = setTimeout(update, 120); };
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
  if (v === 'all') { st.from = dmin; st.to = dmax; }
  else if (v === 'last12') { st.to = dmax; st.from = dmax - 365; }
  else if (v[0] === 'y') { st.from = Math.max(dmin, toDay(v.slice(1) + '-01-01')); st.to = Math.min(dmax, toDay(v.slice(1) + '-12-31')); }
  syncInputs(); schedule();
});
$('sM').addEventListener('change', () => { st.sens[0] = $('sM').checked; schedule(); });
$('sV').addEventListener('change', () => { st.sens[1] = $('sV').checked; schedule(); });
$('conf').addEventListener('change', () => { st.conf = +$('conf').value; schedule(); });
$('units').addEventListener('click', ev => { const u = ev.target.dataset.u; if (!u || (u === 'harm' && !H.hasWeq)) return; st.units = u; if (u === 'harm') st.view = 'one'; syncInputs(); schedule(); });
$('view').addEventListener('click', ev => { const v = ev.target.dataset.v; if (!v || ev.target.disabled) return; st.view = v; syncInputs(); writeHash(); $('tiles').innerHTML = tilesHTML(A); render(); });

// ---------------------------------------------------------------- analytics
const lowerBound = (a, x) => { let lo = 0, hi = a.length; while (lo < hi) { const m = (lo + hi) >> 1; if (a[m] < x) lo = m + 1; else hi = m; } return lo; };
const sampledIn = (j, from, to) => lowerBound(SAMPLED[j], to + 1) - lowerBound(SAMPLED[j], from);

function analyse() {
  const harm = st.units === 'harm', on = st.sens.slice();
  const wt = i => H.sensor[i] === 0 ? 1 : harm ? H.weq[i] : H.pix[i];
  const idx = [], allM = new Map(), inM = new Map();
  const per = SENS.map((name, j) => ({ name, on: on[j], hotspots: 0, count: 0, frp: [], night: 0, days: sampledIn(j, st.from, st.to) }));
  for (let i = 0; i < H.n; i++) {
    const j = H.sensor[i]; if (!on[j] || H.conf[i] < st.conf) continue;
    const la = H.lat[i], lo = H.lon[i]; if (la < st.s || la > st.n || lo < st.w || lo > st.e) continue;
    const w = wt(i), m = H.mon[i];
    let c = allM.get(m); if (!c) allM.set(m, c = [0, 0]); c[j] += w;
    const d = H.day[i]; if (d < st.from || d > st.to) continue;
    idx.push(i); const p = per[j]; p.hotspots++; p.count += w; p.frp.push(H.frp[i]); if (H.night[i]) p.night += w;
    let q = inM.get(m); if (!q) inM.set(m, q = [0, 0]); q[j] += w;
  }
  for (const p of per) { p.rate = p.days ? p.count / p.days : null; p.frp.sort((a, b) => a - b); p.frpMedian = p.frp.length ? p.frp[p.frp.length >> 1] : null; }
  // harmonized activity = mean of the sensors (in MODIS units) that were sampled that month
  const comb = r => { if (!harm) return null; const v = r.filter((x, j) => on[j] && x != null); return v.length ? v.reduce((a, b) => a + b) / v.length : null; };
  const m0 = monOf(st.from), m1 = monOf(st.to), clip = [0, 1].map(j => { const mm = new Map(), a = SAMPLED[j];
    for (let k = lowerBound(a, st.from); k < a.length && a[k] <= st.to; k++) { const x = monOf(a[k]); mm.set(x, (mm.get(x) || 0) + 1); } return mm; });
  const series = [];
  for (let m = m0; m <= m1; m++) {
    const c = inM.get(m) || [0, 0], days = [clip[0].get(m) || 0, clip[1].get(m) || 0], rate = [0, 1].map(j => on[j] && days[j] ? c[j] / days[j] : null);
    series.push({ m, t: 2000 + (m + .5) / 12, count: c, days, rate, comb: comb(rate), fires: [] });
  }
  const wavg = rows => { let s = 0, d = 0; for (const r of rows) { const dd = Math.max(...r.days); if (r.comb != null && dd) { s += r.comb * dd; d += dd; } } return d ? s / d : null; };
  const allRows = [], ks = SMON.flatMap(x => [...x.keys()]), lo = Math.min(...ks), hi = Math.max(...ks);
  for (let m = lo; m <= hi; m++) { const c = allM.get(m) || [0, 0], days = [SMON[0].get(m) || 0, SMON[1].get(m) || 0];
    allRows.push({ days, comb: comb([0, 1].map(j => on[j] && days[j] ? c[j] / days[j] : null)) }); }
  const fires = FIRES.filter(f => f.b >= st.from && f.a <= st.to && f.lat >= st.s && f.lat <= st.n && f.lon >= st.w && f.lon <= st.e).sort((a, b) => b.acres - a.acres);
  for (const f of fires) for (const s of series) if (s.m >= f.ma && s.m <= f.mb) s.fires.push(f);
  const season = MON.map((name, k) => { const rows = series.filter(s => s.m % 12 === k);
    return { name, comb: wavg(rows), rate: [0, 1].map(j => { let c = 0, d = 0; for (const s of rows) { c += s.count[j]; d += s.days[j]; } return on[j] && d ? c / d : null; }) }; });
  let am = 0, av = 0; for (const s of series) if (s.days[0] && s.days[1]) { am += s.count[0] / s.days[0]; av += s.count[1] / s.days[1]; }
  const span = Math.max(st.e - st.w, st.n - st.s), mc = span > 6 ? 0.1 : span > 2 ? 0.05 : span > 0.6 ? 0.02 : 0.01, grid = new Map();
  for (const i of idx) { const k = Math.floor(H.lon[i] / mc) + ',' + Math.floor(H.lat[i] / mc); let g = grid.get(k); if (!g) grid.set(k, g = [0, 0]); g[H.sensor[i]] += wt(i); }
  const rateOf = c => { const r = [0, 1].filter(j => on[j] && per[j].days).map(j => c[j] / per[j].days); return r.length ? r.reduce((a, b) => a + b) / r.length : 0; };
  const cells = [...grid].map(([k, c]) => { const [x, y] = k.split(',').map(Number); return [x * mc, y * mc, rateOf(c)]; }).filter(c => c[2] > 0);
  return { idx, per, series, overall: wavg(series), longTerm: wavg(allRows), fires, season, cells, mc, harm, agree: am ? av / am : null,
    unit: harm ? 'MODIS-equivalent hotspots' : 'raw detections' };
}
const lines = a => st.view === 'one' && a.harm ? [{ key: 'comb', color: 'heat', name: 'Fire activity (MODIS units)' }] :
  [0, 1].filter(j => a.per[j].on).map(j => ({ key: j, color: j, name: SENS[j] + (a.harm && j ? ' (MODIS units)' : '') }));
const pick = (s, key) => key === 'comb' ? s.comb : s.rate[key];
const shown = (a, s) => { const v = lines(a).map(l => pick(s, l.key)).filter(x => x != null); return v.length ? Math.max(...v) : null; };
const peakOf = a => a.series.reduce((b, s) => shown(a, s) != null && (!b || shown(a, s) > shown(a, b)) ? s : b, null);

// ---------------------------------------------------------------- charts (explicit colours, so exports are self-contained)
const NS = 'http://www.w3.org/2000/svg';
const el = (n, a, txt) => { const e = document.createElementNS(NS, n); for (const k in a) e.setAttribute(k, a[k]); if (txt != null) e.textContent = txt; return e; };
const pal = () => ({ bg: css('--surface'), grid: css('--grid'), text: css('--text'), text2: css('--text-2'), text3: css('--text-3'), s: [css('--modis'), css('--viirs')], heat: css('--heat') });
const LIGHT = { bg: '#fcfcfb', grid: '#e9e7e2', text: '#0b0b0b', text2: '#52514e', text3: '#85837c', s: ['#2a78d6', '#eb6834'], heat: '#b3261e' };
const colorOf = (P, c) => c === 'heat' ? P.heat : P.s[c], cssOf = c => c === 'heat' ? '--heat' : SKEY[c];
const niceStep = max => [1, 2, 5].map(k => k * 10 ** Math.floor(Math.log10(max / 4))).concat(10 ** Math.ceil(Math.log10(max / 4))).find(s => max / s <= 5);
function axisY(svg, P, m, W, max, y) {
  const step = niceStep(max || 1), top = Math.ceil((max || 1) / step) * step;
  for (let v = 0; v <= top + 1e-9; v += step) {
    svg.append(el('line', { x1: m.l, x2: W - m.r, y1: y(v, top), y2: y(v, top), stroke: P.grid, 'stroke-width': 1 }));
    svg.append(el('text', { x: m.l - 7, y: y(v, top) + 4, 'text-anchor': 'end', fill: P.text3, 'font-size': 11 }, (+v.toFixed(6)).toLocaleString(undefined, { maximumFractionDigits: 3 })));
  }
  return top;
}
const svgBase = (W, Hh, title) => { const s = el('svg', { xmlns: NS, viewBox: `0 0 ${W} ${Hh}`, width: W, height: Hh, 'font-family': 'system-ui, -apple-system, Segoe UI, sans-serif', role: 'img', 'aria-label': title }); s.append(el('title', {}, title)); return s; };

function monthlyChart(a, W, P, interactive) {
  const Hh = 300, m = { l: 48, r: 14, t: 28, b: 26 }, S = a.series, L = lines(a);
  const svg = svgBase(W, Hh, `Fire activity per sampled day by month (${a.unit})`);
  const t0 = S[0].t - 1 / 24, t1 = S[S.length - 1].t + 1 / 24;
  const x = t => m.l + (t - t0) / Math.max(t1 - t0, 1e-9) * (W - m.l - m.r), y = (v, top) => Hh - m.b - v / top * (Hh - m.t - m.b);
  const top = axisY(svg, P, m, W, Math.max(0, ...S.map(s => shown(a, s) ?? 0)), y);
  const yrs = t1 - t0, every = yrs > 16 ? 4 : yrs > 8 ? 2 : 1;
  if (yrs > 2) for (let yr = Math.ceil(t0); yr <= t1; yr += every) svg.append(el('text', { x: x(yr), y: Hh - 8, 'text-anchor': 'middle', fill: P.text3, 'font-size': 11 }, yr));
  else S.forEach((s, i) => { if (S.length <= 12 || i % 3 === 0) svg.append(el('text', { x: x(s.t), y: Hh - 8, 'text-anchor': 'middle', fill: P.text3, 'font-size': 11 }, MON[s.m % 12] + (s.m % 12 === 0 || i === 0 ? ' ' + (2000 + Math.floor(s.m / 12)) : ''))); });
  for (const l of L) {
    const c = colorOf(P, l.color); let d = '', area = '', run = [];
    const flush = () => { if (run.length > 1 && L.length === 1) area += 'M' + run.map(([px, py]) => px.toFixed(1) + ' ' + py.toFixed(1)).join('L') + `L${run[run.length - 1][0].toFixed(1)} ${Hh - m.b}L${run[0][0].toFixed(1)} ${Hh - m.b}Z`; run = []; };
    S.forEach(s => { const v = pick(s, l.key); if (v == null) { flush(); d += ''; return; } const px = x(s.t), py = y(v, top); d += (run.length ? 'L' : 'M') + px.toFixed(1) + ' ' + py.toFixed(1); run.push([px, py]); });
    flush();
    if (area) svg.append(el('path', { d: area, fill: c, 'fill-opacity': .1 }));
    svg.append(el('path', { d, fill: 'none', stroke: c, 'stroke-width': 2, 'stroke-linejoin': 'round', 'stroke-linecap': 'round' }));
    if (S.length <= 36) S.forEach(s => { const v = pick(s, l.key); if (v != null) svg.append(el('circle', { cx: x(s.t), cy: y(v, top), r: 3, fill: c, stroke: P.bg, 'stroke-width': 1.5 })); });
  }
  // label the three largest recorded fires at their busiest month
  const used = new Set(), notes = [];
  for (const f of a.fires) {
    if (notes.length >= 3) break;
    const rows = S.filter(s => s.m >= f.ma && s.m <= f.mb && shown(a, s) != null); if (!rows.length) continue;
    const s = rows.reduce((b, r) => shown(a, r) > shown(a, b) ? r : b); if (used.has(s.m)) continue;
    used.add(s.m); notes.push({ s, f, v: shown(a, s) });
  }
  notes.sort((p, q) => p.s.t - q.s.t).forEach((nt, i, arr) => {
    const px = x(nt.s.t), py = y(nt.v, top), lift = i > 0 && px - x(arr[i - 1].s.t) < 170 && i % 2 ? 14 : 0, ly = Math.max(12, py - 10 - lift);
    const anchor = px > W - 150 ? 'end' : px < m.l + 90 ? 'start' : 'middle';
    svg.append(el('circle', { cx: px, cy: py, r: 4.5, fill: 'none', stroke: P.text, 'stroke-width': 1.5 }));
    svg.append(el('text', { x: px, y: ly, 'text-anchor': anchor, fill: P.text, 'font-size': 11.5, 'font-weight': 600, stroke: P.bg, 'stroke-width': 3, 'paint-order': 'stroke' }, `${nt.f.name} Fire · ${acresS(nt.f.acres)} ac`));
  });
  if (!S.some(s => shown(a, s) != null)) svg.append(el('text', { x: W / 2, y: Hh / 2, 'text-anchor': 'middle', fill: P.text3, 'font-size': 13 }, 'No sampled days in this range'));
  if (interactive) hoverLine(svg, a, L, x, y, top, m, Hh, W);
  return svg;
}
function hoverLine(svg, a, L, x, y, top, m, Hh, W) {
  const g = el('g'); svg.append(g);
  const hit = el('rect', { x: m.l, y: m.t, width: W - m.l - m.r, height: Hh - m.t - m.b, fill: 'transparent' }); svg.append(hit);
  const tip = $('monthly').querySelector('.tip');
  hit.addEventListener('mousemove', ev => {
    const r = svg.getBoundingClientRect(), px = (ev.clientX - r.left) * W / r.width;
    const s = a.series.reduce((b, q) => Math.abs(x(q.t) - px) < Math.abs(x(b.t) - px) ? q : b);
    g.innerHTML = ''; g.append(el('line', { x1: x(s.t), x2: x(s.t), y1: m.t, y2: Hh - m.b, stroke: css('--text-3'), 'stroke-width': 1 }));
    L.forEach(l => { const v = pick(s, l.key); if (v != null) g.append(el('circle', { cx: x(s.t), cy: y(v, top), r: 4.5, fill: css(cssOf(l.color)), stroke: css('--surface'), 'stroke-width': 2 })); });
    const rows = L.map(l => { const v = pick(s, l.key);
      return `<i class="sw" style="background:var(${cssOf(l.color)})"></i>${l.key === 'comb' ? 'Fire activity' : SENS[l.key]} ${v == null ? '<span class="m">no sampled days</span>' : fmt(v) + ' / day'}`; });
    const fs = s.fires.slice(0, 2).map(f => `${esc(f.name)} <span class="m">${acresS(f.acres)} ac</span>`);
    tip.innerHTML = `<b>${monLabel(s.m)}</b>${rows.join('<br>')}${fs.length ? '<br><span class="m">Burning:</span> ' + fs.join(', ') : ''}`;
    tip.style.display = 'block'; const sx = x(s.t) * r.width / W;
    tip.style.left = (sx + 14 + tip.offsetWidth > r.width ? sx - tip.offsetWidth - 14 : sx + 14) + 'px'; tip.style.top = '30px';
  });
  hit.addEventListener('mouseleave', () => { g.innerHTML = ''; tip.style.display = 'none'; });
}
const barPath = (x, base, w, h, r) => h <= 0 ? '' : `M${x},${base}V${base - h + r}Q${x},${base - h} ${x + r},${base - h}H${x + w - r}Q${x + w},${base - h} ${x + w},${base - h + r}V${base}Z`;
function seasonChart(a, W, P, interactive) {
  const Hh = 240, m = { l: 44, r: 8, t: 10, b: 26 }, L = lines(a), rows = a.season;
  const svg = svgBase(W, Hh, 'Fire season: activity per sampled day by calendar month');
  const y = (v, top) => Hh - m.b - v / top * (Hh - m.t - m.b);
  const top = axisY(svg, P, m, W, Math.max(0, ...rows.flatMap(r => L.map(l => pick(r, l.key) ?? 0))), y);
  const bw = (W - m.l - m.r) / 12, gap = 2, inner = Math.min(L.length === 1 ? 22 : 14, (bw - 6 - gap * (L.length - 1)) / L.length);
  const tip = interactive ? $('season').querySelector('.tip') : null;
  rows.forEach((r, i) => {
    const cx = m.l + bw * (i + .5), x0 = cx - (inner * L.length + gap * (L.length - 1)) / 2;
    svg.append(el('text', { x: cx, y: Hh - 8, 'text-anchor': 'middle', fill: P.text3, 'font-size': 10.5 }, r.name));
    L.forEach((l, k) => { const v = pick(r, l.key); if (v == null) return; const h = Math.max(0, Hh - m.b - y(v, top));
      svg.append(el('path', { d: barPath(x0 + k * (inner + gap), Hh - m.b, inner, h, Math.min(3, inner / 2, h)), fill: colorOf(P, l.color) })); });
    if (!interactive) return;
    const hit = el('rect', { x: m.l + bw * i, y: m.t, width: bw, height: Hh - m.t - m.b, fill: 'transparent' }); svg.append(hit);
    hit.addEventListener('mousemove', ev => {
      const rr = svg.getBoundingClientRect();
      tip.innerHTML = `<b>${r.name}</b>` + L.map(l => `${l.key === 'comb' ? 'Fire activity' : SENS[l.key]} ${fmt(pick(r, l.key))} / day`).join('<br>');
      tip.style.display = 'block'; const sx = ev.clientX - rr.left;
      tip.style.left = (sx + 14 + tip.offsetWidth > rr.width ? sx - tip.offsetWidth - 14 : sx + 14) + 'px'; tip.style.top = '8px'; hit.setAttribute('fill', css('--band'));
    });
    hit.addEventListener('mouseleave', () => { tip.style.display = 'none'; hit.setAttribute('fill', 'transparent'); });
  });
  return svg;
}
function mapSVG(a, W, P) {   // static map for the report (no basemap)
  const asp = Math.cos((st.s + st.n) / 2 * Math.PI / 180), dw = (st.e - st.w) * asp, dh = st.n - st.s;
  const m = { l: 56, r: 10, t: 10, b: 26 }, iw = W - m.l - m.r, sx = Math.min(iw / dw, 520 / dh), ih = dh * sx, Hh = ih + m.t + m.b;
  const svg = svgBase(W, Hh, 'Hotspot density and recorded fires'), x = lon => m.l + (lon - st.w) * asp * sx, y = lat => m.t + (st.n - lat) * sx;
  svg.append(el('rect', { x: m.l, y: m.t, width: dw * sx, height: ih, fill: P.bg, stroke: P.text3, 'stroke-width': 1 }));
  const clip = el('clipPath', { id: 'mapclip' }); clip.append(el('rect', { x: m.l, y: m.t, width: dw * sx, height: ih })); svg.append(clip);
  const g = el('g', { 'clip-path': 'url(#mapclip)' }); svg.append(g);
  const max = Math.max(1e-9, ...a.cells.map(c => c[2]));
  for (const [lo, la, v] of a.cells) g.append(el('rect', { x: x(lo), y: y(la + a.mc), width: Math.max(1, a.mc * asp * sx), height: Math.max(1, a.mc * sx), fill: P.heat, 'fill-opacity': (.12 + .8 * Math.log1p(v) / Math.log1p(max)).toFixed(2) }));
  a.fires.slice(0, 8).forEach(f => { g.append(el('circle', { cx: x(f.lon), cy: y(f.lat), r: 6, fill: 'none', stroke: P.text, 'stroke-width': 1.5 }));
    g.append(el('text', { x: x(f.lon) + 9, y: y(f.lat) + 4, fill: P.text, 'font-size': 11, stroke: P.bg, 'stroke-width': 3, 'paint-order': 'stroke' }, f.name)); });
  const tick = (lo, hi) => { const s = niceStep((hi - lo) * 1.2), out = []; for (let v = Math.ceil(lo / s) * s; v <= hi; v += s) out.push(+v.toFixed(3)); return out; };
  tick(st.w, st.e).forEach(v => svg.append(el('text', { x: x(v), y: Hh - 8, 'text-anchor': 'middle', fill: P.text3, 'font-size': 11 }, lonS(v).replace('.00', ''))));
  tick(st.s, st.n).forEach(v => svg.append(el('text', { x: m.l - 6, y: y(v) + 4, 'text-anchor': 'end', fill: P.text3, 'font-size': 11 }, latS(v).replace('.00', ''))));
  return svg;
}

// ---------------------------------------------------------------- map
let map, cellLayer, boxLayer, fireLayer, drawing = false;
function initMap() {
  if (map) return;
  map = L.map('map', { renderer: L.canvas(), zoomSnap: .25 });
  map.scrollWheelZoom.disable(); map.on('click', () => map.scrollWheelZoom.enable()); map.on('mouseout', () => map.scrollWheelZoom.disable());
  L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Base/MapServer/tile/{z}/{y}/{x}',
    { attribution: 'Basemap © Esri · Fire data: NASA FIRMS, USFS FPA FOD', maxZoom: 13, className: 'base' }).addTo(map);
  L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Reference/MapServer/tile/{z}/{y}/{x}', { maxZoom: 13, className: 'base' }).addTo(map);
  cellLayer = L.layerGroup().addTo(map); fireLayer = L.layerGroup().addTo(map); boxLayer = L.layerGroup().addTo(map);
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
function drawBox() { boxLayer.clearLayers(); L.rectangle([[st.s, st.w], [st.n, st.e]], { color: css('--text'), weight: 1.5, dashArray: '5 4', fill: false, interactive: false }).addTo(boxLayer); }
function fitRegion() { if (map) map.fitBounds([[st.s, st.w], [st.n, st.e]], { padding: [16, 16] }); }
const wikiLink = f => f.wiki_url ? `<a href="${esc(f.wiki_url)}" target="_blank" rel="noopener">Wikipedia</a>` :
  `<a href="https://en.wikipedia.org/w/index.php?search=${encodeURIComponent(f.name + ' Fire ' + f.start.slice(0, 4))}" target="_blank" rel="noopener">Search</a>`;
function drawMap(a) {
  cellLayer.clearLayers(); fireLayer.clearLayers();
  const max = Math.max(1e-9, ...a.cells.map(c => c[2])), col = css('--heat');
  for (const [lo, la, v] of a.cells) L.rectangle([[la, lo], [la + a.mc, lo + a.mc]], { stroke: false, fillColor: col, interactive: false, fillOpacity: .12 + .8 * Math.log1p(v) / Math.log1p(max) }).addTo(cellLayer);
  a.fires.slice(0, 25).forEach(f => L.circleMarker([f.lat, f.lon], { radius: Math.max(5, Math.min(15, Math.sqrt(f.acres) / 50)), color: css('--text'), weight: 1.5, fillOpacity: 0 })
    .bindTooltip(`<b>${esc(f.name)} Fire</b><br>${f.acres.toLocaleString()} acres · from ${f.start}`)
    .bindPopup(`<b>${esc(f.name)} Fire</b><br>${f.acres.toLocaleString()} acres · ${f.start} to ${f.end}<br>${esc(causeS(f.cause))} · ${wikiLink(f)}`).addTo(fireLayer));
  drawBox();
}

// ---------------------------------------------------------------- render
const tile = (k, v, n) => `<div class="tile"><div class="k">${k}</div><div class="v">${v}</div><div class="n">${n}</div></div>`;
function tilesHTML(a) {
  const on = a.per.filter(p => p.on); if (!on.length) return tile('No sensor selected', '–', 'Turn one on under More options');
  const out = [];
  if (a.harm && a.overall != null) {
    const r = a.longTerm ? a.overall / a.longTerm : null;
    out.push(tile('Fire activity', `${fmt(a.overall)} <small>hotspots / day</small>`,
      r == null ? 'MODIS units, per sampled day' : `${r >= 1.05 ? r.toFixed(1) + '× the' : r > .95 ? 'about the' : Math.round(r * 100) + '% of the'} long-term average here`));
  } else out.push(tile('Fire activity', on.map(p => fmt(p.rate)).join(' · ') + ' <small>/ day</small>', on.map(p => p.name).join(' · ') + (a.harm ? '' : ', raw detections')));
  const pk = peakOf(a);
  out.push(tile('Peak month', pk ? monLabel(pk.m) : '–', pk ? `${fmt(shown(a, pk))} / day` + (pk.fires[0] ? ` · ${esc(pk.fires[0].name)} Fire` : '') : 'no sampled days'));
  const f = a.fires[0];
  out.push(f ? tile('Largest recorded fire', `${esc(f.name)} <small>${f.acres.toLocaleString()} acres</small>`, `${f.start.slice(0, 7)} · ${esc(causeS(f.cause))} · ${wikiLink(f)}`) :
    tile('Largest recorded fire', '–', st.from > toDay('2020-12-31') ? 'Fire records run to 2020' : 'None matched in this selection'));
  out.push(tile('Hotspots', on.reduce((s, p) => s + p.hotspots, 0).toLocaleString(), `in ${Math.max(...on.map(p => p.days))} sampled days`));
  return out.join('');
}
function firesHTML(a, max = 8) {
  if (!FIRES.length) return '<tbody><tr><td class="empty">No fire records loaded.</td></tr></tbody>';
  if (!a.fires.length) return `<tbody><tr><td class="empty">No recorded fires of 100+ acres matched here${st.to > toDay('2020-12-31') ? ' (records run to 2020)' : ''}.</td></tr></tbody>`;
  return '<thead><tr><th>Fire</th><th>Started</th><th>Acres</th><th class="l">Cause</th><th>Read more</th></tr></thead><tbody>' +
    a.fires.slice(0, max).map(f => `<tr><td><b>${esc(f.name)}</b></td><td>${f.start}</td><td>${f.acres.toLocaleString()}</td><td class="l">${esc(causeS(f.cause))}</td><td>${wikiLink(f)}</td></tr>`).join('') + '</tbody>';
}
const chartW = id => Math.max(320, Math.round($(id).clientWidth || 600));
function render() {
  const P = pal(), a = A;
  for (const id of ['monthly', 'season']) $(id).innerHTML = '<div class="tip"></div>';
  $('monthly').prepend(monthlyChart(a, chartW('monthly'), P, true));
  $('season').prepend(seasonChart(a, chartW('season'), P, true));
  $('leg1').innerHTML = lines(a).map(l => `<span><i class="sw" style="background:var(${cssOf(l.color)})"></i>${l.name}</span>`).join('');
  $('chnote').textContent = st.view === 'one' && a.harm ? 'Hotspots per sampled day in MODIS units (MODIS and harmonized VIIRS averaged where both flew). Labels mark the largest recorded fires.' :
    `Each sensor separately${a.harm ? ', VIIRS in MODIS units' : ', VIIRS as raw detections'}.` + (a.agree && st.sens[0] && st.sens[1] ? ` Over months both were sampled, VIIRS ÷ MODIS = ${a.agree.toFixed(2)}×${a.harm ? ' (1× = perfect agreement)' : ''}.` : '');
}
function update() {
  if (!H) return;
  A = analyse(); writeHash();
  $('tiles').innerHTML = tilesHTML(A);
  render(); drawMap(A);
  $('mapnote').textContent = 'Darker = more fire per sampled day. Rings = recorded fires; click one for details.';
  $('fires').innerHTML = firesHTML(A);
  $('firecount').textContent = A.fires.length ? `${A.fires.length} fire${A.fires.length > 1 ? 's' : ''} of 100+ acres` : '';
  $('firenote').textContent = A.fires.length > 8 ? `The 8 largest are shown; the CSV download has all ${A.fires.length}.` :
    FIRES.length ? 'US Forest Service records (to 2020) matched to fire spikes in the satellite data.' : '';
}
addEventListener('hashchange', () => { if (H && readHash()) { syncInputs(); fitRegion(); update(); } });
let rz; addEventListener('resize', () => { clearTimeout(rz); rz = setTimeout(() => A && render(), 150); });

// ---------------------------------------------------------------- downloads
function save(name, data, type) {
  const url = URL.createObjectURL(data instanceof Blob ? data : new Blob([data], { type }));
  const a = document.createElement('a'); a.href = url; a.download = name; document.body.append(a); a.click(); a.remove(); setTimeout(() => URL.revokeObjectURL(url), 2000);
}
const stem = () => `fires_${(META.name || 'region').toLowerCase().replace(/\W+/g, '-')}_${dayStr(st.from)}_${dayStr(st.to)}`;
const svgText = svg => '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(svg);
const TITLES = { monthly: 'Fire activity per sampled day', season: 'Fire season (per sampled day, by calendar month)' };
const subtitle = () => `${META.name} · ${boxS()} · ${dayStr(st.from)} to ${dayStr(st.to)} · ${A.unit}`;
function framed(svg, P, title, sub) {
  const W = +svg.getAttribute('width'), H0 = +svg.getAttribute('height'), pad = 46, out = svgBase(W, H0 + pad + 18, title);
  out.append(el('rect', { width: W, height: H0 + pad + 18, fill: P.bg }));
  out.append(el('text', { x: 12, y: 20, fill: P.text, 'font-size': 14, 'font-weight': 650 }, title));
  out.append(el('text', { x: 12, y: 37, fill: P.text2, 'font-size': 11.5 }, sub));
  const g = el('g', { transform: `translate(0 ${pad})` }); [...svg.childNodes].forEach(n => n.nodeName !== 'title' && g.append(n)); out.append(g);
  let lx = 12; lines(A).forEach(l => { out.append(el('rect', { x: lx, y: H0 + pad + 6, width: 14, height: 3, rx: 1.5, fill: colorOf(P, l.color) }));
    out.append(el('text', { x: lx + 19, y: H0 + pad + 11, fill: P.text2, 'font-size': 11.5 }, l.name)); lx += 34 + l.name.length * 6.3; });
  return out;
}
const chartFor = (name, W, P) => name === 'monthly' ? monthlyChart(A, W, P, false) : seasonChart(A, W, P, false);
document.querySelectorAll('.dl').forEach(box => box.addEventListener('click', async ev => {
  const f = ev.target.dataset.f; if (!f || !A) return;
  const name = box.dataset.chart, svg = framed(chartFor(name, 900, LIGHT), LIGHT, TITLES[name], subtitle());
  if (f === 'svg') return save(`${stem()}_${name}.svg`, svgText(svg), 'image/svg+xml');
  const img = new Image(), url = URL.createObjectURL(new Blob([svgText(svg)], { type: 'image/svg+xml' }));
  await new Promise((ok, no) => { img.onload = ok; img.onerror = no; img.src = url; });
  const c = document.createElement('canvas'); c.width = svg.getAttribute('width') * 2; c.height = svg.getAttribute('height') * 2;
  const g = c.getContext('2d'); g.scale(2, 2); g.drawImage(img, 0, 0); URL.revokeObjectURL(url); c.toBlob(b => save(`${stem()}_${name}.png`, b), 'image/png');
}));
const hotspotRows = () => A.idx.map(i => [dayStr(H.day[i]), String(H.time[i]).padStart(4, '0'), H.lat[i].toFixed(4), H.lon[i].toFixed(4), SENS[H.sensor[i]], H.sat[i],
  H.night[i] ? 'N' : 'D', ['low', 'nominal', 'high'][H.conf[i]], H.frp[i].toFixed(1), H.pix[i], H.sensor[i] ? +H.weq[i].toFixed(4) : 1]);
const csvCell = v => /[",\n]/.test(String(v)) ? `"${String(v).replace(/"/g, '""')}"` : v;
function summaryJSON() {
  const pk = peakOf(A);
  return {
    dataset: { name: META.name, source: META.source, sampling: META.sampling, harmonization: META.harmonization, fires_source: FIRESRC && FIRESRC.source },
    selection: { bbox_wsen: [st.w, st.s, st.e, st.n], from: dayStr(st.from), to: dayStr(st.to), sensors: SENS.filter((_, j) => st.sens[j]), units: A.unit,
      min_confidence: ['all', 'nominal', 'high'][st.conf], link: location.href },
    fire_activity_per_sampled_day: A.overall, long_term_average_here: A.longTerm, peak_month: pk && monKey(pk.m), viirs_over_modis_same_months: A.agree,
    sensors: A.per.filter(p => p.on).map(p => ({ sensor: p.name, hotspots: p.hotspots, count: +p.count.toFixed(2), sampled_days: p.days, per_sampled_day: p.rate,
      frp_median_mw: p.frpMedian, night_share: p.count ? +(p.night / p.count).toFixed(3) : null })),
    monthly: A.series.map(s => ({ month: monKey(s.m), fire_activity: s.comb, modis: s.rate[0], viirs: s.rate[1], modis_days: s.days[0], viirs_days: s.days[1], fires_burning: s.fires.map(f => f.name) })),
    season: A.season.map(s => ({ month: s.name, fire_activity: s.comb, modis: s.rate[0], viirs: s.rate[1] })),
    recorded_fires: A.fires.map(({ name, start, end, acres, cause, state, lat, lon, wiki_url }) => ({ name, start, end, acres, cause, state, lat, lon, wiki_url })),
  };
}
function reportHTML() {
  const P = LIGHT, W = 860, sj = summaryJSON(), hm = META.harmonization;
  const sv = name => framed(chartFor(name, W, P), P, TITLES[name], subtitle()).outerHTML;
  return `<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Fire activity report — ${esc(META.name)}</title>
<style>body{margin:0;background:#f4f3f0;color:#0b0b0b;font:14px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif}main{max-width:900px;margin:0 auto;padding:28px 18px 48px;background:#fcfcfb}
h1{font-size:22px;margin:0 0 4px}h2{font-size:15px;margin:26px 0 8px}p,li{color:#52514e}table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums;font-size:13px}
th,td{padding:5px 8px;text-align:right;border-bottom:1px solid #e9e7e2}th:first-child,td:first-child{text-align:left}td.l,th.l{text-align:left}thead th{color:#52514e;border-bottom:1px solid #e2e0da}svg{max-width:100%;height:auto;display:block}
.k{color:#85837c;font-size:12.5px}.tiles{display:grid;grid-template-columns:repeat(4,1fr);gap:10px}.tile{border:1px solid #e2e0da;border-radius:9px;padding:10px 12px}.tile .k{color:#52514e;font-size:12px}.tile .v{font-size:19px;font-weight:650}.tile .v small{font-size:12px;font-weight:500;color:#52514e}.tile .n{color:#85837c;font-size:12px}
@media(max-width:700px){.tiles{grid-template-columns:1fr 1fr}}@media print{body{background:#fff}main{padding:0}h2{break-after:avoid}svg{break-inside:avoid}}</style></head><body><main>
<h1>Fire activity — ${esc(META.name)}</h1>
<p>${boxS()} · ${dayStr(st.from)} to ${dayStr(st.to)} · ${A.unit}${st.sens[0] && st.sens[1] ? '' : ' · ' + SENS.filter((_, j) => st.sens[j]).join('') + ' only'}</p>
<p class="k">Generated ${new Date().toISOString().slice(0, 16).replace('T', ' ')} UTC. Reopen this view: <a href="${esc(location.href)}">${esc(location.href)}</a></p>
<div class="tiles">${tilesHTML(A)}</div>
<h2>Fire activity over time</h2>${sv('monthly')}
<h2>Fire season</h2>${sv('season')}
<h2>Where</h2>${mapSVG(A, W, P).outerHTML}<p class="k">Darker cells = more fire per sampled day; rings = the largest recorded fires. No basemap.</p>
<h2>Recorded fires</h2><table>${firesHTML(A, 25)}</table>${A.fires.length > 25 ? `<p class="k">25 largest of ${A.fires.length}.</p>` : ''}
<h2>Method</h2><ul>
<li>NASA FIRMS active-fire detections${META.sampling ? `, sampled ${esc(META.sampling)}; counts cover sampled days only and rates are per sampled day` : ''}.</li>
${hm ? `<li>Harmonized to MODIS units: VIIRS pixels merged to ${hm.merge_km} km per overpass, then calibrated (${esc(hm.model)}, typical monthly error ×/÷ ${hm.typical_error}); fire activity averages MODIS and harmonized VIIRS where both flew. ${esc(hm.note)}</li>` : ''}
${FIRESRC ? `<li>Recorded fires: ${esc(FIRESRC.source)}; ${esc(FIRESRC.location)}.</li>` : ''}
<li>Long-term average for this area: ${fmt(sj.long_term_average_here)} hotspots per sampled day over the whole record.</li></ul>
</main></body></html>`;
}
document.querySelector('.downloads').addEventListener('click', ev => {
  const b = ev.target.closest('button'); if (!b || !A) return;
  const k = b.dataset.dl, s = stem();
  if (k === 'csv') save(s + '_hotspots.csv', 'date,time,lat,lon,sensor,satellite,daynight,confidence,frp,pixels,modis_equiv\n' + hotspotRows().map(r => r.join(',')).join('\n') + '\n', 'text/csv');
  else if (k === 'geojson') save(s + '_hotspots.geojson', JSON.stringify({ type: 'FeatureCollection', features: hotspotRows().map(r => ({ type: 'Feature', geometry: { type: 'Point', coordinates: [+r[3], +r[2]] },
    properties: { date: r[0], time: r[1], sensor: r[4], satellite: r[5], daynight: r[6], confidence: r[7], frp: +r[8], pixels: r[9], modis_equiv: r[10] } })) }), 'application/geo+json');
  else if (k === 'monthly') save(s + '_monthly.csv', 'month,fire_activity_per_day,modis_per_day,viirs_per_day,modis_sampled_days,viirs_sampled_days,units,fires_burning\n' +
    A.series.map(m => [monKey(m.m), m.comb ?? '', m.rate[0] ?? '', m.rate[1] ?? '', m.days[0], m.days[1], A.harm ? 'modis_equiv' : 'raw', csvCell(m.fires.map(f => f.name).join('; '))].join(',')).join('\n') + '\n', 'text/csv');
  else if (k === 'fires') save(s + '_recorded_fires.csv', 'name,start,end,acres,cause,state,lat,lon,wikipedia\n' +
    A.fires.map(f => [f.name, f.start, f.end, f.acres, f.cause, f.state, f.lat, f.lon, f.wiki_url || ''].map(csvCell).join(',')).join('\n') + '\n', 'text/csv');
  else if (k === 'json') save(s + '_analytics.json', JSON.stringify(summaryJSON(), null, 1), 'application/json');
  else if (k === 'report') save(s + '_report.html', reportHTML(), 'text/html');
});

boot();
