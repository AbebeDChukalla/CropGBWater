/* util.js — plain JS helpers + data loader. Loads before any JSX. */

window.CGBW = {
  data: null,
  paths: {
    summary:    "data/summary.json",
    crops:      "data/crops.json",
    crop_groups:"data/crop_groups.json",
    countries:  "data/countries.json",
    continents: "data/continents.json",
    meta:       "data/meta.json",
    country:    (iso) => `data/country_detail/${iso}.json`,
    world_topo: "https://unpkg.com/world-atlas@2/countries-110m.json",
  },
};

// ─── Loading ─────────────────────────────────────────────────────────
CGBW.loadAll = async function () {
  const fetchJSON = (url) => fetch(url).then((r) => {
    if (!r.ok) throw new Error(`${url} → ${r.status}`);
    return r.json();
  });
  const fetchOpt = (url) => fetch(url).then((r) => r.ok ? r.json() : null).catch(() => null);
  const [summary, crops, groups, countries, continents, meta] = await Promise.all([
    fetchJSON(CGBW.paths.summary),
    fetchJSON(CGBW.paths.crops),
    fetchJSON(CGBW.paths.crop_groups),
    fetchJSON(CGBW.paths.countries),
    fetchOpt(CGBW.paths.continents),
    fetchJSON(CGBW.paths.meta),
  ]);
  CGBW.data = { summary, crops, groups, countries, continents, meta };
  return CGBW.data;
};

CGBW.loadCountry = async function (iso) {
  const r = await fetch(CGBW.paths.country(iso));
  if (!r.ok) return null;
  return r.json();
};

// ─── Formatting ──────────────────────────────────────────────────────
CGBW.fmt = {
  // 7720.31 -> "7,720"  (no decimals for big km³)
  km3: (v) => v == null ? "—" : Math.round(v).toLocaleString("en-US"),
  // 7720.31 -> "7,720.3"
  km3p: (v) => v == null ? "—" : v.toLocaleString("en-US", { maximumFractionDigits: 1 }),
  // 9.7 -> "+9.7%"
  pct: (v) => v == null ? "—" : `${v > 0 ? "+" : ""}${v.toFixed(1)}%`,
  // 7720 -> "7,720"
  num: (v, digits = 0) => v == null ? "—"
    : v.toLocaleString("en-US", { minimumFractionDigits: digits, maximumFractionDigits: digits }),
  // 1455 Mha
  area: (v) => v == null ? "—" : v.toLocaleString("en-US", { maximumFractionDigits: 0 }) + " Mha",
  // 1.14 kg/m³
  prod: (v) => v == null ? "—" : v.toFixed(2),
};

// ─── Color ramps ─────────────────────────────────────────────────────
// Sequential ramp from light paper2 to deep ink (used per-metric)
CGBW.ramps = {
  total: ["#f3efe6","#e7d8a8","#cb9e6f","#a96940","#7c3a23","#3c1306"],
  green: ["#f3efe6","#dfe6cc","#b9cf99","#7ea36b","#4b7d49","#2d5a3d"],
  blue:  ["#f3efe6","#cfdbe9","#9eb9d4","#6fa1cf","#3a719f","#1e4a7a"],
};

CGBW.colorForValue = function (v, breaks, ramp) {
  if (v == null || isNaN(v)) return "var(--paper2)";
  for (let i = breaks.length - 1; i >= 0; i--) {
    if (v >= breaks[i]) return ramp[i + 1];
  }
  return ramp[0];
};

// Compute quantile breaks from an array of numeric values
CGBW.quantileBreaks = function (values, n = 5) {
  const sorted = values.filter((v) => v != null && !isNaN(v)).sort((a, b) => a - b);
  if (sorted.length === 0) return [];
  const breaks = [];
  for (let i = 1; i <= n; i++) {
    const idx = Math.floor((i / (n + 1)) * sorted.length);
    breaks.push(sorted[Math.min(idx, sorted.length - 1)]);
  }
  return breaks;
};

// ─── ISO3 normalisation (UN names → Natural Earth ISO codes) ─────────
// Names in the source data don't always match the Natural Earth feature
// properties exactly. We key by ISO3 instead, which the world-atlas TopoJSON
// exposes as `id` (numeric ISO 3166-1) — so we need a numeric→alpha-3 map.
// We use a small lookup table generated at runtime from the country bundle.
CGBW.iso3ByName = function () {
  if (!CGBW.data) return {};
  return Object.fromEntries(CGBW.data.countries.countries.map((c) => [c.name, c.iso3]));
};

// ─── CountUp helper (no library) ─────────────────────────────────────
// Used inline in JSX to animate from 0 → target on mount.
// Returns a string with the correct precision via the formatter you pass.
// Usage:
//   const v = CGBW.useCountUp(7720, 900, (n) => Math.round(n).toLocaleString());
CGBW.useCountUp = function (target, durationMs = 800, format = (n) => Math.round(n).toString()) {
  const [val, setVal] = React.useState(format(0));
  React.useEffect(() => {
    if (target == null || isNaN(target)) { setVal("—"); return; }
    let raf, start;
    const animate = (t) => {
      if (!start) start = t;
      const p = Math.min(1, (t - start) / durationMs);
      const eased = 1 - Math.pow(1 - p, 3); // ease-out-cubic
      setVal(format(target * eased));
      if (p < 1) raf = requestAnimationFrame(animate);
      else setVal(format(target));   // ensure exact final
    };
    raf = requestAnimationFrame(animate);
    return () => cancelAnimationFrame(raf);
  }, [target, durationMs]);
  return val;
};

// world-atlas/countries-110m uses ISO 3166-1 numeric codes as the feature id.
// Map numeric → alpha-3 for matching against our country bundles.
CGBW.numericToAlpha3 = {
  "004":"AFG","008":"ALB","012":"DZA","024":"AGO","032":"ARG","036":"AUS","040":"AUT","044":"BHS","050":"BGD","051":"ARM","052":"BRB","056":"BEL","064":"BTN","068":"BOL","070":"BIH","072":"BWA","076":"BRA","084":"BLZ","096":"BRN","100":"BGR","104":"MMR","108":"BDI","112":"BLR","116":"KHM","120":"CMR","124":"CAN","140":"CAF","144":"LKA","148":"TCD","152":"CHL","156":"CHN","158":"TWN","170":"COL","178":"COG","180":"COD","188":"CRI","191":"HRV","192":"CUB","196":"CYP","203":"CZE","204":"BEN","208":"DNK","214":"DOM","218":"ECU","222":"SLV","226":"GNQ","231":"ETH","232":"ERI","233":"EST","242":"FJI","246":"FIN","250":"FRA","254":"GUF","258":"PYF","262":"DJI","266":"GAB","268":"GEO","270":"GMB","275":"PSE","276":"DEU","288":"GHA","300":"GRC","320":"GTM","324":"GIN","328":"GUY","332":"HTI","336":"VAT","340":"HND","344":"HKG","348":"HUN","352":"ISL","356":"IND","360":"IDN","364":"IRN","368":"IRQ","372":"IRL","376":"ISR","380":"ITA","384":"CIV","388":"JAM","392":"JPN","398":"KAZ","400":"JOR","404":"KEN","408":"PRK","410":"KOR","414":"KWT","417":"KGZ","418":"LAO","422":"LBN","426":"LSO","428":"LVA","430":"LBR","434":"LBY","438":"LIE","440":"LTU","442":"LUX","446":"MAC","450":"MDG","454":"MWI","458":"MYS","462":"MDV","466":"MLI","470":"MLT","478":"MRT","484":"MEX","496":"MNG","498":"MDA","499":"MNE","504":"MAR","508":"MOZ","512":"OMN","516":"NAM","524":"NPL","528":"NLD","540":"NCL","548":"VUT","554":"NZL","558":"NIC","562":"NER","566":"NGA","578":"NOR","586":"PAK","591":"PAN","598":"PNG","600":"PRY","604":"PER","608":"PHL","616":"POL","620":"PRT","624":"GNB","626":"TLS","630":"PRI","634":"QAT","642":"ROU","643":"RUS","646":"RWA","662":"LCA","678":"STP","682":"SAU","686":"SEN","688":"SRB","694":"SLE","702":"SGP","703":"SVK","704":"VNM","705":"SVN","706":"SOM","710":"ZAF","716":"ZWE","724":"ESP","728":"SSD","729":"SDN","748":"SWZ","752":"SWE","756":"CHE","760":"SYR","762":"TJK","764":"THA","768":"TGO","780":"TTO","784":"ARE","788":"TUN","792":"TUR","795":"TKM","800":"UGA","804":"UKR","807":"MKD","818":"EGY","826":"GBR","834":"TZA","840":"USA","854":"BFA","858":"URY","860":"UZB","862":"VEN","882":"WSM","887":"YEM","894":"ZMB"
};
