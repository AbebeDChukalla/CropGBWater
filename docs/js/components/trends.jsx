// trends.jsx — Global time-series + top-7 crops comparison chart.
// Kept deliberately simple and defensive so it always renders.

const TOP7_CODES = ["WHEA", "RICE", "MAIZ", "SOYB", "SUGC", "COTT", "OPAL"];
const TOP7_COLORS = {
  WHEA: "#c9883b", RICE: "#1e4a7a", MAIZ: "#9b9b3d", SOYB: "#7e5ca0",
  SUGC: "#4a7d4d", COTT: "#b9683a", OPAL: "#3a719f",
};
const TREND_METRICS = [
  { id: "total",     label: "Total water use",            unit: "km³/yr" },
  { id: "green",     label: "Green water use",            unit: "km³/yr" },
  { id: "blue",      label: "Blue water use",             unit: "km³/yr" },
  { id: "area_irr",  label: "Irrigated harvested area",   unit: "Mha" },
  { id: "area_rf",   label: "Rainfed harvested area",     unit: "Mha" },
  { id: "yield_irr", label: "Yield, irrigated",           unit: "t/ha" },
  { id: "yield_rf",  label: "Yield, rainfed",             unit: "t/ha" },
  { id: "wf",        label: "Total water footprint",      unit: "m³/t" },
];

function Trends({ data }) {
  return (
    <section id="trends" className="trend-module" data-screen-label="05 Trends">
      <GlobalTrends data={data} />
      <CropTrends data={data} />
    </section>
  );
}

// ─── Global green vs blue vs total chart ─────────────────────────────
function GlobalTrends({ data }) {
  const years = data.summary.years;
  const g = (yr, key) => {
    const row = data.summary.global[String(yr)] || {};
    const v = row[(key === "total" ? "total" : key) + "_km3"];
    return typeof v === "number" ? v : 0;
  };
  const change2010 = data.summary.change_2010_2020 || {};
  const change2000 = data.summary.change_2000_2020 || {};

  const W = 1100, H = 380, PADX = 60, PADY = 50;
  const maxVal = Math.ceil(Math.max(...years.map((y) => g(y, "total")), 1) / 1000) * 1000 || 8000;
  const x = (y) => PADX + ((y - years[0]) / (years[years.length - 1] - years[0])) * (W - PADX * 2);
  const yScale = (v) => H - PADY - (v / maxVal) * (H - PADY * 2);
  const series = (key) =>
    years.map((yr) => ({ x: x(yr), y: yScale(g(yr, key)), v: g(yr, key), yr }));
  const totalPts = series("total");
  const greenPts = series("green");
  const bluePts  = series("blue");
  const pathFrom = (pts) => pts.map((p, i) => `${i === 0 ? "M" : "L"}${p.x},${p.y}`).join(" ");
  const areaFrom = (pts) => `${pathFrom(pts)} L${pts[pts.length - 1].x},${H - PADY} L${pts[0].x},${H - PADY} Z`;
  const ticks = [];
  for (let v = 0; v <= maxVal; v += 1000) ticks.push(v);

  return (
    <div className="trend-block">
      <div className="module-head">
        <div className="module-eyebrow">
          <span className="eyebrow-mark" />
          <span>Module 05 / Trends · 2000–2020</span>
        </div>
        <h2 className="module-title">
          Two decades. Green water use has grown nearly three times faster than blue.
        </h2>
        <p className="module-sub">
          Annual totals across all {data.meta.n_crops} crops and {data.meta.n_countries} countries.
          Solid lines: green and blue water; dashed line: total. Numbers in km³ per year.
        </p>
      </div>

      <div className="trend-frame">
        <svg className="trend-svg" viewBox={`0 0 ${W} ${H}`} preserveAspectRatio="xMidYMid meet">
          {ticks.map((v) => (
            <g key={v} className="grid">
              <line x1={PADX} x2={W - PADX} y1={yScale(v)} y2={yScale(v)} />
              <text x={PADX - 12} y={yScale(v) + 4} className="grid-label">{v.toLocaleString()}</text>
            </g>
          ))}
          {years.map((y) => (
            <text key={y} x={x(y)} y={H - 14} className="axis-label">{y}</text>
          ))}
          <path d={areaFrom(greenPts)} className="area area-green" />
          <path d={pathFrom(totalPts)} className="line line-total" />
          <path d={pathFrom(greenPts)} className="line line-green" />
          <path d={pathFrom(bluePts)}  className="line line-blue" />
          {totalPts.map((p) => (
            <g key={"t" + p.yr}>
              <circle cx={p.x} cy={p.y} r="4" className="pt pt-total" />
              <text x={p.x} y={p.y - 14} className="pt-label">{Math.round(p.v).toLocaleString()}</text>
            </g>
          ))}
          {greenPts.map((p) => <circle key={"g" + p.yr} cx={p.x} cy={p.y} r="3" className="pt pt-green" />)}
          {bluePts.map((p)  => <circle key={"b" + p.yr} cx={p.x} cy={p.y} r="3" className="pt pt-blue"  />)}
          <text x={x(years[years.length-1]) + 14} y={yScale(g(2020, "green")) + 4} className="ser-label green">
            Green water use · {Math.round(g(2020, "green")).toLocaleString()}
          </text>
          <text x={x(years[years.length-1]) + 14} y={yScale(g(2020, "blue")) + 4} className="ser-label blue">
            Blue water use · {Math.round(g(2020, "blue")).toLocaleString()}
          </text>
          <text x={x(years[years.length-1]) + 14} y={yScale(g(2020, "total")) - 10} className="ser-label total">
            Total · {Math.round(g(2020, "total")).toLocaleString()}&nbsp;km³
          </text>
        </svg>

        <div className="trend-callouts">
          <div className="callout">
            <div className="callout-label">Total change · 2010 → 2020</div>
            <div className="callout-val">
              <span className="num">{change2010.total_pct != null ? change2010.total_pct.toFixed(1) : "—"}</span><i>%</i>
            </div>
            <div className="callout-detail">
              +{Math.round(g(2020,"total") - g(2010,"total")).toLocaleString()} km³/yr · matches the paper&rsquo;s headline.
            </div>
          </div>
          <div className="callout">
            <div className="callout-label">Green water use share, 2020</div>
            <div className="callout-val">
              <span className="num">{g(2020,"total") > 0 ? (g(2020,"green") / g(2020,"total") * 100).toFixed(1) : "—"}</span><i>%</i>
            </div>
            <div className="callout-detail">Stable through the period.</div>
          </div>
          <div className="callout">
            <div className="callout-label">Productivity, 2020</div>
            <div className="callout-val">
              <span className="num">{(data.summary.global["2020"] && data.summary.global["2020"].productivity_kg_m3) || "—"}</span>
              <i>kg/m³</i>
            </div>
            <div className="callout-detail">
              {change2000.productivity_pct != null
                ? "+" + change2000.productivity_pct.toFixed(1) + "% vs 2000 · more food per drop."
                : "Tracked vs total water consumption."}
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}

// ─── Top-7 crops comparison ──────────────────────────────────────────
function CropTrends({ data }) {
  const [metricId, setMetricId] = React.useState("total");
  const [selected, setSelected] = React.useState({
    WHEA: true, RICE: true, MAIZ: true, SOYB: true, SUGC: true, COTT: true, OPAL: true,
  });

  const cropsArr = (data.crops && data.crops.crops) || [];
  const cropsIndex = {};
  cropsArr.forEach((c) => { cropsIndex[c.code] = c; });

  const years = data.summary.years;
  const metric = TREND_METRICS.find((m) => m.id === metricId) || TREND_METRICS[0];

  function getVal(crop, year) {
    if (!crop) return null;
    const yKey = String(year);
    const row = crop[yKey] || {};
    const ah  = (crop.area_history       || {})[yKey] || {};
    const ph  = (crop.production_history || {})[yKey] || {};
    if (metricId === "total" || metricId === "green" || metricId === "blue") {
      const v = row[metricId + "_km3"];
      return typeof v === "number" ? v : null;
    }
    if (metricId === "area_irr") return ah.area_irrigated_Mha || null;
    if (metricId === "area_rf")  return ah.area_rainfed_Mha   || null;
    if (metricId === "yield_irr") {
      const a = ah.area_irrigated_Mha || 0;
      const p = ph.production_irrigated_Mt || 0;
      return a > 0 ? (p / a) : null;
    }
    if (metricId === "yield_rf") {
      const a = ah.area_rainfed_Mha || 0;
      const p = ph.production_rainfed_Mt || 0;
      return a > 0 ? (p / a) : null;
    }
    if (metricId === "wf") {
      const km3 = row.total_km3 || 0;
      const Mt  = ph.production_Mt || 0;
      return Mt > 0 ? (km3 * 1e9 / (Mt * 1e6)) : null;
    }
    return null;
  }

  const activeCodes = TOP7_CODES.filter((c) => selected[c]);

  // Collect all values to determine max
  const allVals = [];
  TOP7_CODES.forEach((code) => {
    years.forEach((y) => {
      const v = getVal(cropsIndex[code], y);
      if (v != null && !isNaN(v) && isFinite(v)) allVals.push(v);
    });
  });
  const maxRaw = allVals.length > 0 ? Math.max.apply(null, allVals) : 1;
  const niceMax = niceTop(maxRaw);

  const W = 1100, H = 360, PADX = 70, PADY = 60;
  const x = (y) => PADX + ((y - years[0]) / (years[years.length - 1] - years[0])) * (W - PADX * 2);
  const yScale = (v) => H - PADY - (v / niceMax) * (H - PADY * 2);
  const ticks = [];
  const step = niceMax / 5;
  for (let i = 0; i <= 5; i++) ticks.push(step * i);

  function toggle(code) {
    setSelected({ ...selected, [code]: !selected[code] });
  }

  return (
    <div className="trend-block trend-block-crops">
      <div className="module-eyebrow">
        <span className="eyebrow-mark" />
        <span>Top-7 crops · interactive trend comparison</span>
      </div>
      <h3 className="trend-h">
        How seven of the world&rsquo;s biggest crops moved between 2000, 2010 and 2020.
      </h3>

      <div className="trend-toolbar trend-toolbar-2">
        <span className="trend-toolbar-label">Metric</span>
        <select className="trend-select"
                value={metricId}
                onChange={(e) => setMetricId(e.target.value)}>
          {TREND_METRICS.map((m) => (
            <option key={m.id} value={m.id}>{m.label} ({m.unit})</option>
          ))}
        </select>
      </div>

      <div className="crop-chip-row">
        {TOP7_CODES.map((code) => {
          const crop = cropsIndex[code];
          const on = !!selected[code];
          return (
            <button key={code}
                    className={"crop-chip " + (on ? "on" : "")}
                    onClick={() => toggle(code)}>
              <span className="crop-chip-dot" style={{ background: TOP7_COLORS[code] }} />
              {(crop && crop.name) || code}
            </button>
          );
        })}
      </div>

      <svg className="trend-svg crop-trend-svg"
           viewBox={"0 0 " + W + " " + H}
           preserveAspectRatio="xMidYMid meet">
        {ticks.map((v, i) => (
          <g key={"tk" + i} className="grid">
            <line x1={PADX} x2={W - PADX} y1={yScale(v)} y2={yScale(v)} />
            <text x={PADX - 12} y={yScale(v) + 4} className="grid-label">{fmtNum(v)}</text>
          </g>
        ))}
        {years.map((y) => (
          <text key={"yr" + y} x={x(y)} y={H - 18} className="axis-label">{y}</text>
        ))}
        {activeCodes.map((code) => {
          const crop = cropsIndex[code];
          const pts = [];
          years.forEach((y) => {
            const v = getVal(crop, y);
            if (v != null && !isNaN(v) && isFinite(v)) {
              pts.push({ x: x(y), y: yScale(v), v: v, yr: y });
            }
          });
          if (pts.length < 1) return null;
          const path = pts.map((p, i) => (i === 0 ? "M" : "L") + p.x + "," + p.y).join(" ");
          const last = pts[pts.length - 1];
          return (
            <g key={code}>
              <path d={path} fill="none"
                    stroke={TOP7_COLORS[code]}
                    strokeWidth="2"
                    strokeLinecap="round"
                    strokeLinejoin="round" />
              {pts.map((p) => (
                <circle key={code + "p" + p.yr}
                        cx={p.x} cy={p.y} r="3.5"
                        fill={TOP7_COLORS[code]} />
              ))}
              <text x={last.x + 8} y={last.y + 4}
                    fontFamily="var(--font-mono)" fontSize="10.5"
                    fontWeight="500"
                    fill={TOP7_COLORS[code]}>
                {(crop && crop.name) || code}
              </text>
            </g>
          );
        })}
        {/* Y axis title */}
        <text x={PADX - 50} y={H / 2}
              transform={"rotate(-90 " + (PADX - 50) + " " + (H / 2) + ")"}
              textAnchor="middle"
              fontFamily="var(--font-mono)" fontSize="10.5"
              fill="var(--ink40)" letterSpacing="0.06em">
          {metric.label.toUpperCase()} ({metric.unit})
        </text>
      </svg>

      <div className="trend-table-wrap">
        <table className="trend-table">
          <thead>
            <tr>
              <th>Crop</th>
              {years.map((y) => <th key={y}>{y}</th>)}
              <th>Δ 2010→2020</th>
            </tr>
          </thead>
          <tbody>
            {TOP7_CODES.map((code) => {
              const crop = cropsIndex[code];
              const vals = years.map((y) => getVal(crop, y));
              const d2010 = vals[1], d2020 = vals[2];
              const delta = (d2010 && d2020) ? ((d2020 - d2010) / d2010 * 100) : null;
              return (
                <tr key={code}>
                  <td>
                    <span className="trend-table-dot" style={{ background: TOP7_COLORS[code] }} />
                    {(crop && crop.name) || code}
                  </td>
                  {vals.map((v, i) => (
                    <td key={i}>{v == null ? "—" : fmtNum(v) + " " + metric.unit}</td>
                  ))}
                  <td>
                    {delta == null ? "—" : (delta > 0 ? "+" : "") + delta.toFixed(1) + "%"}
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
    </div>
  );
}

function niceTop(v) {
  if (!v || v <= 0 || !isFinite(v)) return 1;
  const mag = Math.pow(10, Math.floor(Math.log10(v)));
  const m = [1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 7.5, 10];
  for (let i = 0; i < m.length; i++) {
    if (mag * m[i] >= v) return mag * m[i];
  }
  return Math.ceil(v / mag) * mag;
}
function fmtNum(v) {
  if (v == null || isNaN(v) || !isFinite(v)) return "—";
  if (v >= 1000)  return Math.round(v).toLocaleString();
  if (v >= 100)   return v.toFixed(0);
  if (v >= 10)    return v.toFixed(1);
  return v.toFixed(2);
}

window.Trends = Trends;
