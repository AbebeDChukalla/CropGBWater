// trends.jsx — Global green / blue / total water-use chart (2000-2020)
// + a Top-5 crops yield + water-productivity comparison.
// Deliberately minimal and table-based so it always renders.

const TOP5_CODES = ["RICE", "MAIZ", "WHEA", "SOYB", "SUGC"];
const TOP5_DOTS  = {
  RICE: "#1e4a7a", MAIZ: "#9b9b3d", WHEA: "#c9883b",
  SOYB: "#7e5ca0", SUGC: "#4a7d4d",
};

function Trends({ data }) {
  return (
    <React.Fragment>
      <GlobalTrend data={data} />
      <Top5CropComparison data={data} />
    </React.Fragment>
  );
}

function GlobalTrend({ data }) {
  const years = (data.summary && data.summary.years) || [2000, 2010, 2020];

  const g = (yr, key) => {
    const row = (data.summary && data.summary.global && data.summary.global[String(yr)]) || {};
    const v = row[key + "_km3"];
    return typeof v === "number" && isFinite(v) ? v : 0;
  };
  const change2010 = (data.summary && data.summary.change_2010_2020) || {};
  const change2000 = (data.summary && data.summary.change_2000_2020) || {};

  const W = 1100, H = 380, PADX = 60, PADY = 50;
  const totals = years.map((y) => g(y, "total"));
  const maxRaw = Math.max.apply(null, totals.concat([1]));
  const maxVal = Math.ceil(maxRaw / 1000) * 1000 || 8000;

  const x = (y) => PADX + ((y - years[0]) / (years[years.length - 1] - years[0])) * (W - PADX * 2);
  const yScale = (v) => H - PADY - (v / maxVal) * (H - PADY * 2);

  const series = (key) => years.map((yr) => ({
    x: x(yr), y: yScale(g(yr, key)), v: g(yr, key), yr,
  }));
  const totalPts = series("total");
  const greenPts = series("green");
  const bluePts  = series("blue");

  const pathFrom = (pts) => pts.map((p, i) => (i === 0 ? "M" : "L") + p.x + "," + p.y).join(" ");
  const areaFrom = (pts) =>
    pathFrom(pts) + " L" + pts[pts.length - 1].x + "," + (H - PADY) +
                    " L" + pts[0].x + "," + (H - PADY) + " Z";

  const ticks = [];
  for (let v = 0; v <= maxVal; v += 1000) ticks.push(v);

  const fmtInt = (v) => Math.round(v).toLocaleString();
  const fmtPct = (v) => (v == null || isNaN(v)) ? "—" : (v > 0 ? "+" : "") + v.toFixed(1) + "%";

  const greenShare2020 = g(2020, "total") > 0
    ? (g(2020, "green") / g(2020, "total") * 100).toFixed(1)
    : "—";

  return (
    <section id="trends" className="trend-module" data-screen-label="05 Trends">
      <div className="module-head">
        <div className="module-eyebrow">
          <span className="eyebrow-mark" />
          <span>Module 05 / Trends · 2000–2020</span>
        </div>
        <h2 className="module-title">
          Two decades. Green water use has grown nearly three times faster than blue.
        </h2>
        <p className="module-sub">
          Annual totals across all {data.meta && data.meta.n_crops}&nbsp;crops and {data.meta && data.meta.n_countries}&nbsp;countries.
          Solid lines: green and blue water; dashed line: total. Numbers in km³ per year.
        </p>
      </div>

      <div className="trend-frame">
        <svg className="trend-svg" viewBox={"0 0 " + W + " " + H} preserveAspectRatio="xMidYMid meet">
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
              <text x={p.x} y={p.y - 14} className="pt-label">{fmtInt(p.v)}</text>
            </g>
          ))}
          {greenPts.map((p) => <circle key={"g" + p.yr} cx={p.x} cy={p.y} r="3" className="pt pt-green" />)}
          {bluePts.map((p)  => <circle key={"b" + p.yr} cx={p.x} cy={p.y} r="3" className="pt pt-blue"  />)}
          <text x={x(years[years.length-1]) + 14} y={yScale(g(2020, "green")) + 4} className="ser-label green">
            Green water use · {fmtInt(g(2020, "green"))}
          </text>
          <text x={x(years[years.length-1]) + 14} y={yScale(g(2020, "blue")) + 4} className="ser-label blue">
            Blue water use · {fmtInt(g(2020, "blue"))}
          </text>
          <text x={x(years[years.length-1]) + 14} y={yScale(g(2020, "total")) - 10} className="ser-label total">
            Total · {fmtInt(g(2020, "total"))}&nbsp;km³
          </text>
        </svg>

        <div className="trend-callouts">
          <div className="callout">
            <div className="callout-label">Total change · 2010 → 2020</div>
            <div className="callout-val">
              <span className="num">{change2010.total_pct != null ? change2010.total_pct.toFixed(1) : "—"}</span>
              <i>%</i>
            </div>
            <div className="callout-detail">
              +{fmtInt(g(2020,"total") - g(2010,"total"))} km³/yr.
            </div>
          </div>
          <div className="callout">
            <div className="callout-label">Green water use share, 2020</div>
            <div className="callout-val">
              <span className="num">{greenShare2020}</span><i>%</i>
            </div>
            <div className="callout-detail">Stable through the period.</div>
          </div>
          <div className="callout">
            <div className="callout-label">Productivity, 2020</div>
            <div className="callout-val">
              <span className="num">
                {(data.summary && data.summary.global && data.summary.global["2020"]
                  && data.summary.global["2020"].productivity_kg_m3) || "—"}
              </span>
              <i>kg/m³</i>
            </div>
            <div className="callout-detail">
              {change2000.productivity_pct != null
                ? fmtPct(change2000.productivity_pct) + " vs 2000 · more food per drop."
                : "Tracked vs total water consumption."}
            </div>
          </div>
        </div>
      </div>
    </section>
  );
}

// ─── Top-5 crops: yield (RF/IRR) + water productivity (WP) ────────────
function Top5CropComparison({ data }) {
  const years = (data.summary && data.summary.years) || [2000, 2010, 2020];
  const cropsArr = (data.crops && data.crops.crops) || [];
  const cropsIdx = {};
  cropsArr.forEach((c) => { if (c && c.code) cropsIdx[c.code] = c; });

  // For one crop + one year, derive RF yield, IRR yield, total WP.
  function getMetrics(crop, year) {
    if (!crop) return { yr: null, yi: null, wp: null };
    const yKey = String(year);
    const ah = (crop.area_history       || {})[yKey] || {};
    const ph = (crop.production_history || {})[yKey] || {};
    const wat = crop[yKey] || {};

    const aRf  = num(ah.area_rainfed_Mha);
    const aIrr = num(ah.area_irrigated_Mha);
    const pRf  = num(ph.production_rainfed_Mt);
    const pIrr = num(ph.production_irrigated_Mt);
    const prod = num(ph.production_Mt);
    const km3  = num(wat.total_km3);

    // Mt / Mha == t/ha
    const yr = aRf  > 0 ? pRf  / aRf  : null;
    const yi = aIrr > 0 ? pIrr / aIrr : null;
    // production_Mt × 1e9 kg  /  total_km3 × 1e9 m³  = kg/m³
    const wp = km3 > 0 ? prod / km3 : null;
    return { yr, yi, wp };
  }

  // Build a small table for one metric.
  function MetricTable({ title, metricKey, unit, dp }) {
    const fmt = (v) => v == null || isNaN(v) ? "—" : v.toFixed(dp);
    return (
      <div className="cmp-block">
        <h4 className="cmp-title">{title} <span className="cmp-unit">({unit})</span></h4>
        <table className="cmp-table">
          <thead>
            <tr>
              <th>Crop</th>
              {years.map((y) => <th key={y}>{y}</th>)}
              <th>Δ 2010→2020</th>
            </tr>
          </thead>
          <tbody>
            {TOP5_CODES.map((code) => {
              const crop = cropsIdx[code];
              const vals = years.map((y) => getMetrics(crop, y)[metricKey]);
              const d10 = vals[1], d20 = vals[2];
              const delta = (d10 != null && d20 != null && d10 !== 0)
                ? (d20 - d10) / d10 * 100 : null;
              return (
                <tr key={code}>
                  <td>
                    <span className="cmp-dot" style={{ background: TOP5_DOTS[code] }} />
                    {(crop && crop.name) || code}
                  </td>
                  {vals.map((v, i) => <td key={i}>{fmt(v)}</td>)}
                  <td className="cmp-delta">
                    {delta == null ? "—" : (delta > 0 ? "+" : "") + delta.toFixed(1) + "%"}
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
    );
  }

  return (
    <section className="cmp-section" data-screen-label="05b Top-5 crops comparison">
      <div className="module-eyebrow" style={{ marginBottom: 6 }}>
        <span className="eyebrow-mark" />
        <span>Top-5 crops · yield + water productivity</span>
      </div>
      <h2 className="module-title" style={{ marginTop: 8, fontSize: "clamp(22px, 2.4vw, 32px)" }}>
        Five major crops, three decades — what changed.
      </h2>
      <p className="module-sub" style={{ marginBottom: 18 }}>
        Yield is computed as production ÷ harvested area, separately for rainfed and
        irrigated systems. Water productivity (WP) is production ÷ total water
        consumption (kg per m³).
      </p>

      <div className="cmp-grid">
        <MetricTable title="Yield, rainfed"   metricKey="yr" unit="t/ha"  dp={2} />
        <MetricTable title="Yield, irrigated" metricKey="yi" unit="t/ha"  dp={2} />
        <MetricTable title="Water productivity" metricKey="wp" unit="kg/m³" dp={2} />
      </div>
    </section>
  );
}

function num(v) {
  const n = typeof v === "number" ? v : parseFloat(v);
  return Number.isFinite(n) ? n : 0;
}

window.Trends = Trends;
