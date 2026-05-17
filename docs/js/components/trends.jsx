// trends.jsx — Global green / blue / total water-use chart (2000-2020).
// Deliberately minimal and matches the v1.5 working version, with only
// guards added against missing numeric fields so it always renders.

function Trends({ data }) {
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
              +{fmtInt(g(2020,"total") - g(2010,"total"))} km³/yr · matches the paper&rsquo;s headline.
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

window.Trends = Trends;
