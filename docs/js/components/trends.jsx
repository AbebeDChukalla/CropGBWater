// trends.jsx — Time-series chart of green vs blue vs total water, 2000-2020.

function Trends({ data }) {
  const years = data.summary.years;
  const g = (yr, key) => data.summary.global[String(yr)][`${key}_km3`];
  const change2010 = data.summary.change_2010_2020;
  const change2000 = data.summary.change_2000_2020;

  const W = 1100, H = 380, PADX = 60, PADY = 50;
  // Auto y-scale to a tidy round number above the max total.
  const maxVal = Math.ceil(Math.max(...years.map((y) => g(y, "total"))) / 1000) * 1000;
  const yRange = maxVal;
  const x = (y) => PADX + ((y - years[0]) / (years[years.length - 1] - years[0])) * (W - PADX * 2);
  const yScale = (v) => H - PADY - (v / yRange) * (H - PADY * 2);

  const series = (key) =>
    years.map((yr) => ({ x: x(yr), y: yScale(g(yr, key)), v: g(yr, key), yr }));
  const totalPts = series("total");
  const greenPts = series("green");
  const bluePts  = series("blue");

  const pathFrom = (pts) => pts.map((p, i) => `${i === 0 ? "M" : "L"}${p.x},${p.y}`).join(" ");
  const areaFrom = (pts) => `${pathFrom(pts)} L${pts[pts.length - 1].x},${H - PADY} L${pts[0].x},${H - PADY} Z`;

  // y-axis ticks at every 1000
  const ticks = [];
  for (let v = 0; v <= maxVal; v += 1000) ticks.push(v);

  return (
    <section id="trends" className="trend-module" data-screen-label="05 Trends">
      <div className="module-head">
        <div className="module-eyebrow">
          <span className="eyebrow-mark" />
          <span>Module 05 / Trends</span>
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
            <g key={`t${p.yr}`}>
              <circle cx={p.x} cy={p.y} r="4" className="pt pt-total" />
              <text x={p.x} y={p.y - 14} className="pt-label">{CGBW.fmt.km3(p.v)}</text>
            </g>
          ))}
          {greenPts.map((p) => <circle key={`g${p.yr}`} cx={p.x} cy={p.y} r="3" className="pt pt-green" />)}
          {bluePts.map((p)  => <circle key={`b${p.yr}`} cx={p.x} cy={p.y} r="3" className="pt pt-blue"  />)}

          <text x={x(years[years.length-1]) + 14} y={yScale(g(2020, "green")) + 4} className="ser-label green">
            Green water use · {CGBW.fmt.km3(g(2020, "green"))}
          </text>
          <text x={x(years[years.length-1]) + 14} y={yScale(g(2020, "blue")) + 4} className="ser-label blue">
            Blue water use · {CGBW.fmt.km3(g(2020, "blue"))}
          </text>
          <text x={x(years[years.length-1]) + 14} y={yScale(g(2020, "total")) - 10} className="ser-label total">
            Total · {CGBW.fmt.km3(g(2020, "total"))}&nbsp;km³
          </text>
        </svg>

        <div className="trend-callouts">
          <div className="callout">
            <div className="callout-label">Total change · 2010 → 2020</div>
            <div className="callout-val"><span className="num">{CGBW.fmt.pct(change2010.total_pct).replace("%","")}</span><i>%</i></div>
            <div className="callout-detail">
              +{CGBW.fmt.km3(g(2020,"total") - g(2010,"total"))} km³/yr · matches the paper&rsquo;s headline trend.
            </div>
          </div>
          <div className="callout">
            <div className="callout-label">Green water use share, 2020</div>
            <div className="callout-val">
              <span className="num">{(g(2020,"green") / g(2020,"total") * 100).toFixed(1)}</span><i>%</i>
            </div>
            <div className="callout-detail">
              Stable through the period — irrigation expansion offset by rainfed-area gains.
            </div>
          </div>
          <div className="callout">
            <div className="callout-label">Productivity, 2020</div>
            <div className="callout-val">
              <span className="num">{data.summary.global["2020"].productivity_kg_m3 ?? "—"}</span>
              <i>kg/m³</i>
            </div>
            <div className="callout-detail">
              {change2000.productivity_pct != null ?
                `${CGBW.fmt.pct(change2000.productivity_pct)} vs 2000 · more food per drop, but absolute use still climbs.`
                : "Tracked vs total water consumption."}
            </div>
          </div>
        </div>
      </div>
    </section>
  );
}

window.Trends = Trends;
