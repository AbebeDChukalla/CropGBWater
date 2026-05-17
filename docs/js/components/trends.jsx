// trends.jsx — Time-series chart of green vs blue vs total water (global),
// PLUS a top-7 crop trend comparison with metric switcher (water / area /
// yield / WF per ton), PLUS callouts. Interactive, animated, downloadable.

const TOP7_CODES = ["WHEA", "RICE", "MAIZ", "SOYB", "SUGC", "COTT", "OPAL"];
const TOP7_COLORS = {
  WHEA: "#c9883b", RICE: "#1e4a7a", MAIZ: "#9b9b3d", SOYB: "#7e5ca0",
  SUGC: "#4a7d4d", COTT: "#b9683a", OPAL: "#3a719f",
};

const TREND_METRICS = [
  { id: "total",    label: "Total water use (km³/yr)",     unit: "km³" },
  { id: "green",    label: "Green water use (km³/yr)",     unit: "km³" },
  { id: "blue",     label: "Blue water use (km³/yr)",      unit: "km³" },
  { id: "area_irr", label: "Irrigated harvested area (Mha)", unit: "Mha" },
  { id: "area_rf",  label: "Rainfed harvested area (Mha)",   unit: "Mha" },
  { id: "yield_irr",label: "Yield, irrigated (t/ha)",        unit: "t/ha" },
  { id: "yield_rf", label: "Yield, rainfed (t/ha)",          unit: "t/ha" },
  { id: "wf",       label: "Total water footprint (m³/t)",   unit: "m³/t" },
];

function Trends({ data }) {
  return (
    <section id="trends" className="trend-module" data-screen-label="05 Trends">
      <GlobalTrends data={data} />
      <CropTrends data={data} />
    </section>
  );
}

// ─── 1. Global green vs blue vs total chart (original) ────────────────
function GlobalTrends({ data }) {
  const years = data.summary.years;
  const g = (yr, key) => data.summary.global[String(yr)][`${key}_km3`];
  const change2010 = data.summary.change_2010_2020 || {};
  const change2000 = data.summary.change_2000_2020 || {};

  const W = 1100, H = 380, PADX = 60, PADY = 50;
  const maxVal = Math.ceil(Math.max(...years.map((y) => g(y, "total"))) / 1000) * 1000;
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

  const downloadCSV = () => {
    const rows = [["year","green_km3","blue_km3","total_km3"]];
    years.forEach((y) => rows.push([y, g(y,"green"), g(y,"blue"), g(y,"total")]));
    triggerDownload("cropgbwater_global_trends.csv", rows.map((r) => r.join(",")).join("\n"));
  };

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

      <div className="trend-toolbar">
        <button className="pill" onClick={downloadCSV} title="Download underlying values">⬇ CSV</button>
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
              +{CGBW.fmt.km3(g(2020,"total") - g(2010,"total"))} km³/yr · matches the paper&rsquo;s headline.
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
              {change2000.productivity_pct != null
                ? `${CGBW.fmt.pct(change2000.productivity_pct)} vs 2000 · more food per drop, but absolute use still climbs.`
                : "Tracked vs total water consumption."}
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}

// ─── 2. Top-7 crops: multi-metric trend comparison ───────────────────
function CropTrends({ data }) {
  const [metricId, setMetricId] = React.useState("total");
  const [selected, setSelected] = React.useState(() => new Set(TOP7_CODES));
  const [hover, setHover] = React.useState(null);

  const cropsIndex = React.useMemo(() => {
    const m = {};
    (data.crops.crops || []).forEach((c) => { m[c.code] = c; });
    return m;
  }, [data]);

  const metric = TREND_METRICS.find((m) => m.id === metricId) || TREND_METRICS[0];
  const years = data.summary.years;

  // Compute a per-year value for a given crop and metric.
  const getVal = (crop, year) => {
    if (!crop) return null;
    const yKey = String(year);
    switch (metricId) {
      case "total":
      case "green":
      case "blue":  return crop[yKey]?.[`${metricId}_km3`] ?? null;
      case "area_irr":  return crop.area_history?.[yKey]?.area_irrigated_Mha ?? null;
      case "area_rf":   return crop.area_history?.[yKey]?.area_rainfed_Mha ?? null;
      case "yield_irr": {
        const a = crop.area_history?.[yKey]?.area_irrigated_Mha ?? 0;
        const p = crop.production_history?.[yKey]?.production_irrigated_Mt ?? 0;
        return a > 0 ? (p / a) : null;     // Mt/Mha == t/ha
      }
      case "yield_rf": {
        const a = crop.area_history?.[yKey]?.area_rainfed_Mha ?? 0;
        const p = crop.production_history?.[yKey]?.production_rainfed_Mt ?? 0;
        return a > 0 ? (p / a) : null;
      }
      case "wf": {
        const km3 = crop[yKey]?.total_km3 ?? 0;
        const Mt  = crop.production_history?.[yKey]?.production_Mt ?? 0;
        return Mt > 0 ? (km3 * 1e9 / (Mt * 1e6)) : null;  // m³ per tonne
      }
      default: return null;
    }
  };

  // Chart dims
  const W = 1100, H = 360, PADX = 60, PADY = 60;
  const allVals = TOP7_CODES.flatMap((code) =>
    years.map((y) => getVal(cropsIndex[code], y))).filter((v) => v != null);
  const maxVal = allVals.length ? Math.max(...allVals) : 1;
  const niceMax = niceTop(maxVal);
  const x = (y) => PADX + ((y - years[0]) / (years[years.length - 1] - years[0])) * (W - PADX * 2);
  const yScale = (v) => H - PADY - (v / niceMax) * (H - PADY * 2);
  const ticks = niceTicks(0, niceMax, 5);

  const downloadCSV = () => {
    const rows = [["crop", ...years]];
    TOP7_CODES.forEach((code) => {
      const crop = cropsIndex[code];
      rows.push([crop?.name || code, ...years.map((y) => {
        const v = getVal(crop, y); return v != null ? v.toFixed(3) : "";
      })]);
    });
    triggerDownload(
      `cropgbwater_top7_${metricId}.csv`,
      rows.map((r) => r.join(",")).join("\n")
    );
  };

  const toggle = (code) => {
    const s = new Set(selected);
    if (s.has(code)) s.delete(code); else s.add(code);
    if (s.size === 0) s.add(code);   // never leave empty
    setSelected(s);
  };

  return (
    <div className="trend-block trend-block-crops">
      <div className="module-eyebrow">
        <span className="eyebrow-mark" />
        <span>Top-7 crops · interactive trend comparison</span>
      </div>
      <h3 className="trend-h">
        How {selected.size} of the world&rsquo;s biggest crops moved between 2000, 2010 and 2020.
      </h3>

      <div className="trend-toolbar trend-toolbar-2">
        <span className="trend-toolbar-label">Metric</span>
        <select className="trend-select" value={metricId} onChange={(e) => setMetricId(e.target.value)}>
          {TREND_METRICS.map((m) => <option key={m.id} value={m.id}>{m.label}</option>)}
        </select>
        <span style={{ flex: 1 }} />
        <button className="pill" onClick={downloadCSV}>⬇ CSV</button>
      </div>

      <div className="crop-chip-row">
        {TOP7_CODES.map((code) => {
          const c = cropsIndex[code];
          const on = selected.has(code);
          return (
            <button key={code}
                    className={`crop-chip ${on ? "on" : ""}`}
                    style={{ "--c": TOP7_COLORS[code] }}
                    onClick={() => toggle(code)}>
              <span className="crop-chip-dot" />
              {c?.name || code}
            </button>
          );
        })}
      </div>

      <div className="trend-frame trend-frame-crops">
        <svg className="trend-svg" viewBox={`0 0 ${W} ${H}`} preserveAspectRatio="xMidYMid meet"
             onMouseLeave={() => setHover(null)}>
          {/* gridlines + Y labels */}
          {ticks.map((v) => (
            <g key={v} className="grid">
              <line x1={PADX} x2={W - PADX} y1={yScale(v)} y2={yScale(v)} />
              <text x={PADX - 12} y={yScale(v) + 4} className="grid-label">
                {fmtNum(v)}
              </text>
            </g>
          ))}
          {/* X axis: years + invisible hover columns */}
          {years.map((y, i) => (
            <g key={y}>
              <text x={x(y)} y={H - 18} className="axis-label">{y}</text>
              <rect x={x(y) - (W - PADX * 2) / years.length / 2}
                    y={PADY / 2}
                    width={(W - PADX * 2) / years.length}
                    height={H - PADY}
                    fill="transparent"
                    onMouseEnter={() => setHover(y)} />
            </g>
          ))}

          {/* Series for each selected crop */}
          {TOP7_CODES.filter((c) => selected.has(c)).map((code) => {
            const crop = cropsIndex[code];
            const pts = years.map((y) => {
              const v = getVal(crop, y);
              return v == null ? null : { x: x(y), y: yScale(v), v, yr: y };
            });
            const filtered = pts.filter(Boolean);
            if (filtered.length < 2) return null;
            const path = filtered.map((p, i) => `${i === 0 ? "M" : "L"}${p.x},${p.y}`).join(" ");
            return (
              <g key={code}>
                <path d={path} fill="none" stroke={TOP7_COLORS[code]} strokeWidth="2.2"
                      className="crop-line" strokeLinecap="round" strokeLinejoin="round" />
                {filtered.map((p) => (
                  <circle key={p.yr} cx={p.x} cy={p.y} r="4"
                          fill={TOP7_COLORS[code]} className="crop-pt" />
                ))}
                {/* End-of-line label */}
                <text x={filtered[filtered.length - 1].x + 8}
                      y={filtered[filtered.length - 1].y + 4}
                      className="crop-line-label"
                      style={{ fill: TOP7_COLORS[code] }}>
                  {crop?.name || code}
                </text>
              </g>
            );
          })}

          {/* Hover guide line + tooltip */}
          {hover != null && (
            <g>
              <line x1={x(hover)} x2={x(hover)} y1={PADY / 2} y2={H - PADY}
                    stroke="rgba(20,20,15,0.18)" strokeDasharray="2 3" />
              <text x={x(hover)} y={PADY / 2 - 6} textAnchor="middle"
                    fontFamily="var(--font-mono)" fontSize="10.5" fill="var(--ink70)">
                {hover}
              </text>
            </g>
          )}
        </svg>

        {/* hover stats panel */}
        <div className="trend-hover-panel">
          <div className="trend-hover-title">
            {metric.label}
            {hover != null && <span className="trend-hover-year"> · {hover}</span>}
          </div>
          <div className="trend-hover-list">
            {TOP7_CODES.filter((c) => selected.has(c)).map((code) => {
              const crop = cropsIndex[code];
              const v = hover != null ? getVal(crop, hover) : getVal(crop, 2020);
              return (
                <div key={code} className="trend-hover-row">
                  <span className="trend-hover-swatch" style={{ background: TOP7_COLORS[code] }} />
                  <span className="trend-hover-name">{crop?.name || code}</span>
                  <span className="trend-hover-val">
                    {v == null ? "—" : `${fmtNum(v)} ${metric.unit}`}
                  </span>
                </div>
              );
            })}
          </div>
          {hover == null && (
            <div className="trend-hover-hint">Hover the chart to compare a specific year.</div>
          )}
        </div>
      </div>
    </div>
  );
}

// ─── helpers ─────────────────────────────────────────────────────────
function niceTop(v) {
  if (v <= 0) return 1;
  const mag = Math.pow(10, Math.floor(Math.log10(v)));
  for (const m of [1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 7.5, 10]) {
    if (mag * m >= v) return mag * m;
  }
  return Math.ceil(v / mag) * mag;
}
function niceTicks(min, max, n) {
  const step = (max - min) / n;
  const out = [];
  for (let i = 0; i <= n; i++) out.push(min + i * step);
  return out;
}
function fmtNum(v) {
  if (v == null) return "—";
  if (v >= 1000)  return Math.round(v).toLocaleString();
  if (v >= 100)   return v.toFixed(0);
  if (v >= 10)    return v.toFixed(1);
  return v.toFixed(2);
}
function triggerDownload(filename, text) {
  const blob = new Blob([text], { type: "text/csv;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = filename;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1500);
}

window.Trends = Trends;
