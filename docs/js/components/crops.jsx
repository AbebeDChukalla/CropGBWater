// crops.jsx — group donut + production/area/yield strip + per-crop grid.

// Pleasant categorical palette ranging from green → blue → warm
const GROUP_COLORS = [
  "#2d5a3d", "#4a7d4d", "#7ea36b",
  "#9b9b3d", "#c9883b", "#b9683a",
  "#7e5ca0", "#3a719f", "#1e4a7a", "#4e4e4e",
];

function CropExplorer({ data }) {
  const [mode, setMode] = React.useState("groups");
  const [activeCrop, setActiveCrop] = React.useState(null);
  const [hoverSlice, setHoverSlice] = React.useState(null);

  return (
    <section className="crops">
      <div className="module-head">
        <div className="module-eyebrow">
          <span className="eyebrow-mark" />
          <span>Module 04 / Crop composition · production · area · yield</span>
        </div>
        <h2 className="module-title">
          Cereals dominate water use. Vegetables and fruits drink the most blue water per kilo.
        </h2>
        <p className="module-sub">
          Toggle between {data.meta.n_groups} crop groups and the full list of {data.meta.n_crops} crops.
          Numbers are 2020 globals — green / blue split below each tile,
          production · harvested area · yield · irrigated share at top.
        </p>
      </div>

      {/* Production / Area / Yield / Irr-share strip */}
      <ProductionStrip data={data} />

      <div className="crop-mode" role="tablist">
        <button onClick={() => setMode("groups")} className={mode === "groups" ? "on" : ""}>
          By group ({data.meta.n_groups})
        </button>
        <button onClick={() => setMode("crops")} className={mode === "crops" ? "on" : ""}>
          By crop ({data.meta.n_crops})
        </button>
      </div>

      <div className="crops-layout">
        <DonutPanel groups={data.groups.groups}
                    hover={hoverSlice}
                    onHover={setHoverSlice}
                    total={data.summary.global["2020"].total_km3} />

        {mode === "groups" ? (
          <div className="crops-grid">
            {data.groups.groups.map((g, i) => (
              <article key={g.id} className="crop-card"
                       onMouseEnter={() => setHoverSlice(i)}
                       onMouseLeave={() => setHoverSlice(null)}
                       data-active={hoverSlice === i ? "1" : "0"}
                       style={{ animationDelay: `${i * 35}ms` }}>
                <div className="crop-card-head">
                  <span className="crop-rank">{String(i + 1).padStart(2, "0")}</span>
                  <span className="crop-name">{g.name}</span>
                </div>
                <div className="crop-share">
                  <span className="num">{g.share_pct.toFixed(1)}</span><i>% of global</i>
                </div>
                <div className="crop-split">
                  <span className="crop-split-bar">
                    <i className="seg-green" style={{ flexBasis: `${g.green_pct}%` }} />
                    <i className="seg-blue"  style={{ flexBasis: `${g.blue_pct}%` }} />
                  </span>
                  <div className="crop-split-labels">
                    <span><span className="dot dot-green" />{g.green_pct}% green use</span>
                    <span><span className="dot dot-blue"  />{g.blue_pct}% blue use</span>
                  </div>
                </div>
                <div className="crop-meta">
                  <span>{CGBW.fmt.km3p(g.total_km3)} km³/yr</span>
                </div>
              </article>
            ))}
          </div>
        ) : (
          <CropChipGrid data={data} active={activeCrop} setActive={setActiveCrop} />
        )}
      </div>
    </section>
  );
}

// ─── Production / Area / Yield / Irrigated share strip ───────────────
function ProductionStrip({ data }) {
  const g = data.summary.global["2020"];
  const c = data.summary.change_2010_2020;
  return (
    <div className="prod-strip">
      <div className="prod-cell">
        <span className="prod-cell-label">Crop production, 2020</span>
        <span className="prod-cell-val"><span className="num">{CGBW.fmt.km3(g.production_Mt)}</span><i>Mt</i></span>
        <span className="prod-cell-meta">{CGBW.fmt.pct(c.production_pct)} vs 2010</span>
      </div>
      <div className="prod-cell">
        <span className="prod-cell-label">Harvested area, 2020</span>
        <span className="prod-cell-val"><span className="num">{CGBW.fmt.km3(g.area_Mha)}</span><i>Mha</i></span>
        <span className="prod-cell-meta">{CGBW.fmt.pct(c.area_pct)} vs 2010</span>
      </div>
      <div className="prod-cell">
        <span className="prod-cell-label">Average yield, 2020</span>
        <span className="prod-cell-val"><span className="num">{g.yield_ton_ha?.toFixed(2)}</span><i>t/ha</i></span>
        <span className="prod-cell-meta">{CGBW.fmt.pct(c.yield_pct)} vs 2010</span>
      </div>
      <div className="prod-cell" style={{ borderRight: 0 }}>
        <span className="prod-cell-label">Irrigated cropland share</span>
        <span className="prod-cell-val"><span className="num">{g.irrigated_share_area_pct?.toFixed(1)}</span><i>%</i></span>
        <div className="irrf-strip" aria-label="Rainfed vs irrigated area share">
          <i className="irrf-rainfed" style={{ width: `${100 - (g.irrigated_share_area_pct || 0)}%` }} />
          <i className="irrf-irrig"   style={{ width: `${g.irrigated_share_area_pct || 0}%` }} />
        </div>
        <div className="irrf-legend">
          <span>{(100 - (g.irrigated_share_area_pct || 0)).toFixed(1)}% rainfed</span>
          <span>{(g.irrigated_share_area_pct || 0).toFixed(1)}% irrigated</span>
        </div>
      </div>
    </div>
  );
}

// ─── Donut chart (pure SVG, no library) ──────────────────────────────
function DonutPanel({ groups, hover, onHover, total }) {
  // Compute slice geometry
  const sum = groups.reduce((s, g) => s + (g.total_km3 || 0), 0);
  let cumAngle = -Math.PI / 2; // start at 12 o'clock
  const slices = groups.map((g, i) => {
    const frac = (g.total_km3 || 0) / sum;
    const a0 = cumAngle;
    const a1 = cumAngle + frac * Math.PI * 2;
    cumAngle = a1;
    return { ...g, a0, a1, frac, color: GROUP_COLORS[i % GROUP_COLORS.length] };
  });

  const cx = 100, cy = 100, rOuter = 86, rInner = 54;

  const arcPath = (a0, a1, ro, ri) => {
    const large = a1 - a0 > Math.PI ? 1 : 0;
    const x0 = cx + ro * Math.cos(a0), y0 = cy + ro * Math.sin(a0);
    const x1 = cx + ro * Math.cos(a1), y1 = cy + ro * Math.sin(a1);
    const x2 = cx + ri * Math.cos(a1), y2 = cy + ri * Math.sin(a1);
    const x3 = cx + ri * Math.cos(a0), y3 = cy + ri * Math.sin(a0);
    return `M ${x0} ${y0} A ${ro} ${ro} 0 ${large} 1 ${x1} ${y1} L ${x2} ${y2} A ${ri} ${ri} 0 ${large} 0 ${x3} ${y3} Z`;
  };

  const hoveredGroup = hover != null ? slices[hover] : null;

  return (
    <div className="donut-panel">
      <div className="donut-panel-title">Crop-group share of total water use · 2020</div>

      <svg className="donut-svg" viewBox="0 0 200 200" role="img" aria-label="Donut chart of crop group share">
        {slices.map((s, i) => (
          <path key={s.id}
                d={arcPath(s.a0, s.a1, rOuter, rInner)}
                fill={s.color}
                className="donut-slice"
                style={{ animationDelay: `${i * 60}ms`,
                         transform: hover === i ? "scale(1.04)" : "scale(1)" }}
                onMouseEnter={() => onHover(i)}
                onMouseLeave={() => onHover(null)}>
            <title>{s.name}: {s.share_pct}%</title>
          </path>
        ))}

        {/* center label */}
        <text x={cx} y={cy - 4} textAnchor="middle" className="donut-center-num">
          {hoveredGroup ? hoveredGroup.share_pct.toFixed(1) : Math.round(total).toLocaleString("en-US")}
        </text>
        <text x={cx} y={cy + 14} textAnchor="middle" className="donut-center-lab">
          {hoveredGroup ? `${hoveredGroup.name} · %` : "km³/yr total"}
        </text>
      </svg>

      <div className="donut-legend">
        {slices.map((s, i) => (
          <div key={s.id} className="donut-legend-row"
               style={{ animationDelay: `${300 + i * 30}ms` }}
               onMouseEnter={() => onHover(i)}
               onMouseLeave={() => onHover(null)}>
            <span className="donut-legend-swatch" style={{ background: s.color }} />
            <span className="donut-legend-name">{s.name}</span>
            <span className="donut-legend-val">{CGBW.fmt.km3p(s.total_km3)}</span>
            <span className="donut-legend-pct">{s.share_pct.toFixed(1)}%</span>
          </div>
        ))}
      </div>
    </div>
  );
}

function CropChipGrid({ data, active, setActive }) {
  const crops = data.crops.crops;
  const featured = crops.slice(0, 16);
  return (
    <div className="crops-grid">
      {featured.map((c, i) => (
        <article key={c.code}
                 className="crop-card"
                 data-active={active === c.code ? "1" : "0"}
                 style={{ animationDelay: `${i * 30}ms` }}
                 onMouseEnter={() => setActive(c.code)}
                 onMouseLeave={() => setActive(null)}>
          <div className="crop-card-head">
            <span className="crop-rank">{String(i + 1).padStart(2, "0")}</span>
            <span className="crop-name">{c.name}</span>
          </div>
          <div className="crop-share">
            <span className="num">{(c.share_2020_pct || 0).toFixed(1)}</span><i>% of global</i>
          </div>
          <div className="crop-split">
            <span className="crop-split-bar">
              <i className="seg-green" style={{ flexBasis: `${c.green_share_pct || 0}%` }} />
              <i className="seg-blue"  style={{ flexBasis: `${c.blue_share_pct  || 0}%` }} />
            </span>
            <div className="crop-split-labels">
              <span><span className="dot dot-green" />{Math.round(c.green_share_pct || 0)}% green</span>
              <span><span className="dot dot-blue"  />{Math.round(c.blue_share_pct  || 0)}% blue</span>
            </div>
          </div>
          <div className="crop-meta">
            <span>{c.group}</span>
            <span>{c.area_2020_Mha != null ? `${c.area_2020_Mha.toFixed(1)} Mha` : "—"}</span>
          </div>
          <div className="crop-meta" style={{ borderTop: 0, paddingTop: 0 }}>
            <span>{c.production_2020_Mt ? `${CGBW.fmt.km3p(c.production_2020_Mt)} Mt prod` : "—"}</span>
            <span>{c.yield_2020_ton_ha ? `${c.yield_2020_ton_ha} t/ha` : "—"}</span>
          </div>
        </article>
      ))}
    </div>
  );
}

window.CropExplorer = CropExplorer;
