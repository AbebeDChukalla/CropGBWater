// crops.jsx — Crop explorer. Switches between group view and individual-crop chip grid.

function CropExplorer({ data }) {
  const [mode, setMode] = React.useState("groups"); // "groups" | "crops"
  const [activeCrop, setActiveCrop] = React.useState(null);

  return (
    <section id="crops" className="crops-module" data-screen-label="04 Crop Breakdown">
      <div className="module-head">
        <div className="module-eyebrow">
          <span className="eyebrow-mark" />
          <span>Module 04 / Crop composition</span>
        </div>
        <h2 className="module-title">
          Cereals alone account for nearly half of all crop water.
          Vegetables and fruits draw the most blue water per kilo.
        </h2>
        <p className="module-sub">
          Toggle between {data.meta.n_groups} crop groups and the full list of {data.meta.n_crops} crops.
          Numbers are 2020 global totals — green / blue split shown below each tile.
        </p>
      </div>

      <div className="crop-mode" role="tablist">
        <button onClick={() => setMode("groups")} className={mode === "groups" ? "on" : ""}>
          By group ({data.meta.n_groups})
        </button>
        <button onClick={() => setMode("crops")} className={mode === "crops" ? "on" : ""}>
          By crop ({data.meta.n_crops})
        </button>
      </div>

      {mode === "groups" ? (
        <div className="crops-grid">
          {data.groups.groups.map((g, i) => (
            <article key={g.id} className="crop-card" data-active="0">
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
                  <span><span className="dot dot-green" />{g.green_pct}% green</span>
                  <span><span className="dot dot-blue"  />{g.blue_pct}% blue</span>
                </div>
              </div>
              <div className="crop-meta">
                <span>{CGBW.fmt.km3p(g.total_km3)} km³/yr</span>
                <span>{CGBW.fmt.km3p(g.green_km3)} + {CGBW.fmt.km3p(g.blue_km3)} </span>
              </div>
            </article>
          ))}
        </div>
      ) : (
        <CropChipGrid data={data} active={activeCrop} setActive={setActiveCrop} />
      )}
    </section>
  );
}

function CropChipGrid({ data, active, setActive }) {
  const crops = data.crops.crops;
  // Top 12 by 2020 total become the highlight cards; rest as chips.
  const featured = crops.slice(0, 12);
  const rest = crops.slice(12);
  return (
    <>
      <div className="crops-grid">
        {featured.map((c, i) => (
          <article key={c.code}
                   className="crop-card"
                   data-active={active === c.code ? "1" : "0"}
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
              <span>{CGBW.fmt.pct(c.trend_2010_2020_pct)} since 2010</span>
            </div>
          </article>
        ))}
      </div>
      <div style={{ marginTop: 36 }}>
        <span style={{ fontFamily: "var(--font-mono)", fontSize: 11, color: "var(--ink40)",
                       letterSpacing: ".06em", textTransform: "uppercase" }}>
          Remaining {rest.length} crops · sorted by 2020 total
        </span>
        <div className="ctrl-body" style={{ marginTop: 14, flexWrap: "wrap", gap: 6 }}>
          {rest.map((c) => (
            <button key={c.code}
                    className="chip"
                    title={`${c.name}: ${CGBW.fmt.km3p(c["2020"].total_km3)} km³/yr · ${c.group}`}>
              <span className="chip-dot" style={{ background:
                (c.green_share_pct || 0) >= 80 ? "var(--green)"
              : (c.blue_share_pct  || 0) >= 30 ? "var(--blue)"
              : "var(--accent)" }} />
              {c.name}
              <span style={{ fontFamily: "var(--font-mono)", fontSize: 10.5, color: "var(--ink40)", marginLeft: 4 }}>
                {CGBW.fmt.km3p(c["2020"].total_km3)}
              </span>
            </button>
          ))}
        </div>
      </div>
    </>
  );
}

window.CropExplorer = CropExplorer;
