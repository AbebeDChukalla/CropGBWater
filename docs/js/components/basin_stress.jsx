// basin_stress.jsx — Comparison figure: Blue-water stress in 14 major river
// basins, 2010 vs 2020. Sourced from Chukalla et al. 2025, Nature Food,
// Figure 5 (Supplementary Source Data, MOESM5_ESM.xlsx).

function BasinStress() {
  const [data, setData] = React.useState(null);
  const [sort, setSort] = React.useState("delta"); // delta | ws2020 | ws2010

  React.useEffect(() => {
    fetch("data/basin_stress.json")
      .then((r) => r.json())
      .then((d) => setData(d))
      .catch((e) => console.warn("basin_stress.json failed", e));
  }, []);

  if (!data) return null;

  const sorted = React.useMemo(() => {
    const arr = [...data.basins];
    if (sort === "delta")  arr.sort((a, b) => (b.delta_pct ?? -999) - (a.delta_pct ?? -999));
    if (sort === "ws2020") arr.sort((a, b) => b.ws_2020 - a.ws_2020);
    if (sort === "ws2010") arr.sort((a, b) => b.ws_2010 - a.ws_2010);
    return arr;
  }, [data, sort]);

  // Scale bars to the max across both years for visual consistency.
  const maxVal = Math.max(...data.basins.flatMap((b) => [b.ws_2010, b.ws_2020]));

  const deltaTone = (d) => {
    if (d == null) return "low";
    if (d >= 50)   return "very-high";
    if (d >= 15)   return "high";
    if (d >= 0)    return "mid";
    return "neg";
  };

  return (
    <section className="basin-stress" aria-labelledby="bs-h">
      <div className="module-eyebrow">
        <span className="eyebrow-mark" />
        <span>Comparison · river-basin blue water stress</span>
      </div>
      <h2 id="bs-h" className="module-title">
        Where blue water stress jumped most between 2010 and 2020.
      </h2>
      <p className="module-sub">
        Blue water consumption-to-availability stress in 14 major river basins
        from the underlying paper (Figure 5). Higher values mean a larger share
        of the basin&rsquo;s renewable blue water is consumed by agriculture.
        Niger, Murray and Mekong saw the steepest increases; only Mississippi
        and Orange saw small decreases.
      </p>

      <div className="bs-toolbar">
        <span className="bs-toolbar-label">Sort by</span>
        <button className={`pill ${sort==="delta"  ? "on" : ""}`} onClick={() => setSort("delta")}>Δ 2010 → 2020</button>
        <button className={`pill ${sort==="ws2020" ? "on" : ""}`} onClick={() => setSort("ws2020")}>2020 stress</button>
        <button className={`pill ${sort==="ws2010" ? "on" : ""}`} onClick={() => setSort("ws2010")}>2010 stress</button>
        <span style={{ flex: 1 }} />
        <span className="bs-legend">
          <span><span className="dot" style={{ background: "var(--blueSoft)" }} /> 2010</span>
          <span><span className="dot" style={{ background: "var(--blue)" }} /> 2020</span>
        </span>
      </div>

      <ol className="bs-list">
        {sorted.map((b, i) => (
          <li key={b.basin} className="bs-row" style={{ animationDelay: `${i * 40}ms` }}>
            <span className="bs-pos">{String(i + 1).padStart(2, "0")}</span>
            <span className="bs-name">{b.basin}</span>
            <div className="bs-bars">
              <div className="bs-bar bs-bar-2010" title={`2010: ${b.ws_2010}`}>
                <i style={{ width: `${(b.ws_2010 / maxVal) * 100}%` }} />
                <span className="bs-bar-label">{b.ws_2010.toFixed(1)}</span>
              </div>
              <div className="bs-bar bs-bar-2020" title={`2020: ${b.ws_2020}`}>
                <i style={{ width: `${(b.ws_2020 / maxVal) * 100}%` }} />
                <span className="bs-bar-label">{b.ws_2020.toFixed(1)}</span>
              </div>
            </div>
            <span className={`bs-delta tone-${deltaTone(b.delta_pct)}`}>
              {b.delta_pct == null ? "—" : `${b.delta_pct > 0 ? "+" : ""}${b.delta_pct.toFixed(1)}%`}
            </span>
          </li>
        ))}
      </ol>

      <div className="bs-source">
        Source: {data.source}. Stress values are unitless ratios scaled ×100.
      </div>
    </section>
  );
}

window.BasinStress = BasinStress;
