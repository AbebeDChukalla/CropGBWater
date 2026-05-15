// ranking.jsx — Continental rollup strip + country ranking table.

function CountryRanking({ data, onSelectCountry }) {
  const [sortBy, setSortBy] = React.useState("total");
  const [continent, setContinent] = React.useState("all");
  const [limit, setLimit] = React.useState(20);
  const countries = data.countries.countries;
  const continents = (data.continents && data.continents.continents) || [];
  const max = Math.max(...countries.map((c) => c.total_km3 || 0));

  const filtered = React.useMemo(() => {
    return continent === "all" ? countries : countries.filter((c) => c.continent === continent);
  }, [countries, continent]);

  const sorted = React.useMemo(() => {
    const arr = [...filtered];
    if (sortBy === "trend") arr.sort((a, b) => (b.trend_2010_2020_pct || 0) - (a.trend_2010_2020_pct || 0));
    else if (sortBy === "blue")  arr.sort((a, b) => (b.blue_pct || 0) - (a.blue_pct || 0));
    else if (sortBy === "green") arr.sort((a, b) => (b.green_pct || 0) - (a.green_pct || 0));
    else arr.sort((a, b) => (b.total_km3 || 0) - (a.total_km3 || 0));
    return arr.slice(0, limit);
  }, [filtered, sortBy, limit]);

  return (
    <section className="ranking">
      <div className="module-head">
        <div className="module-eyebrow">
          <span className="eyebrow-mark" />
          <span>Module 05 / Continental &amp; country ranking</span>
        </div>
        <h2 className="module-title">
          India and China each draw nearly a thousand km³ a year. Together with the US, they use a third of global crop water.
        </h2>
        <p className="module-sub">
          Top strip shows continental rollups. The list below ranks countries by 2020 total
          consumption — bars show the green / blue split, trend is the change since 2010.
          Click any row for the full per-country crop breakdown.
        </p>
      </div>

      {/* Continental strip */}
      <div className="continent-strip">
        {continents.map((c, i) => (
          <button key={c.name}
                  className="continent-card"
                  onClick={() => setContinent(continent === c.name ? "all" : c.name)}
                  style={{
                    animationDelay: `${i * 60}ms`,
                    background: continent === c.name ? "var(--paper2)" : "var(--paper)",
                    cursor: "pointer",
                  }}>
            <span className="continent-card-label">{c.name} · {c.n_countries}</span>
            <span className="continent-card-val">
              <span className="num">{CGBW.fmt.km3p(c.total_km3)}</span><i>km³</i>
            </span>
            <span className="continent-card-bar">
              <i style={{ background: "var(--green)", width: `${c.green_pct}%`, height: "100%" }} />
              <i style={{ background: "var(--blue)",  width: `${c.blue_pct }%`, height: "100%" }} />
            </span>
            <span className="continent-card-meta">
              <span>{c.green_pct}%g · {c.blue_pct}%b</span>
              <span>{c.share_pct}% world</span>
            </span>
          </button>
        ))}
      </div>

      <div className="rank-toolbar">
        <span className="rank-toolbar-label">Continent</span>
        <SortPill active={continent === "all"} onClick={() => setContinent("all")}>All</SortPill>
        {continents.map((c) => (
          <SortPill key={c.name} active={continent === c.name} onClick={() => setContinent(c.name)}>
            {c.name}
          </SortPill>
        ))}
        <span style={{ flex: 1, minWidth: 12 }} />
        <span className="rank-toolbar-label">Sort by</span>
        <SortPill active={sortBy === "total"} onClick={() => setSortBy("total")}>Total water</SortPill>
        <SortPill active={sortBy === "green"} onClick={() => setSortBy("green")}>Green share</SortPill>
        <SortPill active={sortBy === "blue"}  onClick={() => setSortBy("blue")}>Blue share</SortPill>
        <SortPill active={sortBy === "trend"} onClick={() => setSortBy("trend")}>Δ 2010→2020</SortPill>
        <span style={{ flex: 1, minWidth: 12 }} />
        <span className="rank-toolbar-label">Show</span>
        {[20, 50, 172].map((n) => (
          <SortPill key={n} active={limit === n} onClick={() => setLimit(n)}>
            {n === 172 ? "All" : `Top ${n}`}
          </SortPill>
        ))}
      </div>

      <div className="rank-wrap">
        <ol className="rank-list">
          {sorted.map((c, i) => (
            <li key={c.iso3}
                className="rank-row"
                style={{ animationDelay: `${i * 20}ms` }}
                onClick={() => onSelectCountry && onSelectCountry(c.iso3)}
                role="button"
                tabIndex={0}>
              <span className="rank-pos">{String(i + 1).padStart(2, "0")}</span>
              <span className="rank-flag">{c.iso3}</span>
              <span className="rank-name">{c.name}</span>
              <div className="rank-bar">
                <div className="rank-bar-track" style={{ width: `${((c.total_km3 || 0) / max) * 100}%` }}>
                  <span className="rank-seg rank-green" style={{ flexBasis: `${c.green_pct}%` }} />
                  <span className="rank-seg rank-blue"  style={{ flexBasis: `${c.blue_pct}%` }} />
                </div>
                <div className="rank-bar-meta">
                  <span className="rank-share">{c.green_pct}% green use · {c.blue_pct}% blue use</span>
                </div>
              </div>
              <span className="rank-total"><span className="num">{CGBW.fmt.km3p(c.total_km3)}</span><i>km³</i></span>
              <span className={`rank-trend ${trendClass(c.trend_2010_2020_pct)}`}>
                {c.trend_2010_2020_pct == null ? "—" : CGBW.fmt.pct(c.trend_2010_2020_pct)}
              </span>
            </li>
          ))}
        </ol>
      </div>
    </section>
  );
}

function trendClass(v) {
  if (v == null) return "low";
  if (v < 0) return "neg";
  if (v >= 15) return "high";
  if (v >= 5) return "mid";
  return "low";
}

function SortPill({ active, onClick, children }) {
  return <button className={active ? "pill on" : "pill"} onClick={onClick}>{children}</button>;
}

window.CountryRanking = CountryRanking;
