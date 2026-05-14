// ranking.jsx — Top countries by water consumption, sortable by metric.

function CountryRanking({ data, onSelectCountry }) {
  const [sortBy, setSortBy] = React.useState("total");
  const [limit, setLimit] = React.useState(20);
  const countries = data.countries.countries;
  const max = Math.max(...countries.map((c) => c.total_km3 || 0));

  const sorted = React.useMemo(() => {
    const arr = [...countries];
    if (sortBy === "trend") arr.sort((a, b) => (b.trend_2010_2020_pct || 0) - (a.trend_2010_2020_pct || 0));
    else if (sortBy === "blue")  arr.sort((a, b) => (b.blue_pct || 0) - (a.blue_pct || 0));
    else if (sortBy === "green") arr.sort((a, b) => (b.green_pct || 0) - (a.green_pct || 0));
    else arr.sort((a, b) => (b.total_km3 || 0) - (a.total_km3 || 0));
    return arr.slice(0, limit);
  }, [countries, sortBy, limit]);

  return (
    <section id="ranking" className="rank-module" data-screen-label="03 Country Ranking">
      <div className="module-head">
        <div className="module-eyebrow">
          <span className="eyebrow-mark" />
          <span>Module 03 / Country ranking</span>
        </div>
        <h2 className="module-title">
          Two countries — India and China — each draw nearly a thousand km³ a year.
        </h2>
        <p className="module-sub">
          Ranked by 2020 total consumption. Bars show the green / blue composition;
          trend is the change since 2010. Click any row for full breakdown.
        </p>
      </div>

      <div className="rank-toolbar">
        <span className="rank-toolbar-label">Sort by</span>
        <SortPill active={sortBy === "total"} onClick={() => setSortBy("total")}>Total water</SortPill>
        <SortPill active={sortBy === "green"} onClick={() => setSortBy("green")}>Green share</SortPill>
        <SortPill active={sortBy === "blue"}  onClick={() => setSortBy("blue")}>Blue share</SortPill>
        <SortPill active={sortBy === "trend"} onClick={() => setSortBy("trend")}>Δ 2010→2020</SortPill>
        <span style={{ flex: 1 }} />
        <span className="rank-toolbar-label">Show</span>
        {[20, 50, 172].map((n) => (
          <SortPill key={n} active={limit === n} onClick={() => setLimit(n)}>
            {n === 172 ? "All" : `Top ${n}`}
          </SortPill>
        ))}
      </div>

      <ol className="rank-list">
        {sorted.map((c, i) => (
          <li key={c.iso3}
              className="rank-row"
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
                <span className="rank-share">{c.green_pct}% green · {c.blue_pct}% blue</span>
              </div>
            </div>
            <span className="rank-total"><span className="num">{CGBW.fmt.km3p(c.total_km3)}</span><i>km³</i></span>
            <span className={`rank-trend ${trendClass(c.trend_2010_2020_pct)}`}>
              {c.trend_2010_2020_pct == null ? "—" : CGBW.fmt.pct(c.trend_2010_2020_pct)}
            </span>
          </li>
        ))}
      </ol>
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
