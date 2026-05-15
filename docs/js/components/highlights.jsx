// highlights.jsx — auto-rotating ribbon under the hero. Cycles through the
// key finding of each dashboard page; click jumps to that page.

function Highlights({ data, onNavigate }) {
  const slides = React.useMemo(() => buildSlides(data), [data]);
  const [idx, setIdx] = React.useState(0);
  const [paused, setPaused] = React.useState(false);

  React.useEffect(() => {
    if (paused) return;
    const t = setInterval(() => setIdx((i) => (i + 1) % slides.length), 5500);
    return () => clearInterval(t);
  }, [paused, slides.length]);

  const cur = slides[idx];

  return (
    <div className="highlights"
         onMouseEnter={() => setPaused(true)}
         onMouseLeave={() => setPaused(false)}
         aria-roledescription="carousel">
      <div className="highlights-card"
           key={idx}
           style={{ background: cur.bg, color: cur.fg }}
           onClick={() => onNavigate && onNavigate(cur.target)}
           role="button"
           tabIndex={0}>
        <div className="highlights-eyebrow">
          <span className="highlights-dot" style={{ background: cur.fg, opacity: 0.6 }} />
          <span>{cur.module}</span>
        </div>
        <div className="highlights-stat">
          <span className="highlights-num">{cur.stat}</span>
          <span className="highlights-unit">{cur.unit}</span>
        </div>
        <div className="highlights-line">{cur.line}</div>
        <div className="highlights-cta">
          Explore {cur.target.charAt(0).toUpperCase() + cur.target.slice(1)} →
        </div>
      </div>

      <div className="highlights-dots" role="tablist" aria-label="Highlight slides">
        {slides.map((s, i) => (
          <button key={i}
                  role="tab"
                  aria-selected={i === idx}
                  aria-label={s.module}
                  className={`highlights-tick ${i === idx ? "on" : ""}`}
                  onClick={() => { setIdx(i); setPaused(true); }} />
        ))}
        <button className="highlights-pause"
                onClick={() => setPaused((p) => !p)}
                aria-label={paused ? "Resume" : "Pause"}>
          {paused ? "▶" : "❚❚"}
        </button>
      </div>
    </div>
  );
}

function buildSlides(d) {
  const g20 = d.summary.global["2020"];
  const change = d.summary.change_2010_2020;
  const top3countries = d.countries.countries.slice(0, 3);
  const top3sum = top3countries.reduce((s, c) => s + (c.total_km3 || 0), 0);
  const top3pct = Math.round(top3sum / g20.total_km3 * 100);
  const cereals = d.groups.groups.find((g) => g.id === "cereals");
  const asia = d.continents && d.continents.continents.find((c) => c.name === "Asia");
  const greenPct = (g20.green_km3 / g20.total_km3 * 100);

  return [
    {
      target: "overview",
      module: "Overview · headline",
      stat: `+${Math.round(change.total_pct)}%`,
      unit: "in a decade",
      line: `Total crop water use rose from ${Math.round(d.summary.global["2010"].total_km3).toLocaleString()} km³ in 2010 to ${Math.round(g20.total_km3).toLocaleString()} km³ in 2020 — driven almost entirely by rainfall, not irrigation.`,
      bg: "linear-gradient(135deg, #1e4a7a 0%, #2d5a3d 100%)",
      fg: "#f3efe6",
    },
    {
      target: "atlas",
      module: "Atlas · where",
      stat: `${asia ? asia.share_pct : 52}%`,
      unit: "of global crop water is consumed in Asia",
      line: `Asia draws ${asia ? Math.round(asia.total_km3).toLocaleString() : "≈3,500"} km³/yr — more than Africa, Europe, and the Americas combined.`,
      bg: "linear-gradient(135deg, #2d5a3d 0%, #4a7d4d 100%)",
      fg: "#f3efe6",
    },
    {
      target: "trends",
      module: "Trends · how it changed",
      stat: `${greenPct.toFixed(1)}%`,
      unit: "is green water (rainfall)",
      line: `Green water consumption grew ${change.green_pct}% from 2010 to 2020; blue water (irrigation) grew only ${change.blue_pct}%. The green share has stayed stable for two decades.`,
      bg: "linear-gradient(135deg, #2d5a3d 0%, #7ea36b 100%)",
      fg: "#f3efe6",
    },
    {
      target: "crops",
      module: "Crops · what's grown",
      stat: cereals ? `${cereals.share_pct}%` : "47%",
      unit: "of crop water is for cereals",
      line: `Cereals (wheat, maize, rice…) dominate global crop water use. Rice alone accounts for over 1,000 km³/yr — about 15% of all crop water.`,
      bg: "linear-gradient(135deg, #c9883b 0%, #8b5e2b 100%)",
      fg: "#f3efe6",
    },
    {
      target: "ranking",
      module: "Countries · who",
      stat: `${top3pct}%`,
      unit: `from ${top3countries.map((c) => c.name).join(", ")}`,
      line: `Just three countries — ${top3countries.map((c) => c.name).join(", ")} — account for over a third of global crop water consumption, each drawing 400–1,000 km³/yr.`,
      bg: "linear-gradient(135deg, #1e4a7a 0%, #3a719f 100%)",
      fg: "#f3efe6",
    },
  ];
}

window.Highlights = Highlights;
