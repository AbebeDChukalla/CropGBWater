// ai_insights.jsx — "AI Insights" page.
//
// Honest framing: this is a DEMO panel. No live model is called. Insight cards
// are pre-computed JSON authored against the same numbers as the rest of the
// atlas. The UX *feels* AI-assisted (generate button, animated reveal,
// filter chips) but every number traces back to data/summary.json.
//
// Loads data/ai_insights.json on mount. If the file is missing, the panel
// renders an honest empty state instead of crashing.

function AIInsights({ data, onNavigate }) {
  const [doc, setDoc] = React.useState(null);
  const [err, setErr] = React.useState(null);
  const [filter, setFilter] = React.useState("all");
  const [revealed, setRevealed] = React.useState(false);   // false = "press Generate"
  const [animKey, setAnimKey] = React.useState(0);         // bumps to retrigger reveal

  React.useEffect(() => {
    fetch("data/ai_insights.json")
      .then((r) => { if (!r.ok) throw new Error("ai_insights.json " + r.status); return r.json(); })
      .then((d) => setDoc(d))
      .catch((e) => setErr(e.message));
  }, []);

  const visibleCards = React.useMemo(() => {
    if (!doc) return [];
    if (filter === "all") return doc.cards;
    return doc.cards.filter((c) => c.category === filter);
  }, [doc, filter]);

  function handleAction(actionId) {
    const map = { generate: "all", hotspots: "hotspot", compare: "compare", explain: "trend" };
    setFilter(map[actionId] || "all");
    setRevealed(true);
    setAnimKey((k) => k + 1);
  }

  if (err) {
    return (
      <div className="ai-page">
        <AIHeader onNavigate={onNavigate} />
        <p className="ai-empty">Couldn&rsquo;t load <code>data/ai_insights.json</code>: {err}</p>
      </div>
    );
  }
  if (!doc) return <div className="ai-page"><div className="ai-skel" /></div>;

  return (
    <div className="ai-page">
      <AIHeader onNavigate={onNavigate} />

      <div className="ai-actions" role="group" aria-label="Insight actions">
        {doc.suggested_actions.map((a) => (
          <button key={a.id}
                  className="ai-action"
                  onClick={() => handleAction(a.id)}
                  title={a.hint}>
            <ActionIcon id={a.id} />
            <span>{a.label}</span>
          </button>
        ))}
      </div>

      <div className="ai-filterbar" role="tablist" aria-label="Filter insights">
        {doc.categories.map((c) => (
          <button key={c.id}
                  role="tab"
                  aria-selected={filter === c.id}
                  className={`ai-chip ${filter === c.id ? "on" : ""}`}
                  onClick={() => { setFilter(c.id); setRevealed(true); setAnimKey((k) => k + 1); }}>
            {c.label}
          </button>
        ))}
      </div>

      {!revealed ? (
        <div className="ai-prompt-cta">
          <div className="ai-prompt-orb" aria-hidden="true"><span/><span/><span/></div>
          <p>Press <b>Generate insights</b> to surface the {doc.cards.length} pre-computed cards for this release.</p>
          <p className="ai-prompt-foot">This is a demo. No live AI model is called — see the Roadmap for the live-AI plan.</p>
        </div>
      ) : (
        <div className="ai-grid" key={animKey}>
          {visibleCards.map((card, i) => (
            <InsightCard key={card.id}
                         card={card}
                         delay={i * 80}
                         onNavigate={onNavigate} />
          ))}
          {visibleCards.length === 0 && (
            <p className="ai-empty">No insights match this filter.</p>
          )}
        </div>
      )}

      <p className="ai-disclaimer">
        <span className="ai-badge">Demo</span>
        Insight text is authored against the same numbers as the rest of the atlas
        (see <button className="linklike" onClick={() => onNavigate && onNavigate("method")}>Method</button>).
        Live LLM integration is on the <button className="linklike" onClick={() => onNavigate && onNavigate("roadmap")}>Roadmap</button>.
      </p>
    </div>
  );
}

function AIHeader({ onNavigate }) {
  return (
    <header className="ai-header">
      <span className="ai-eyebrow"><span className="ai-eyebrow-dot" /> AI Insights · Demo</span>
      <h1 className="ai-title">Analytical highlights from the CropGBWater atlas</h1>
      <p className="ai-sub">
        Auto-curated cards that summarise the headline findings, hotspots and
        anomalies in the dataset. Use the actions below to filter by intent —
        or open the Copilot (bottom-right) to ask in your own words.
      </p>
    </header>
  );
}

function InsightCard({ card, delay = 0, onNavigate }) {
  const toneClass = `tone-${card.tone || "blue"}`;
  return (
    <article className={`ai-card ${toneClass} ${card.needs_attention ? "needs-attention" : ""}`}
             style={{ animationDelay: delay + "ms" }}>
      <header className="ai-card-head">
        <span className="ai-card-cat">{card.category}</span>
        {card.needs_attention && <span className="ai-card-flag">⚠ flagged</span>}
      </header>
      <div className="ai-card-stat">
        <span className="num">{card.stat}</span>
        <span className="ai-card-stat-label">{card.stat_label}</span>
      </div>
      <h3 className="ai-card-title">{card.title}</h3>
      <p className="ai-card-body">{card.body}</p>
      {card.evidence && card.evidence.length > 0 && (
        <ul className="ai-card-evidence">
          {card.evidence.map((e, i) => <li key={i}>{e}</li>)}
        </ul>
      )}
      {card.navigates_to && (
        <button className="ai-card-cta"
                onClick={() => onNavigate && onNavigate(card.navigates_to)}>
          Open {card.navigates_to} →
        </button>
      )}
    </article>
  );
}

function ActionIcon({ id }) {
  // Tiny inline SVGs — no external icon dependency.
  const common = { width: 16, height: 16, viewBox: "0 0 24 24", fill: "none",
                   stroke: "currentColor", strokeWidth: 1.6,
                   strokeLinecap: "round", strokeLinejoin: "round",
                   "aria-hidden": true };
  switch (id) {
    case "generate":
      return (
        <svg {...common}><path d="M12 3l1.8 4.6L18 9l-4.2 1.4L12 15l-1.8-4.6L6 9l4.2-1.4L12 3z"/><path d="M19 14l.8 2 2 .8-2 .8-.8 2-.8-2-2-.8 2-.8.8-2z"/></svg>
      );
    case "hotspots":
      return (
        <svg {...common}><path d="M12 21s-7-6-7-12a7 7 0 1114 0c0 6-7 12-7 12z"/><circle cx="12" cy="9" r="2.5"/></svg>
      );
    case "compare":
      return (
        <svg {...common}><path d="M4 6h7"/><path d="M4 18h7"/><path d="M13 6h7"/><path d="M13 18h7"/><path d="M7 6v12"/><path d="M17 6v12"/></svg>
      );
    case "explain":
      return (
        <svg {...common}><path d="M3 17l5-5 4 4 7-8"/><path d="M14 8h5v5"/></svg>
      );
    default:
      return <span aria-hidden="true">·</span>;
  }
}

window.AIInsights = AIInsights;
