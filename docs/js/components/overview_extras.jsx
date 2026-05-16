// overview_extras.jsx — Overview-page enhancements that live BELOW the Hero:
//   1. CTA strip (Explore Atlas / View Trends / Explore Crops / Methodology)
//   2. "Why It Matters" educational cards (green vs blue water concepts)
//   3. Animated conceptual flow diagram (Climate -> Soil -> ET -> Green/Blue -> Crop -> Country)
//   4. Use Cases cards (irrigation planning, water accounting, food security, ...)

function OverviewExtras({ data, onNavigate }) {
  return (
    <>
      <CtaStrip onNavigate={onNavigate} />
      <WhyItMatters />
      <ConceptualFlow />
      <UseCases onNavigate={onNavigate} />
    </>
  );
}

// ─── CTA strip ──────────────────────────────────────────────────────────
function CtaStrip({ onNavigate }) {
  const items = [
    { id: "atlas",   label: "Explore Atlas",   sub: "Interactive country map", icon: "🌐" },
    { id: "trends",  label: "View Trends",     sub: "2000 → 2010 → 2020",      icon: "📈" },
    { id: "crops",   label: "Explore Crops",   sub: "46 crops, 10 groups",     icon: "🌾" },
    { id: "method",  label: "Methodology",     sub: "Citation & sources",      icon: "📚" },
  ];
  return (
    <nav className="cta-strip" aria-label="Jump to a section">
      {items.map((it, i) => (
        <button key={it.id}
                className="cta-card"
                style={{ animationDelay: `${i * 80}ms` }}
                onClick={() => onNavigate(it.id)}>
          <span className="cta-icon" aria-hidden="true">{it.icon}</span>
          <span className="cta-body">
            <span className="cta-label">{it.label}</span>
            <span className="cta-sub">{it.sub}</span>
          </span>
          <span className="cta-arrow" aria-hidden="true">→</span>
        </button>
      ))}
    </nav>
  );
}

// ─── Why It Matters ────────────────────────────────────────────────────
function WhyItMatters() {
  const cards = [
    {
      icon: "🌧️",
      title: "Green water = rainfall",
      body: "The portion of soil moisture coming from rainfall and stored in the root zone. It feeds rainfed crops at no extraction cost from rivers or aquifers."
    },
    {
      icon: "💧",
      title: "Blue water = irrigation",
      body: "Surface water and groundwater withdrawn to irrigate crops. It's a limited, contested resource — measurable abstraction with hydrological consequences."
    },
    {
      icon: "🍚",
      title: "Food production hinges on both",
      body: "84% of global crop water is green — but the blue 16% disproportionately supports calorie-dense staples like rice and wheat in arid regions."
    },
    {
      icon: "⚠️",
      title: "Why it's getting harder",
      body: "Population growth, dietary shifts, climate variability, and aquifer depletion are pushing irrigation demand higher just as variable rainfall stresses rainfed systems."
    },
  ];
  return (
    <section className="why-it-matters" aria-labelledby="why-h">
      <div className="module-eyebrow"><span className="eyebrow-mark" /><span>Why it matters</span></div>
      <h2 id="why-h" className="module-title">
        The two flavours of water that grow your food.
      </h2>

      <div className="why-grid">
        {cards.map((c, i) => (
          <article key={i} className="why-card" style={{ animationDelay: `${i * 80}ms` }}>
            <div className="why-icon" aria-hidden="true">{c.icon}</div>
            <h3 className="why-title">{c.title}</h3>
            <p className="why-body">{c.body}</p>
          </article>
        ))}
      </div>
    </section>
  );
}

// ─── Animated conceptual flow diagram ──────────────────────────────────
function ConceptualFlow() {
  // 6 stages — each animates in with a stagger; arrows pulse between them.
  const stages = [
    { id: "climate", icon: "☀️", label: "Climate",          sub: "Precipitation · ET₀ · temperature" },
    { id: "soil",    icon: "🟤", label: "Soil moisture",     sub: "Root-zone storage · CN · AWC" },
    { id: "eta",     icon: "🌫️", label: "ET actual",        sub: "Crop water requirement" },
    { id: "split",   icon: "🔀", label: "Green / Blue split", sub: "Rainfed vs irrigation share" },
    { id: "crop",   icon: "🌾", label: "Crop water use",     sub: "Per crop · m³ per ton" },
    { id: "agg",    icon: "🗺️", label: "Country & global",  sub: "Aggregated totals · trends" },
  ];
  return (
    <section className="flow-section" aria-labelledby="flow-h">
      <div className="module-eyebrow"><span className="eyebrow-mark" /><span>From model to atlas</span></div>
      <h2 id="flow-h" className="module-title">
        How a raindrop becomes a country-level number.
      </h2>
      <p className="module-sub">
        The chain that turns climate inputs into the green and blue water figures
        you see across the dashboard. Hover each stage for what it represents.
      </p>

      <div className="flow-diagram" role="img" aria-label="Conceptual flow diagram">
        {stages.map((s, i) => (
          <React.Fragment key={s.id}>
            <div className="flow-stage" style={{ animationDelay: `${i * 110}ms` }} tabIndex={0}>
              <div className="flow-icon" aria-hidden="true">{s.icon}</div>
              <div className="flow-stage-body">
                <div className="flow-stage-label">{s.label}</div>
                <div className="flow-stage-sub">{s.sub}</div>
              </div>
            </div>
            {i < stages.length - 1 && (
              <div className="flow-arrow" aria-hidden="true" style={{ animationDelay: `${i * 110 + 55}ms` }}>
                <svg viewBox="0 0 24 12" width="28" height="14">
                  <path d="M0 6 H 22" stroke="currentColor" strokeWidth="1.4" />
                  <path d="M18 2 L 22 6 L 18 10" stroke="currentColor" strokeWidth="1.4" fill="none" strokeLinecap="round" strokeLinejoin="round" />
                </svg>
              </div>
            )}
          </React.Fragment>
        ))}
      </div>
    </section>
  );
}

// ─── Use Cases ─────────────────────────────────────────────────────────
function UseCases({ onNavigate }) {
  const cases = [
    { icon: "🚿", title: "Irrigation planning",
      body: "Identify where blue water demand is concentrated and where irrigation expansion or efficiency gains have the biggest leverage." },
    { icon: "📊", title: "National water accounting",
      body: "Country-level green/blue/total volumes consistent with global modelling outputs, broken down per crop and technology." },
    { icon: "🍞", title: "Food security",
      body: "Trace which staple crops depend most on irrigation and where rainfall variability poses the largest production risk." },
    { icon: "🌡️", title: "Climate adaptation",
      body: "Trends 2000–2020 reveal where green water has expanded vs where blue water is offsetting precipitation deficits." },
    { icon: "🎓", title: "Teaching & education",
      body: "A scientifically credible, fully sourced atlas suitable for university, high school, and outreach use, with downloadable data." },
    { icon: "🤝", title: "Policy & sustainable management",
      body: "Underpinning evidence for SDG 6.4 (water-use efficiency), national water plans, and basin-scale management dialogues." },
  ];
  return (
    <section className="usecases" aria-labelledby="uc-h">
      <div className="module-eyebrow"><span className="eyebrow-mark" /><span>Use cases</span></div>
      <h2 id="uc-h" className="module-title">
        Built for the people deciding where the next drop goes.
      </h2>

      <div className="usecases-grid">
        {cases.map((c, i) => (
          <article key={i} className="usecase-card" style={{ animationDelay: `${i * 60}ms` }}>
            <div className="usecase-icon" aria-hidden="true">{c.icon}</div>
            <h3 className="usecase-title">{c.title}</h3>
            <p className="usecase-body">{c.body}</p>
          </article>
        ))}
      </div>
    </section>
  );
}

window.OverviewExtras = OverviewExtras;
