// roadmap.jsx — Digital Twin / future roadmap page.
//
// Narrative + SVG architecture diagram + four "concept" scenario mockup
// cards. Every forward-looking element is explicitly labelled as concept
// or roadmap so visitors don't mistake mock-ups for live tools.

function Roadmap({ data, onNavigate }) {
  return (
    <div className="rm-page">
      <header className="rm-header">
        <span className="rm-eyebrow"><span className="rm-eyebrow-dot" /> Roadmap · concept</span>
        <h1 className="rm-title">Toward an AI-assisted Digital Twin for agricultural water use</h1>
        <p className="rm-lede">
          CropGBWater is evolving from a static atlas into an intelligent layer
          for agricultural water analytics — integrating geospatial data, crop
          water consumption, water productivity, forecasting and scenario
          exploration to support sustainable decision-making.
        </p>
        <p className="rm-sub">
          What you see today is the <b>retrospective atlas</b> (2000–2020) plus
          a <b>demo AI layer</b>. The components below are concept mock-ups for
          the next phase — they are not live tools.
        </p>
      </header>

      <section className="rm-section">
        <h2 className="rm-h2">Architecture</h2>
        <p className="rm-p">
          The diagram shows how data, models, AI and users connect.
          <span className="rm-legend"><span className="rm-leg rm-leg-live">■ live today</span> <span className="rm-leg rm-leg-demo">■ demo</span> <span className="rm-leg rm-leg-concept">■ concept</span></span>
        </p>
        <ArchitectureDiagram />
      </section>

      <section className="rm-section">
        <h2 className="rm-h2">Scenario explorer · concept</h2>
        <p className="rm-p">
          Four scenario modules are planned for the next phase. Each one will
          let users perturb the baseline atlas with a small set of physical
          parameters and re-render the consequence on the map.
        </p>
        <div className="rm-scenarios">
          <ScenarioCard
            tone="amber"
            tag="Climate"
            title="Climate change impact"
            body="Apply a +1.5°C / +2°C / +3°C warming pattern to the gridded crop calendar and re-compute green/blue water demand per crop. Surface the basins where blue-water demand grows fastest under each scenario."
            inputs={["Warming level: 1.5 / 2 / 3°C", "Precipitation change: ± %", "CO₂ fertilisation: on / off"]}
          />
          <ScenarioCard
            tone="blue"
            tag="Expansion"
            title="Irrigation expansion"
            body="Slider for additional irrigated area (0–100 Mha) distributed across rainfed cropland with highest water-stress reduction potential. Shows the marginal blue-water cost per additional ton of cereal production."
            inputs={["Additional irrigated area (Mha)", "Allocation rule (yield gap / poverty)", "Efficiency assumption"]}
          />
          <ScenarioCard
            tone="green"
            tag="Savings"
            title="Water-saving interventions"
            body="Toggle interventions — drip conversion, mulching, deficit irrigation, crop switching — and project the basin-level blue-water savings. Pairs with the productivity dashboard to show the yield trade-off."
            inputs={["Drip adoption (% irrigated area)", "Mulching adoption (% rainfed)", "Crop-switch rules"]}
          />
          <ScenarioCard
            tone="ink"
            tag="Forecast"
            title="2030 / 2050 demand projection"
            body="Forward extrapolation of the 2000–2020 trajectory under SSP1–5 socio-economic scenarios. Outputs: total km³, green/blue share, top-20 country trajectories, and the basins crossing planetary-boundary thresholds."
            inputs={["SSP pathway", "Year (2030 / 2050)", "Diet shift (optional)"]}
          />
        </div>
      </section>

      <section className="rm-section">
        <h2 className="rm-h2">Why no live LLM yet?</h2>
        <p className="rm-p rm-p-narrow">
          The site is hosted on GitHub Pages — a static host that cannot keep
          secrets. A live LLM would either expose an API key in the browser or
          require a paid serverless proxy. For an open, citable, free-to-run
          research atlas we chose pre-computed insights and a keyword-matched
          assistant instead. A serverless proxy with rate limits is the
          natural next step, but only once the cost model is clear.
        </p>
        <div className="rm-callout">
          <b>You can still try the demo AI layer now:</b>
          <div className="rm-cta-row">
            <button className="rm-cta" onClick={() => onNavigate && onNavigate("ai")}>Open AI Insights</button>
            <button className="rm-cta rm-cta-ghost"
                    onClick={() => window.dispatchEvent(new CustomEvent("cgbw:open-copilot"))}>
              Try the Copilot
            </button>
          </div>
        </div>
      </section>
    </div>
  );
}

// ─── Architecture diagram (SVG, no external deps) ─────────────────
function ArchitectureDiagram() {
  // Nodes are positioned in a 3-row flow: data sources → engines → users.
  return (
    <div className="rm-arch-wrap">
      <svg viewBox="0 0 880 460"
           role="img"
           aria-label="Architecture: data sources flow into modelling and AI engines, then surface as dashboards and a copilot for users"
           className="rm-arch">
        <defs>
          <linearGradient id="arch-blue" x1="0" y1="0" x2="0" y2="1">
            <stop offset="0%" stopColor="#3a719f"/>
            <stop offset="100%" stopColor="#1e4a7a"/>
          </linearGradient>
          <linearGradient id="arch-green" x1="0" y1="0" x2="0" y2="1">
            <stop offset="0%" stopColor="#7ea36b"/>
            <stop offset="100%" stopColor="#2d5a3d"/>
          </linearGradient>
          <linearGradient id="arch-amber" x1="0" y1="0" x2="0" y2="1">
            <stop offset="0%" stopColor="#d8a05a"/>
            <stop offset="100%" stopColor="#a96940"/>
          </linearGradient>
          <marker id="arch-arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto">
            <path d="M0 0 L10 5 L0 10 z" fill="#4a4842"/>
          </marker>
          <filter id="arch-shadow" x="-20%" y="-20%" width="140%" height="140%">
            <feGaussianBlur in="SourceAlpha" stdDeviation="1.4"/>
            <feOffset dx="0" dy="1.6" result="off"/>
            <feComponentTransfer><feFuncA type="linear" slope="0.22"/></feComponentTransfer>
            <feMerge><feMergeNode/><feMergeNode in="SourceGraphic"/></feMerge>
          </filter>
        </defs>

        {/* Row labels */}
        <text x="20" y="40"  className="rm-arch-rowlbl">Data</text>
        <text x="20" y="220" className="rm-arch-rowlbl">Engines</text>
        <text x="20" y="400" className="rm-arch-rowlbl">Users</text>

        {/* Subtle row guides */}
        <line x1="120" y1="60"  x2="860" y2="60"  className="rm-arch-guide"/>
        <line x1="120" y1="240" x2="860" y2="240" className="rm-arch-guide"/>
        <line x1="120" y1="420" x2="860" y2="420" className="rm-arch-guide"/>

        {/* DATA ROW (top) — three sources */}
        <ArchNode x={150} y={20}  w={170} h={80} fill="url(#arch-green)" status="live"
                  title="Remote sensing"     sub="MODIS · SPAM · GAEZ"/>
        <ArchNode x={355} y={20}  w={170} h={80} fill="url(#arch-blue)" status="live"
                  title="CropGBWater data"   sub="46 crops · 2000–2020"/>
        <ArchNode x={560} y={20}  w={170} h={80} fill="url(#arch-amber)" status="concept"
                  title="Real-time streams"  sub="weather · irrigation IoT"/>

        {/* ENGINES ROW (middle) — three engines */}
        <ArchNode x={150} y={200} w={170} h={80} fill="url(#arch-green)" status="live"
                  title="Crop water model"   sub="PCR-GLOBWB · GAEZ"/>
        <ArchNode x={355} y={200} w={170} h={80} fill="url(#arch-blue)" status="demo"
                  title="AI Agent"           sub="pre-computed insights"/>
        <ArchNode x={560} y={200} w={170} h={80} fill="url(#arch-amber)" status="concept"
                  title="Digital Twin core"  sub="scenario · forecast"/>

        {/* USERS ROW (bottom) — three surfaces */}
        <ArchNode x={150} y={380} w={170} h={60} fill="#f3efe6" status="live"
                  title="Dashboards"         sub="maps · charts · rankings" dark/>
        <ArchNode x={355} y={380} w={170} h={60} fill="#f3efe6" status="demo"
                  title="Copilot"            sub="natural language" dark/>
        <ArchNode x={560} y={380} w={170} h={60} fill="#f3efe6" status="concept"
                  title="Scenario explorer"  sub="what-if sliders" dark/>

        {/* Edges: data → engines */}
        <ArchEdge x1={235} y1={100} x2={235} y2={200}/>
        <ArchEdge x1={440} y1={100} x2={440} y2={200}/>
        <ArchEdge x1={645} y1={100} x2={645} y2={200}/>
        {/* cross-feeds in engines row */}
        <ArchEdge x1={320} y1={240} x2={355} y2={240} arrow={false}/>
        <ArchEdge x1={525} y1={240} x2={560} y2={240} arrow={false}/>
        {/* Engines → users */}
        <ArchEdge x1={235} y1={280} x2={235} y2={380}/>
        <ArchEdge x1={440} y1={280} x2={440} y2={380}/>
        <ArchEdge x1={645} y1={280} x2={645} y2={380}/>

        {/* End-users node on the right */}
        <g transform="translate(770, 200)">
          <circle r="42" fill="#14140f"/>
          <text textAnchor="middle" y="-2" className="rm-arch-userlbl">Users</text>
          <text textAnchor="middle" y="16" className="rm-arch-usersub">researchers · policy · farmers</text>
        </g>
        <ArchEdge x1={730} y1={240} x2={770} y2={240} arrow={false}/>
      </svg>
      <p className="rm-arch-cap">Data → modelling → AI analysis → user-facing surfaces. The dashed connectors are bidirectional: users can ask the Copilot, which queries the AI Agent, which reads the same numbers the dashboards render.</p>
    </div>
  );
}

function ArchNode({ x, y, w, h, fill, title, sub, status, dark }) {
  const statusColor = { live: "#2d5a3d", demo: "#1e4a7a", concept: "#a96940" }[status] || "#8a877d";
  return (
    <g filter="url(#arch-shadow)" transform={`translate(${x}, ${y})`}>
      <rect width={w} height={h} rx="10" fill={fill} stroke="rgba(0,0,0,0.08)"/>
      <text x={w / 2} y={h / 2 - 4} textAnchor="middle"
            className={dark ? "rm-arch-title rm-arch-title-dark" : "rm-arch-title"}>{title}</text>
      <text x={w / 2} y={h / 2 + 16} textAnchor="middle"
            className={dark ? "rm-arch-sub rm-arch-sub-dark" : "rm-arch-sub"}>{sub}</text>
      <rect x="8" y={h - 18} width="46" height="13" rx="6.5" fill="rgba(255,255,255,0.9)" stroke={statusColor}/>
      <text x="31" y={h - 8} textAnchor="middle" className="rm-arch-status" fill={statusColor}>{status}</text>
    </g>
  );
}

function ArchEdge({ x1, y1, x2, y2, arrow = true }) {
  return (
    <line x1={x1} y1={y1} x2={x2} y2={y2}
          className="rm-arch-edge"
          markerEnd={arrow ? "url(#arch-arrow)" : undefined}/>
  );
}

function ScenarioCard({ tag, title, body, inputs, tone }) {
  return (
    <article className={`rm-scn tone-${tone}`}>
      <header className="rm-scn-head">
        <span className="rm-scn-tag">{tag}</span>
        <span className="rm-scn-badge">concept</span>
      </header>
      <h3 className="rm-scn-title">{title}</h3>
      <p className="rm-scn-body">{body}</p>
      <div className="rm-scn-inputs">
        <span className="rm-scn-inputs-lbl">Planned inputs</span>
        <ul>{inputs.map((i, k) => <li key={k}>{i}</li>)}</ul>
      </div>
      <div className="rm-scn-mock">
        <div className="rm-scn-mock-bar" style={{ width: "32%" }} />
        <div className="rm-scn-mock-bar" style={{ width: "58%" }} />
        <div className="rm-scn-mock-bar" style={{ width: "74%" }} />
        <span className="rm-scn-mock-cap">illustrative</span>
      </div>
    </article>
  );
}

window.Roadmap = Roadmap;
