// app.jsx — tab-paged dashboard with one module visible at a time.
// Modules stay mounted (display:none) so React + Leaflet keep their state.

const PAGES = [
  { id: "overview", label: "Overview" },
  { id: "atlas",    label: "Atlas" },
  { id: "trends",   label: "Trends" },
  { id: "crops",    label: "Crops" },
  { id: "ranking",  label: "Countries" },
  { id: "ai",       label: "AI Insights" },
  { id: "roadmap",  label: "Roadmap" },
  { id: "method",   label: "Method" },
  { id: "contact",  label: "Contact" },
];

function App() {
  const [data, setData] = React.useState(null);
  const [err, setErr] = React.useState(null);
  const [activeISO, setActiveISO] = React.useState(null);

  // URL hash drives the active page so deep-links and back/forward work.
  const [page, setPage] = React.useState(() => {
    const h = (location.hash || "#overview").replace("#", "");
    return PAGES.some((p) => p.id === h) ? h : "overview";
  });

  React.useEffect(() => {
    CGBW.loadAll()
      .then((d) => { setData(d); document.body.dataset.loading = "0"; })
      .catch((e) => { console.error(e); setErr(e.message); document.body.dataset.loading = "0"; });
  }, []);

  // Sync hash → state and state → hash
  React.useEffect(() => {
    const onHash = () => {
      const h = (location.hash || "#overview").replace("#", "");
      if (PAGES.some((p) => p.id === h)) setPage(h);
    };
    window.addEventListener("hashchange", onHash);
    return () => window.removeEventListener("hashchange", onHash);
  }, []);

  React.useEffect(() => {
    if (location.hash.replace("#", "") !== page) {
      history.replaceState(null, "", `#${page}`);
    }
    // Fire a custom event so each module can react (Atlas uses it to invalidate size)
    window.dispatchEvent(new CustomEvent("cgbw:page", { detail: page }));
    // Always reset scroll to top on page change
    window.scrollTo({ top: 0, behavior: "instant" });
  }, [page]);

  // Keyboard arrows ← → to navigate between pages
  React.useEffect(() => {
    const onKey = (e) => {
      if (e.target.tagName === "INPUT" || e.target.tagName === "TEXTAREA") return;
      if (activeISO) return;
      const idx = PAGES.findIndex((p) => p.id === page);
      if (e.key === "ArrowRight" && idx < PAGES.length - 1) setPage(PAGES[idx + 1].id);
      if (e.key === "ArrowLeft"  && idx > 0)                setPage(PAGES[idx - 1].id);
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [page, activeISO]);

  if (err) {
    return (
      <div className="boot-error">
        <h1>Couldn&rsquo;t load data</h1>
        <p>Serve the page through a static server, not <code>file://</code>.</p>
        <pre>{err}</pre>
      </div>
    );
  }
  if (!data) return null;

  return (
    <div className="page">
      <TopBar meta={data.meta} pages={PAGES} active={page} onNavigate={setPage} />

      <main className="stage">
        <Module id="overview" active={page === "overview"}>
          <Hero data={data} />
          <Highlights data={data} onNavigate={setPage} />
          <OverviewExtras data={data} onNavigate={setPage} />
        </Module>
        <Module id="atlas" active={page === "atlas"}>
          <Atlas data={data} onSelectCountry={setActiveISO} />
        </Module>
        <Module id="trends" active={page === "trends"}>
          <Trends data={data} />
          {/* BasinStress temporarily removed to isolate render error */}
        </Module>
        <Module id="crops" active={page === "crops"}>
          <CropExplorer data={data} />
        </Module>
        <Module id="ranking" active={page === "ranking"}>
          <CountryRanking data={data} onSelectCountry={setActiveISO} />
        </Module>
        <Module id="ai" active={page === "ai"}>
          <AIInsights data={data} onNavigate={setPage} />
        </Module>
        <Module id="roadmap" active={page === "roadmap"}>
          <Roadmap data={data} onNavigate={setPage} />
        </Module>
        <Module id="method" active={page === "method"}>
          <MethodFooter data={data} />
        </Module>
        <Module id="contact" active={page === "contact"}>
          <FeedbackPage data={data} />
        </Module>
      </main>

      <PageNav pages={PAGES} active={page} onNavigate={setPage} />

      {activeISO && <CountrySheet iso={activeISO} onClose={() => setActiveISO(null)} />}

      <TweaksPanel />
      <Copilot onNavigate={setPage} />
    </div>
  );
}

// One mounted shell per module. We use display:none for inactive ones so
// Leaflet keeps its state across switches; the .active class fires the page-in
// animation each time it's re-applied via the keyed wrapper inside.
function Module({ id, active, children }) {
  return (
    <section className={`module ${active ? "active" : ""}`} data-page={id} aria-hidden={!active}>
      {active && <div className="module-animate" key={id}>{children}</div>}
    </section>
  );
}

// Bottom navigation: previous / next + step indicator. Visible on every page.
function PageNav({ pages, active, onNavigate }) {
  const idx = pages.findIndex((p) => p.id === active);
  const prev = idx > 0 ? pages[idx - 1] : null;
  const next = idx < pages.length - 1 ? pages[idx + 1] : null;
  return (
    <nav className="pagenav" aria-label="Page navigation">
      <button className="pagenav-btn"
              disabled={!prev}
              onClick={() => prev && onNavigate(prev.id)}>
        <span className="pagenav-dir">← Previous</span>
        <span className="pagenav-label">{prev ? prev.label : "—"}</span>
      </button>

      <div className="pagenav-steps" role="tablist">
        {pages.map((p, i) => (
          <button key={p.id}
                  role="tab"
                  aria-selected={p.id === active}
                  className={`pagenav-step ${p.id === active ? "on" : ""}`}
                  title={`${i + 1}. ${p.label}`}
                  onClick={() => onNavigate(p.id)}>
            <span className="pagenav-step-num">{String(i + 1).padStart(2, "0")}</span>
            <span className="pagenav-step-label">{p.label}</span>
          </button>
        ))}
      </div>

      <button className="pagenav-btn pagenav-btn-next"
              disabled={!next}
              onClick={() => next && onNavigate(next.id)}>
        <span className="pagenav-dir">Next →</span>
        <span className="pagenav-label">{next ? next.label : "—"}</span>
      </button>
    </nav>
  );
}

const root = ReactDOM.createRoot(document.getElementById("root"));
root.render(<App />);
