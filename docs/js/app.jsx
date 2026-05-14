// app.jsx — root composition. Loads all data, mounts every module.

function App() {
  const [data, setData] = React.useState(null);
  const [err, setErr] = React.useState(null);
  const [activeISO, setActiveISO] = React.useState(null);

  React.useEffect(() => {
    CGBW.loadAll()
      .then((d) => {
        setData(d);
        document.body.dataset.loading = "0";
      })
      .catch((e) => {
        console.error(e);
        setErr(e.message);
        document.body.dataset.loading = "0";
      });
  }, []);

  if (err) {
    return (
      <div style={{
        padding: 80, maxWidth: 720, margin: "0 auto",
        fontFamily: "var(--font-body)", color: "var(--ink)"
      }}>
        <h1 style={{ fontFamily: "var(--font-display)", fontSize: 36 }}>Couldn&rsquo;t load data</h1>
        <p style={{ color: "var(--ink70)" }}>
          The static JSON in <code>webapp/data/</code> failed to load.
          Make sure you&rsquo;re serving the page from a local web server
          (e.g. <code>python -m http.server</code>) — opening the file directly
          via <code>file://</code> blocks <code>fetch()</code>.
        </p>
        <pre style={{ fontFamily: "var(--font-mono)", color: "var(--accent)" }}>{err}</pre>
      </div>
    );
  }
  if (!data) return null;

  return (
    <div className="page">
      <TopBar meta={data.meta} />
      <Hero data={data} />
      <Atlas data={data} onSelectCountry={setActiveISO} />
      <CropExplorer data={data} />
      <Trends data={data} />
      <CountryRanking data={data} onSelectCountry={setActiveISO} />
      <MethodFooter data={data} />

      {activeISO && <CountrySheet iso={activeISO} onClose={() => setActiveISO(null)} />}

      <TweaksPanel />
    </div>
  );
}

const root = ReactDOM.createRoot(document.getElementById("root"));
root.render(<App />);
