// tweaks.jsx — floating theme + density controls (palette, dark mode, fonts).

const PALETTES = {
  earth: {
    paper:  "#f3efe6", paper2: "#ebe6d8",
    ink:    "#14140f", ink70:  "#4a4842", ink40: "#8a877d", ink20: "#b9b6ac",
    line:   "rgba(20,20,15,0.10)",
    green:  "#2d5a3d", greenSoft: "#7ea36b",
    blue:   "#1e4a7a", blueSoft:  "#6fa1cf",
    accent: "#c9883b",
  },
  oceanic: {
    paper:  "#f4f4f1", paper2: "#e9ebe6",
    ink:    "#0d1117", ink70: "#3b424c", ink40: "#7f8693", ink20: "#b8bcc4",
    line:   "rgba(13,17,23,0.10)",
    green:  "#1e6f5c", greenSoft: "#6fb29d",
    blue:   "#0e3a64", blueSoft:  "#5a96c8",
    accent: "#d68c45",
  },
  monochrome: {
    paper:  "#fafaf7", paper2: "#f1f0eb",
    ink:    "#1a1a1a", ink70: "#454545", ink40: "#8a8a8a", ink20: "#c0c0c0",
    line:   "rgba(0,0,0,0.10)",
    green:  "#4a4a4a", greenSoft: "#9b9b9b",
    blue:   "#1a1a1a", blueSoft:  "#6b6b6b",
    accent: "#1a1a1a",
  },
  bold: {
    paper:  "#f6f4ee", paper2: "#ede8dc",
    ink:    "#101820", ink70: "#3a4452", ink40: "#7e8896", ink20: "#b6bcc4",
    line:   "rgba(16,24,32,0.10)",
    green:  "#1f8a5b", greenSoft: "#74c39d",
    blue:   "#1d4ed8", blueSoft:  "#7aa4ef",
    accent: "#e8893b",
  },
};

const DARK = {
  paper: "#0f110d", paper2: "#15171f",
  ink: "#f1eee5", ink70: "#bcb9af", ink40: "#83807a", ink20: "#4d4b46",
  line: "rgba(241,238,229,0.10)",
};

const FONTS = {
  instrument: { display: "'Instrument Serif', serif",       body: "'Geist', system-ui, sans-serif",          mono: "'JetBrains Mono', ui-monospace, monospace" },
  newsreader: { display: "'Newsreader', 'Georgia', serif",  body: "'IBM Plex Sans', system-ui, sans-serif",  mono: "'IBM Plex Mono', ui-monospace, monospace" },
  schibsted:  { display: "'Schibsted Grotesk', sans-serif", body: "'Schibsted Grotesk', system-ui, sans-serif", mono: "'JetBrains Mono', ui-monospace, monospace" },
};

function TweaksPanel() {
  const [open, setOpen] = React.useState(false);
  const [palette, setPalette] = React.useState("earth");
  const [dark, setDark] = React.useState(false);
  const [fontPair, setFontPair] = React.useState("instrument");
  const [density, setDensity] = React.useState("spacious");

  React.useEffect(() => {
    const r = document.documentElement.style;
    const base = PALETTES[palette];
    const merged = dark ? { ...base, ...DARK } : base;
    Object.entries(merged).forEach(([k, v]) => r.setProperty(`--${k}`, v));
    const f = FONTS[fontPair];
    r.setProperty("--font-display", f.display);
    r.setProperty("--font-body",    f.body);
    r.setProperty("--font-mono",    f.mono);
    document.body.dataset.density = density;
    document.body.dataset.dark = dark ? "1" : "0";
  }, [palette, dark, fontPair, density]);

  return (
    <div className="tweaks">
      {open && (
        <div className="tweaks-panel on">
          <div className="tweaks-section">
            <span className="tweaks-section-label">Palette</span>
            <div className="ctrl-body">
              {Object.keys(PALETTES).map((p) => (
                <button key={p} className={`chip ${palette === p ? "on" : ""}`}
                        onClick={() => setPalette(p)}>{p}</button>
              ))}
            </div>
          </div>
          <div className="tweaks-section">
            <span className="tweaks-section-label">Dark mode</span>
            <div className="ctrl-body">
              <button className={`chip ${!dark ? "on" : ""}`} onClick={() => setDark(false)}>Light</button>
              <button className={`chip ${dark ? "on" : ""}`}  onClick={() => setDark(true)}>Dark</button>
            </div>
          </div>
          <div className="tweaks-section">
            <span className="tweaks-section-label">Type</span>
            <div className="ctrl-body">
              <button className={`chip ${fontPair === "instrument" ? "on" : ""}`} onClick={() => setFontPair("instrument")}>Serif</button>
              <button className={`chip ${fontPair === "newsreader" ? "on" : ""}`} onClick={() => setFontPair("newsreader")}>Editorial</button>
              <button className={`chip ${fontPair === "schibsted" ? "on" : ""}`}  onClick={() => setFontPair("schibsted")}>Sans</button>
            </div>
          </div>
          <div className="tweaks-section">
            <span className="tweaks-section-label">Density</span>
            <div className="ctrl-body">
              <button className={`chip ${density === "compact" ? "on" : ""}`} onClick={() => setDensity("compact")}>compact</button>
              <button className={`chip ${density === "spacious" ? "on" : ""}`} onClick={() => setDensity("spacious")}>spacious</button>
            </div>
          </div>
        </div>
      )}
      <button className="tweaks-toggle" onClick={() => setOpen((v) => !v)}>
        {open ? "Close tweaks" : "Tweaks"}
      </button>
    </div>
  );
}

window.TweaksPanel = TweaksPanel;
