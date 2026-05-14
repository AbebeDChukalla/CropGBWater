// topbar.jsx — sticky top bar with scrollspy nav.

const NAV_SECTIONS = [
  { id: "overview", label: "Overview" },
  { id: "atlas",    label: "Atlas" },
  { id: "crops",    label: "Crops" },
  { id: "trends",   label: "Trends" },
  { id: "ranking",  label: "Countries" },
  { id: "method",   label: "Method" },
];

function TopBar({ meta }) {
  const [active, setActive] = React.useState("overview");

  React.useEffect(() => {
    const targets = NAV_SECTIONS
      .map((s) => document.getElementById(s.id))
      .filter(Boolean);
    if (targets.length === 0) return;
    const io = new IntersectionObserver(
      (entries) => {
        const visible = entries.filter((e) => e.isIntersecting);
        if (visible.length === 0) return;
        const top = visible.sort((a, b) => a.boundingClientRect.top - b.boundingClientRect.top)[0];
        setActive(top.target.id);
      },
      { rootMargin: "-30% 0px -60% 0px", threshold: 0 }
    );
    targets.forEach((t) => io.observe(t));
    return () => io.disconnect();
  }, []);

  return (
    <header className="topbar">
      <a className="brand" href="#overview">
        <span className="brand-mark">
          <svg viewBox="0 0 20 20" width="20" height="20" aria-hidden="true">
            <circle cx="7" cy="10" r="5.5" fill="var(--green)" />
            <circle cx="13" cy="10" r="5.5" fill="var(--blue)" style={{ mixBlendMode: "multiply", opacity: .92 }} />
          </svg>
        </span>
        <span className="brand-word">
          Crops<span className="brand-g">Green</span><span className="brand-b">Blue</span>Water
        </span>
      </a>
      <nav className="topnav">
        {NAV_SECTIONS.map((s) => (
          <a key={s.id} href={`#${s.id}`} className={active === s.id ? "on" : ""}>
            {s.label}
          </a>
        ))}
      </nav>
      <div className="topbar-actions">
        <span className="topbar-meta">v{meta?.version || "1.0"} · {meta?.release || "2026"} release</span>
        <a className="topbar-cta" href="#method">Cite & download</a>
      </div>
    </header>
  );
}

window.TopBar = TopBar;
