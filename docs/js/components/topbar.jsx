// topbar.jsx — sticky top bar with sharper CropGBWater logo + tab nav.

function BrandMark({ size = 44 }) {
  // Sharper: drop-shadow filter on circles, thicker plant strokes, glow on water-drop.
  return (
    <svg viewBox="0 0 100 100" width={size} height={size} aria-hidden="true" className="brand-mark-svg">
      <defs>
        {/* Soft drop shadow under the two circles */}
        <filter id="cgbw-shadow" x="-20%" y="-20%" width="140%" height="140%">
          <feGaussianBlur in="SourceAlpha" stdDeviation="0.9" />
          <feOffset dx="0" dy="0.9" result="off" />
          <feComponentTransfer result="dark">
            <feFuncA type="linear" slope="0.32" />
          </feComponentTransfer>
          <feMerge>
            <feMergeNode in="dark" />
            <feMergeNode in="SourceGraphic" />
          </feMerge>
        </filter>
        {/* Highlight glow for the plant + drop */}
        <filter id="cgbw-glow" x="-30%" y="-30%" width="160%" height="160%">
          <feGaussianBlur in="SourceGraphic" stdDeviation="0.7" result="b" />
          <feMerge>
            <feMergeNode in="b" />
            <feMergeNode in="SourceGraphic" />
          </feMerge>
        </filter>
        <clipPath id="cgbw-overlap-clip">
          <circle cx="36" cy="50" r="30" />
        </clipPath>
        {/* Radial sheen on each circle for a "shining" look */}
        <radialGradient id="cgbw-green-sheen" cx="35%" cy="35%" r="65%">
          <stop offset="0%"  stopColor="#5a8a6e" />
          <stop offset="60%" stopColor="#2d5a3d" />
          <stop offset="100%" stopColor="#234a32" />
        </radialGradient>
        <radialGradient id="cgbw-blue-sheen" cx="35%" cy="35%" r="65%">
          <stop offset="0%"  stopColor="#3f73a8" />
          <stop offset="60%" stopColor="#1e4a7a" />
          <stop offset="100%" stopColor="#163a62" />
        </radialGradient>
        <radialGradient id="cgbw-grey-sheen" cx="50%" cy="40%" r="65%">
          <stop offset="0%"  stopColor="#f0eee6" />
          <stop offset="100%" stopColor="#c8c6bd" />
        </radialGradient>
      </defs>

      <g filter="url(#cgbw-shadow)">
        <circle cx="36" cy="50" r="30" fill="url(#cgbw-green-sheen)" />
        <circle cx="64" cy="50" r="30" fill="url(#cgbw-blue-sheen)" />
        {/* Grey intersection — pure CSS would clip the green circle within the blue */}
        <g clipPath="url(#cgbw-overlap-clip)">
          <circle cx="64" cy="50" r="30" fill="url(#cgbw-grey-sheen)" />
        </g>
      </g>

      {/* Plant — sharper, with stronger strokes + glow */}
      <g transform="translate(50, 56)" filter="url(#cgbw-glow)">
        {/* Stem */}
        <path d="M0 12 Q0 0 0 -12" stroke="#0d3320" strokeWidth="2.2" strokeLinecap="round" fill="none" />
        {/* Left leaf */}
        <path d="M0 -1 Q-7 -3 -9 -10 Q-2 -10 0 -3 Z"
              fill="#7ea36b" stroke="#0d3320" strokeWidth="1.4" strokeLinejoin="round" />
        {/* Right leaf (higher) */}
        <path d="M0 -6 Q7 -8 9 -15 Q2 -15 0 -8 Z"
              fill="#7ea36b" stroke="#0d3320" strokeWidth="1.4" strokeLinejoin="round" />
        {/* Water drop on top with subtle highlight */}
        <ellipse cx="0" cy="-14.5" rx="1.9" ry="2.4" fill="#6fa1cf" stroke="#163a62" strokeWidth="0.9" />
        <ellipse cx="-0.5" cy="-15.2" rx="0.6" ry="0.9" fill="#e4eef8" opacity="0.85" />
      </g>
    </svg>
  );
}

function Logo({ stacked = true }) {
  return (
    <span className={`logo ${stacked ? "logo-stacked" : ""}`}>
      <span className="logo-crop">Crop</span>
      <span className="logo-letter logo-G">
        G
        {stacked && <span className="logo-letter-stack logo-G-stack">Green</span>}
      </span>
      <span className="logo-letter logo-B">
        B
        {stacked && <span className="logo-letter-stack logo-B-stack">Blue</span>}
      </span>
      <span className="logo-water">Water</span>
    </span>
  );
}

function TopBar({ meta, pages, active, onNavigate }) {
  const [navOpen, setNavOpen] = React.useState(false);
  return (
    <header className="topbar">
      <button className="brand" onClick={() => onNavigate("overview")} aria-label="CropGBWater home">
        <BrandMark size={44} />
        <Logo stacked={true} />
      </button>

      <button className="topnav-toggle"
              aria-label="Toggle navigation"
              aria-expanded={navOpen}
              onClick={() => setNavOpen((v) => !v)}>
        <span /><span /><span />
      </button>

      <nav className={`topnav ${navOpen ? "open" : ""}`} role="tablist">
        {pages.map((p) => (
          <button key={p.id}
                  role="tab"
                  aria-selected={active === p.id}
                  className={active === p.id ? "on" : ""}
                  onClick={() => { onNavigate(p.id); setNavOpen(false); }}>
            {p.label}
          </button>
        ))}
      </nav>

      <div className="topbar-actions">
        <span className="topbar-meta">v{meta?.version || "1.0"} · {meta?.release || "2026"}</span>
        <button className="topbar-cta" onClick={() => onNavigate("method")}>Cite &amp; download</button>
      </div>
    </header>
  );
}

window.TopBar = TopBar;
window.BrandMark = BrandMark;
window.Logo = Logo;
