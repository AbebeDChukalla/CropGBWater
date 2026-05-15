// topbar.jsx — sticky top bar with CropGBWater logo:
//
//   Brand mark: two overlapping glossy globes (green left, blue right) with
//   world-map grid lines, drop-shape at the centre containing a shiny plant.
//
//   Wordmark:   Crop[G][B]Water — G in green, B in blue, rest in charcoal.
//               Small horizontal labels "Green" / "Blue" above the G / B,
//               in their respective colours.

function BrandMark({ size = 56 }) {
  // viewBox is wider than tall because the two globes sit side by side
  const VBW = 120, VBH = 80;
  return (
    <svg viewBox={`0 0 ${VBW} ${VBH}`}
         width={Math.round(size * VBW / VBH)}
         height={size}
         aria-hidden="true"
         className="brand-mark-svg">
      <defs>
        {/* Globe gradients — glossy 3D feel with brighter highlight upper-left */}
        <radialGradient id="cgbw-globe-green" cx="32%" cy="28%" r="78%">
          <stop offset="0%"  stopColor="#bce39f" />
          <stop offset="35%" stopColor="#52a04a" />
          <stop offset="70%" stopColor="#2d5a3d" />
          <stop offset="100%" stopColor="#173d24" />
        </radialGradient>
        <radialGradient id="cgbw-globe-blue" cx="32%" cy="28%" r="78%">
          <stop offset="0%"  stopColor="#a6d4ee" />
          <stop offset="35%" stopColor="#3f7bb7" />
          <stop offset="70%" stopColor="#1e4a7a" />
          <stop offset="100%" stopColor="#0d2a4a" />
        </radialGradient>
        {/* Specular shine for both globes */}
        <radialGradient id="cgbw-shine" cx="30%" cy="22%" r="32%">
          <stop offset="0%"   stopColor="#ffffff" stopOpacity="0.75" />
          <stop offset="60%"  stopColor="#ffffff" stopOpacity="0.10" />
          <stop offset="100%" stopColor="#ffffff" stopOpacity="0" />
        </radialGradient>
        {/* White drop fill — soft inner gradient */}
        <linearGradient id="cgbw-drop" x1="0%" y1="0%" x2="0%" y2="100%">
          <stop offset="0%"  stopColor="#ffffff" />
          <stop offset="100%" stopColor="#e8eaec" />
        </linearGradient>
        {/* Leaf gradients (glossy bright green) */}
        <radialGradient id="cgbw-leaf" cx="30%" cy="25%" r="80%">
          <stop offset="0%"  stopColor="#d4f5a8" />
          <stop offset="40%" stopColor="#7ec74a" />
          <stop offset="100%" stopColor="#2c7a14" />
        </radialGradient>
        {/* Stem gradient */}
        <linearGradient id="cgbw-stem" x1="0%" y1="100%" x2="0%" y2="0%">
          <stop offset="0%"  stopColor="#2c5a18" />
          <stop offset="100%" stopColor="#5a9d2f" />
        </linearGradient>

        {/* Clip-paths for the globe grid lines */}
        <clipPath id="cgbw-clip-green"><circle cx="40" cy="40" r="30" /></clipPath>
        <clipPath id="cgbw-clip-blue"><circle cx="80" cy="40" r="30" /></clipPath>

        {/* Soft drop shadow */}
        <filter id="cgbw-shadow" x="-20%" y="-20%" width="140%" height="140%">
          <feGaussianBlur in="SourceAlpha" stdDeviation="1.2" />
          <feOffset dx="0" dy="1.4" result="off" />
          <feComponentTransfer result="dark">
            <feFuncA type="linear" slope="0.35" />
          </feComponentTransfer>
          <feMerge>
            <feMergeNode in="dark" />
            <feMergeNode in="SourceGraphic" />
          </feMerge>
        </filter>
      </defs>

      {/* GREEN GLOBE ------------------------------------------------- */}
      <g filter="url(#cgbw-shadow)">
        <circle cx="40" cy="40" r="30" fill="url(#cgbw-globe-green)" />
        {/* Grid: meridians + parallels, clipped inside the circle */}
        <g clipPath="url(#cgbw-clip-green)"
           stroke="rgba(255,255,255,0.35)" strokeWidth="0.5" fill="none">
          {/* parallels (horizontal arcs) */}
          <line x1="10" y1="40" x2="70" y2="40" />
          <path d="M 11 28 Q 40 26 69 28" />
          <path d="M 11 52 Q 40 54 69 52" />
          <path d="M 14 18 Q 40 14 66 18" />
          <path d="M 14 62 Q 40 66 66 62" />
          {/* meridians (vertical ellipses) */}
          <ellipse cx="40" cy="40" rx="8"  ry="30" />
          <ellipse cx="40" cy="40" rx="18" ry="30" />
          <line x1="40" y1="10" x2="40" y2="70" />
        </g>
        {/* Specular shine */}
        <circle cx="40" cy="40" r="30" fill="url(#cgbw-shine)" pointerEvents="none" />
      </g>

      {/* BLUE GLOBE -------------------------------------------------- */}
      <g filter="url(#cgbw-shadow)">
        <circle cx="80" cy="40" r="30" fill="url(#cgbw-globe-blue)" />
        <g clipPath="url(#cgbw-clip-blue)"
           stroke="rgba(255,255,255,0.35)" strokeWidth="0.5" fill="none">
          <line x1="50" y1="40" x2="110" y2="40" />
          <path d="M 51 28 Q 80 26 109 28" />
          <path d="M 51 52 Q 80 54 109 52" />
          <path d="M 54 18 Q 80 14 106 18" />
          <path d="M 54 62 Q 80 66 106 62" />
          <ellipse cx="80" cy="40" rx="8"  ry="30" />
          <ellipse cx="80" cy="40" rx="18" ry="30" />
          <line x1="80" y1="10" x2="80" y2="70" />
        </g>
        <circle cx="80" cy="40" r="30" fill="url(#cgbw-shine)" pointerEvents="none" />
      </g>

      {/* WHITE DROP / LEAF SHAPE AT CENTRE -------------------------- */}
      <g>
        {/* Drop shape: a vertical leaf/teardrop centred on (60, 40) */}
        <path d="M 60 14
                 C 47 26, 45 46, 60 66
                 C 75 46, 73 26, 60 14 Z"
              fill="url(#cgbw-drop)"
              stroke="rgba(0,0,0,0.15)" strokeWidth="0.6" />
        {/* Drop highlight (top-left wet sheen) */}
        <path d="M 54 22 Q 50 32, 53 40"
              stroke="rgba(255,255,255,0.85)" strokeWidth="1.4"
              strokeLinecap="round" fill="none" />
      </g>

      {/* SHINY PLANT INSIDE DROP ----------------------------------- */}
      <g transform="translate(60, 48)">
        {/* Stem */}
        <path d="M 0 10 Q 0 0 0 -8"
              stroke="url(#cgbw-stem)" strokeWidth="1.8"
              strokeLinecap="round" fill="none" />
        {/* Left leaf (rounded teardrop, tilted up-left) */}
        <path d="M 0 -2 C -3 -4 -7 -6 -9 -10 C -10 -14 -5 -15 -2 -12 C 0 -8 0 -5 0 -2 Z"
              fill="url(#cgbw-leaf)"
              stroke="#1f4413" strokeWidth="0.55" strokeLinejoin="round" />
        {/* Left leaf glossy highlight */}
        <path d="M -7 -10 Q -5 -12 -3 -11"
              stroke="#e8ffce" strokeWidth="1.2" strokeLinecap="round"
              fill="none" opacity="0.85" />
        {/* Right leaf (higher, tilted up-right) */}
        <path d="M 0 -6 C 3 -8 7 -10 9 -14 C 10 -18 5 -19 2 -16 C 0 -12 0 -9 0 -6 Z"
              fill="url(#cgbw-leaf)"
              stroke="#1f4413" strokeWidth="0.55" strokeLinejoin="round" />
        {/* Right leaf glossy highlight */}
        <path d="M 3 -16 Q 6 -17 8 -15"
              stroke="#e8ffce" strokeWidth="1.2" strokeLinecap="round"
              fill="none" opacity="0.85" />
        {/* Tiny specular dot on the top leaf */}
        <ellipse cx="6" cy="-16" rx="1" ry="0.5" fill="#ffffff" opacity="0.7" />
      </g>
    </svg>
  );
}

function Logo() {
  // Wordmark: Crop + [G + "Green" label] + [B + "Blue" label] + Water
  return (
    <span className="logo">
      <span className="logo-crop">Crop</span>
      <span className="logo-letter logo-G">
        <span className="logo-letter-label">Green</span>
        <span className="logo-letter-main">G</span>
      </span>
      <span className="logo-letter logo-B">
        <span className="logo-letter-label">Blue</span>
        <span className="logo-letter-main">B</span>
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
        <BrandMark size={50} />
        <Logo />
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
