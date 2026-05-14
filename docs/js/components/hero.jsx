// hero.jsx — editorial hero anchored on the headline +9% trend story.

function Hero({ data }) {
  const g = data.summary.global["2020"];
  const change = data.summary.change_2010_2020;
  const total = g.total_km3;
  const greenPct = (g.green_km3 / total) * 100;
  const bluePct  = (g.blue_km3  / total) * 100;

  return (
    <section id="overview" className="hero" data-screen-label="01 Hero">
      <div className="hero-grid">
        <div className="eyebrow eyebrow-lg">
          <span className="eyebrow-mark" />
          <span>Atlas of agricultural green- and blue water use · 2016 · {data.meta.n_crops} crops</span>
        </div>

        <h1 className="hero-headline">
          The world&rsquo;s farms<br/>
          drink <em className="hl-blue">{Math.round(change.total_pct)}% </em>
          more water<br/>
          than they did <span className="hl-faint">a decade ago.</span>
        </h1>

        <p className="hero-lede">
          Across <Num>{data.meta.n_crops}</Num> crops and <Num>{CGBW.fmt.km3p(g.area_Mha)}</Num> million hectares,
          total crop water consumption rose from <Num>{CGBW.fmt.km3(data.summary.global["2010"].total_km3)}</Num> km³ in 2010
          to <Num>{CGBW.fmt.km3(total)}</Num> km³ in 2020 — driven almost entirely by green water (rainfall absorbed by crops),
          while irrigation barely grew.
        </p>

        <div className="split">
          <div className="split-cols" aria-label="Green vs blue water share of 2020 total">
            <div className="split-col split-green" style={{ flexBasis: `${greenPct}%` }}>
              <div className="split-meta">
                <span className="dot dot-green" />
                <span>Green water</span>
                <span className="split-meta-sub">soil moisture from rainfall</span>
              </div>
              <div className="split-fig">
                <span className="split-pct"><Num>{greenPct.toFixed(1)}</Num><i>%</i></span>
                <span className="split-val"><Num>{CGBW.fmt.km3(g.green_km3)}</Num> km³/yr</span>
              </div>
            </div>
            <div className="split-col split-blue" style={{ flexBasis: `${bluePct}%` }}>
              <div className="split-meta">
                <span className="dot dot-blue" />
                <span>Blue water</span>
                <span className="split-meta-sub">surface + groundwater withdrawal</span>
              </div>
              <div className="split-fig">
                <span className="split-pct"><Num>{bluePct.toFixed(1)}</Num><i>%</i></span>
                <span className="split-val"><Num>{CGBW.fmt.km3(g.blue_km3)}</Num> km³/yr</span>
              </div>
            </div>
          </div>
          <div className="split-axis">
            <span>0</span>
            <span>{CGBW.fmt.km3(total * 0.25)}</span>
            <span>{CGBW.fmt.km3(total * 0.50)}</span>
            <span>{CGBW.fmt.km3(total * 0.75)}</span>
            <span>{CGBW.fmt.km3(total)}&nbsp;km³</span>
          </div>
        </div>

        <div className="kpi-row">
          <Kpi label="Total crop water, 2020"
               value={CGBW.fmt.km3(total)}
               unit="km³/yr"
               delta={`${CGBW.fmt.pct(change.total_pct)} vs 2010`} tone="ink" />
          <Kpi label="Green water growth"
               value={CGBW.fmt.pct(change.green_pct).replace("%","")}
               unit="%"
               delta="2010 → 2020" tone="green" />
          <Kpi label="Blue water growth"
               value={CGBW.fmt.pct(change.blue_pct).replace("%","")}
               unit="%"
               delta="2010 → 2020" tone="blue" />
          <Kpi label="Harvested area, 2020"
               value={CGBW.fmt.km3(g.area_Mha)}
               unit="Mha"
               delta={`${CGBW.fmt.pct(change.area_pct)} vs 2010`} tone="ink" />
        </div>

        <p className="hero-source">
          Source: {data.meta.citation} · Method: {data.meta.model} · Coverage: {data.meta.n_crops} crops, {data.meta.resolution} · Built {data.meta.built_at.slice(0,10)}.
        </p>

        <a className="hero-scroll" href="#atlas">
          <span>Scroll to explore the atlas</span>
          <svg width="14" height="14" viewBox="0 0 14 14" fill="none">
            <path d="M7 2v9m0 0L3 7m4 4l4-4" stroke="currentColor" strokeWidth="1.2" strokeLinecap="round"/>
          </svg>
        </a>
      </div>
    </section>
  );
}

const Num = ({ children }) => <span className="num">{children}</span>;

function Kpi({ label, value, unit, delta, tone = "ink" }) {
  return (
    <div className={`kpi kpi-${tone}`}>
      <div className="kpi-label">{label}</div>
      <div className="kpi-value">
        <span className="num">{value}</span>
        <span className="kpi-unit">{unit}</span>
      </div>
      <div className="kpi-delta">{delta}</div>
    </div>
  );
}

window.Hero = Hero;
window.Num = Num;
window.Kpi = Kpi;
