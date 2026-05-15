// hero.jsx — overview page with count-ups and animated key indicators.

function Hero({ data }) {
  const g = data.summary.global["2020"];
  const change = data.summary.change_2010_2020;
  const total = g.total_km3;
  const greenPct = (g.green_km3 / total) * 100;
  const bluePct  = (g.blue_km3  / total) * 100;

  const headlinePct = CGBW.useCountUp(Math.round(change.total_pct), 900,
    (n) => `${Math.round(n)}%`);
  const totalCount = CGBW.useCountUp(total, 1100,
    (n) => Math.round(n).toLocaleString("en-US"));
  const greenPctCount = CGBW.useCountUp(greenPct, 1100, (n) => n.toFixed(1));
  const bluePctCount  = CGBW.useCountUp(bluePct,  1100, (n) => n.toFixed(1));
  const greenVal = CGBW.useCountUp(g.green_km3, 1100, (n) => Math.round(n).toLocaleString("en-US"));
  const blueVal  = CGBW.useCountUp(g.blue_km3,  1100, (n) => Math.round(n).toLocaleString("en-US"));

  return (
    <section className="hero">
      <div className="eyebrow eyebrow-lg">
        <span className="eyebrow-mark" />
        <span>Atlas of agricultural green- and blue water use · 2026 · {data.meta.n_crops} crops</span>
      </div>

      <h1 className="hero-headline">
        The world&rsquo;s farms drink <em className="hl-blue highlight-on">{headlinePct}</em> more water than they did <span className="hl-faint">a decade ago.</span>
      </h1>

      <p className="hero-lede">
        Across <Num>{data.meta.n_crops}</Num> crops and <Num>{CGBW.fmt.km3p(g.area_Mha)}</Num> Mha,
        total crop water consumption rose from <Num>{CGBW.fmt.km3(data.summary.global["2010"].total_km3)}</Num> km³
        in 2010 to <Num>{CGBW.fmt.km3(total)}</Num> km³ in 2020 — driven almost entirely by
        <em> green water consumption</em> (rainfall), not blue water (irrigation).
      </p>

      <div className="split">
        <div className="split-cols" aria-label="Green water use vs blue water use, 2020">
          <div className="split-col split-green" style={{ flexBasis: `${greenPct}%` }}>
            <div className="split-meta">
              <span className="dot dot-green" />
              <span>Green water use</span>
              <span className="split-meta-sub">soil moisture from rainfall</span>
            </div>
            <div className="split-fig">
              <span className="split-pct"><Num>{greenPctCount}</Num><i>%</i></span>
              <span className="split-val"><Num>{greenVal}</Num> km³/yr</span>
            </div>
          </div>
          <div className="split-col split-blue" style={{ flexBasis: `${bluePct}%` }}>
            <div className="split-meta">
              <span className="dot dot-blue" />
              <span>Blue water use</span>
              <span className="split-meta-sub">surface + groundwater</span>
            </div>
            <div className="split-fig">
              <span className="split-pct"><Num>{bluePctCount}</Num><i>%</i></span>
              <span className="split-val"><Num>{blueVal}</Num> km³/yr</span>
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
             value={totalCount} unit="km³/yr"
             delta={`${CGBW.fmt.pct(change.total_pct)} vs 2010`} tone="ink" />
        <Kpi label="Green water use growth"
             value={CGBW.fmt.pct(change.green_pct).replace("%","")}
             unit="%" delta="2010 → 2020" tone="green" />
        <Kpi label="Blue water use growth"
             value={CGBW.fmt.pct(change.blue_pct).replace("%","")}
             unit="%" delta="2010 → 2020" tone="blue" />
        <Kpi label="Harvested area, 2020"
             value={CGBW.fmt.km3(g.area_Mha)} unit="Mha"
             delta={`${CGBW.fmt.pct(change.area_pct)} vs 2010`} tone="ink" />
      </div>

      <p className="hero-source">
        Source: {data.meta.paper_authors} ({data.meta.paper_year}) · {data.meta.paper_journal} · doi.org/{data.meta.paper_doi}.
        Method: {data.meta.model}. Resolution: {data.meta.resolution}.
      </p>
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
