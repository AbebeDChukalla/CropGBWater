// footer.jsx — Method + full Citation & Data Sources + Country detail sheet.

function MethodFooter({ data }) {
  const m = data.meta;
  const paperCitation =
    `${m.paper_authors} (${m.paper_year}). ${m.paper_title}. ${m.paper_journal}, ${m.paper_volume}. https://doi.org/${m.paper_doi}`;
  const dataCitation =
    `${m.paper_authors} (${m.paper_year}). ${m.data_title} [Data set]. Zenodo. https://doi.org/${m.data_doi}`;

  const bibtex = `@article{chukalla2025cropwater,
  author  = {Chukalla, A. D. and Mekonnen, M. M. and Gunathilake, D. and Wolkeba, F. T. and Gunasekara, B. and Vanham, D.},
  title   = {Global spatially explicit crop water consumption shows an overall increase of 9% for 46 agricultural crops from 2010 to 2020},
  journal = {Nature Food},
  year    = {2025},
  volume  = {6},
  doi     = {${m.paper_doi}},
  url     = {${m.paper_url}}
}`;

  return (
    <section className="method">
      <div className="module-head">
        <div className="module-eyebrow">
          <span className="eyebrow-mark" />
          <span>Module 06 / Methodology · Citation · Downloads</span>
        </div>
        <h2 className="module-title">Citation, data sources &amp; methodology.</h2>
        <p className="module-sub">
          All dashboard figures derive from the dataset cited below.
          Both the paper and the data archive are openly available under CC&nbsp;BY&nbsp;4.0 — please cite both when using the data.
        </p>
      </div>

      <div className="method-grid">
        {/* ── Citations column ── */}
        <div className="method-col">
          <div className="citation-block">
            <h4>Paper citation</h4>
            <p className="citation-text">
              {paperCitation.split(m.paper_url)[0]}
              <a href={m.paper_url} target="_blank" rel="noopener">{m.paper_url}</a>
            </p>
            <div className="citation-actions">
              <a href={m.paper_url} target="_blank" rel="noopener">↗ Open in Nature Food</a>
              <a href={`https://doi.org/${m.paper_doi}`} target="_blank" rel="noopener">DOI</a>
              <a href={`data:text/plain;charset=utf-8,${encodeURIComponent(bibtex)}`}
                 download="cropgbwater.bib">Download BibTeX</a>
              <a href={`data:text/plain;charset=utf-8,${encodeURIComponent(paperCitation)}`}
                 download="cropgbwater-paper-citation.txt">Plain-text citation</a>
            </div>
          </div>

          <div className="citation-block" style={{ marginTop: 14 }}>
            <h4>Data &amp; software citation (Zenodo)</h4>
            <p className="citation-text">
              {dataCitation.split(m.data_url)[0]}
              <a href={m.data_url} target="_blank" rel="noopener">{m.data_url}</a>
            </p>
            <div className="citation-actions">
              <a href={m.data_url} target="_blank" rel="noopener">↗ Open in Zenodo</a>
              <a href={`https://doi.org/${m.data_doi}`} target="_blank" rel="noopener">DOI</a>
              <a href={`data:text/plain;charset=utf-8,${encodeURIComponent(dataCitation)}`}
                 download="cropgbwater-data-citation.txt">Plain-text citation</a>
            </div>
          </div>

          <div className="citation-block" style={{ marginTop: 14 }}>
            <h4>Secondary dataset — SPAM IFPRI</h4>
            <p className="citation-text">
              {m.spam_publisher} ({m.spam_year}). {m.spam_title}.
              {" "}<a href={m.spam_url} target="_blank" rel="noopener">{m.spam_url}</a>
              {" "}Licence: {m.spam_license}.
            </p>
            <div className="citation-actions">
              <a href={m.spam_url} target="_blank" rel="noopener">↗ Open in Harvard Dataverse</a>
              <a href={`https://doi.org/${m.spam_doi}`} target="_blank" rel="noopener">DOI</a>
            </div>
          </div>

          <div className="citation-block" style={{ marginTop: 14 }}>
            <h4>Methodology &amp; metadata</h4>
            <p className="citation-text" style={{ fontSize: 13.5, lineHeight: 1.6 }}>
              Crop water consumption is estimated at a <strong>{m.resolution}</strong> for 46 crops
              across three reference years (2000, 2010, 2020). Green water consumption is the
              portion of precipitation absorbed from the root zone by rainfed and irrigated crops;
              blue water consumption is the portion withdrawn from surface and groundwater to
              irrigate crops. Pixel-level estimates are aggregated to crop, country, continent and
              global totals using SPAM-derived cropland masks.
              Full methodology, input datasets, variable definitions, equations, and uncertainty
              treatment are documented in the paper&rsquo;s Methods + Supplementary Information,
              and the dataset README on Zenodo.
            </p>
            <div className="citation-actions">
              <a href={`${m.paper_url}#Sec`} target="_blank" rel="noopener">↗ Supplementary information</a>
              <a href={`${m.data_url}#metadata`} target="_blank" rel="noopener">↗ Metadata (Zenodo)</a>
            </div>
          </div>
        </div>

        {/* ── Right column: dataset metadata + downloads ── */}
        <div className="method-col">
          <div className="method-meta">
            <Meta label="Coverage"     value={`${m.n_crops} crops · ${m.n_countries} countries · ${m.n_continents} continents`} />
            <Meta label="Resolution"   value={m.resolution} />
            <Meta label="Time range"   value={m.years.join(", ")} />
            <Meta label="License"      value={m.license} />
            <Meta label="Paper DOI"    value={<a href={m.paper_url} target="_blank" rel="noopener">{m.paper_doi}</a>} />
            <Meta label="Dataset DOI"  value={<a href={m.data_url}  target="_blank" rel="noopener">{m.data_doi}</a>} />
            <Meta label="SPAM DOI"     value={<a href={m.spam_url}  target="_blank" rel="noopener">{m.spam_doi}</a>} />
            <Meta label="Last built"   value={m.built_at.slice(0, 10)} />
          </div>

          <div style={{ marginTop: 18 }}>
            <div className="method-eyebrow">Download — this dashboard&rsquo;s feed</div>
            <ul className="action-list">
              <li>
                <span className="action-key">JSON</span>
                <span>Country totals (172 countries)</span>
                <a className="action-size" href="data/countries.json" download>23 KB</a>
              </li>
              <li>
                <span className="action-key">JSON</span>
                <span>Continental rollup (6 continents)</span>
                <a className="action-size" href="data/continents.json" download>1 KB</a>
              </li>
              <li>
                <span className="action-key">JSON</span>
                <span>Per-crop totals · area · production · yield</span>
                <a className="action-size" href="data/crops.json" download>16 KB</a>
              </li>
              <li>
                <span className="action-key">JSON</span>
                <span>Global summary (per year)</span>
                <a className="action-size" href="data/summary.json" download>12 KB</a>
              </li>
            </ul>

            <div className="method-eyebrow" style={{ marginTop: 14 }}>Full dataset — Zenodo</div>
            <ul className="action-list">
              <li>
                <span className="action-key">CSV</span>
                <span>Country × crop water footprint (per year)</span>
                <a className="action-size" href={m.data_url} target="_blank" rel="noopener">↗ Zenodo</a>
              </li>
              <li>
                <span className="action-key">NetCDF</span>
                <span>Gridded outputs · all metrics · 0.083°</span>
                <a className="action-size" href={m.data_url} target="_blank" rel="noopener">↗ Zenodo</a>
              </li>
              <li>
                <span className="action-key">PDF</span>
                <span>Paper (Nature Food)</span>
                <a className="action-size" href={m.paper_url} target="_blank" rel="noopener">↗ Open</a>
              </li>
            </ul>
          </div>
        </div>
      </div>
    </section>
  );
}

function Meta({ label, value }) {
  return (
    <div className="meta-row">
      <span className="meta-label">{label}</span>
      <span className="meta-value">{value}</span>
    </div>
  );
}

window.MethodFooter = MethodFooter;

// ─── Country detail sheet ────────────────────────────────────────────
function CountrySheet({ iso, onClose }) {
  const [detail, setDetail] = React.useState(null);
  const [year, setYear] = React.useState("2020");        // 2010 ↔ 2020 toggle
  const [metric, setMetric] = React.useState("water");   // water | area | yield | wp

  React.useEffect(() => {
    if (!iso) return;
    setDetail(null);
    setYear("2020");
    setMetric("water");
    CGBW.loadCountry(iso).then(setDetail).catch((e) => console.warn(e));
  }, [iso]);

  React.useEffect(() => {
    const onKey = (e) => e.key === "Escape" && onClose();
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [onClose]);

  if (!iso) return null;

  // Year-specific country totals (for the KPI strip)
  const yearTotal = detail && detail.years && detail.years[year] ? detail.years[year] : null;
  const yt = yearTotal || {};
  const ytGreenPct = (yt.total_km3 || 0) > 0 ? Math.round((yt.green_km3 || 0) / yt.total_km3 * 100) : 0;
  const ytBluePct  = (yt.total_km3 || 0) > 0 ? Math.round((yt.blue_km3  || 0) / yt.total_km3 * 100) : 0;

  // Year-specific per-crop list
  const cropsForYear = (detail && detail.crops_by_year && detail.crops_by_year[year])
    || (detail && detail.crops) || [];

  return (
    <div className="country-sheet" onClick={onClose}>
      <div className="country-sheet-body" onClick={(e) => e.stopPropagation()}>
        <button className="cs-close" onClick={onClose}>Close (Esc)</button>
        {!detail ? (
          <div style={{ color: "var(--ink40)", fontFamily: "var(--font-mono)", fontSize: 13 }}>
            Loading {iso}…
          </div>
        ) : (
          <>
            <div className="cs-head">
              <span className="cs-eyebrow">Country detail · {year}</span>
              <h3 className="cs-name">{detail.name}</h3>
              <span style={{ fontFamily: "var(--font-mono)", fontSize: 11.5, color: "var(--ink40)" }}>
                {detail.iso3} · {detail.continent || "—"}
              </span>
            </div>

            {/* Year toggle — drives the whole sheet */}
            <div className="cs-year-toggle" role="tablist" aria-label="Select year">
              {["2010", "2020"].map((y) => (
                <button key={y}
                        role="tab"
                        aria-selected={year === y}
                        className={year === y ? "on" : ""}
                        onClick={() => setYear(y)}>
                  {y}
                </button>
              ))}
            </div>

            <div className="cs-row">
              <div className="cs-cell">
                <span className="cs-cell-label">Total water use</span>
                <span className="cs-cell-val"><span className="num">{CGBW.fmt.km3p(yt.total_km3)}</span><i>km³</i></span>
              </div>
              <div className="cs-cell">
                <span className="cs-cell-label">Green water use</span>
                <span className="cs-cell-val" style={{ color: "var(--green)" }}>
                  <span className="num">{CGBW.fmt.km3p(yt.green_km3)}</span><i>km³ · {ytGreenPct}%</i>
                </span>
              </div>
              <div className="cs-cell">
                <span className="cs-cell-label">Blue water use</span>
                <span className="cs-cell-val" style={{ color: "var(--blue)" }}>
                  <span className="num">{CGBW.fmt.km3p(yt.blue_km3)}</span><i>km³ · {ytBluePct}%</i>
                </span>
              </div>
            </div>

            {/* Only render the harvested-area row if at least one value is non-zero */}
            {(detail.area_2020_Mha > 0 || detail.area_2020_rainfed_Mha > 0 || detail.area_2020_irrigated_Mha > 0) && (
              <div className="cs-row">
                <div className="cs-cell">
                  <span className="cs-cell-label">Harvested area</span>
                  <span className="cs-cell-val"><span className="num">{CGBW.fmt.km3p(detail.area_2020_Mha)}</span><i>Mha</i></span>
                </div>
                <div className="cs-cell">
                  <span className="cs-cell-label">Rainfed</span>
                  <span className="cs-cell-val" style={{ color: "var(--green)" }}>
                    <span className="num">{CGBW.fmt.km3p(detail.area_2020_rainfed_Mha)}</span><i>Mha</i>
                  </span>
                </div>
                <div className="cs-cell">
                  <span className="cs-cell-label">Irrigated</span>
                  <span className="cs-cell-val" style={{ color: "var(--blue)" }}>
                    <span className="num">{CGBW.fmt.km3p(detail.area_2020_irrigated_Mha)}</span><i>Mha{detail.irrigated_share_pct != null && detail.irrigated_share_pct > 0 ? ` · ${detail.irrigated_share_pct}%` : ""}</i>
                  </span>
                </div>
              </div>
            )}

            <div>
              <div className="cs-section-label">Trend · 2000 → 2010 → 2020 · km³/yr</div>
              <div className="cs-trend">
                {Object.entries(detail.years).map(([yr, v]) => (
                  <div key={yr} className="cs-trend-cell">
                    <span className="cs-trend-label">{yr}</span>
                    <span className="cs-trend-val">
                      <span className="num">{CGBW.fmt.km3p(v.total_km3)}</span><i>km³</i>
                    </span>
                  </div>
                ))}
              </div>
              {detail.trend_2010_2020_pct != null && (
                <div style={{ marginTop: 10, fontFamily: "var(--font-mono)", fontSize: 12,
                              color: "var(--accent)" }}>
                  Δ 2010 → 2020: {CGBW.fmt.pct(detail.trend_2010_2020_pct)}
                </div>
              )}
            </div>

            <CropDetailPanel crops={cropsForYear} year={year} metric={metric} setMetric={setMetric} />
            </>
        )}
      </div>
    </div>
  );
}

// ─── Per-crop detail panel ───────────────────────────────────────────
// Shows a metric-tab strip (Water · Area · Yield · WP) and a clean
// per-crop list for one metric at a time, driven by the selected year.
function CropDetailPanel({ crops, year, metric, setMetric }) {
  const METRICS = [
    { id: "water", label: "Water",  unit: "km³/yr",  hint: "% green vs blue" },
    { id: "area",  label: "Area",   unit: "Mha",     hint: "rainfed vs irrigated" },
    { id: "yield", label: "Yield",  unit: "t/ha",    hint: "rainfed vs irrigated" },
    { id: "wp",    label: "WP",     unit: "kg/m³",   hint: "total · rainfed / irrigated" },
  ];
  const max = crops.length > 0 ? (crops[0].total_km3 || 1) : 1;

  return (
    <div className="cs-detail">
      <div className="cs-section-label">
        Top crops · {year} — pick a metric below
      </div>

      <div className="cs-metric-tabs" role="tablist" aria-label="Metric">
        {METRICS.map((m) => (
          <button key={m.id}
                  role="tab"
                  aria-selected={metric === m.id}
                  className={metric === m.id ? "on" : ""}
                  onClick={() => setMetric(m.id)}>
            <span className="cs-metric-tab-label">{m.label}</span>
            <span className="cs-metric-tab-unit">{m.unit}</span>
          </button>
        ))}
      </div>

      <div className="cs-crop-cards">
        {crops.slice(0, 14).map((c) => (
          <CropMetricRow key={c.code} c={c} metric={metric} max={max} />
        ))}
        {crops.length === 0 && (
          <div className="cs-no-data">No crop records for {year}.</div>
        )}
      </div>

      <div className="cs-crop-legend">
        <span><i className="dot dot-green" /> rainfed</span>
        <span><i className="dot dot-blue"  /> irrigated</span>
        <span>n/d = not reported in source</span>
      </div>
    </div>
  );
}

function CropMetricRow({ c, metric, max }) {
  // Pull values per metric
  let bar = 0, primary = "", subParts = null;
  if (metric === "water") {
    bar = (c.total_km3 / max) * 100;
    primary = CGBW.fmt.km3p(c.total_km3) + " km³";
    subParts = (
      <>
        <span style={{ color: "var(--green)" }}>{c.green_pct}% green</span>
        <span style={{ color: "var(--ink40)" }}> · </span>
        <span style={{ color: "var(--blue)" }}>{c.blue_pct}% blue</span>
      </>
    );
  } else if (metric === "area") {
    const rf  = c.area_rainfed_Mha,   ir = c.area_irrigated_Mha;
    const tot = (rf || 0) + (ir || 0);
    bar = max > 0 ? (tot / max) * 50 : 0;     // independent scale, capped at 50%
    primary = (rf != null || ir != null) ? (tot.toFixed(2) + " Mha") : "n/d";
    subParts = (rf != null || ir != null) ? (
      <>
        <span style={{ color: "var(--green)" }}>{rf != null ? rf.toFixed(2) : "—"} rf</span>
        <span style={{ color: "var(--ink40)" }}> · </span>
        <span style={{ color: "var(--blue)" }}>{ir != null ? ir.toFixed(2) : "—"} irr</span>
      </>
    ) : null;
  } else if (metric === "yield") {
    const rf = c.yield_rainfed_ton_ha, ir = c.yield_irrigated_ton_ha;
    primary = (rf != null && ir != null) ? `${rf.toFixed(2)} / ${ir.toFixed(2)} t/ha`
            : (rf != null ? `${rf.toFixed(2)} rf · — irr` : (ir != null ? `— rf · ${ir.toFixed(2)} irr` : "n/d"));
    subParts = (rf != null || ir != null) ? (
      <>
        <span style={{ color: "var(--green)" }}>{rf != null ? rf.toFixed(2) : "—"} rainfed</span>
        <span style={{ color: "var(--ink40)" }}> · </span>
        <span style={{ color: "var(--blue)" }}>{ir != null ? ir.toFixed(2) : "—"} irrigated</span>
      </>
    ) : null;
    // bar uses max yield among top crops as scale
    bar = (rf || ir) ? ((Math.max(rf || 0, ir || 0)) / 15) * 100 : 0;  // assume top yield ~15 t/ha
  } else { // wp
    const t = c.wp_total_kg_m3, rf = c.wp_rainfed_kg_m3, ir = c.wp_irrigated_kg_m3;
    primary = t != null ? `${t.toFixed(2)} kg/m³` : "n/d";
    subParts = (rf != null || ir != null) ? (
      <>
        <span style={{ color: "var(--green)" }}>{rf != null ? rf.toFixed(2) : "—"} rainfed</span>
        <span style={{ color: "var(--ink40)" }}> · </span>
        <span style={{ color: "var(--blue)" }}>{ir != null ? ir.toFixed(2) : "—"} irrigated</span>
      </>
    ) : null;
    bar = t != null ? Math.min(100, (t / 6) * 100) : 0;   // scale 0-6 kg/m³
  }

  return (
    <div className="cs-crop-card">
      <div className="cs-crop-card-head">
        <span className="cs-crop-card-name" title={c.code}>{c.code}</span>
        <span className="cs-crop-card-primary">{primary}</span>
      </div>
      {bar > 0 && (
        <div className="cs-crop-card-bar">
          <i style={{ width: `${Math.max(2, bar)}%`,
                      background: metric === "water"
                        ? "var(--ink)"
                        : metric === "area" ? "var(--ink40)"
                        : metric === "yield" ? "var(--greenSoft)"
                        : "var(--accent)" }} />
        </div>
      )}
      {subParts && (
        <div className="cs-crop-card-sub">{subParts}</div>
      )}
    </div>
  );
}

window.CountrySheet = CountrySheet;
