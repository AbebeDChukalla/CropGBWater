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
            <h4>Methodology &amp; metadata</h4>
            <p className="citation-text" style={{ fontSize: 13.5, lineHeight: 1.6 }}>
              Crop water consumption is estimated at a <strong>{m.resolution}</strong> using a coupled
              crop-growth and hydrological model (<strong>{m.model}</strong>). Green water consumption is
              partitioned from precipitation reaching the root zone; blue water consumption is the
              irrigation fraction sourced from surface and groundwater withdrawals. Grid cells are
              aggregated to country totals via SPAM-derived cropland masks.
              Full methodology, variable definitions, and uncertainty bounds are documented in the
              paper Supplementary Information; metadata for the gridded NetCDF outputs is in the
              Zenodo archive&rsquo;s <code>README</code>.
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
            <Meta label="Model"        value={m.model} />
            <Meta label="License"      value={m.license} />
            <Meta label="Paper DOI"    value={<a href={m.paper_url} target="_blank" rel="noopener">{m.paper_doi}</a>} />
            <Meta label="Dataset DOI"  value={<a href={m.data_url}  target="_blank" rel="noopener">{m.data_doi}</a>} />
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
  React.useEffect(() => {
    if (!iso) return;
    setDetail(null);
    CGBW.loadCountry(iso).then(setDetail).catch((e) => console.warn(e));
  }, [iso]);

  React.useEffect(() => {
    const onKey = (e) => e.key === "Escape" && onClose();
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [onClose]);

  if (!iso) return null;

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
              <span className="cs-eyebrow">Country detail · 2020</span>
              <h3 className="cs-name">{detail.name}</h3>
              <span style={{ fontFamily: "var(--font-mono)", fontSize: 11.5, color: "var(--ink40)" }}>
                {detail.iso3} · {detail.continent || "—"}
              </span>
            </div>

            <div className="cs-row">
              <div className="cs-cell">
                <span className="cs-cell-label">Total water use</span>
                <span className="cs-cell-val"><span className="num">{CGBW.fmt.km3p(detail.total_km3)}</span><i>km³</i></span>
              </div>
              <div className="cs-cell">
                <span className="cs-cell-label">Green water use</span>
                <span className="cs-cell-val" style={{ color: "var(--green)" }}>
                  <span className="num">{CGBW.fmt.km3p(detail.green_km3)}</span><i>km³ · {detail.green_pct}%</i>
                </span>
              </div>
              <div className="cs-cell">
                <span className="cs-cell-label">Blue water use</span>
                <span className="cs-cell-val" style={{ color: "var(--blue)" }}>
                  <span className="num">{CGBW.fmt.km3p(detail.blue_km3)}</span><i>km³ · {detail.blue_pct}%</i>
                </span>
              </div>
            </div>

            {detail.area_2020_Mha != null && (
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
                    <span className="num">{CGBW.fmt.km3p(detail.area_2020_irrigated_Mha)}</span><i>Mha · {detail.irrigated_share_pct ?? "—"}%</i>
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

            <div>
              <div className="cs-section-label">Top crops by water consumption · 2020</div>
              <div className="cs-crop-list">
                {detail.crops.slice(0, 10).map((c) => {
                  const max = detail.crops[0].total_km3 || 1;
                  return (
                    <div key={c.code} className="cs-crop-row">
                      <span className="cs-crop-name">{c.code}</span>
                      <div className="cs-crop-bar" style={{ width: `${(c.total_km3 / max) * 100}%` }}>
                        <i className="seg-green" style={{ flexBasis: `${c.green_pct}%` }} />
                        <i className="seg-blue"  style={{ flexBasis: `${c.blue_pct }%` }} />
                      </div>
                      <span className="cs-crop-val">{CGBW.fmt.km3p(c.total_km3)} km³</span>
                    </div>
                  );
                })}
              </div>
            </div>
          </>
        )}
      </div>
    </div>
  );
}

window.CountrySheet = CountrySheet;
