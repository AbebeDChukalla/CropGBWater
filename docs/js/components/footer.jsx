// footer.jsx — methodology + download + bottom strip.

function MethodFooter({ data }) {
  const m = data.meta;
  return (
    <footer id="method" className="method-footer" data-screen-label="06 Methodology">
      <div className="method-grid">
        <div className="method-col">
          <div className="method-eyebrow">Methodology</div>
          <h3 className="method-title">How the numbers are made.</h3>
          <p className="method-body">
            Crop water consumption is estimated at a {m.resolution} using a coupled
            crop-growth and hydrological model. Green water is partitioned from
            precipitation reaching the root zone; blue water is the irrigation
            fraction sourced from surface and groundwater withdrawals. Grid cells
            are aggregated to country totals via SPAM-derived cropland masks. All
            figures on this dashboard come from those aggregations — see the
            paper for the full chain.
          </p>
        </div>

        <div className="method-col method-meta">
          <Meta label="Coverage"     value={`${m.n_crops} crops · ${m.n_countries} countries`} />
          <Meta label="Resolution"   value={m.resolution} />
          <Meta label="Time range"   value={m.years.join(", ")} />
          <Meta label="Model"        value={m.model} />
          <Meta label="Citation"     value={m.citation.length > 60 ? m.citation.slice(0, 60) + "…" : m.citation} />
          <Meta label="License"      value={m.license} />
          <Meta label="Last built"   value={m.built_at.slice(0, 10)} />
          <Meta label="DOI"          value={m.doi} />
        </div>

        <div className="method-col method-actions">
          <div className="method-eyebrow">Download & cite</div>
          <ul className="action-list">
            <li>
              <span className="action-key">JSON</span>
              <span>Country totals (this dashboard&rsquo;s feed)</span>
              <a className="action-size" href="data/countries.json" download>23 KB</a>
            </li>
            <li>
              <span className="action-key">JSON</span>
              <span>Crops × 3 years (per-crop bundle)</span>
              <a className="action-size" href="data/crops.json" download>16 KB</a>
            </li>
            <li>
              <span className="action-key">CSV</span>
              <span>Country × crop water footprint (2020)</span>
              <a className="action-size" href="../CountryCropsCSVs/Y2020_country_WaterFootPrint.csv" download>0.9 MB</a>
            </li>
            <li>
              <span className="action-key">PDF</span>
              <span>Original methodology paper</span>
              <a className="action-size" href="../s43016-025-01231-x.pdf" download>4 MB</a>
            </li>
            <li>
              <span className="action-key">BibTeX</span>
              <span>Cite this dataset</span>
              <span className="action-size">0.4 KB</span>
            </li>
          </ul>
        </div>
      </div>

      <div className="method-bottom">
        <span>CropsGreenBlueWater · {m.release}</span>
        <span>·</span>
        <span>Data: {m.license}</span>
        <span>·</span>
        <span>Built for researchers, policymakers &amp; the curious public</span>
        <span className="method-bottom-spacer" />
        <span>v{m.version}</span>
      </div>
    </footer>
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
              <span style={{ fontFamily: "var(--font-mono)", fontSize: 11, color: "var(--ink40)" }}>
                {detail.iso3}
              </span>
            </div>

            <div className="cs-row">
              <div className="cs-cell">
                <span className="cs-cell-label">Total water</span>
                <span className="cs-cell-val"><span className="num">{CGBW.fmt.km3p(detail.total_km3)}</span><i>km³</i></span>
              </div>
              <div className="cs-cell">
                <span className="cs-cell-label">Green</span>
                <span className="cs-cell-val" style={{ color: "var(--green)" }}>
                  <span className="num">{CGBW.fmt.km3p(detail.green_km3)}</span><i>km³ · {detail.green_pct}%</i>
                </span>
              </div>
              <div className="cs-cell">
                <span className="cs-cell-label">Blue</span>
                <span className="cs-cell-val" style={{ color: "var(--blue)" }}>
                  <span className="num">{CGBW.fmt.km3p(detail.blue_km3)}</span><i>km³ · {detail.blue_pct}%</i>
                </span>
              </div>
            </div>

            <div>
              <div className="cs-section-label">Trend · 2000 → 2010 → 2020</div>
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
                <div style={{ marginTop: 12, fontFamily: "var(--font-mono)", fontSize: 12,
                              color: "var(--accent)" }}>
                  Δ 2010 → 2020: {CGBW.fmt.pct(detail.trend_2010_2020_pct)}
                </div>
              )}
            </div>

            <div>
              <div className="cs-section-label">Top crops by water consumption · 2020</div>
              <div className="cs-crop-list">
                {detail.crops.slice(0, 12).map((c) => {
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
              {detail.crops.length > 12 && (
                <div style={{ marginTop: 10, fontFamily: "var(--font-mono)", fontSize: 11,
                              color: "var(--ink40)" }}>
                  Showing 12 of {detail.crops.length} crops.
                </div>
              )}
            </div>
          </>
        )}
      </div>
    </div>
  );
}

window.CountrySheet = CountrySheet;
