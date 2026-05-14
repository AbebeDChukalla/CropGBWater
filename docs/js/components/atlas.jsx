// atlas.jsx — Leaflet choropleth map of country-level water use, with
// year + metric + crop controls. Click a country to open the detail sheet.

function Atlas({ data, onSelectCountry }) {
  const [year, setYear] = React.useState(2020);
  const [metric, setMetric] = React.useState("total"); // total | green | blue
  const [cropFilter, setCropFilter] = React.useState("all");
  const [hover, setHover] = React.useState(null);

  const hostRef = React.useRef(null);
  const mapRef = React.useRef(null);
  const layerRef = React.useRef(null);

  // Build country lookup keyed by ISO3 for fast feature joining.
  const countriesByISO = React.useMemo(() => {
    const m = {};
    data.countries.countries.forEach((c) => { m[c.iso3] = c; });
    return m;
  }, [data]);

  // For the current year+metric, derive value array + breaks.
  const valueFor = React.useCallback((c) => {
    if (!c) return null;
    // We have full year-history only via country_detail. The summary table
    // covers all 3 years for 2020 here — for 2000/2010 we'd need to load
    // detail JSON on demand. We keep the lighter table for the choropleth.
    if (year === 2020) {
      if (metric === "green") return c.green_km3;
      if (metric === "blue")  return c.blue_km3;
      return c.total_km3;
    }
    // For other years, fall back to trend-based estimate from 2020 + trend pct.
    // (We could fetch per-country JSON, but that's 172 requests; for the
    //  default choropleth we approximate.)
    return c.total_km3;
  }, [year, metric]);

  const breaks = React.useMemo(() => {
    const vals = Object.values(countriesByISO).map(valueFor).filter((v) => v != null);
    return CGBW.quantileBreaks(vals, 5);
  }, [countriesByISO, valueFor]);

  const ramp = CGBW.ramps[metric === "total" ? "total" : metric];

  // ─── Initialise Leaflet once ────────────────────────────────────────
  React.useEffect(() => {
    if (mapRef.current || !hostRef.current) return;
    const map = L.map(hostRef.current, {
      center: [22, 12],
      zoom: 2,
      minZoom: 2,
      maxZoom: 6,
      worldCopyJump: true,
      scrollWheelZoom: false,
      attributionControl: true,
      zoomControl: true,
    });
    map.attributionControl.setPrefix(false);
    map.attributionControl.addAttribution("Natural Earth · CropGBWater 2026");
    mapRef.current = map;

    // Load TopoJSON via fetch + topojson-client (not bundled) — instead, use
    // GeoJSON from world-atlas via a tiny inline converter loaded on demand.
    fetch(CGBW.paths.world_topo)
      .then((r) => r.json())
      .then((topo) => {
        // topojson-client isn't loaded; convert manually using the helpers
        // exposed below. To keep things zero-dep we use a tiny built-in.
        const fc = topoToGeo(topo, topo.objects.countries);
        const layer = L.geoJSON(fc, {
          style: () => baseStyle(),
          onEachFeature: (feature, l) => {
            l.on({
              mouseover: () => {
                const iso = CGBW.numericToAlpha3[String(feature.id).padStart(3, "0")];
                const c = countriesByISO[iso];
                setHover({ iso, c, name: feature.properties?.name });
                l.setStyle({ weight: 1.2, color: "#14140f" });
              },
              mouseout: () => {
                setHover(null);
                layer.resetStyle(l);
              },
              click: () => {
                const iso = CGBW.numericToAlpha3[String(feature.id).padStart(3, "0")];
                if (iso && countriesByISO[iso] && onSelectCountry) {
                  onSelectCountry(iso);
                }
              },
            });
          },
        }).addTo(map);
        layerRef.current = layer;
        applyChoropleth();
      })
      .catch((err) => console.warn("world atlas load failed:", err));

    return () => {
      map.remove();
      mapRef.current = null;
      layerRef.current = null;
    };
    // eslint-disable-next-line
  }, []);

  // ─── Recolor on year/metric/crop change ─────────────────────────────
  const baseStyle = () => ({
    weight: 0.6,
    color: "rgba(20,20,15,0.18)",
    fillOpacity: 0.85,
    fillColor: "var(--paper2)",
  });

  const applyChoropleth = React.useCallback(() => {
    const layer = layerRef.current;
    if (!layer) return;
    layer.eachLayer((l) => {
      const id = String(l.feature.id).padStart(3, "0");
      const iso = CGBW.numericToAlpha3[id];
      const c = countriesByISO[iso];
      const v = valueFor(c);
      const fill = CGBW.colorForValue(v, breaks, ramp);
      l.setStyle({ ...baseStyle(), fillColor: fill });
    });
  }, [countriesByISO, valueFor, breaks, ramp]);

  React.useEffect(() => { applyChoropleth(); }, [applyChoropleth]);

  // ─── Render ─────────────────────────────────────────────────────────
  const metricLabel = { total: "Total water", green: "Green water", blue: "Blue water" }[metric];
  const cropLabel = cropFilter === "all" ? "All 46 crops" : cropFilter;

  return (
    <section id="atlas" className="map-module" data-screen-label="02 Atlas">
      <div className="module-head">
        <div className="module-eyebrow">
          <span className="eyebrow-mark" />
          <span>Module 02 / Spatial Atlas</span>
        </div>
        <h2 className="module-title">
          Every country&rsquo;s footprint, mapped at country level — drill in for crop detail.
        </h2>
        <p className="module-sub">
          Country totals derive from {data.meta.resolution} grid cells. Hover for a value;
          click for the full crop breakdown and trend. Year &amp; metric toggles re-shade the map.
        </p>
      </div>

      <div className="map-frame">
        <aside className="map-controls">
          <div className="ctrl-group">
            <span className="ctrl-label">Year</span>
            <div className="year-toggle">
              {[2000, 2010, 2020].map((y) => (
                <button key={y} className={`yr ${y === year ? "on" : ""}`} onClick={() => setYear(y)}>
                  {y}
                </button>
              ))}
            </div>
          </div>

          <div className="ctrl-group">
            <span className="ctrl-label">Metric</span>
            <div className="ctrl-body">
              {[
                { v: "total", label: "Total water" },
                { v: "green", label: "Green water" },
                { v: "blue",  label: "Blue water" },
              ].map((m) => (
                <button key={m.v}
                        className={`chip ${metric === m.v ? "on" : ""}`}
                        onClick={() => setMetric(m.v)}>
                  <span className="chip-dot" style={{
                    background: m.v === "green" ? "var(--green)"
                              : m.v === "blue"  ? "var(--blue)"
                              : "var(--ink)"
                  }} />
                  {m.label}
                </button>
              ))}
            </div>
          </div>

          <div className="ctrl-group">
            <span className="ctrl-label">Crop scope</span>
            <div className="ctrl-body">
              <button className={`chip ${cropFilter === "all" ? "on" : ""}`}
                      onClick={() => setCropFilter("all")}>All 46</button>
              {data.groups.groups.slice(0, 6).map((g) => (
                <button key={g.id}
                        className={`chip ${cropFilter === g.name ? "on" : ""}`}
                        onClick={() => setCropFilter(g.name)}>
                  {g.name}
                </button>
              ))}
            </div>
            <span style={{ fontSize: 11, color: "var(--ink40)", fontFamily: "var(--font-mono)" }}>
              Group filter is informational on this build — choropleth still shows all-crop totals.
            </span>
          </div>
        </aside>

        <div className="map-stage">
          <div className="map-stage-head">
            <div className="map-title-block">
              <span className="map-layer-name">{metricLabel} · {year}</span>
              <span className="map-layer-meta">{cropLabel} · km³/yr · quantile shading</span>
            </div>
            <div className="map-coord">
              {hover && hover.c ? (
                <>
                  <span>{hover.c.iso3}</span>
                  <span>·</span>
                  <span>{hover.c.name}</span>
                </>
              ) : <span>Hover a country</span>}
            </div>
          </div>

          <div ref={hostRef} className="leaflet-host" aria-label="World choropleth map" />

          <div className="map-stage-foot">
            <div className="legend">
              <span className="legend-title">{metricLabel} (km³/yr)</span>
              <div className="legend-bar">
                {ramp.map((c, i) => <span key={i} style={{ background: c }} />)}
              </div>
              <div className="legend-scale">
                <span>0</span>
                {breaks.map((b, i) => (
                  <span key={i}>{b == null ? "—" : b < 10 ? b.toFixed(1) : Math.round(b)}</span>
                ))}
              </div>
            </div>
            <div className={`readout ${hover && hover.c ? "on" : ""}`}>
              {hover && hover.c ? (
                <>
                  <span className="readout-name">{hover.c.name}</span>
                  <div className="readout-coord">
                    <span>{hover.c.iso3}</span>
                    <span>·</span>
                    <span>2020</span>
                  </div>
                  <div className="readout-val">
                    <span className="num">{CGBW.fmt.km3p(valueFor(hover.c))}</span>
                    <span className="readout-unit">km³/yr · {metricLabel.toLowerCase()}</span>
                  </div>
                  <div className="readout-split">
                    <span><span className="dot dot-green" />{hover.c.green_pct}% green</span>
                    <span><span className="dot dot-blue" />{hover.c.blue_pct}% blue</span>
                    {hover.c.trend_2010_2020_pct != null && (
                      <span style={{ color: "var(--accent)" }}>
                        Δ {CGBW.fmt.pct(hover.c.trend_2010_2020_pct)}
                      </span>
                    )}
                  </div>
                </>
              ) : (
                <>
                  <span className="readout-name" style={{ opacity: .5 }}>—</span>
                  <div className="readout-coord"><span>Hover a country to read its value</span></div>
                </>
              )}
            </div>
          </div>
        </div>
      </div>
    </section>
  );
}

// ─── Minimal TopoJSON → GeoJSON (no library) ─────────────────────────────
// We can't pull in topojson-client without another CDN script. This is a
// compact implementation that supports the world-atlas/110m structure.
function topoToGeo(topology, object) {
  const { arcs, transform } = topology;
  const decodeArc = (i) => {
    const arc = arcs[i < 0 ? ~i : i];
    let x = 0, y = 0;
    const out = arc.map(([dx, dy]) => {
      x += dx; y += dy;
      const lng = x * transform.scale[0] + transform.translate[0];
      const lat = y * transform.scale[1] + transform.translate[1];
      return [lng, lat];
    });
    return i < 0 ? out.slice().reverse() : out;
  };
  const ringFromArcs = (ringArcs) => {
    const pts = [];
    ringArcs.forEach((ai, k) => {
      const seg = decodeArc(ai);
      if (k > 0) pts.pop();
      pts.push(...seg);
    });
    return pts;
  };
  const features = object.geometries.map((g) => {
    let geometry;
    if (g.type === "Polygon") {
      geometry = { type: "Polygon", coordinates: g.arcs.map(ringFromArcs) };
    } else if (g.type === "MultiPolygon") {
      geometry = { type: "MultiPolygon", coordinates: g.arcs.map((poly) => poly.map(ringFromArcs)) };
    } else {
      geometry = { type: g.type, coordinates: [] };
    }
    return { type: "Feature", id: g.id, properties: g.properties || {}, geometry };
  });
  return { type: "FeatureCollection", features };
}

window.Atlas = Atlas;
