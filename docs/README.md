# CropsGreenBlueWater — dashboard

Static, dependency-free dashboard built on the aggregated CSV outputs in
`Combined_WF_Outputs/`. Real numbers for 46 crops × 172 countries × 3 years
(2000, 2010, 2020). React + Leaflet via CDN — no build step required.

## Run

The page uses `fetch()` to load JSON, which is blocked under `file://`.
Serve it through any static server:

```
# from this folder
python -m http.server 7733

# then open
http://127.0.0.1:7733/
```

## Rebuild data

If the source CSVs change, regenerate the static JSON bundles:

```
python webapp/scripts/build_static_data.py
```

Outputs into `webapp/data/`:

- `summary.json`         — global totals per year, plus deltas
- `crops.json`           — per-crop totals, group, trend
- `crop_groups.json`     — group-level breakdown for the Crops module
- `countries.json`       — country ranking with green/blue split + trend
- `country_detail/<ISO3>.json` — one file per country with top crops + 3-year history
- `meta.json`            — dataset metadata for the footer

## Files

```
index.html                       page shell, loads React + Leaflet + Babel from CDN
css/styles.css                   design tokens + module styles
js/util.js                       data loader, formatters, color ramps, ISO3 lookups
js/components/
  topbar.jsx                     sticky nav with scrollspy
  hero.jsx                       headline + 4-KPI block + green/blue split
  atlas.jsx                      Leaflet choropleth map (year/metric/scope controls)
  ranking.jsx                    sortable country table
  crops.jsx                      crop-group grid + 46-crop detail
  trends.jsx                     2000-2020 line chart with callouts
  footer.jsx                     methodology + downloads + country detail sheet
  tweaks.jsx                     palette / dark mode / font / density switcher
js/app.jsx                       composition root
scripts/build_static_data.py     CSV → JSON aggregator
```

## What is and isn't here

- Country totals and per-country crop breakdowns are real, from the
  `CountryCropsCSVs/` outputs.
- The atlas is a country-level choropleth using Natural Earth boundaries
  (loaded from `unpkg.com/world-atlas`). Per-grid-cell rasters are not
  served here — see the M4 milestone in the build brief for the tile-server
  path required to ship them.
- Year/metric/crop-scope toggles re-shade the choropleth. The "crop scope"
  filter is currently informational; we'd need per-crop country totals
  baked into a denser JSON to make it interactive — straightforward to add.
- The "Sign in" / "Cite & download" CTA on the top bar links into the
  Methodology section; no auth in v1.
