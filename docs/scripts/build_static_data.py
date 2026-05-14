"""Aggregate the CropGBWater CSV outputs into the static JSON bundles the
webapp consumes. Run from repo root:

    python webapp/scripts/build_static_data.py

Reads from Combined_WF_Outputs/:
  Cross_year_summary_WC.csv
  Cross_year_summary_area.csv
  Cross_year_summary_production.csv
  Crop_summary_{2000,2010,2020}.csv
  CountryCropsCSVs/Y{2000,2010,2020}_country_WaterFootPrint.csv
  CountryCropsCSVs/Y{2000,2010,2020}_country_HarvestedArea.csv
  CountryCropsCSVs/Y{2000,2010,2020}_country_CropYield.csv

Writes:
  webapp/data/summary.json
  webapp/data/crops.json
  webapp/data/crop_groups.json
  webapp/data/countries.json
  webapp/data/country_detail/{ISO3}.json
  webapp/data/meta.json
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT  = ROOT / "webapp" / "data"
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "country_detail").mkdir(parents=True, exist_ok=True)

YEARS = [2000, 2010, 2020]

# Crop -> group mapping (from Crop_summary_*.csv `group` column, kept here so
# we can classify by code without re-reading every time).
CROP_GROUPS = {
    "Cereals":        ["BARL","MAIZ","OCER","PMIL","RICE","SMIL","SORG","WHEA"],
    "Oil crops":      ["GROU","OOIL","OPAL","RAPE","SESA","SOYB","SUNF"],
    "Roots & tubers": ["CASS","OROO","POTA","SWPO","YAMS"],
    "Pulses":         ["BEAN","CHIC","COWP","LENT","OPUL","PIGE"],
    "Fibres":         ["COTT","OFIB","RUBB"],
    "Sugar crops":    ["SUGB","SUGC"],
    "Fruits":         ["BANA","CITR","PLNT","TEMF","TROF"],
    "Vegetables":     ["ONIO","TOMA","VEGE"],
    "Beverages":      ["ACOF","RCOF","COCO","TEAS"],
    "Other":          ["CNUT","OCRO","TOBA"],
}

CROP_GROUP_LOOKUP: dict[str, str] = {}
for g, codes in CROP_GROUPS.items():
    for c in codes:
        CROP_GROUP_LOOKUP[c] = g

ISO3_TO_NAME: dict[str, str] = {}


def _safe(v):
    """Return None for NaN / ±inf, else round-friendly float."""
    if v is None:
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    if math.isnan(f) or math.isinf(f):
        return None
    return f


def _round(v, digits=2):
    f = _safe(v)
    return None if f is None else round(f, digits)


import re

_NAME_OVERRIDES = {
    "UnitedStates": "United States",
    "UnitedKingdom": "United Kingdom",
    "RussianFederation": "Russia",
    "IranIslamicRepublicof": "Iran",
    "VietNam": "Vietnam",
    "SyrianArabRepublic": "Syria",
    "RepublicofKorea": "South Korea",
    "DemocraticPeoplesRepublicofKorea": "North Korea",
    "DemocraticRepublicoftheCongo": "DR Congo",
    "Congo": "Congo",
    "CotedIvoire": "Côte d'Ivoire",
    "CoteDivoire": "Côte d'Ivoire",
    "Czechia": "Czechia",
    "CzechRepublic": "Czechia",
    "UnitedRepublicofTanzania": "Tanzania",
    "BoliviaPlurinationalStateof": "Bolivia",
    "VenezuelaBolivarianRepublicof": "Venezuela",
    "TaiwanProvinceofChina": "Taiwan",
    "LaoPeoplesDemocraticRepublic": "Laos",
    "RepublicofMoldova": "Moldova",
    "PapuaNewGuinea": "Papua New Guinea",
    "SouthAfrica": "South Africa",
    "SouthSudan": "South Sudan",
    "SaudiArabia": "Saudi Arabia",
    "DominicanRepublic": "Dominican Republic",
    "CostaRica": "Costa Rica",
    "PuertoRico": "Puerto Rico",
    "ElSalvador": "El Salvador",
    "SriLanka": "Sri Lanka",
    "BurkinaFaso": "Burkina Faso",
    "SierraLeone": "Sierra Leone",
    "NewZealand": "New Zealand",
    "EquatorialGuinea": "Equatorial Guinea",
    "WesternSahara": "Western Sahara",
    "CentralAfricanRepublic": "Central African Republic",
    "CapeVerde": "Cape Verde",
    "NorthMacedonia": "North Macedonia",
    "BosniaandHerzegovina": "Bosnia and Herzegovina",
    "Turkiye": "Türkiye",
    "Turkey": "Türkiye",
    "ChinaHongKongSAR": "Hong Kong",
    "ChinaTaiwanProvinceof": "Taiwan",
    "ChinamainlandPartoftheChinaPRC": "China",
    "AntiguaandBarbuda": "Antigua and Barbuda",
    "BruneiDarussalam": "Brunei",
    "TrinidadandTobago": "Trinidad and Tobago",
    "SaintVincentandtheGrenadines": "St. Vincent",
    "SaintLucia": "St. Lucia",
    "SaintKittsandNevis": "St. Kitts",
}

def clean_country_name(raw: str) -> str:
    """Map known camelcase names + insert spaces before caps as fallback."""
    if not isinstance(raw, str):
        return raw
    raw = raw.strip()
    if raw in _NAME_OVERRIDES:
        return _NAME_OVERRIDES[raw]
    # Generic camelcase → "Camel Case"
    return re.sub(r"(?<=[a-z])(?=[A-Z])", " ", raw)


# ─── Cross-year global totals ────────────────────────────────────────────────

def _wide_yearly(path: Path, metric_order: list[str]) -> tuple[pd.DataFrame, dict[int, dict[str, str]]]:
    """Read a Cross_year_summary_*.csv whose year header doesn't forward-fill.

    Returns (data_frame_without_unit_row, {year: {metric_name: column_name}}).
    """
    df = pd.read_csv(path, header=1, low_memory=False)  # use 2nd row as header
    # First column is Crop_Code, 2nd Crop_Name, 3rd Unit
    df.columns = list(df.columns)
    # Drop the unit row if present (the cell under "Crop_Code" is NaN there)
    df = df[df.iloc[:, 0].notna()]
    # Build per-year column maps (3 fixed columns then N metrics per year)
    cols_per_year = len(metric_order)
    year_cols: dict[int, dict[str, str]] = {}
    # Skip header offset: original columns 0,1,2 are Crop_Code, Crop_Name, Unit;
    # the file we just parsed has its first data-row as numbers — actual column
    # names use the metric names plus suffixes like ".1", ".2" for repeats.
    cols = list(df.columns)
    base = cols[3:3 + cols_per_year]
    for yi, yr in enumerate(YEARS):
        start = 3 + yi * cols_per_year
        block = cols[start:start + cols_per_year]
        year_cols[yr] = dict(zip(metric_order, block))
    return df, year_cols


def build_summary() -> dict:
    """Global green/blue/total per year, plus area and production."""
    wc_metrics = ["WCgn_irr_m3", "WCgn_rf_m3", "WC_gn_m3", "WCbl_m3", "WCpaddy_m3", "WC_Total_m3"]
    df, year_cols = _wide_yearly(ROOT / "Cross_year_summary_WC.csv", wc_metrics)
    crop_codes = df.iloc[:, 0].astype(str).tolist()
    crop_names = df.iloc[:, 1].astype(str).tolist()

    out: dict = {"years": YEARS, "global": {}, "by_crop": {}}

    for yr in YEARS:
        gn = pd.to_numeric(df[year_cols[yr]["WC_gn_m3"]],   errors="coerce").sum()
        bl = pd.to_numeric(df[year_cols[yr]["WCbl_m3"]],    errors="coerce").sum()
        tot = pd.to_numeric(df[year_cols[yr]["WC_Total_m3"]], errors="coerce").sum()
        # m3/yr -> km3/yr
        out["global"][str(yr)] = {
            "green_km3": _round(gn / 1e9, 1),
            "blue_km3":  _round(bl / 1e9, 1),
            "total_km3": _round(tot / 1e9, 1),
        }

    # Per-crop totals across all years
    for i, code in enumerate(crop_codes):
        if not isinstance(code, str) or code == "nan":
            continue
        crop_row = {"code": code, "name": crop_names[i],
                    "group": CROP_GROUP_LOOKUP.get(code, "Other")}
        for yr in YEARS:
            gn = pd.to_numeric(df[year_cols[yr]["WC_gn_m3"]],   errors="coerce").iloc[i]
            bl = pd.to_numeric(df[year_cols[yr]["WCbl_m3"]],    errors="coerce").iloc[i]
            tot = pd.to_numeric(df[year_cols[yr]["WC_Total_m3"]], errors="coerce").iloc[i]
            crop_row[str(yr)] = {
                "green_km3": _round((gn or 0) / 1e9, 2),
                "blue_km3":  _round((bl or 0) / 1e9, 2),
                "total_km3": _round((tot or 0) / 1e9, 2),
            }
        out["by_crop"][code] = crop_row

    # Area + production
    area_df, area_cols = _wide_yearly(ROOT / "Cross_year_summary_area.csv", ["RF_ha", "IRR_ha", "Total_ha"])
    prod_df, prod_cols = _wide_yearly(ROOT / "Cross_year_summary_production.csv",
                                       ["production_rf_ton", "production_irr_ton", "production_total_ton"])

    for yr in YEARS:
        ha_total = pd.to_numeric(area_df[area_cols[yr]["Total_ha"]], errors="coerce").sum()
        prod_total = pd.to_numeric(prod_df[prod_cols[yr]["production_total_ton"]], errors="coerce").sum()
        out["global"][str(yr)]["area_Mha"] = _round(ha_total / 1e6, 1)
        out["global"][str(yr)]["production_Mt"] = _round(prod_total / 1e6, 1)
        kg = prod_total * 1000
        m3 = pd.to_numeric(df[year_cols[yr]["WC_Total_m3"]], errors="coerce").sum()
        out["global"][str(yr)]["productivity_kg_m3"] = _round(kg / m3, 2) if m3 else None

    # Trend deltas
    g00, g10, g20 = (out["global"][str(y)] for y in YEARS)
    out["change_2010_2020"] = {
        "green_pct": _round((g20["green_km3"] - g10["green_km3"]) / g10["green_km3"] * 100, 1),
        "blue_pct":  _round((g20["blue_km3"]  - g10["blue_km3"])  / g10["blue_km3"]  * 100, 1),
        "total_pct": _round((g20["total_km3"] - g10["total_km3"]) / g10["total_km3"] * 100, 1),
        "area_pct":  _round((g20["area_Mha"]  - g10["area_Mha"])  / g10["area_Mha"]  * 100, 1),
    }
    out["change_2000_2020"] = {
        "green_pct": _round((g20["green_km3"] - g00["green_km3"]) / g00["green_km3"] * 100, 1),
        "blue_pct":  _round((g20["blue_km3"]  - g00["blue_km3"])  / g00["blue_km3"]  * 100, 1),
        "total_pct": _round((g20["total_km3"] - g00["total_km3"]) / g00["total_km3"] * 100, 1),
        "productivity_pct": _round((g20["productivity_kg_m3"] - g00["productivity_kg_m3"])
                                   / g00["productivity_kg_m3"] * 100, 1)
                            if g00.get("productivity_kg_m3") else None,
    }
    return out


# ─── Per-crop summary with area + group share ────────────────────────────────

def build_crops(summary: dict) -> dict:
    """Per-crop record: name, group, by-year totals + 2020 area + green/blue share."""
    crops = []
    by_crop = summary["by_crop"]
    g2020_total = summary["global"]["2020"]["total_km3"]

    # Pull area for each crop from Crop_summary_2020.csv. Columns:
    #   group, Crop_Code, Crop_Name, WCgn_rf_Mm3, WCgn_irr_Mm3, WC_gn_Mm3, WCbl_Mm3,
    #   of_which_Paddy_Mm3, WC_blgn_Mm3, rainfed_Mha, irrigated_Mha, total_Mha
    cs20 = pd.read_csv(ROOT / "Crop_summary_2020.csv", low_memory=False)
    # Drop unit row
    cs20 = cs20[cs20["Crop_Code"].notna()]
    area_lookup: dict[str, float] = {}
    for _, row in cs20.iterrows():
        code = str(row["Crop_Code"]).strip()
        ha = pd.to_numeric(row.get("total_Mha"), errors="coerce")
        area_lookup[code] = _safe(ha)

    for code, c in by_crop.items():
        rec = dict(c)
        rec["share_2020_pct"] = _round(c["2020"]["total_km3"] / g2020_total * 100, 2) if g2020_total else None
        tot20 = c["2020"]["total_km3"] or 0
        rec["green_share_pct"] = _round((c["2020"]["green_km3"] or 0) / tot20 * 100, 1) if tot20 else None
        rec["blue_share_pct"]  = _round((c["2020"]["blue_km3"]  or 0) / tot20 * 100, 1) if tot20 else None
        # Trend
        t10 = c["2010"]["total_km3"]
        t20 = c["2020"]["total_km3"]
        rec["trend_2010_2020_pct"] = _round((t20 - t10) / t10 * 100, 1) if t10 else None
        # Area
        rec["area_2020_Mha"] = area_lookup.get(code)
        crops.append(rec)

    crops.sort(key=lambda r: r["2020"]["total_km3"] or 0, reverse=True)
    return {"crops": crops, "groups": list(CROP_GROUPS.keys())}


def build_crop_groups(summary: dict) -> dict:
    """Aggregate by crop group for the crop-breakdown grid."""
    g2020_total = summary["global"]["2020"]["total_km3"]
    groups: dict[str, dict] = {g: {"name": g, "id": g.lower().replace(" ", "_").replace("&", "and"),
                                   "green_km3": 0, "blue_km3": 0, "total_km3": 0}
                                for g in CROP_GROUPS}
    for code, c in summary["by_crop"].items():
        g = CROP_GROUP_LOOKUP.get(code, "Other")
        groups[g]["green_km3"] += (c["2020"]["green_km3"] or 0)
        groups[g]["blue_km3"]  += (c["2020"]["blue_km3"]  or 0)
        groups[g]["total_km3"] += (c["2020"]["total_km3"] or 0)

    out = []
    for g, rec in groups.items():
        tot = rec["total_km3"]
        rec["green_pct"] = round(rec["green_km3"] / tot * 100) if tot else 0
        rec["blue_pct"]  = round(rec["blue_km3"]  / tot * 100) if tot else 0
        rec["share_pct"] = _round(tot / g2020_total * 100, 1) if g2020_total else 0
        rec["green_km3"] = _round(rec["green_km3"], 1)
        rec["blue_km3"]  = _round(rec["blue_km3"], 1)
        rec["total_km3"] = _round(rec["total_km3"], 1)
        out.append(rec)
    out.sort(key=lambda r: r["total_km3"] or 0, reverse=True)
    return {"groups": out}


# ─── Country aggregates ───────────────────────────────────────────────────────

def _read_country_wf(year: int) -> pd.DataFrame:
    """Load country WF CSV, drop unit row, coerce numerics."""
    path = ROOT / "CountryCropsCSVs" / f"Y{year}_country_WaterFootPrint.csv"
    df = pd.read_csv(path, low_memory=False)
    # Row 0 contains units in the same shape as columns; drop it
    df = df.iloc[1:].reset_index(drop=True)
    # Coerce numeric columns
    for col in df.columns:
        if col not in ("ISO3", "Country"):
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def _country_totals(df: pd.DataFrame) -> pd.DataFrame:
    """Compute green/blue/total per country by summing across crops."""
    # Green = sum of *_RF_WFgn_m3_yr + *_IRR_WFgn_m3_yr
    # Blue  = sum of *_IRR_WFbl_m3_yr
    gn_cols = [c for c in df.columns if c.endswith("_WFgn_m3_yr")]
    bl_cols = [c for c in df.columns if c.endswith("_WFbl_m3_yr")]
    out = pd.DataFrame({
        "ISO3": df["ISO3"],
        "Country": df["Country"],
        "green_m3": df[gn_cols].sum(axis=1, min_count=1),
        "blue_m3":  df[bl_cols].sum(axis=1, min_count=1),
    })
    out["total_m3"] = out["green_m3"].fillna(0) + out["blue_m3"].fillna(0)
    return out


def build_countries() -> dict:
    """Per-country totals for 2000/2010/2020 + 2020 ranking with green/blue split."""
    totals: dict[int, pd.DataFrame] = {}
    for yr in YEARS:
        df = _read_country_wf(yr)
        totals[yr] = _country_totals(df)

    # Build the master country list off 2020 (most complete)
    master = totals[2020].copy()
    master = master.dropna(subset=["ISO3"])
    master = master[master["total_m3"] > 0]
    master = master.sort_values("total_m3", ascending=False).reset_index(drop=True)

    # Update global ISO3 → name lookup for country detail files
    for _, row in master.iterrows():
        ISO3_TO_NAME[row["ISO3"]] = row["Country"]

    # Compute trend 2010 → 2020 per country
    t10 = totals[2010].set_index("ISO3")["total_m3"]
    t20 = totals[2020].set_index("ISO3")["total_m3"]
    trend = ((t20 - t10) / t10 * 100).replace([float("inf"), -float("inf")], None)

    out_rows = []
    for _, row in master.iterrows():
        iso = row["ISO3"]
        total = row["total_m3"]
        green = row["green_m3"] or 0
        blue  = row["blue_m3"] or 0
        if total <= 0 or not isinstance(iso, str):
            continue
        out_rows.append({
            "iso3": iso,
            "name": clean_country_name(row["Country"]),
            "green_km3": _round(green / 1e9, 2),
            "blue_km3":  _round(blue / 1e9, 2),
            "total_km3": _round(total / 1e9, 2),
            "green_pct": round(green / total * 100) if total else None,
            "blue_pct":  round(blue  / total * 100) if total else None,
            "trend_2010_2020_pct": _round(trend.get(iso), 1),
        })

    return {"countries": out_rows}


def build_country_detail(countries: dict) -> None:
    """One JSON per country with per-crop breakdown + year trend."""
    # Reload all three years' country WF files once
    dfs: dict[int, pd.DataFrame] = {yr: _read_country_wf(yr) for yr in YEARS}

    crop_codes = sorted({c.split("_")[0] for c in dfs[2020].columns
                         if any(c.endswith(suf) for suf in
                                ("_RF_WFgn_m3_yr","_IRR_WFgn_m3_yr","_IRR_WFbl_m3_yr"))})

    for country in countries["countries"]:
        iso = country["iso3"]
        detail = dict(country)
        detail["years"] = {}

        for yr in YEARS:
            df = dfs[yr]
            row = df[df["ISO3"] == iso]
            if row.empty:
                continue
            row = row.iloc[0]
            gn_cols = [c for c in df.columns if c.endswith("_WFgn_m3_yr")]
            bl_cols = [c for c in df.columns if c.endswith("_WFbl_m3_yr")]
            g = pd.to_numeric(row[gn_cols], errors="coerce").sum()
            b = pd.to_numeric(row[bl_cols], errors="coerce").sum()
            detail["years"][str(yr)] = {
                "green_km3": _round(g / 1e9, 2),
                "blue_km3":  _round(b / 1e9, 2),
                "total_km3": _round((g + b) / 1e9, 2),
            }

        # Per-crop breakdown for 2020
        df = dfs[2020]
        row = df[df["ISO3"] == iso]
        crops = []
        if not row.empty:
            row = row.iloc[0]
            for code in crop_codes:
                gn_rf  = pd.to_numeric(row.get(f"{code}_RF_WFgn_m3_yr"),  errors="coerce")
                gn_irr = pd.to_numeric(row.get(f"{code}_IRR_WFgn_m3_yr"), errors="coerce")
                bl     = pd.to_numeric(row.get(f"{code}_IRR_WFbl_m3_yr"), errors="coerce")
                g = (gn_rf if not pd.isna(gn_rf) else 0) + (gn_irr if not pd.isna(gn_irr) else 0)
                b = bl if not pd.isna(bl) else 0
                tot = g + b
                if tot <= 0:
                    continue
                crops.append({
                    "code": code,
                    "green_km3": _round(g / 1e9, 3),
                    "blue_km3":  _round(b / 1e9, 3),
                    "total_km3": _round(tot / 1e9, 3),
                    "green_pct": round(g / tot * 100) if tot else 0,
                    "blue_pct":  round(b / tot * 100) if tot else 0,
                })
            crops.sort(key=lambda r: r["total_km3"] or 0, reverse=True)
        detail["crops"] = crops[:20]  # top 20 crops per country

        with (OUT / "country_detail" / f"{iso}.json").open("w", encoding="utf-8") as f:
            json.dump(detail, f, separators=(",", ":"))


# ─── Meta ────────────────────────────────────────────────────────────────────

def build_meta(summary: dict, crops: dict, countries: dict) -> dict:
    return {
        "version": "1.0.0",
        "release": "2026",
        "license": "CC BY 4.0",
        "n_crops":      len(summary["by_crop"]),
        "n_groups":     len(CROP_GROUPS),
        "n_countries":  len(countries["countries"]),
        "years":        YEARS,
        "resolution":   "0.083° (~10 km) grid",
        "model":        "PCR-GLOBWB + GAEZ crop calendar",
        "doi":          "10.xxxx/cropsgbw.2026",
        "citation":     "Global spatially explicit crop water consumption shows an overall "
                        "increase of 9% for 46 agricultural crops from 2010 to 2020.",
        "built_at":     pd.Timestamp.utcnow().isoformat(timespec="seconds") + "Z",
    }


# ─── Entrypoint ──────────────────────────────────────────────────────────────

def main() -> None:
    print("[1/6] summary.json")
    summary = build_summary()
    (OUT / "summary.json").write_text(json.dumps(summary, separators=(",", ":")))

    print("[2/6] crops.json")
    crops = build_crops(summary)
    (OUT / "crops.json").write_text(json.dumps(crops, separators=(",", ":")))

    print("[3/6] crop_groups.json")
    groups = build_crop_groups(summary)
    (OUT / "crop_groups.json").write_text(json.dumps(groups, separators=(",", ":")))

    print("[4/6] countries.json")
    countries = build_countries()
    (OUT / "countries.json").write_text(json.dumps(countries, separators=(",", ":")))

    print(f"[5/6] country_detail/* ({len(countries['countries'])} files)")
    build_country_detail(countries)

    print("[6/6] meta.json")
    meta = build_meta(summary, crops, countries)
    (OUT / "meta.json").write_text(json.dumps(meta, indent=2))

    # Summary of what was written
    print()
    print(f"  Years:       {summary['years']}")
    for y in YEARS:
        g = summary["global"][str(y)]
        print(f"    {y}: green {g['green_km3']:>8}  blue {g['blue_km3']:>7}  total {g['total_km3']:>8} km3/yr  (area {g['area_Mha']} Mha)")
    print(f"  Crops:       {len(summary['by_crop'])}")
    print(f"  Countries:   {len(countries['countries'])}")
    c = summary['change_2010_2020']
    print(f"  d 2010-2020: total {c['total_pct']}%  green {c['green_pct']}%  blue {c['blue_pct']}%")
    print(f"\nDone. Bundles written to {OUT}")


if __name__ == "__main__":
    main()
