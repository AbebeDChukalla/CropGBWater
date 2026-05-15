"""Aggregate the CropGBWater CSV outputs into the static JSON bundles the
webapp consumes.

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
  webapp/data/continents.json
  webapp/data/country_detail/{ISO3}.json
  webapp/data/meta.json
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT  = ROOT / "webapp" / "data"
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "country_detail").mkdir(parents=True, exist_ok=True)

YEARS = [2000, 2010, 2020]

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
CROP_GROUP_LOOKUP: dict[str, str] = {c: g for g, codes in CROP_GROUPS.items() for c in codes}

# ─── ISO3 → continent ────────────────────────────────────────────────
ISO3_CONTINENT: dict[str, str] = {}
def _seed(continent: str, *codes: str):
    for c in codes:
        ISO3_CONTINENT[c] = continent

_seed("Africa", "DZA","AGO","BEN","BWA","BFA","BDI","CMR","CPV","CAF","TCD","COM","COG","COD","CIV","DJI","EGY","GNQ","ERI","SWZ","ETH","GAB","GMB","GHA","GIN","GNB","KEN","LSO","LBR","LBY","MDG","MWI","MLI","MRT","MUS","MAR","MOZ","NAM","NER","NGA","RWA","STP","SEN","SLE","SOM","ZAF","SSD","SDN","TZA","TGO","TUN","UGA","ZMB","ZWE","ESH","REU","MYT","SHN","SYC")
_seed("Asia", "AFG","ARM","AZE","BHR","BGD","BTN","BRN","KHM","CHN","CYP","GEO","IND","IDN","IRN","IRQ","ISR","JPN","JOR","KAZ","KWT","KGZ","LAO","LBN","MYS","MDV","MNG","MMR","NPL","PRK","OMN","PAK","PSE","PHL","QAT","SAU","SGP","KOR","LKA","SYR","TWN","TJK","THA","TLS","TUR","TKM","ARE","UZB","VNM","YEM","HKG","MAC")
_seed("Europe", "ALB","AND","AUT","BLR","BEL","BIH","BGR","HRV","CZE","DNK","EST","FIN","FRA","DEU","GRC","HUN","ISL","IRL","ITA","LVA","LIE","LTU","LUX","MLT","MDA","MCO","MNE","NLD","MKD","NOR","POL","PRT","ROU","RUS","SMR","SRB","SVK","SVN","ESP","SWE","CHE","UKR","GBR","VAT","FRO","GIB")
_seed("North America", "CAN","USA","MEX","GTM","BLZ","SLV","HND","NIC","CRI","PAN","CUB","DOM","HTI","JAM","BHS","BRB","ATG","TTO","GRD","DMA","LCA","VCT","KNA","PRI","BMU")
_seed("South America", "ARG","BOL","BRA","CHL","COL","ECU","GUY","PRY","PER","SUR","URY","VEN","GUF","FLK")
_seed("Oceania", "AUS","NZL","FJI","PNG","SLB","VUT","WSM","TON","KIR","MHL","FSM","NRU","PLW","TUV","NCL","PYF")

ISO3_TO_NAME: dict[str, str] = {}


def _safe(v):
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
}

def clean_country_name(raw: str) -> str:
    if not isinstance(raw, str):
        return raw
    raw = raw.strip()
    if raw in _NAME_OVERRIDES:
        return _NAME_OVERRIDES[raw]
    return re.sub(r"(?<=[a-z])(?=[A-Z])", " ", raw)


# ─── Wide-yearly CSV reader ──────────────────────────────────────────
def _wide_yearly(path: Path, metric_order: list[str]):
    df = pd.read_csv(path, header=1, low_memory=False)
    df = df[df.iloc[:, 0].notna()]
    cols = list(df.columns)
    cols_per_year = len(metric_order)
    year_cols: dict[int, dict[str, str]] = {}
    for yi, yr in enumerate(YEARS):
        start = 3 + yi * cols_per_year
        block = cols[start:start + cols_per_year]
        year_cols[yr] = dict(zip(metric_order, block))
    return df, year_cols


# ─── Global summary ──────────────────────────────────────────────────
def build_summary() -> dict:
    wc_metrics = ["WCgn_irr_m3", "WCgn_rf_m3", "WC_gn_m3", "WCbl_m3", "WCpaddy_m3", "WC_Total_m3"]
    df, year_cols = _wide_yearly(ROOT / "Cross_year_summary_WC.csv", wc_metrics)
    crop_codes = df.iloc[:, 0].astype(str).tolist()
    crop_names = df.iloc[:, 1].astype(str).tolist()

    out: dict = {"years": YEARS, "global": {}, "by_crop": {}}

    for yr in YEARS:
        gn = pd.to_numeric(df[year_cols[yr]["WC_gn_m3"]],   errors="coerce").sum()
        bl = pd.to_numeric(df[year_cols[yr]["WCbl_m3"]],    errors="coerce").sum()
        tot = pd.to_numeric(df[year_cols[yr]["WC_Total_m3"]], errors="coerce").sum()
        gn_rf  = pd.to_numeric(df[year_cols[yr]["WCgn_rf_m3"]],  errors="coerce").sum()
        gn_irr = pd.to_numeric(df[year_cols[yr]["WCgn_irr_m3"]], errors="coerce").sum()
        out["global"][str(yr)] = {
            "green_km3": _round(gn / 1e9, 1),
            "blue_km3":  _round(bl / 1e9, 1),
            "total_km3": _round(tot / 1e9, 1),
            "green_rainfed_km3":   _round(gn_rf  / 1e9, 1),
            "green_irrigated_km3": _round(gn_irr / 1e9, 1),
        }

    # Per-crop
    for i, code in enumerate(crop_codes):
        if not isinstance(code, str) or code == "nan":
            continue
        crop_row = {"code": code, "name": crop_names[i],
                    "group": CROP_GROUP_LOOKUP.get(code, "Other")}
        for yr in YEARS:
            gn = pd.to_numeric(df[year_cols[yr]["WC_gn_m3"]],   errors="coerce").iloc[i]
            bl = pd.to_numeric(df[year_cols[yr]["WCbl_m3"]],    errors="coerce").iloc[i]
            tot = pd.to_numeric(df[year_cols[yr]["WC_Total_m3"]], errors="coerce").iloc[i]
            gn_rf  = pd.to_numeric(df[year_cols[yr]["WCgn_rf_m3"]],  errors="coerce").iloc[i]
            gn_irr = pd.to_numeric(df[year_cols[yr]["WCgn_irr_m3"]], errors="coerce").iloc[i]
            crop_row[str(yr)] = {
                "green_km3": _round((gn or 0) / 1e9, 2),
                "blue_km3":  _round((bl or 0) / 1e9, 2),
                "total_km3": _round((tot or 0) / 1e9, 2),
                "green_rainfed_km3":   _round((gn_rf  or 0) / 1e9, 2),
                "green_irrigated_km3": _round((gn_irr or 0) / 1e9, 2),
            }
        out["by_crop"][code] = crop_row

    # Area + production + yield
    area_df, area_cols = _wide_yearly(ROOT / "Cross_year_summary_area.csv",
                                       ["RF_ha", "IRR_ha", "Total_ha"])
    prod_df, prod_cols = _wide_yearly(ROOT / "Cross_year_summary_production.csv",
                                       ["production_rf_ton", "production_irr_ton", "production_total_ton"])

    for yr in YEARS:
        rf_ha   = pd.to_numeric(area_df[area_cols[yr]["RF_ha"]],   errors="coerce").sum()
        irr_ha  = pd.to_numeric(area_df[area_cols[yr]["IRR_ha"]],  errors="coerce").sum()
        tot_ha  = pd.to_numeric(area_df[area_cols[yr]["Total_ha"]], errors="coerce").sum()
        rf_t    = pd.to_numeric(prod_df[prod_cols[yr]["production_rf_ton"]],    errors="coerce").sum()
        irr_t   = pd.to_numeric(prod_df[prod_cols[yr]["production_irr_ton"]],   errors="coerce").sum()
        tot_t   = pd.to_numeric(prod_df[prod_cols[yr]["production_total_ton"]], errors="coerce").sum()
        out["global"][str(yr)]["area_Mha"]          = _round(tot_ha / 1e6, 1)
        out["global"][str(yr)]["area_rainfed_Mha"]   = _round(rf_ha  / 1e6, 1)
        out["global"][str(yr)]["area_irrigated_Mha"] = _round(irr_ha / 1e6, 1)
        out["global"][str(yr)]["production_Mt"]       = _round(tot_t / 1e6, 1)
        out["global"][str(yr)]["production_rainfed_Mt"]   = _round(rf_t  / 1e6, 1)
        out["global"][str(yr)]["production_irrigated_Mt"] = _round(irr_t / 1e6, 1)
        out["global"][str(yr)]["yield_ton_ha"] = _round(tot_t / tot_ha, 2) if tot_ha else None
        out["global"][str(yr)]["irrigated_share_area_pct"] = _round(irr_ha / tot_ha * 100, 1) if tot_ha else None
        out["global"][str(yr)]["irrigated_share_prod_pct"] = _round(irr_t  / tot_t  * 100, 1) if tot_t  else None
        kg = tot_t * 1000
        m3 = pd.to_numeric(df[year_cols[yr]["WC_Total_m3"]], errors="coerce").sum()
        out["global"][str(yr)]["productivity_kg_m3"] = _round(kg / m3, 2) if m3 else None

    # Trend deltas
    g00, g10, g20 = (out["global"][str(y)] for y in YEARS)
    def _delta(a, b):
        return _round((b - a) / a * 100, 1) if a else None
    out["change_2010_2020"] = {
        "green_pct": _delta(g10["green_km3"], g20["green_km3"]),
        "blue_pct":  _delta(g10["blue_km3"],  g20["blue_km3"]),
        "total_pct": _delta(g10["total_km3"], g20["total_km3"]),
        "area_pct":  _delta(g10["area_Mha"],  g20["area_Mha"]),
        "production_pct": _delta(g10["production_Mt"], g20["production_Mt"]),
        "yield_pct": _delta(g10["yield_ton_ha"], g20["yield_ton_ha"]),
    }
    out["change_2000_2020"] = {
        "green_pct": _delta(g00["green_km3"], g20["green_km3"]),
        "blue_pct":  _delta(g00["blue_km3"],  g20["blue_km3"]),
        "total_pct": _delta(g00["total_km3"], g20["total_km3"]),
        "productivity_pct": _delta(g00["productivity_kg_m3"], g20["productivity_kg_m3"]),
    }
    return out


# ─── Per-crop with area, production, yield, irr/rf share ─────────────
def build_crops(summary: dict) -> dict:
    crops = []
    by_crop = summary["by_crop"]
    g2020_total = summary["global"]["2020"]["total_km3"]

    # Per-crop area and irr/rf split from Crop_summary_2020.csv
    cs20 = pd.read_csv(ROOT / "Crop_summary_2020.csv", low_memory=False)
    cs20 = cs20[cs20["Crop_Code"].notna()]
    area_lookup: dict[str, dict] = {}
    for _, row in cs20.iterrows():
        code = str(row["Crop_Code"]).strip()
        rf_ha  = _safe(pd.to_numeric(row.get("rainfed_Mha"), errors="coerce"))
        irr_ha = _safe(pd.to_numeric(row.get("irrigated_Mha"), errors="coerce"))
        tot_ha = _safe(pd.to_numeric(row.get("total_Mha"), errors="coerce"))
        area_lookup[code] = {
            "area_2020_Mha":           tot_ha,
            "area_2020_rainfed_Mha":   rf_ha,
            "area_2020_irrigated_Mha": irr_ha,
            "irrigated_share_pct":     round(irr_ha / tot_ha * 100) if (tot_ha and irr_ha is not None) else None,
        }

    # Per-crop production from Cross_year_summary_production.csv
    prod_df, prod_cols = _wide_yearly(
        ROOT / "Cross_year_summary_production.csv",
        ["production_rf_ton", "production_irr_ton", "production_total_ton"],
    )
    prod_lookup: dict[str, dict] = {}
    for _, row in prod_df.iterrows():
        code = str(row.iloc[0]).strip()
        per_year = {}
        for yr in YEARS:
            tot = _safe(pd.to_numeric(row[prod_cols[yr]["production_total_ton"]], errors="coerce"))
            rf  = _safe(pd.to_numeric(row[prod_cols[yr]["production_rf_ton"]],    errors="coerce"))
            irr = _safe(pd.to_numeric(row[prod_cols[yr]["production_irr_ton"]],   errors="coerce"))
            per_year[str(yr)] = {
                "production_Mt":           _round(tot / 1e6, 2) if tot else None,
                "production_rainfed_Mt":   _round(rf  / 1e6, 2) if rf  else None,
                "production_irrigated_Mt": _round(irr / 1e6, 2) if irr else None,
            }
        prod_lookup[code] = per_year

    # Per-crop harvested area history
    area_yr_df, area_yr_cols = _wide_yearly(
        ROOT / "Cross_year_summary_area.csv", ["RF_ha", "IRR_ha", "Total_ha"]
    )
    area_yr_lookup: dict[str, dict] = {}
    for _, row in area_yr_df.iterrows():
        code = str(row.iloc[0]).strip()
        per_year = {}
        for yr in YEARS:
            tot = _safe(pd.to_numeric(row[area_yr_cols[yr]["Total_ha"]], errors="coerce"))
            rf  = _safe(pd.to_numeric(row[area_yr_cols[yr]["RF_ha"]],    errors="coerce"))
            irr = _safe(pd.to_numeric(row[area_yr_cols[yr]["IRR_ha"]],   errors="coerce"))
            per_year[str(yr)] = {
                "area_Mha":           _round(tot / 1e6, 2) if tot else None,
                "area_rainfed_Mha":   _round(rf  / 1e6, 2) if rf  else None,
                "area_irrigated_Mha": _round(irr / 1e6, 2) if irr else None,
            }
        area_yr_lookup[code] = per_year

    for code, c in by_crop.items():
        rec = dict(c)
        rec["share_2020_pct"] = _round(c["2020"]["total_km3"] / g2020_total * 100, 2) if g2020_total else None
        tot20 = c["2020"]["total_km3"] or 0
        rec["green_share_pct"] = _round((c["2020"]["green_km3"] or 0) / tot20 * 100, 1) if tot20 else None
        rec["blue_share_pct"]  = _round((c["2020"]["blue_km3"]  or 0) / tot20 * 100, 1) if tot20 else None
        t10 = c["2010"]["total_km3"]
        t20 = c["2020"]["total_km3"]
        rec["trend_2010_2020_pct"] = _round((t20 - t10) / t10 * 100, 1) if t10 else None

        # area / production / yield
        ar = area_lookup.get(code, {})
        rec["area_2020_Mha"]           = ar.get("area_2020_Mha")
        rec["area_2020_rainfed_Mha"]   = ar.get("area_2020_rainfed_Mha")
        rec["area_2020_irrigated_Mha"] = ar.get("area_2020_irrigated_Mha")
        rec["irrigated_share_pct"]     = ar.get("irrigated_share_pct")
        rec["production_history"] = prod_lookup.get(code, {})
        rec["area_history"]       = area_yr_lookup.get(code, {})

        prod_2020 = rec["production_history"].get("2020", {})
        prod_t = prod_2020.get("production_Mt")
        area_ha = rec["area_2020_Mha"]
        rec["production_2020_Mt"] = prod_t
        # yield: t / (Mha * 1e6 ha/Mha) → t/ha
        if prod_t and area_ha:
            rec["yield_2020_ton_ha"] = _round((prod_t * 1e6) / (area_ha * 1e6), 2)
        else:
            rec["yield_2020_ton_ha"] = None

        crops.append(rec)

    crops.sort(key=lambda r: r["2020"]["total_km3"] or 0, reverse=True)
    return {"crops": crops, "groups": list(CROP_GROUPS.keys())}


def build_crop_groups(summary: dict) -> dict:
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


# ─── Country tables ──────────────────────────────────────────────────
def _read_country_wf(year: int) -> pd.DataFrame:
    path = ROOT / "CountryCropsCSVs" / f"Y{year}_country_WaterFootPrint.csv"
    df = pd.read_csv(path, low_memory=False)
    df = df.iloc[1:].reset_index(drop=True)
    for col in df.columns:
        if col not in ("ISO3", "Country"):
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def _read_country_area(year: int) -> pd.DataFrame | None:
    p = ROOT / "CountryCropsCSVs" / f"Y{year}_country_HarvestedArea.csv"
    if not p.exists(): return None
    df = pd.read_csv(p, low_memory=False)
    df = df.iloc[1:].reset_index(drop=True)
    for col in df.columns:
        if col not in ("ISO3", "Country"):
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def _country_totals(df: pd.DataFrame) -> pd.DataFrame:
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


def _country_area_totals(df: pd.DataFrame) -> pd.DataFrame:
    """Sum harvested area across all crops, splitting RF vs IRR by column suffix."""
    rf_cols  = [c for c in df.columns if c.endswith("_RF_ha")]
    irr_cols = [c for c in df.columns if c.endswith("_IRR_ha")]
    out = pd.DataFrame({
        "ISO3": df["ISO3"],
        "area_rf_ha":  df[rf_cols].sum(axis=1, min_count=1) if rf_cols else 0,
        "area_irr_ha": df[irr_cols].sum(axis=1, min_count=1) if irr_cols else 0,
    })
    out["area_ha"] = out["area_rf_ha"].fillna(0) + out["area_irr_ha"].fillna(0)
    return out


def build_countries() -> dict:
    totals: dict[int, pd.DataFrame] = {yr: _country_totals(_read_country_wf(yr)) for yr in YEARS}
    area_2020 = _read_country_area(2020)
    area_tots = _country_area_totals(area_2020) if area_2020 is not None else None

    master = totals[2020].copy()
    master = master.dropna(subset=["ISO3"])
    master = master[master["total_m3"] > 0]
    master = master.sort_values("total_m3", ascending=False).reset_index(drop=True)

    for _, row in master.iterrows():
        ISO3_TO_NAME[row["ISO3"]] = row["Country"]

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
        rec = {
            "iso3": iso,
            "name": clean_country_name(row["Country"]),
            "continent": ISO3_CONTINENT.get(iso, "Other"),
            "green_km3": _round(green / 1e9, 2),
            "blue_km3":  _round(blue / 1e9, 2),
            "total_km3": _round(total / 1e9, 2),
            "green_pct": round(green / total * 100) if total else None,
            "blue_pct":  round(blue  / total * 100) if total else None,
            "trend_2010_2020_pct": _round(trend.get(iso), 1),
        }
        if area_tots is not None:
            a = area_tots[area_tots["ISO3"] == iso]
            if not a.empty:
                ha   = _safe(a.iloc[0]["area_ha"])
                rf   = _safe(a.iloc[0]["area_rf_ha"])
                irr  = _safe(a.iloc[0]["area_irr_ha"])
                rec["area_2020_Mha"]           = _round((ha or 0) / 1e6, 2)
                rec["area_2020_rainfed_Mha"]   = _round((rf or 0) / 1e6, 2)
                rec["area_2020_irrigated_Mha"] = _round((irr or 0) / 1e6, 2)
                rec["irrigated_share_pct"] = round(irr / ha * 100) if (ha and irr is not None and ha > 0) else None
        out_rows.append(rec)

    return {"countries": out_rows}


def build_continents(countries: dict) -> dict:
    """Continent-level rollup from the country bundle."""
    agg: dict[str, dict] = {}
    for c in countries["countries"]:
        k = c["continent"]
        bucket = agg.setdefault(k, {
            "name": k, "n_countries": 0,
            "green_km3": 0, "blue_km3": 0, "total_km3": 0,
            "area_Mha": 0, "area_rainfed_Mha": 0, "area_irrigated_Mha": 0,
        })
        bucket["n_countries"] += 1
        for fld in ("green_km3", "blue_km3", "total_km3",
                    "area_Mha", "area_rainfed_Mha", "area_irrigated_Mha"):
            bucket[fld] += (c.get(fld) or 0)

    out = []
    grand_total = sum(b["total_km3"] for b in agg.values()) or 1
    for k, b in agg.items():
        tot = b["total_km3"]
        b["green_pct"] = round(b["green_km3"] / tot * 100) if tot else 0
        b["blue_pct"]  = round(b["blue_km3"]  / tot * 100) if tot else 0
        b["share_pct"] = _round(tot / grand_total * 100, 1)
        a = b["area_Mha"]
        b["irrigated_share_pct"] = round(b["area_irrigated_Mha"] / a * 100) if a else None
        for fld in ("green_km3", "blue_km3", "total_km3",
                    "area_Mha", "area_rainfed_Mha", "area_irrigated_Mha"):
            b[fld] = _round(b[fld], 1)
        out.append(b)
    out.sort(key=lambda r: r["total_km3"] or 0, reverse=True)
    return {"continents": out}


def build_country_detail(countries: dict) -> None:
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
        detail["crops"] = crops[:20]

        with (OUT / "country_detail" / f"{iso}.json").open("w", encoding="utf-8") as f:
            json.dump(detail, f, separators=(",", ":"))


def build_meta(summary: dict, crops: dict, countries: dict) -> dict:
    return {
        "version": "1.1.0",
        "release": "2026",
        "license": "CC BY 4.0",
        "n_crops":      len(summary["by_crop"]),
        "n_groups":     len(CROP_GROUPS),
        "n_countries":  len(countries["countries"]),
        "n_continents": 6,
        "years":        YEARS,
        "resolution":   "0.083° (~10 km) grid",
        "model":        "PCR-GLOBWB + GAEZ crop calendar",
        # Paper
        "paper_authors":   "Chukalla, A.D., Mekonnen, M.M., Gunathilake, D., Wolkeba, F.T., Gunasekara, B., Vanham, D.",
        "paper_year":      2025,
        "paper_title":     "Global spatially explicit crop water consumption shows an overall increase of 9% for 46 agricultural crops from 2010 to 2020",
        "paper_journal":   "Nature Food",
        "paper_volume":    "Volume 6",
        "paper_doi":       "10.1038/s43016-025-01231-x",
        "paper_url":       "https://doi.org/10.1038/s43016-025-01231-x",
        # Dataset (Zenodo)
        "data_doi":        "10.5281/zenodo.17059989",
        "data_url":        "https://doi.org/10.5281/zenodo.17059989",
        "data_title":      "Global spatially explicit crop water consumption shows an overall increase of 9% for 46 agricultural crops from 2010 to 2020: Data and software",
        # Headline citation string used in the hero source line
        "citation":     "Chukalla et al. (2025), Nature Food. doi.org/10.1038/s43016-025-01231-x",
        "built_at":     pd.Timestamp.utcnow().isoformat(timespec="seconds") + "Z",
    }


def main() -> None:
    print("[1/7] summary.json")
    summary = build_summary()
    (OUT / "summary.json").write_text(json.dumps(summary, separators=(",", ":")))

    print("[2/7] crops.json")
    crops = build_crops(summary)
    (OUT / "crops.json").write_text(json.dumps(crops, separators=(",", ":")))

    print("[3/7] crop_groups.json")
    groups = build_crop_groups(summary)
    (OUT / "crop_groups.json").write_text(json.dumps(groups, separators=(",", ":")))

    print("[4/7] countries.json")
    countries = build_countries()
    (OUT / "countries.json").write_text(json.dumps(countries, separators=(",", ":")))

    print("[5/7] continents.json")
    continents = build_continents(countries)
    (OUT / "continents.json").write_text(json.dumps(continents, separators=(",", ":")))

    print(f"[6/7] country_detail/* ({len(countries['countries'])} files)")
    build_country_detail(countries)

    print("[7/7] meta.json")
    meta = build_meta(summary, crops, countries)
    (OUT / "meta.json").write_text(json.dumps(meta, indent=2))

    print()
    print(f"  Years:       {summary['years']}")
    for y in YEARS:
        g = summary["global"][str(y)]
        print(f"    {y}: WC green {g['green_km3']:>7}  blue {g['blue_km3']:>6}  total {g['total_km3']:>7} km3  area {g['area_Mha']} Mha  prod {g['production_Mt']} Mt  yield {g['yield_ton_ha']} t/ha  irr-area {g['irrigated_share_area_pct']}%")
    print(f"  Crops:       {len(summary['by_crop'])}")
    print(f"  Countries:   {len(countries['countries'])}")
    print(f"  Continents:  {len(continents['continents'])}")
    c = summary['change_2010_2020']
    print(f"  d 2010-2020: total {c['total_pct']}%  green {c['green_pct']}%  blue {c['blue_pct']}%  area {c['area_pct']}%  prod {c['production_pct']}%")
    print(f"\nDone. Bundles written to {OUT}")


if __name__ == "__main__":
    main()
