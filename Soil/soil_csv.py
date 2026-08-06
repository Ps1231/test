#!/usr/bin/env python3
"""
Per-row soil hydraulic-parameter CSV builder
=============================================
Reads a KSRSAC soil shapefile and, for EVERY polygon/row, computes:
    Field Capacity (FC), Wilting Point (WP), Bulk Density (BD),
    Saturation (SAT), Plant-Available Water (PAW, v/v and mm),
using the Saxton & Rawls (2006) pedotransfer functions from
sand / clay / organic-matter. Writes one CSV row per input row,
keeping the original attributes plus the derived columns.

Usage:
    python3 make_soil_csv.py               # auto-detect *.shp here
    python3 make_soil_csv.py Mandya Mysuru # specific districts
Output: ./output/<District>_soil_parameters.csv
"""
import os, re, sys, glob, warnings
import numpy as np, pandas as pd, geopandas as gpd
warnings.filterwarnings("ignore")

TARGET_EPSG = 32643          # WGS84 / UTM 43N (from .prj)
PARTICLE_DENSITY = 2.65      # g/cm3, mineral soil


def parse_depth(s):
    """'100-150(145)' -> 145.0 ; else class midpoint ; else NaN."""
    if not isinstance(s, str):
        return np.nan
    m = re.search(r"\(\s*([\d.]+)\s*\)", s)
    if m:
        return float(m.group(1))
    nums = re.findall(r"[\d.]+", s)
    if len(nums) >= 2:
        return (float(nums[0]) + float(nums[1])) / 2.0
    return float(nums[0]) if nums else np.nan


def saxton_rawls(sand, clay, om):
    """
    Saxton & Rawls (2006). Vectorized.
    sand, clay in %(0-100); om = organic matter %.
    Returns FC, WP, SAT (v/v) and BD (g/cm3), each a numpy array.
    """
    S = np.clip(np.asarray(sand, float), 0, 100) / 100.0
    C = np.clip(np.asarray(clay, float), 0, 100) / 100.0
    OM = np.clip(np.asarray(om, float), 0, 8)

    # Wilting point (1500 kPa)
    t1500 = (-0.024*S + 0.487*C + 0.006*OM
             + 0.005*S*OM - 0.013*C*OM + 0.068*S*C + 0.031)
    WP = t1500 + (0.14*t1500 - 0.02)

    # Field capacity (33 kPa)
    t33 = (-0.251*S + 0.195*C + 0.011*OM
           + 0.006*S*OM - 0.027*C*OM + 0.452*S*C + 0.299)
    FC = t33 + (1.283*t33**2 - 0.374*t33 - 0.015)

    # Saturation and bulk density
    ts33 = (0.278*S + 0.034*C + 0.022*OM
            - 0.018*S*OM - 0.027*C*OM - 0.584*S*C + 0.078)
    ts33 = ts33 + (0.636*ts33 - 0.107)
    SAT = FC + ts33 - 0.097*S + 0.043
    BD = (1.0 - SAT) * PARTICLE_DENSITY

    WP = np.clip(WP, 0, 1)
    FC = np.clip(FC, 0, 1)
    SAT = np.clip(SAT, 0, 1)
    return FC, WP, SAT, BD


def build(name, outdir):
    gdf = gpd.read_file(f"{name}.shp")
    if gdf.crs is None:
        gdf.set_crs(epsg=TARGET_EPSG, inplace=True, allow_override=True)

    sand = pd.to_numeric(gdf["Sand_Perce"], errors="coerce")
    silt = pd.to_numeric(gdf["Silt_Perce"], errors="coerce")
    clay = pd.to_numeric(gdf["Clay_Perce"], errors="coerce")
    oc = pd.to_numeric(gdf["Organic_Ca"], errors="coerce").fillna(0.0)
    om = oc * 1.724                      # organic C -> organic matter
    depth = gdf["Soil_Depth"].apply(parse_depth)

    FC, WP, SAT, BD = saxton_rawls(sand.values, clay.values, om.values)
    paw = FC - WP

    # polygon centroid (lon/lat) so each row is a locatable point
    cen = gdf.geometry.centroid
    ll = gpd.GeoSeries(cen, crs=gdf.crs).to_crs(4326)

    df = pd.DataFrame({
        "FID": range(len(gdf)),
        "SoilSeries": gdf.get("Soil_Serie"),
        "SoilCode": gdf.get("SoilCode"),
        "Texture": gdf["Texture"],
        "Order": gdf.get("Orders"),
        "SubGroup": gdf.get("Sub_Groups"),
        "Centroid_Lon": ll.x.round(5),
        "Centroid_Lat": ll.y.round(5),
        "Depth_cm": depth.round(0),
        "Sand_pct": sand.round(1),
        "Silt_pct": silt.round(1),
        "Clay_pct": clay.round(1),
        "OrganicC_pct": oc.round(2),
        "pH": pd.to_numeric(gdf.get("pH"), errors="coerce").round(2),
        "EC_dS_m": pd.to_numeric(gdf.get("EC_dcpermi"), errors="coerce").round(2),
        # ---- derived hydraulic parameters ----
        "BulkDensity_g_cm3": np.round(BD, 3),
        "FieldCapacity_vv": np.round(FC, 3),
        "WiltingPoint_vv": np.round(WP, 3),
        "Saturation_vv": np.round(SAT, 3),
        "PAW_vv": np.round(paw, 3),
        "PAW_mm_profile": np.round(paw * depth.values * 10, 1),
    })

    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f"{name}_soil_parameters.csv")
    df.to_csv(path, index=False)
    return df, path


def main():
    args = sys.argv[1:] or sorted(
        {os.path.splitext(f)[0] for f in glob.glob("*.shp")
         if all(os.path.exists(os.path.splitext(f)[0] + e)
                for e in (".shp", ".dbf", ".shx"))})
    if not args:
        print("No complete shapefile sets (.shp/.dbf/.shx) found here."); return

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
    for name in args:
        df, path = build(name, outdir)
        print(f"[{name}] {len(df)} rows written -> {path}")
        print(df[["Texture", "Depth_cm", "BulkDensity_g_cm3",
                  "FieldCapacity_vv", "WiltingPoint_vv", "PAW_vv"]]
              .head(10).to_string(index=False))
        print()


if __name__ == "__main__":
    main()
