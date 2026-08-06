#!/usr/bin/env python3
"""
Soil Shapefile Analyzer
=======================
Analyzes district-level soil mapping-unit shapefiles (as produced by
KSRSAC/RRSC for the Kharif Paddy crop-simulation model) and generates:
  - a text/markdown report of soil taxonomy, texture and physico-chemical stats
  - a multi-panel PNG dashboard (choropleth maps + distributions)
  - a CSV export of the attribute table

Usage:
    python3 soil_analysis.py Mandya [Mysuru ...]
If no districts are given, it auto-detects every complete <name>.shp/.dbf/.shx
set in the current directory.
"""

import os
import sys
import glob
import warnings
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

warnings.filterwarnings("ignore")

# CRS of the dataset (WGS84 / UTM Zone 43N) per the project .prj / metadata
TARGET_EPSG = 32643

# Attribute groups
NUMERIC = {
    "Sand_Perce": "Sand %", "Silt_Perce": "Silt %", "Clay_Perce": "Clay %",
    "pH": "pH", "EC_dcpermi": "EC (dS/m)", "Organic_Ca": "Organic C (%)",
    "Cation_Exc": "CEC", "Base_Satur": "Base Saturation %",
    "Excha_Calc": "Exch. Ca", "Excha_Magn": "Exch. Mg",
    "Excha_Sodi": "Exch. Na", "Excha_Pota": "Exch. K",
    "Min_Depth_": "Min depth", "Max_Depth_": "Max depth",
}
CATEGORICAL = {
    "Texture": "Surface texture class", "Orders": "Soil order (USDA)",
    "Sub_Orders": "Sub-order", "Groups": "Great group",
    "Sub_Groups": "Sub-group", "Soil_Serie": "Soil series",
    "Soil_Depth": "Soil depth class", "Family_Tex": "Family texture",
    "Family_Min": "Family mineralogy",
}


def load_district(name):
    """Load one district shapefile, fix CRS + invalid geometry."""
    gdf = gpd.read_file(f"{name}.shp")
    if gdf.crs is None:
        gdf.set_crs(epsg=TARGET_EPSG, inplace=True, allow_override=True)
    # repair invalid geometries (dissolve artifacts)
    gdf["geometry"] = gdf.geometry.buffer(0)
    gdf["area_km2"] = gdf.geometry.area / 1e6
    return gdf


def summarize(gdf, name):
    """Build a markdown report string for one district."""
    L = []
    A = L.append
    A(f"## {name} District — Soil Analysis\n")
    A(f"- **Mapping units (polygons):** {len(gdf)}")
    A(f"- **Attributes:** {len([c for c in gdf.columns if c != 'geometry'])}")
    A(f"- **CRS:** EPSG:{TARGET_EPSG} (WGS84 / UTM 43N)")
    A(f"- **Mapped soil area (sum of units):** {gdf['area_km2'].sum():,.0f} km²")
    A(f"  *(units overlap by soil type; not a simple district partition)*\n")

    # Taxonomy
    A("### Soil Taxonomy (USDA Orders)\n")
    A("| Order | Units | Share |")
    A("|---|---|---|")
    vc = gdf["Orders"].value_counts()
    for k, v in vc.items():
        A(f"| {k} | {v} | {v/len(gdf)*100:.1f}% |")
    A("")

    # Texture
    A("### Surface Texture Classes\n")
    A("| Texture | Units | Share |")
    A("|---|---|---|")
    vc = gdf["Texture"].value_counts()
    for k, v in vc.items():
        A(f"| {k} | {v} | {v/len(gdf)*100:.1f}% |")
    A("")

    # Numeric summary
    A("### Physico-chemical Properties\n")
    A("| Property | Min | Mean | Median | Max | Std |")
    A("|---|---|---|---|---|---|")
    for col, lbl in NUMERIC.items():
        if col in gdf.columns:
            s = pd.to_numeric(gdf[col], errors="coerce")
            A(f"| {lbl} | {s.min():.2f} | {s.mean():.2f} | "
              f"{s.median():.2f} | {s.max():.2f} | {s.std():.2f} |")
    A("")

    # Interpretive notes
    A("### Interpretive Notes\n")
    ph = pd.to_numeric(gdf["pH"], errors="coerce")
    A(f"- **Reaction:** pH ranges {ph.min():.1f}–{ph.max():.1f} "
      f"(mean {ph.mean():.1f}); "
      f"{(ph<6.5).mean()*100:.0f}% of units are acidic (<6.5), "
      f"{((ph>=6.5)&(ph<=7.5)).mean()*100:.0f}% near-neutral, "
      f"{(ph>7.5).mean()*100:.0f}% alkaline (>7.5).")
    ec = pd.to_numeric(gdf["EC_dcpermi"], errors="coerce")
    A(f"- **Salinity:** EC mean {ec.mean():.2f} dS/m; "
      f"{(ec>1).mean()*100:.0f}% of units exceed 1 dS/m.")
    oc = pd.to_numeric(gdf["Organic_Ca"], errors="coerce")
    A(f"- **Organic carbon:** mean {oc.mean():.2f}% "
      f"({'generally low' if oc.mean()<0.5 else 'moderate'}).")
    clay = pd.to_numeric(gdf["Clay_Perce"], errors="coerce")
    A(f"- **Clay content:** mean {clay.mean():.0f}% "
      f"(range {clay.min():.0f}–{clay.max():.0f}%).")
    top_series = gdf["Soil_Serie"].value_counts().head(5)
    A(f"- **Dominant soil series:** " +
      ", ".join(f"{k} ({v})" for k, v in top_series.items()) + ".")
    A("")
    return "\n".join(L)


def make_dashboard(gdf, name, outpath):
    """Multi-panel visualization for one district."""
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.25)
    fig.suptitle(f"{name} District — Soil Mapping Units Dashboard",
                 fontsize=18, fontweight="bold", y=0.98)

    def choropleth(ax, col, title, cmap, categorical=False):
        gdf.plot(column=col, ax=ax, cmap=cmap, legend=True,
                 categorical=categorical, edgecolor="white", linewidth=0.15,
                 legend_kwds={"fontsize": 6, "loc": "lower left"}
                 if categorical else {"shrink": 0.6})
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_axis_off()

    # Map panels
    choropleth(fig.add_subplot(gs[0, 0]), "Orders",
               "Soil Orders (USDA Taxonomy)", "tab10", categorical=True)
    choropleth(fig.add_subplot(gs[0, 1]), "Texture",
               "Surface Texture Class", "Set3", categorical=True)
    choropleth(fig.add_subplot(gs[0, 2]), "pH", "Soil pH", "RdYlGn_r")
    choropleth(fig.add_subplot(gs[1, 0]), "Clay_Perce", "Clay %", "YlOrBr")
    choropleth(fig.add_subplot(gs[1, 1]), "Organic_Ca",
               "Organic Carbon %", "Greens")
    choropleth(fig.add_subplot(gs[1, 2]), "Cation_Exc",
               "Cation Exchange Capacity", "PuBu")

    # Distribution panels
    ax = fig.add_subplot(gs[2, 0])
    gdf["Orders"].value_counts().plot.barh(ax=ax, color="steelblue")
    ax.set_title("Soil Order Frequency", fontsize=11, fontweight="bold")
    ax.invert_yaxis()

    ax = fig.add_subplot(gs[2, 1])
    sand = pd.to_numeric(gdf["Sand_Perce"], errors="coerce")
    silt = pd.to_numeric(gdf["Silt_Perce"], errors="coerce")
    clay = pd.to_numeric(gdf["Clay_Perce"], errors="coerce")
    ax.hist([sand, silt, clay], bins=15,
            label=["Sand", "Silt", "Clay"],
            color=["#d9a066", "#8fb339", "#a0522d"])
    ax.set_title("Particle-size Distribution", fontsize=11, fontweight="bold")
    ax.set_xlabel("% content"); ax.legend(fontsize=8)

    ax = fig.add_subplot(gs[2, 2])
    ph = pd.to_numeric(gdf["pH"], errors="coerce")
    ax.hist(ph, bins=20, color="indianred", edgecolor="white")
    ax.axvspan(6.5, 7.5, alpha=0.15, color="green")
    ax.set_title("pH Distribution (neutral band shaded)",
                 fontsize=11, fontweight="bold")
    ax.set_xlabel("pH")

    fig.savefig(outpath, dpi=130, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    args = sys.argv[1:]
    if not args:
        args = sorted({os.path.splitext(f)[0] for f in glob.glob("*.shp")
                       if all(os.path.exists(os.path.splitext(f)[0] + e)
                              for e in (".shp", ".dbf", ".shx"))})
    if not args:
        print("No complete shapefile sets found."); return

    outdir = "./outputs"
    os.makedirs(outdir, exist_ok=True)
    report = ["# Soil Mapping-Unit Analysis Report\n",
              f"_Districts analyzed: {', '.join(args)}_\n"]

    for name in args:
        if not all(os.path.exists(f"{name}{e}") for e in (".shp", ".dbf", ".shx")):
            print(f"[skip] {name}: incomplete file set")
            report.append(f"## {name}\n_Incomplete file set — skipped._\n")
            continue
        print(f"[proc] {name} ...")
        gdf = load_district(name)
        report.append(summarize(gdf, name))
        img = os.path.join(outdir, f"{name}_soil_dashboard.png")
        make_dashboard(gdf, name, img)
        csv = os.path.join(outdir, f"{name}_soil_attributes.csv")
        gdf.drop(columns="geometry").to_csv(csv, index=False)
        print(f"       -> {img}\n       -> {csv}")

    rpt = os.path.join(outdir, "soil_report.md")
    with open(rpt, "w") as f:
        f.write("\n".join(report))
    print(f"[done] report -> {rpt}")


if __name__ == "__main__":
    main()

