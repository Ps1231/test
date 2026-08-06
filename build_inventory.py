"""
=========================================================
BAH 2026 — PS-6 :: data_inventory.csv builder
=========================================================
FIRST job, before any export.

Answers: what scenes actually exist over Mandya,
for which dates, at what cloud cover. You cannot plan a
compositing window without knowing whether that window has
any usable optical data in it.

Output: data_inventory.csv
        gap_statistics.csv
=========================================================
"""

import csv
from pathlib import Path

import ee

from gee_config import (SEASONS, COLLECTIONS, S2_CLOUD_MAX, LANDSAT_CLOUD_MAX,
                        S1_ORBIT_PASS, S1_RELATIVE_ORBIT)
import gee_collections as gc
from gee_export import make_windows


OUT_DIR = Path(__file__).resolve().parent / "inventory"
OUT_DIR.mkdir(exist_ok=True)


def _scene_rows(col, sensor, props):
    """Pull per-scene metadata into plain dicts."""
    def extract(img):
        d = {"id": img.get("system:index"),
             "date": ee.Date(img.get("system:time_start"))
                       .format("YYYY-MM-dd")}
        for p in props:
            d[p] = img.get(p)
        return ee.Feature(None, d)

    feats = col.map(extract).getInfo()["features"]

    rows = []
    for f in feats:
        p = f["properties"]
        p["sensor"] = sensor
        rows.append(p)
    return rows


def inventory_season(aoi, season):
    """Collect scene-level metadata for every sensor in one season."""
    start, end = SEASONS[season]
    rows = []

    print(f"\n{season}  ({start} -> {end})")

    # --- Sentinel-2 ---
    s2 = (ee.ImageCollection(COLLECTIONS["sentinel2"])
          .filterBounds(aoi).filterDate(start, end)
          .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", S2_CLOUD_MAX)))
    n = s2.size().getInfo()
    print(f"  Sentinel-2      : {n:>4} scenes")
    if n:
        rows += _scene_rows(s2, "Sentinel-2",
                            ["CLOUDY_PIXEL_PERCENTAGE", "MGRS_TILE"])

    # --- Sentinel-1 ---
    s1 = (ee.ImageCollection(COLLECTIONS["sentinel1"])
          .filterBounds(aoi).filterDate(start, end)
          .filter(ee.Filter.eq("instrumentMode", "IW"))
          .filter(ee.Filter.eq("orbitProperties_pass", S1_ORBIT_PASS))
          .filter(ee.Filter.eq("relativeOrbitNumber_start", S1_RELATIVE_ORBIT)))
    n = s1.size().getInfo()
    print(f"  Sentinel-1      : {n:>4} scenes")
    if n:
        rows += _scene_rows(s1, "Sentinel-1",
                            ["orbitProperties_pass",
                             "relativeOrbitNumber_start"])

    # --- Landsat 8/9 ---
    for key, label in [("landsat8", "Landsat-8"), ("landsat9", "Landsat-9")]:
        ls = (ee.ImageCollection(COLLECTIONS[key])
              .filterBounds(aoi).filterDate(start, end)
              .filter(ee.Filter.lt("CLOUD_COVER", LANDSAT_CLOUD_MAX)))
        n = ls.size().getInfo()
        print(f"  {label}       : {n:>4} scenes")
        if n:
            rows += _scene_rows(ls, label, ["CLOUD_COVER", "WRS_PATH", "WRS_ROW"])

    # --- MODIS LST ---
    lst = (ee.ImageCollection(COLLECTIONS["modis_lst"])
           .filterBounds(aoi).filterDate(start, end))
    print(f"  MODIS LST       : {lst.size().getInfo():>4} composites")

    # --- VIIRS ---
    vi = (ee.ImageCollection(COLLECTIONS["viirs_ndvi"])
          .filterBounds(aoi).filterDate(start, end))
    print(f"  VIIRS NDVI      : {vi.size().getInfo():>4} composites")

    for r in rows:
        r["season"] = season

    return rows


def gap_analysis(aoi, season):
    """
    Per-window optical availability.

    This is the number that decides whether SAR-informed gap
    filling is a nice-to-have or the whole ballgame. Kharif
    windows with zero clear optical scenes are expected —
    that is precisely the problem the fusion is solving.
    """
    start, end = SEASONS[season]
    windows = make_windows(start, end)

    s2 = (ee.ImageCollection(COLLECTIONS["sentinel2"])
          .filterBounds(aoi)
          .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", S2_CLOUD_MAX)))
    s1 = (ee.ImageCollection(COLLECTIONS["sentinel1"])
          .filterBounds(aoi)
          .filter(ee.Filter.eq("instrumentMode", "IW"))
          .filter(ee.Filter.eq("orbitProperties_pass", S1_ORBIT_PASS))
          .filter(ee.Filter.eq("relativeOrbitNumber_start", S1_RELATIVE_ORBIT)))

    rows = []
    print(f"\n  window   S2   S1   status")
    print("  " + "-" * 34)

    for t0, t1, idx in windows:
        n_s2 = s2.filterDate(t0, t1).size().getInfo()
        n_s1 = s1.filterDate(t0, t1).size().getInfo()

        status = ("OK" if n_s2 > 0 else
                  "SAR-only" if n_s1 > 0 else "EMPTY")

        print(f"  t{idx+1:02d} {t0[5:]}  {n_s2:>3}  {n_s1:>3}   {status}")

        rows.append({
            "season": season, "window": idx + 1,
            "start": t0, "end": t1,
            "n_sentinel2": n_s2, "n_sentinel1": n_s1,
            "status": status,
        })

    n_empty = sum(r["status"] == "EMPTY" for r in rows)
    n_sar = sum(r["status"] == "SAR-only" for r in rows)
    print(f"\n  {len(rows)} windows | {n_sar} SAR-only | {n_empty} empty")

    return rows


def write_csv(rows, path):
    if not rows:
        return
    keys = sorted({k for r in rows for k in r})
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {path.name}  ({len(rows)} rows)")


def main(seasons=None):
    gc.initialize()
    aoi = gc.load_aoi()

    seasons = seasons or ["rabi_2023_24"]

    all_scenes, all_gaps = [], []
    for s in seasons:
        all_scenes += inventory_season(aoi, s)
        all_gaps += gap_analysis(aoi, s)

    print()
    tag = "_".join(seasons)
    write_csv(all_scenes, OUT_DIR / f"data_inventory_{tag}.csv")
    write_csv(all_gaps, OUT_DIR / f"gap_statistics_{tag}.csv")


if __name__ == "__main__":
    main()
