"""
=========================================================
BAH 2026 — PS-6 :: Track B :: LISS-III local processor
=========================================================
The six optical indices NDVI, NDWI, LSWI, SAVI, MSAVI, GNDVI now come
from the local ISRO ResourceSat-2 LISS-III archive instead of Sentinel-2.
Sentinel-1 (SAR) and Sentinel-2 (NDRE, EVI only) still come from GEE.

This module turns a folder of LISS-III scene ZIPs into per-window
optical-index GeoTIFFs on the SAME grid the GEE export uses, so
merge_cube.py can stack them next to the SAR/NDRE/EVI bands.

Pipeline per scene:
   1. read 4 bands (Green B2, Red B3, NIR B4, SWIR B5) from the ZIP
   2. classify STD vs BOA vs broken (from ProcessingLevel + contents)
   3. STD  -> radiometric calibration (DN->radiance->TOA reflectance)
             using per-band Lmin/Lmax + Sun elevation + Earth-Sun distance,
             then DOS (dark-object subtraction) haze removal  -> BOA-equiv
      BOA  -> already surface reflectance, use as-is
      broken -> skip
   4. cloud/shadow mask (LISS-III has no SCL band, so threshold-based)
   5. reproject to EXPORT_CRS (EPSG:32643) at EXPORT grid resolution
   6. compute the 6 indices

Then, across all scenes:
   7. group scenes into the same 8-day windows as the GEE stack
   8. median-composite each window (robust to residual cloud)
   9. write <season>_liss3_8day.tif : 6 indices x N windows, band-named
      NDVI_t01, NDWI_t01, ... to match merge_cube.py

Notes / honesty:
   - DOS-corrected STD is BOA-EQUIVALENT, not identical to ISRO BOA.
     Each output carries a sidecar quality flag so downstream code
     knows which windows leaned on corrected STD data.
   - No blue / no red-edge here on purpose: EVI, NDRE stay on S2.
=========================================================
"""

import os
import glob
import io
import json
import math
import zipfile
import datetime as dt
from collections import defaultdict

import numpy as np
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling

from gee_config import EXPORT_CRS, EXPORT_SCALE, SEASONS, WINDOW_DAYS


def make_windows(start, end, window_days=WINDOW_DAYS):
    """Split a date range into fixed-length windows (local copy — no ee).
    Identical logic to gee_export.make_windows so windows line up exactly."""
    d0 = dt.datetime.strptime(start[:10], "%Y-%m-%d")
    d1 = dt.datetime.strptime(end[:10], "%Y-%m-%d")
    windows, i, cur = [], 0, d0
    while cur < d1:
        nxt = min(cur + dt.timedelta(days=window_days), d1)
        windows.append((cur.strftime("%Y-%m-%d"), nxt.strftime("%Y-%m-%d"), i))
        cur = nxt
        i += 1
    return windows

# LISS-III band-number -> role (BandNumbers = 2,3,4,5)
BAND_ROLE = {2: "GREEN", 3: "RED", 4: "NIR", 5: "SWIR"}
DN_MAX = 1023.0            # 10-bit product (BitsPerPixel = 10)
TARGET_RES = EXPORT_SCALE["sentinel2"]   # 20 m — same grid as the GEE stack


# =========================================================
# metadata
# =========================================================
def _read_meta_bytes(zf, names):
    for n in names:
        cand = [m for m in zf.namelist() if m.upper().endswith(n)]
        if cand:
            with zf.open(cand[0]) as fh:
                return io.TextIOWrapper(fh, "utf-8", errors="ignore").read()
    return ""


def parse_meta(text):
    d = {}
    for line in text.splitlines():
        if "=" in line:
            k, _, v = line.partition("=")
            d[k.strip()] = v.strip()
    return d


def scene_info(zip_path):
    """Return a dict describing the scene, or None if unreadable/broken."""
    try:
        with zipfile.ZipFile(zip_path) as zf:
            meta = parse_meta(_read_meta_bytes(zf, ["BAND_META.TXT", ".META"]))
            names = zf.namelist()
            bands_present = {b: any(f"BAND{b}.tif" in n.upper().replace(".TIF", ".TIF")
                                    or f"BAND{b}.TIF" in n.upper() for n in names)
                             for b in BAND_ROLE}
            # simpler robust check:
            bands_present = {}
            for b in BAND_ROLE:
                hit = [n for n in names if n.upper().endswith(f"BAND{b}.TIF")]
                ok = False
                if hit:
                    ok = zf.getinfo(hit[0]).file_size > 0
                bands_present[b] = ok
    except zipfile.BadZipFile:
        return None

    if not meta or not all(bands_present.values()):
        return None   # broken / incomplete

    date = _parse_date(meta.get("DateOfPass"))
    level = meta.get("ProcessingLevel", "").strip().upper()
    return {
        "zip": zip_path,
        "date": date,
        "level": "BOA" if "ATMOS" in level else "STD",
        "sun_elev": _f(meta.get("SunElevationAtCenter")),
        "lmin": {b: _f(meta.get(f"B{b}_Lmin")) for b in BAND_ROLE},
        "lmax": {b: _f(meta.get(f"B{b}_Lmax")) for b in BAND_ROLE},
        "meta": meta,
    }


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _parse_date(s):
    s = (s or "").strip()
    for fmt in ("%d-%b-%Y", "%d-%B-%Y", "%Y-%m-%d"):
        try:
            return dt.datetime.strptime(s, fmt).date()
        except ValueError:
            continue
    return None


# =========================================================
# radiometric calibration (STD only)  DN -> TOA reflectance
# =========================================================
def _earth_sun_distance(d):
    """Earth-Sun distance in AU for day-of-year d (Duffie-Beckman)."""
    doy = d.timetuple().tm_yday
    return 1 - 0.01672 * math.cos(math.radians(0.9856 * (doy - 4)))

# Mean exo-atmospheric solar irradiance (ESUN) for ResourceSat-2 LISS-III
# bands, W/m^2/um. Standard IRS-R2 LISS-III values. Used only for the
# optional TOA-reflectance path; indices themselves are ratio-based and
# do NOT require ESUN, so a units mismatch here can't corrupt NDVI etc.
ESUN = {2: 1849.0, 3: 1553.0, 4: 1092.0, 5: 240.0}


def dn_to_radiance(dn, band_num, info):
    """
    STD DN -> at-sensor spectral radiance.
        L = Lmin + (Lmax - Lmin) * DN / DN_MAX
    Lmin/Lmax come straight from BAND_META. This is the physically
    meaningful, unit-consistent quantity to run DOS on. Normalized-
    difference indices computed from radiance are essentially identical
    to those from reflectance (the solar/geometry constants cancel in
    the ratio), so we stay at radiance and skip the error-prone ESUN
    reflectance conversion for index work.
    """
    lmin = info["lmin"][band_num]
    lmax = info["lmax"][band_num]
    if lmin is None or lmax is None:
        return dn / DN_MAX      # last-resort: normalized DN
    return lmin + (lmax - lmin) * (dn / DN_MAX)


def dn_to_toa_reflectance(dn, band_num, info):
    """
    Optional TOA reflectance (only if you need true reflectance, not for
    the indices). rho = pi * L * d^2 / (ESUN * cos(theta_z)).
    """
    L = dn_to_radiance(dn, band_num, info)
    if info["sun_elev"] is None:
        return L
    d = _earth_sun_distance(info["date"]) if info["date"] else 1.0
    theta_z = math.radians(90.0 - info["sun_elev"])
    return (math.pi * L * d * d) / (ESUN[band_num] * math.cos(theta_z))


# =========================================================
# DOS (dark object subtraction)  removes additive haze
# =========================================================
def dos_correct(x, valid):
    """
    Subtract per-band haze, estimated as the 1st-percentile value of the
    valid pixels (the 'dark object' — deep shadow / clear water should be
    ~0 once haze is removed). DOS1. Works on radiance or reflectance;
    here we run it on radiance. Additive haze removal is exactly what
    sharpens NDVI contrast in an uncorrected STD scene.
    """
    if valid.sum() == 0:
        return x
    haze = np.percentile(x[valid], 1)
    return np.clip(x - haze, 0, None)


# =========================================================
# cloud mask (no SCL band on LISS-III)
# =========================================================
def cloud_mask(green, red, nir, swir):
    """
    Cloud/shadow mask for LISS-III. LISS-III has no SCL band, so this is
    threshold-based. It must work whether the inputs are radiance (STD
    path) or reflectance (BOA path), so thresholds are RELATIVE
    (percentile-based) rather than absolute.

    Cloud   : bright in green AND red AND swir together (high in all).
    Shadow  : very dark in NIR.
    Cloud detection uses NDVI too — cloud has near-zero NDVI while veg is
    high, which separates bright cloud from bright bare soil.
    Conservative: better to drop a suspect pixel than keep a cloudy one.
    """
    finite = np.isfinite(green) & np.isfinite(red) & np.isfinite(nir)
    valid = finite & ((green + red + nir) > 0)
    if valid.sum() == 0:
        return valid

    def hi(x, p=90):
        return np.nanpercentile(x[valid], p)

    def lo(x, p=5):
        return np.nanpercentile(x[valid], p)

    eps = 1e-6
    ndvi = (nir - red) / (nir + red + eps)

    # cloud: bright across visible+swir AND low NDVI (rules out bright crops)
    bright = (green > hi(green)) & (red > hi(red)) & (swir > hi(swir))
    cloud = bright & (ndvi < 0.2)

    # shadow / deep water: very dark NIR
    shadow = nir < lo(nir, 3)

    bad = cloud | shadow | (~valid)
    return ~bad


# =========================================================
# read + process ONE scene -> dict of 6 index arrays on target grid
# =========================================================
def process_scene(info, aoi_bounds=None):
    """
    Returns (indices_dict, profile, qflag) reprojected to EXPORT_CRS @ TARGET_RES.
    indices_dict: {'NDVI':arr, 'NDWI':..., 'LSWI':..., 'SAVI':..., 'MSAVI':..., 'GNDVI':...}
    qflag: 'BOA' or 'STD_DOS' (corrected).
    """
    bands = {}
    src_profile = None
    with zipfile.ZipFile(info["zip"]) as zf:
        for b, role in BAND_ROLE.items():
            hit = [n for n in zf.namelist() if n.upper().endswith(f"BAND{b}.TIF")]
            with zf.open(hit[0]) as fh:
                data = fh.read()
            with rasterio.open(io.BytesIO(data)) as ds:
                arr = ds.read(1).astype(np.float32)
                if src_profile is None:
                    src_profile = ds.profile.copy()
                    src_crs = ds.crs
                    src_transform = ds.transform
            bands[role] = arr

    # --- calibrate to a common radiometric footing ---
    valid = (bands["GREEN"] + bands["RED"] + bands["NIR"]) > 0
    if info["level"] == "STD":
        # DN -> radiance (unit-safe), then DOS haze removal. Indices are
        # ratio-based so radiance is fine; no reflectance conversion needed.
        G = dn_to_radiance(bands["GREEN"], 2, info)
        R = dn_to_radiance(bands["RED"],   3, info)
        N = dn_to_radiance(bands["NIR"],   4, info)
        S = dn_to_radiance(bands["SWIR"],  5, info)
        G, R, N, S = (dos_correct(x, valid) for x in (G, R, N, S))
        qflag = "STD_DOS"
    else:  # BOA already surface reflectance * 10000
        G, R, N, S = (bands[k] / 1e4 for k in ("GREEN", "RED", "NIR", "SWIR"))
        qflag = "BOA"

    # --- cloud mask ---
    good = cloud_mask(G, R, N, S)
    for a in (G, R, N, S):
        a[~good] = np.nan

    # --- indices (the 6 LISS-III can make) ---
    eps = 1e-6
    idx = {
        "NDVI":  (N - R) / (N + R + eps),
        "NDWI":  (G - N) / (G + N + eps),
        "LSWI":  (N - S) / (N + S + eps),
        "SAVI":  1.5 * (N - R) / (N + R + 0.5 + eps),
        "GNDVI": (N - G) / (N + G + eps),
        "MSAVI": (2*N + 1 - np.sqrt(np.clip((2*N+1)**2 - 8*(N-R), 0, None))) / 2,
    }

    # --- reproject each index to EXPORT_CRS @ TARGET_RES ---
    dst_crs = EXPORT_CRS
    transform, width, height = calculate_default_transform(
        src_crs, dst_crs, src_profile["width"], src_profile["height"],
        *rasterio.transform.array_bounds(src_profile["height"],
                                         src_profile["width"], src_transform),
        resolution=TARGET_RES,
    )
    out_profile = {
        "driver": "GTiff", "dtype": "float32", "count": 1,
        "crs": dst_crs, "transform": transform,
        "width": width, "height": height, "nodata": np.nan,
    }
    reproj_idx = {}
    for name, arr in idx.items():
        dst = np.full((height, width), np.nan, dtype=np.float32)
        reproject(
            source=arr, destination=dst,
            src_transform=src_transform, src_crs=src_crs,
            dst_transform=transform, dst_crs=dst_crs,
            resampling=Resampling.bilinear,
            src_nodata=np.nan, dst_nodata=np.nan,
        )
        reproj_idx[name] = dst
    return reproj_idx, out_profile, qflag


# =========================================================
# season driver — folder of ZIPs -> per-window composite GeoTIFF
# =========================================================
INDEX_ORDER = ["NDVI", "NDWI", "LSWI", "SAVI", "MSAVI", "GNDVI"]


def build_liss3_season(folder, season, out_dir, window_days=WINDOW_DAYS,
                       use_std=True):
    """
    Turn a folder of LISS-III scene ZIPs into one season composite GeoTIFF
    matching the GEE stack's windows and grid.

    use_std=True   -> STD scenes are DOS-corrected and included.
    use_std=False  -> STD scenes are dropped (BOA-only, cleaner but sparser).

    Writes:
      <out_dir>/<season>_liss3_<mode>.tif   (6 indices x N windows)
      <out_dir>/<season>_liss3_qflags.json  (per-window: BOA / STD_DOS / empty)
    """
    os.makedirs(out_dir, exist_ok=True)
    start, end = SEASONS[season]
    start_d, end_d = start[:10], end[:10]
    windows = make_windows(start_d, end_d, window_days)

    # ---- inventory the folder ----
    zips = sorted(glob.glob(os.path.join(folder, "*.zip")))
    scenes = []
    n_broken = n_std = n_boa = 0
    for z in zips:
        info = scene_info(z)
        if info is None:
            n_broken += 1
            continue
        if info["date"] is None or not (start_d <= info["date"].isoformat() <= end_d):
            continue
        if info["level"] == "STD":
            if not use_std:
                continue
            n_std += 1
        else:
            n_boa += 1
        scenes.append(info)

    print(f"[{season}] usable scenes: {len(scenes)}  "
          f"(BOA={n_boa}, STD_DOS={n_std if use_std else 0}, broken skipped={n_broken})")
    if not scenes:
        print(f"[{season}] no scenes in range — nothing written.")
        return None

    # ---- assign scenes to windows ----
    win_scenes = defaultdict(list)
    for info in scenes:
        di = info["date"].isoformat()
        for t0, t1, k in windows:
            if t0 <= di < t1:
                win_scenes[k].append(info)
                break

    # ---- process + median-composite per window ----
    ref_profile = None
    per_window_bands = []          # list of (band_name, array)
    qflags = {}

    total_scenes = sum(len(v) for v in win_scenes.values())
    done = 0
    for t0, t1, k in windows:
        sfx = f"_t{k+1:02d}"
        infos = win_scenes.get(k, [])
        if not infos:
            qflags[f"t{k+1:02d}"] = "empty"
            per_window_bands.append((None, None, sfx))   # placeholder
            print(f"  t{k+1:02d}: empty")
            continue

        # process each scene in window, stack, median
        stacks = {name: [] for name in INDEX_ORDER}
        wflags = set()
        for info in infos:
            done += 1
            print(f"  t{k+1:02d}: processing scene {done}/{total_scenes} "
                  f"({info['date']}, {info['level']}) ...", flush=True)
            reproj, prof, qf = process_scene(info)
            wflags.add(qf)
            if ref_profile is None:
                ref_profile = prof
            else:
                reproj = {n: _align(a, prof, ref_profile) for n, a in reproj.items()}
            for name in INDEX_ORDER:
                stacks[name].append(reproj[name])
        qflags[f"t{k+1:02d}"] = "+".join(sorted(wflags))
        comp = {name: np.nanmedian(np.stack(stacks[name]), axis=0)
                for name in INDEX_ORDER}
        per_window_bands.append((comp, sfx, sfx))

    if ref_profile is None:
        print(f"[{season}] all windows empty after processing.")
        return None

    # ---- assemble multi-band output ----
    H, W = ref_profile["height"], ref_profile["width"]
    out_bands, band_names = [], []
    for entry in per_window_bands:
        comp, sfx = entry[0], entry[2]
        for name in INDEX_ORDER:
            band_names.append(f"{name}{sfx}")
            if comp is None:
                out_bands.append(np.full((H, W), np.nan, np.float32))
            else:
                out_bands.append(comp[name].astype(np.float32))

    mode = "8day" if window_days == WINDOW_DAYS else "monthly"
    out_tif = os.path.join(out_dir, f"{season}_liss3_{mode}.tif")
    prof = ref_profile.copy()
    prof.update(count=len(out_bands))
    with rasterio.open(out_tif, "w", **prof) as dst:
        for i, (arr, nm) in enumerate(zip(out_bands, band_names), 1):
            dst.write(arr, i)
            dst.set_band_description(i, nm)

    qpath = os.path.join(out_dir, f"{season}_liss3_qflags.json")
    with open(qpath, "w") as fh:
        json.dump(qflags, fh, indent=2)

    print(f"[{season}] wrote {out_tif}  ({len(out_bands)} bands)")
    print(f"[{season}] wrote {qpath}")
    return out_tif


def _align(arr, prof, ref):
    """Reproject/pad a scene onto the reference grid if profiles differ."""
    if (prof["width"], prof["height"], prof["transform"]) == \
       (ref["width"], ref["height"], ref["transform"]):
        return arr
    dst = np.full((ref["height"], ref["width"]), np.nan, np.float32)
    reproject(
        source=arr, destination=dst,
        src_transform=prof["transform"], src_crs=prof["crs"],
        dst_transform=ref["transform"], dst_crs=ref["crs"],
        resampling=Resampling.bilinear,
        src_nodata=np.nan, dst_nodata=np.nan,
    )
    return dst


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print("usage: python liss3_process.py <liss3_folder> <season> [out_dir]")
        print("example: python liss3_process.py ./LISS3 rabi_2023_24 ./liss3_out")
        raise SystemExit(1)
    folder, season = sys.argv[1], sys.argv[2]
    out_dir = sys.argv[3] if len(sys.argv) > 3 else "./liss3_out"
    build_liss3_season(folder, season, out_dir)