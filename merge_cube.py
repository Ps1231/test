"""
=========================================================
BAH 2026 — PS-6 :: Track B :: cube merger
=========================================================
Combine the two halves of the Track B feature set into one 15-feature
8-day cube, band-named exactly like the original single-source stack so
the rest of data_prep is unchanged.

Inputs
------
  GEE stack tiles     : <season>_<mode>_r*c*.tif   (from run_ingest.py stack)
        9 bands/window : NDRE, EVI  +  VV,VH,VH_VV,RVI,VV_contrast,VV_corr
                         +  LST_day
        (auto-split sub-tiles like r0c1_s0.tif are matched too)
  LISS-III season tif : <season>_liss3_<mode>.tif  (from liss3_process.py)
        6 bands/window : NDVI, NDWI, LSWI, SAVI, MSAVI, GNDVI

Output
------
  <season>_<mode>_cube.tif : 15 bands/window, ordered per window as
     NDVI, EVI, NDWI, NDRE, LSWI, SAVI, GNDVI, MSAVI,
     VV, VH, VH_VV, RVI, VV_contrast, VV_corr, LST_day
  matching the original BANDS_PER_WINDOW so downstream code is untouched.

How the merge works
-------------------
Both inputs are in EXPORT_CRS (EPSG:32643). The GEE tiles are mosaicked
into one raster; that mosaic's grid becomes the reference. LISS-III bands
are resampled onto it, and per-window bands are interleaved in the
canonical feature order.
"""

import os
import glob
import numpy as np
import rasterio
from rasterio.merge import merge as rio_merge
from rasterio.warp import reproject, Resampling

# canonical per-window feature order (matches original BANDS_PER_WINDOW)
FEATURE_ORDER = [
    "NDVI", "EVI", "NDWI", "NDRE", "LSWI", "SAVI", "GNDVI", "MSAVI",
    "VV", "VH", "VH_VV", "RVI", "VV_contrast", "VV_corr", "LST_day",
]
# which of those come from the LISS-III tif (the rest come from the GEE mosaic)
LISS_FEATURES = {"NDVI", "NDWI", "LSWI", "SAVI", "GNDVI", "MSAVI"}


def _band_index(path):
    """Return {band_description: (path, band_number)} for a raster."""
    out = {}
    with rasterio.open(path) as ds:
        for i in range(1, ds.count + 1):
            desc = ds.descriptions[i - 1]
            if desc:
                out[desc] = (path, i)
    return out


def _mosaic_tiles(tile_paths, out_path):
    """Mosaic GEE export tiles into one raster, preserving band descriptions."""
    srcs = [rasterio.open(p) for p in tile_paths]
    mosaic, transform = rio_merge(srcs)
    prof = srcs[0].profile.copy()
    descs = srcs[0].descriptions
    prof.update(height=mosaic.shape[1], width=mosaic.shape[2],
                transform=transform, count=mosaic.shape[0])
    with rasterio.open(out_path, "w", **prof) as dst:
        dst.write(mosaic)
        for i, d in enumerate(descs, 1):
            if d:
                dst.set_band_description(i, d)
    for s in srcs:
        s.close()
    return out_path


def _to_ref(arr, ds, ref_transform, ref_crs, H, W):
    """Resample one band array onto the reference grid."""
    dst = np.full((H, W), np.nan, np.float32)
    reproject(
        source=arr, destination=dst,
        src_transform=ds.transform, src_crs=ds.crs,
        dst_transform=ref_transform, dst_crs=ref_crs,
        resampling=Resampling.bilinear,
        src_nodata=ds.nodata, dst_nodata=np.nan,
    )
    return dst


def merge_season_cube(season, gee_dir, liss_tif, out_dir, mode="8day"):
    """Build <season>_<mode>_cube.tif from the GEE tiles + LISS-III tif."""
    os.makedirs(out_dir, exist_ok=True)

    # ---- 1. gather + mosaic GEE tiles (matches r*c* AND auto-split _s* tiles) ----
    tiles = sorted(glob.glob(os.path.join(gee_dir, f"{season}_{mode}_r*c*.tif")))
    if not tiles:
        tiles = sorted(glob.glob(os.path.join(gee_dir, f"{season}_{mode}_full.tif")))
    if not tiles:
        raise FileNotFoundError(
            f"No GEE stack tiles for {season}_{mode}_* in {gee_dir}")
    print(f"[{season}] found {len(tiles)} GEE stack tile(s)")

    if len(tiles) == 1:
        gee_mosaic = tiles[0]
    else:
        gee_mosaic = os.path.join(out_dir, f"{season}_{mode}_gee_mosaic.tif")
        _mosaic_tiles(tiles, gee_mosaic)
        print(f"[{season}] mosaicked -> {os.path.basename(gee_mosaic)}")

    gee_bands = _band_index(gee_mosaic)
    liss_bands = _band_index(liss_tif)

    with rasterio.open(gee_mosaic) as ref:
        ref_prof = ref.profile.copy()
        ref_transform = ref.transform
        ref_crs = ref.crs
        H, W = ref.height, ref.width

    # ---- 2. how many windows? (from GEE band names like NDRE_t01) ----
    win_ids = sorted({name.split("_t")[-1] for name in gee_bands if "_t" in name})
    if not win_ids:
        raise ValueError("No windowed band names (…_tNN) found in the GEE mosaic.")
    print(f"[{season}] windows: {len(win_ids)}  ({win_ids[0]}..{win_ids[-1]})")

    # ---- 3. assemble output bands in canonical order ----
    out_arrays, out_names = [], []
    for wid in win_ids:
        for feat in FEATURE_ORDER:
            bname = f"{feat}_t{wid}"
            src_map = liss_bands if feat in LISS_FEATURES else gee_bands
            if bname not in src_map:
                # missing (e.g. empty window) -> fully masked band
                out_arrays.append(np.full((H, W), np.nan, np.float32))
                out_names.append(bname)
                continue
            path, bidx = src_map[bname]
            with rasterio.open(path) as ds:
                arr = ds.read(bidx).astype(np.float32)
                if (ds.width, ds.height, ds.transform) != (W, H, ref_transform):
                    arr = _to_ref(arr, ds, ref_transform, ref_crs, H, W)
            out_arrays.append(arr)
            out_names.append(bname)

    # ---- 4. write cube ----
    out_tif = os.path.join(out_dir, f"{season}_{mode}_cube.tif")
    prof = ref_prof.copy()
    prof.update(count=len(out_arrays), dtype="float32", nodata=np.nan)
    with rasterio.open(out_tif, "w", **prof) as dst:
        for i, (arr, nm) in enumerate(zip(out_arrays, out_names), 1):
            dst.write(arr.astype(np.float32), i)
            dst.set_band_description(i, nm)

    print(f"[{season}] wrote {out_tif}")
    print(f"           {len(out_arrays)} bands = {len(FEATURE_ORDER)} feat "
          f"x {len(win_ids)} windows")
    return out_tif


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 4:
        print("usage: python merge_cube.py <season> <gee_dir> <liss_tif> [out_dir]")
        print('example: python merge_cube.py rabi_2023_24 '
              '"/media/dell/New Volume/PS6/gee_stacks" '
              '"/media/dell/New Volume/PS6/liss3_out/rabi_2023_24_liss3_8day.tif" '
              '"/media/dell/New Volume/PS6/cubes"')
        raise SystemExit(1)
    season, gee_dir, liss_tif = sys.argv[1], sys.argv[2], sys.argv[3]
    out_dir = sys.argv[4] if len(sys.argv) > 4 else "./cubes"
    merge_season_cube(season, gee_dir, liss_tif, out_dir)