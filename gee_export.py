"""
=========================================================
BAH 2026 — PS-6 :: Compositing + Drive export
=========================================================
Turns the filtered collections into fixed-window composites
and pushes them to Google Drive as GeoTIFFs.

Windowing is parameterised (WINDOW_DAYS = 8 or 7) so the
8-day vs IMD-week question can be settled with ISRO without
rewriting anything.
=========================================================
"""

import ee

from gee_config import (
    SEASONS, WINDOW_DAYS, WINDOW_MODE,
    EXPORT_SCALE, EXPORT_CRS, DRIVE_FOLDER, MAX_PIXELS,
    EXPORT_TILE_ROWS, EXPORT_TILE_COLS,
    ASSET_ROOT,
)
import gee_collections as gc


# =========================================================
# WINDOW GENERATION
# =========================================================

def make_windows(start, end, window_days=WINDOW_DAYS):
    """
    Split a date range into fixed-length windows.

    Returns a list of (start_str, end_str, index) tuples.
    window_days=8  -> standard 8-day composites
    window_days=7  -> IMD standard weeks
    """
    from datetime import datetime, timedelta

    d0 = datetime.strptime(start, "%Y-%m-%d")
    d1 = datetime.strptime(end, "%Y-%m-%d")

    windows = []
    i = 0
    cur = d0
    while cur < d1:
        nxt = min(cur + timedelta(days=window_days), d1)
        windows.append((
            cur.strftime("%Y-%m-%d"),
            nxt.strftime("%Y-%m-%d"),
            i,
        ))
        cur = nxt
        i += 1

    return windows


# =========================================================
# COMPOSITING
# =========================================================

def composite_optical(collection, t0, t1):
    """
    Median composite for a window.

    Median (not mean) because it is robust to residual cloud
    that survived the SCL mask.
    """
    return collection.filterDate(t0, t1).median()


def composite_sar(collection, t0, t1):
    """
    Mean composite in LINEAR power, returned in dB.

    Averaging dB directly is mathematically wrong — dB is a
    log scale, so the mean of decibels is not the decibel of
    the mean. Convert, average, convert back.
    """
    sub = collection.filterDate(t0, t1)

    def to_linear(img):
        vv = ee.Image(10).pow(img.select("VV").divide(10))
        vh = ee.Image(10).pow(img.select("VH").divide(10))
        return (
            vv.rename("VV").addBands(vh.rename("VH"))
            .addBands(img.select(["VH_VV", "RVI"]))
        )

    lin_mean = sub.map(to_linear).mean()

    vv_db = lin_mean.select("VV").log10().multiply(10).rename("VV")
    vh_db = lin_mean.select("VH").log10().multiply(10).rename("VH")

    return ee.Image.cat([
        vv_db, vh_db,
        lin_mean.select("VH_VV"),
        lin_mean.select("RVI"),
    ])

def _safe_composite(col, t0, t1, band_names, reducer="median"):
    """
    Composite one window, returning a fully-masked image with the
    correct band names when the window contains no scenes.

    Empty windows are expected in Kharif — monsoon cloud empties the
    optical stream, and Sentinel-1's 12-day repeat leaves some 8-day
    windows with no acquisition at all. A masked placeholder keeps
    band count and band order identical across every window, which
    the downstream [pixels, T, F] cube requires.

    Masked, not zero-filled: a zero would read as a real NDVI of 0
    and corrupt any interpolation across the gap.
    """
    sub = col.filterDate(t0, t1)

    composite = sub.median() if reducer == "median" else sub.mean()

    empty = (ee.Image.constant([0] * len(band_names))
             .rename(band_names)
             .updateMask(ee.Image.constant(0))
             .float())

    return ee.Image(ee.Algorithms.If(sub.size().gt(0), composite, empty))

def build_season_stack(aoi, season, window_days=WINDOW_DAYS,
                       sensors=("s2", "s1", "lst")):
    """
    Build one multi-band image for an entire season.

    Band naming: <VAR>_t<NN>  e.g. NDVI_t01, VH_VV_t07

    This is the object you either export to Drive, or feed
    straight into sampleRegions() once ground truth arrives.
    """
    start, end = SEASONS[season]
    year = int(start[:4])

    windows = make_windows(start, end, window_days)

    s2 = gc.sentinel2(aoi, start, end) if "s2" in sensors else None
    s1 = gc.sentinel1(aoi, start, end) if "s1" in sensors else None
    lst = gc.modis_lst(aoi, start, end) if "lst" in sensors else None

    bands = []

    for t0, t1, idx in windows:
        sfx = f"_t{idx + 1:02d}"

        if s2 is not None:
            keep = ["NDVI", "EVI", "NDWI", "NDRE", "LSWI", "SAVI"]
            opt = _safe_composite(s2, t0, t1, keep, "median")
            bands.append(
                opt.select(keep).rename([b + sfx for b in keep])
            )

        if s1 is not None:
            keep = ["VV", "VH", "VH_VV", "RVI"]
            sar = ee.Image(ee.Algorithms.If(
                s1.filterDate(t0, t1).size().gt(0),
                composite_sar(s1, t0, t1),
                (ee.Image.constant([0] * len(keep))
                 .rename(keep)
                 .updateMask(ee.Image.constant(0))
                 .float()),
            ))
            bands.append(
                sar.select(keep).rename([b + sfx for b in keep])
            )

        if lst is not None:
            th = _safe_composite(lst, t0, t1, ["LST_day"], "mean")
            bands.append(th.select(["LST_day"]).rename(["LST_day" + sfx]))
    stack = ee.Image.cat(bands)

    mask = gc.cropland_mask(aoi, year)
    return stack.updateMask(mask).clip(aoi).toFloat()


# =========================================================
# EXPORT TILING
# =========================================================

def make_export_tiles(aoi, rows=EXPORT_TILE_ROWS, cols=EXPORT_TILE_COLS):
    """
    Split the AOI bounding box into a rows x cols grid.

    The Mandya + Mysuru AOI is ~11,266 km2. A single export
    task over that extent will hit GEE limits or run for
    hours; tiling keeps each task manageable and lets failed
    tiles be retried individually.
    """
    b = aoi.bounds().coordinates().get(0).getInfo()
    xs = [p[0] for p in b]
    ys = [p[1] for p in b]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)

    dx = (xmax - xmin) / cols
    dy = (ymax - ymin) / rows

    tiles = []
    for r in range(rows):
        for c in range(cols):
            rect = ee.Geometry.Rectangle([
                xmin + c * dx, ymin + r * dy,
                xmin + (c + 1) * dx, ymin + (r + 1) * dy,
            ])
            tiles.append((f"r{r}c{c}", rect.intersection(aoi, ee.ErrorMargin(1))))

    return tiles


# =========================================================
# EXPORTERS
# =========================================================

def export_to_drive(image, description, region, scale,
                    folder=DRIVE_FOLDER, crs=EXPORT_CRS):
    """Start one Drive export task and return it."""
    task = ee.batch.Export.image.toDrive(
        image=image,
        description=description,
        folder=folder,
        fileNamePrefix=description,
        region=region,
        scale=scale,
        crs=crs,
        maxPixels=MAX_PIXELS,
        fileFormat="GeoTIFF",
    )
    task.start()
    print(f"  started :: {description}")
    return task


def export_season_stack(aoi, season, sensors=("s2", "s1", "lst"),
                        scale=None, tiled=True):
    """
    Export a whole season's composite stack to Drive.

    Files land as:
      <season>_<mode>_<tile>.tif   e.g. rabi_2023_24_8day_r0c1.tif
    """
    scale = scale or EXPORT_SCALE["sentinel2"]
    stack = build_season_stack(aoi, season, sensors=sensors)

    print(f"\nExporting {season} @ {scale} m ({WINDOW_MODE})")
    # band count is known from config — don't force a getInfo() on the
    # full stack, it materialises all 253 bands just to count them
    n_win = len(make_windows(*SEASONS[season]))
    print(f"  windows: {n_win}  (~{n_win * 11} bands)")

    tasks = []
    if tiled:
        for name, tile in make_export_tiles(aoi):
            desc = f"{season}_{WINDOW_MODE}_{name}"
            tasks.append(export_to_drive(stack, desc, tile, scale))
    else:
        desc = f"{season}_{WINDOW_MODE}_full"
        tasks.append(export_to_drive(stack, desc, aoi, scale))

    return tasks


def export_to_asset(image, asset_id, region, scale, crs=EXPORT_CRS):
    """
    Export to a GEE asset instead of Drive.

    Much better for iteration: the stack stays server-side,
    so sampleRegions() can be re-run repeatedly without
    recomputing composites or downloading anything.
    """
    task = ee.batch.Export.image.toAsset(
        image=image,
        description=asset_id.split("/")[-1],
        assetId=asset_id,
        region=region,
        scale=scale,
        crs=crs,
        maxPixels=MAX_PIXELS,
    )
    task.start()
    print(f"  started asset :: {asset_id}")
    return task


def export_meteorology(aoi, season):
    """
    Export CHIRPS rainfall and ERA5 weather as window sums /
    means. These are coarse (5-11 km) so a single untiled
    export is fine and fast.
    """
    start, end = SEASONS[season]
    windows = make_windows(start, end)

    rain = gc.chirps(aoi, start, end)
    met = gc.era5(aoi, start, end)

    rain_bands, met_bands = [], []
    for t0, t1, idx in windows:
        sfx = f"_t{idx + 1:02d}"
        rain_bands.append(
            rain.filterDate(t0, t1).sum().rename("rain" + sfx)
        )
        met_bands.append(
            met.filterDate(t0, t1).mean()
            .rename([b + sfx for b in met.first().bandNames().getInfo()])
        )

    tasks = [
        export_to_drive(
            ee.Image.cat(rain_bands), f"{season}_chirps_{WINDOW_MODE}",
            aoi, EXPORT_SCALE["chirps"],
        ),
        export_to_drive(
            ee.Image.cat(met_bands), f"{season}_era5_{WINDOW_MODE}",
            aoi, EXPORT_SCALE["weather"],
        ),
    ]
    return tasks


def export_static(aoi):
    """Terrain and cropland mask — one-off, season-independent."""
    return [
        export_to_drive(
            gc.terrain(aoi), "static_terrain", aoi, EXPORT_SCALE["dem"]
        ),
        export_to_drive(
            gc.cropland_mask(aoi, 2024), "static_cropmask",
            aoi, EXPORT_SCALE["lulc"]
        ),
    ]


# =========================================================
# TASK MONITORING
# =========================================================

def monitor(tasks, interval=30):
    """Poll running tasks until all finish."""
    import time

    while True:
        states = [t.status().get("state") for t in tasks]
        done = sum(s in ("COMPLETED", "FAILED", "CANCELLED") for s in states)

        print(f"  {done}/{len(tasks)} finished :: " +
              ", ".join(f"{s}" for s in set(states)))

        if done == len(tasks):
            break
        time.sleep(interval)

    for t in tasks:
        st = t.status()
        if st.get("state") == "FAILED":
            print(f"  FAILED {st.get('description')}: "
                  f"{st.get('error_message')}")
