"""
=========================================================
BAH 2026 — PS-6 :: Collection builders
=========================================================
One function per sensor. Each returns a filtered, masked
ee.ImageCollection over the AOI for a given date range.

Nothing here downloads. Composite + export live in
gee_export.py so the same collections can be reused for
sampleRegions() once ISRO ground truth arrives.
=========================================================
"""

import ee

from gee_config import (
    AOI_FILE,
    COLLECTIONS,
    S2_BANDS, S2_SCALE_FACTOR, S2_CLOUD_MAX, S2_SCL_DROP,
    S1_BANDS, S1_ORBIT_PASS, S1_RELATIVE_ORBIT,
    LANDSAT_BANDS, LANDSAT_THERMAL, LANDSAT_MULT, LANDSAT_ADD,
    LANDSAT_CLOUD_MAX,
    ERA5_BANDS,
    DW_CROP_CLASSES, DW_NON_CROP_CLASSES,
    GEE_PROJECT,
)


# =========================================================
# INIT + AOI
# =========================================================

def initialize(project=GEE_PROJECT):
    """
    Initialise Earth Engine.

    Priority:
      1. Service-account JSON string in GEE_SERVICE_ACCOUNT_KEY
      2. Service-account key file at GEE_KEY_PATH
      3. Interactive auth (local dev only)
    """
    import os, json
    from pathlib import Path

    email = os.environ.get("GEE_SERVICE_ACCOUNT_EMAIL", "").strip()

    # 1. JSON string directly in env (deployment)
    key_str = os.environ.get("GEE_SERVICE_ACCOUNT_KEY", "").strip()
    if email and key_str:
        key_dict = json.loads(key_str.strip("'\""))
        cred = ee.ServiceAccountCredentials(
            email, key_data=json.dumps(key_dict)
        )
        ee.Initialize(cred, project=project)
        print(f"Earth Engine :: service account ({email})")
        return

    # 2. Key file path
    key_path = os.environ.get("GEE_KEY_PATH", "").strip()
    if email and key_path:
        p = Path(key_path)
        if not p.is_absolute():
            p = Path(__file__).resolve().parent / key_path
        if p.exists():
            cred = ee.ServiceAccountCredentials(
                email, key_data=p.read_text(encoding="utf-8")
            )
            ee.Initialize(cred, project=project)
            print(f"Earth Engine :: key file ({p.name})")
            return

    # 3. Interactive — local only
    try:
        ee.Initialize(project=project)
    except Exception:
        ee.Authenticate()
        ee.Initialize(project=project)
    print(f"Earth Engine :: interactive ({project})")


def load_aoi(path=AOI_FILE):
    """
    Load the dissolved Mandya + Mysuru AOI as an ee.Geometry.

    Reads the GeoJSON locally (no geemap dependency) so this
    works in a bare python env.
    """
    import json

    with open(path) as f:
        gj = json.load(f)

    if gj.get("type") == "FeatureCollection":
        geoms = [ee.Geometry(feat["geometry"]) for feat in gj["features"]]
        geom = ee.Geometry.MultiPolygon(
            ee.List([g.coordinates() for g in geoms]).flatten()
        ) if len(geoms) > 1 else geoms[0]
    elif gj.get("type") == "Feature":
        geom = ee.Geometry(gj["geometry"])
    else:
        geom = ee.Geometry(gj)

    # simplify slightly — district polygons carry a lot of
    # vertices and unsimplified geometry slows every filter
    return geom.simplify(maxError=30)


# =========================================================
# SENTINEL-2  (optical, 10-20 m, 5-day revisit)
# =========================================================

def mask_s2(image):
    """Drop cloud, shadow, cirrus and snow via the SCL band."""
    scl = image.select("SCL")
    mask = ee.Image.constant(1)
    for cls in S2_SCL_DROP:
        mask = mask.And(scl.neq(cls))
    return image.updateMask(mask)


def sentinel2(aoi, start, end):
    """
    Sentinel-2 L2A surface reflectance.

    Scaled to true reflectance (0-1) and cloud-masked.
    Adds NDVI, EVI, NDWI, NDRE, LSWI, SAVI.
    """
    col = (
        ee.ImageCollection(COLLECTIONS["sentinel2"])
        .filterBounds(aoi)
        .filterDate(start, end)
        .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", S2_CLOUD_MAX))
        .map(mask_s2)
        .select(S2_BANDS)
    )

    def add_indices(img):
        img = img.multiply(S2_SCALE_FACTOR).copyProperties(
            img, img.propertyNames()
        )
        img = ee.Image(img)

        ndvi = img.normalizedDifference(["B8", "B4"]).rename("NDVI")
        ndwi = img.normalizedDifference(["B3", "B8"]).rename("NDWI")
        ndre = img.normalizedDifference(["B8", "B5"]).rename("NDRE")
        lswi = img.normalizedDifference(["B8", "B11"]).rename("LSWI")

        evi = img.expression(
            "2.5 * (N - R) / (N + 6*R - 7.5*B + 1)",
            {"N": img.select("B8"), "R": img.select("B4"), "B": img.select("B2")},
        ).rename("EVI")

        savi = img.expression(
            "1.5 * (N - R) / (N + R + 0.5)",
            {"N": img.select("B8"), "R": img.select("B4")},
        ).rename("SAVI")

        # GNDVI — green NDVI, more sensitive to chlorophyll / canopy N.
        # Uses Green (B3) instead of Red; better at high biomass where
        # NDVI saturates (e.g. dense sugarcane, mid-season paddy).
        gndvi = img.normalizedDifference(["B8", "B3"]).rename("GNDVI")

        # MSAVI2 — modified SAVI, self-calibrating soil correction.
        # Handles bare-soil / early-season pixels better than SAVI
        # without needing a manual L factor. Standard MSAVI2 form:
        #   0.5 * (2*NIR + 1 - sqrt((2*NIR+1)^2 - 8*(NIR - Red)))
        msavi = img.expression(
            "0.5 * (2*N + 1 - sqrt((2*N + 1)**2 - 8*(N - R)))",
            {"N": img.select("B8"), "R": img.select("B4")},
        ).rename("MSAVI")

        return img.addBands([ndvi, ndwi, ndre, lswi, evi, savi,
                             gndvi, msavi])

    return col.map(add_indices)


# =========================================================
# SENTINEL-1  (SAR, 10 m, all-weather)
# =========================================================

def sentinel1(aoi, start, end, orbit_pass=S1_ORBIT_PASS,
              relative_orbit=S1_RELATIVE_ORBIT):
    """
    Sentinel-1 GRD, IW mode, dual-pol.

    GEE's S1_GRD is already calibrated to sigma0 (dB),
    thermal-noise-removed and terrain-corrected. Speckle
    filtering is NOT applied by GEE — added here.

    CRITICAL: one orbit pass AND one relative orbit only.
    Mixing look angles corrupts the backscatter time series.
    """
    col = (
        ee.ImageCollection(COLLECTIONS["sentinel1"])
        .filterBounds(aoi)
        .filterDate(start, end)
        .filter(ee.Filter.eq("instrumentMode", "IW"))
        .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VV"))
        .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VH"))
        .filter(ee.Filter.eq("orbitProperties_pass", orbit_pass))
        .select(S1_BANDS)
    )

    if relative_orbit is not None:
        col = col.filter(
            ee.Filter.eq("relativeOrbitNumber_start", relative_orbit)
        )

    def speckle_and_ratio(img):
        # Refined-Lee is heavy in GEE; a focal mean over a
        # 30 m circle is the standard lightweight substitute.
        sm = img.focal_mean(30, "circle", "meters")

        vv = sm.select("VV")
        vh = sm.select("VH")

        # dB -> linear power BEFORE any ratio or averaging
        vv_lin = ee.Image(10).pow(vv.divide(10))
        vh_lin = ee.Image(10).pow(vh.divide(10))

        ratio = vh_lin.divide(vv_lin).rename("VH_VV")
        rvi = vh_lin.multiply(4).divide(
            vv_lin.add(vh_lin)
        ).rename("RVI")

        # GLCM texture on VV backscatter — captures the spatial
        # roughness pattern, not just brightness. Crops differ in
        # canopy structure (flooded paddy = smooth, sugarcane =
        # rough), so texture adds a discriminator the per-pixel
        # backscatter alone misses.
        # glcmTexture needs an integer image; scale dB to int first.
        vv_int = sm.select("VV").multiply(100).toInt32()
        glcm = vv_int.glcmTexture(size=3)
        # keep two of the most informative measures:
        #   contrast  = local intensity variation (roughness)
        #   corr      = linear dependency of grey levels (structure)
        vv_contrast = glcm.select("VV_contrast").rename(
            "VV_contrast").toFloat()
        vv_corr = glcm.select("VV_corr").rename("VV_corr").toFloat()

        return (
            sm.rename(["VV", "VH"])
            .addBands([ratio, rvi, vv_contrast, vv_corr])
            .copyProperties(img, ["system:time_start"])
        )

    return col.map(speckle_and_ratio)


def inspect_s1_orbits(aoi, start, end, orbit_pass=S1_ORBIT_PASS):
    """
    Run this ONCE before anything else.

    Prints which relative orbits cover the AOI and how many
    scenes each has. Pin the best one in gee_config.py.
    """
    col = (
        ee.ImageCollection(COLLECTIONS["sentinel1"])
        .filterBounds(aoi)
        .filterDate(start, end)
        .filter(ee.Filter.eq("instrumentMode", "IW"))
        .filter(ee.Filter.eq("orbitProperties_pass", orbit_pass))
    )

    orbits = col.aggregate_array("relativeOrbitNumber_start").getInfo()

    from collections import Counter
    counts = Counter(orbits)

    print(f"\nSentinel-1 {orbit_pass} orbits over AOI ({start} to {end})")
    print("-" * 50)
    for orbit, n in counts.most_common():
        print(f"  relative orbit {orbit:>3} : {n:>3} scenes")
    print("-" * 50)
    print("Pin the highest-count orbit as S1_RELATIVE_ORBIT.\n")

    return counts


# =========================================================
# LANDSAT 8 / 9  (optical + thermal, 30 m, 16-day)
# =========================================================

def _mask_landsat(img):
    """QA_PIXEL bit 3 = cloud, bit 4 = cloud shadow."""
    qa = img.select("QA_PIXEL")
    mask = (
        qa.bitwiseAnd(1 << 3).eq(0)
        .And(qa.bitwiseAnd(1 << 4).eq(0))
    )
    return img.updateMask(mask)


def landsat(aoi, start, end, include_thermal=True):
    """
    Landsat 8 + 9 Collection-2 Level-2, merged.

    Fills Sentinel-2 gaps and — more importantly — carries a
    thermal band, which is the only 30 m LST source available.
    """
    def prep(col_id):
        col = (
            ee.ImageCollection(col_id)
            .filterBounds(aoi)
            .filterDate(start, end)
            .filter(ee.Filter.lt("CLOUD_COVER", LANDSAT_CLOUD_MAX))
            .map(_mask_landsat)
        )
        return col

    merged = prep(COLLECTIONS["landsat8"]).merge(prep(COLLECTIONS["landsat9"]))

    def scale_and_index(img):
        opt = (
            img.select(LANDSAT_BANDS)
            .multiply(LANDSAT_MULT)
            .add(LANDSAT_ADD)
        )
        ndvi = opt.normalizedDifference(["SR_B5", "SR_B4"]).rename("NDVI")
        ndwi = opt.normalizedDifference(["SR_B3", "SR_B5"]).rename("NDWI")

        out = opt.addBands([ndvi, ndwi])

        if include_thermal:
            lst = (
                img.select(LANDSAT_THERMAL)
                .multiply(0.00341802).add(149.0).subtract(273.15)
                .rename("LST")
            )
            out = out.addBands(lst)

        return out.copyProperties(img, ["system:time_start"])

    return merged.map(scale_and_index)


# =========================================================
# MODIS  (baseline only — sensor is end-of-life)
# =========================================================

def modis_ndvi(aoi, start, end):
    """
    MOD13Q1 16-day NDVI, 250 m.

    Panel flagged MODIS as end-of-life. Kept ONLY as a
    long-baseline reference for multi-year VCI min/max,
    where its 2000-onwards archive is genuinely useful.
    Use VIIRS for anything forward-looking.
    """
    return (
        ee.ImageCollection(COLLECTIONS["modis_ndvi"])
        .filterBounds(aoi)
        .filterDate(start, end)
        .select(["NDVI", "EVI"])
        .map(lambda i: i.multiply(1e-4)
             .copyProperties(i, ["system:time_start"]))
    )


def viirs_ndvi(aoi, start, end):
    """
    VNP13A1 16-day vegetation indices, 500 m.
    The designated MODIS successor — panel explicitly
    asked for this to be explored.
    """
    return (
        ee.ImageCollection(COLLECTIONS["viirs_ndvi"])
        .filterBounds(aoi)
        .filterDate(start, end)
        .select(["NDVI", "EVI"])
        .map(lambda i: i.multiply(1e-4)
             .copyProperties(i, ["system:time_start"]))
    )


def modis_lst(aoi, start, end):
    """
    MOD11A2 8-day land surface temperature, 1 km.

    Panel called LST out specifically as missing from most
    submissions and valuable for moisture-stress work.
    Native 8-day cadence matches the compositing window.
    """
    return (
        ee.ImageCollection(COLLECTIONS["modis_lst"])
        .filterBounds(aoi)
        .filterDate(start, end)
        .select(["LST_Day_1km", "LST_Night_1km"])
        .map(lambda i: i.multiply(0.02).subtract(273.15)
             .rename(["LST_day", "LST_night"])
             .copyProperties(i, ["system:time_start"]))
    )


# =========================================================
# METEOROLOGY  (drives the FAO-56 engine)
# =========================================================

def chirps(aoi, start, end):
    """CHIRPS daily rainfall, ~5.5 km, mm/day."""
    return (
        ee.ImageCollection(COLLECTIONS["chirps"])
        .filterBounds(aoi)
        .filterDate(start, end)
        .select("precipitation")
    )


def era5(aoi, start, end):
    """
    ERA5-Land daily aggregates, ~11 km.

    Carries every Penman-Monteith input in one collection —
    no need to hunt separate temperature / wind / radiation
    sources. Unit conversions are applied downstream in the
    FAO-56 module, not here, so raw values stay inspectable.

      temperature      K        -> subtract 273.15
      wind u/v         m/s @10m -> u2 = u10 * 0.748
      solar radiation  J/m2     -> divide by 1e6 for MJ/m2/day
      precipitation    m        -> multiply by 1000 for mm
    """
    return (
        ee.ImageCollection(COLLECTIONS["era5"])
        .filterBounds(aoi)
        .filterDate(start, end)
        .select(ERA5_BANDS)
    )


# =========================================================
# TERRAIN + LAND COVER
# =========================================================

def terrain(aoi):
    """SRTM elevation with derived slope and aspect."""
    dem = ee.Image(COLLECTIONS["srtm"]).clip(aoi)
    return ee.Image.cat([
        dem.rename("elevation"),
        ee.Terrain.slope(dem).rename("slope"),
        ee.Terrain.aspect(dem).rename("aspect"),
    ]).toFloat()          # <-- SRTM is Int16, slope/aspect Float32


def cropland_mask(aoi, year):
    """
    Binary cropland mask from Dynamic World annual mode.

    1 = agricultural, 0 = everything else. Applied to every
    exported stack so advisory pixels never land on rooftops,
    water bodies or forest.
    """
    dw = (
        ee.ImageCollection(COLLECTIONS["dynamicworld"])
        .filterBounds(aoi)
        .filterDate(f"{year}-01-01", f"{year + 1}-01-01")
        .select("label")
        .mode()
    )
    non_ag = DW_NON_CROP_CLASSES
    ag = [c for c in range(9) if c not in non_ag]
    return dw.remap(
        non_ag + ag,
        [0] * len(non_ag) + [1] * len(ag)
    ).rename("crop_mask")
