"""
=========================================================
BAH 2026 — PS-6 :: GEE Track A Configuration
=========================================================
Study area : Mandya + Mysuru districts, Karnataka
             (KRS / Visvesvaraya Canal command + rainfed)

All GEE-side constants live here. Nothing else in the
pipeline should hardcode a date, a scale, or a band name.
=========================================================
"""

from pathlib import Path

# =========================================================
# PROJECT
# =========================================================

GEE_PROJECT = "orbital-valor-500812-v4"          # <-- your Cloud project id

ROOT_DIR = Path(__file__).resolve().parent
DATA_DIR = ROOT_DIR / "data"
BOUNDARY_DIR = DATA_DIR / "boundaries"

AOI_FILE = BOUNDARY_DIR / "mandya_mysuru_aoi.geojson"

# Drive folder that every Export task writes into
DRIVE_FOLDER = "BAH2026_GEE"

# GEE asset root for intermediate stacks (faster to re-sample
# than re-downloading from Drive)
ASSET_ROOT = f"projects/{GEE_PROJECT}/assets/bah2026"

MAX_PIXELS = 1e13


# =========================================================
# COORDINATE REFERENCE
# =========================================================
# EPSG:32643 (UTM 43N)  -> area / distance / export grid
# EPSG:4326             -> GEE filtering, web display

CRS_METRIC = "EPSG:32643"
CRS_GEOGRAPHIC = "EPSG:4326"

EXPORT_CRS = CRS_METRIC

# AOI bounding box (EPSG:4326) from the dissolved district shapefile
AOI_BBOX = [75.91, 11.74, 77.33, 13.06]


# =========================================================
# TIME PERIODS
# =========================================================
# Panel protocol: train on 2023 + 2024, TEST on 2025.
# 2026 is excluded — abnormal year, crop area very low.
# Start with RABI (clear skies, fewer optical gaps).

SEASONS = {
    # ---- 2023-24 ----
    "kharif_2023":    ("2023-06-01T00:00:00Z", "2023-10-31T23:59:59Z"),
    "rabi_2023_24":   ("2023-11-01T00:00:00Z", "2024-02-28T23:59:59Z"),
    "summer_2024":    ("2024-03-01T00:00:00Z", "2024-05-31T23:59:59Z"),
    # ---- 2024-25 ----
    "kharif_2024":    ("2024-06-01T00:00:00Z", "2024-10-31T23:59:59Z"),
    "rabi_2024_25":   ("2024-11-01T00:00:00Z", "2025-02-28T23:59:59Z"),
    "summer_2025":    ("2025-03-01T00:00:00Z", "2025-05-31T23:59:59Z"),
    # ---- 2025 (test) ----
    "kharif_2025":    ("2025-06-01T00:00:00Z", "2025-11-30T23:59:59Z"),
    "rabi_2025_26":   ("2025-12-01T00:00:00Z", "2026-02-28T23:59:59Z"),
    "summer_2026":    ("2026-03-01T00:00:00Z", "2026-05-31T23:59:59Z")
}

DEFAULT_SEASON = "rabi_2023_24"

TRAIN_SEASONS = ["rabi_2023_24", "kharif_2024", "rabi_2024_25"]
TEST_SEASONS = ["kharif_2025"]

DEFAULT_SEASON = "rabi_2023_24"


# =========================================================
# TEMPORAL COMPOSITING
# =========================================================
# Panel suggested IMD standard weeks (7 d) over 8-day windows
# for India-specificity. EOS-06 NDVI is natively 8-day.
# Keep this parameterised, generate both, decide with ISRO.

WINDOW_DAYS = 8          # set to 7 for IMD weeks
WINDOW_MODE = "8day"     # "8day" | "imdweek"  (used in filenames)


# =========================================================
# EXPORT SCALES (metres)
# =========================================================
# 10 m over the full 11,266 km2 district AOI is ~113 billion
# cells — do NOT export that. Use 20 m for full-AOI rasters,
# drop to 10 m only after clipping to the canal command core.

EXPORT_SCALE = {
    "sentinel2": 20,
    "sentinel1": 20,
    "landsat":   30,
    "modis":    250,
    "viirs":    500,
    "lst":     1000,
    "weather": 11132,   # ERA5-Land native ~0.1 deg
    "chirps":  5566,    # CHIRPS native ~0.05 deg
    "dem":       30,
    "lulc":      10,
}

# Split the AOI into a grid of export tiles so no single
# Drive task exceeds GEE limits.
EXPORT_TILE_ROWS = 3
EXPORT_TILE_COLS = 3


# =========================================================
# EARTH ENGINE COLLECTIONS
# =========================================================

COLLECTIONS = {
    # --- Optical ---
    "sentinel2":  "COPERNICUS/S2_SR_HARMONIZED",
    "landsat8":   "LANDSAT/LC08/C02/T1_L2",
    "landsat9":   "LANDSAT/LC09/C02/T1_L2",

    # --- SAR ---
    "sentinel1":  "COPERNICUS/S1_GRD",

    # --- Vegetation baseline ---
    "modis_ndvi": "MODIS/061/MOD13Q1",      # deprecated, baseline only
    "viirs_ndvi": "NASA/VIIRS/002/VNP13A1",  # MODIS replacement

    # --- Thermal ---
    "modis_lst":  "MODIS/061/MOD11A2",

    # --- Meteorology ---
    "chirps":     "UCSB-CHG/CHIRPS/DAILY",
    "era5":       "ECMWF/ERA5_LAND/DAILY_AGGR",

    # --- Terrain ---
    "srtm":       "USGS/SRTMGL1_003",

    # --- Land cover ---
    "dynamicworld": "GOOGLE/DYNAMICWORLD/V1",
    "worldcover":   "ESA/WorldCover/v200",

    # --- Soil (backup until ISRO provides soil maps) ---
    "soilgrids_sand": "projects/soilgrids-isric/sand_mean",
    "soilgrids_clay": "projects/soilgrids-isric/clay_mean",
}


# =========================================================
# FILTERS
# =========================================================

S2_CLOUD_MAX = 60          # CLOUDY_PIXEL_PERCENTAGE
LANDSAT_CLOUD_MAX = 60     # CLOUD_COVER

# Sentinel-1: pick ONE orbit pass and ONE relative orbit.
# Mixing orbits changes the look angle and corrupts any
# backscatter time series. Run inspect_s1_orbits() first,
# then pin S1_RELATIVE_ORBIT to the best-covered number.
S1_ORBIT_PASS = "DESCENDING"
S1_RELATIVE_ORBIT = 165  # e.g. 63 — set after inspection

# Sentinel-2 Scene Classification Layer classes to drop
S2_SCL_DROP = [3, 8, 9, 10, 11]
# 3 cloud shadow | 8 cloud med | 9 cloud high | 10 cirrus | 11 snow

# Dynamic World: classes counted as agricultural land
DW_CROP_CLASSES = [4]              # 4 = crops
DW_NON_CROP_CLASSES = [0, 1, 6, 8]  # water, trees, built, snow/ice


# =========================================================
# BANDS
# =========================================================

S2_BANDS = ["B2", "B3", "B4", "B5", "B8", "B11", "B12"]
S2_SCALE_FACTOR = 1e-4

S1_BANDS = ["VV", "VH"]

LANDSAT_BANDS = ["SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7"]
LANDSAT_THERMAL = "ST_B10"
LANDSAT_MULT = 2.75e-5
LANDSAT_ADD = -0.2

ERA5_BANDS = [
    "temperature_2m",
    "temperature_2m_max",
    "temperature_2m_min",
    "dewpoint_temperature_2m",
    "u_component_of_wind_10m",
    "v_component_of_wind_10m",
    "surface_solar_radiation_downwards_sum",
    "surface_pressure",
    "total_precipitation_sum",
]


# =========================================================
# CROP CLASSES  (Mandya / Mysuru — corrected from Punjab)
# =========================================================

CROP_CLASSES = {
    1: "Sugarcane",
    2: "Paddy",
    3: "Ragi",
    4: "Maize",
    5: "Pulses",
    6: "Plantation",   # banana / mulberry / coconut
    7: "Other",
    0: "Fallow",
}
