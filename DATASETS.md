# Track A — GEE Datasets Reference

**Study area:** Mandya + Mysuru districts, Karnataka
**Extent:** ~11,266 km² · bbox 75.91–77.33°E, 11.74–13.06°N
**Master CRS:** EPSG:32643 (UTM 43N) · EPSG:4326 for GEE filtering

---

## 1. Why these datasets

Three objectives, three different data needs.

| Objective | What it needs | Primary datasets |
|---|---|---|
| **O1** Crop type classification | Dense multi-temporal signal, cloud-resilient | Sentinel-2, Sentinel-1 |
| **O2** Moisture stress by growth stage | Vegetation + water + thermal | Sentinel-2, Sentinel-1, MODIS LST, Landsat thermal |
| **O3** Water deficit + advisory | Meteorology for Penman-Monteith | ERA5-Land, CHIRPS, SRTM |

Everything below is on GEE. Nothing needs a download, an account whitelisting, or ISRO permission — it runs today.

---

## 2. Dataset profiles

### Sentinel-2 L2A — `COPERNICUS/S2_SR_HARMONIZED`

| Property | Value |
|---|---|
| Resolution | 10 m (VNIR), 20 m (red-edge, SWIR) |
| Revisit | 5 days (2 satellites) |
| Archive | 2017 → present |
| Correction | Surface reflectance (Sen2Cor applied) |
| Scale factor | ÷10,000 |

**Why:** Primary optical source. Only sensor here with **red-edge bands** (B5, B6, B7), which give NDRE — the most sensitive pre-visual stress index available and a strong discriminator between crops with similar NDVI.

**Filters:** `CLOUDY_PIXEL_PERCENTAGE < 60`, SCL mask drops classes 3/8/9/10/11.

**Downstream:** NDVI, EVI, NDWI, NDRE, LSWI, SAVI → phenology curves (SOS/POS/EOS) → VCI → NDVI→Kc transform for the FAO-56 engine.

**What it looks like:** Over Rabi, an unbroken rise-plateau-decline per pixel. Over Kharif, a broken series — monsoon cloud will empty entire windows. That gap is the reason SAR exists in this pipeline.

**Caveat:** LSWI needs B11 at 20 m; if you export at 10 m, B11 is resampled, not measured.

---

### Sentinel-1 GRD — `COPERNICUS/S1_GRD`

| Property | Value |
|---|---|
| Resolution | 10 m |
| Revisit | 6–12 days |
| Band | C-band, dual-pol VV + VH |
| Preprocessing | Already calibrated to σ⁰ (dB), thermal-noise-removed, terrain-corrected |
| Not applied | Speckle filtering — we add it |

**Why:** The all-weather backbone. Cloud is transparent to C-band, so Kharif continuity depends entirely on this. Also carries structural and soil-moisture information that optical cannot see.

**Filters:** `instrumentMode = IW`, VV **and** VH present, **one orbit pass, one relative orbit**.

> **The single most important filter in this pipeline.** Different relative orbits image the same field at different incidence angles. Backscatter changes with incidence angle. Mix two orbits and your time series contains a step change that has nothing to do with the crop. Run `run_ingest.py orbits` first and pin the winner.

**dB handling:** σ⁰ is logarithmic. Never average dB directly.

```python
linear = 10 ** (db / 10)
mean_linear = linear.mean()
mean_db = 10 * log10(mean_linear)
```

The `composite_sar()` function does this. So does the VH/VV ratio.

**Downstream:** VV, VH, VH/VV, RVI → SAR relative SMI via change detection → gap filling when optical is missing.

**What it looks like:** Paddy gives the diagnostic **V-shape** — VH collapses to −22…−18 dB at flooding (specular reflection off standing water), then rises sharply through tillering. That single feature is the strongest paddy discriminator available and should be your first-pass rule.

---

### Landsat 8/9 — `LANDSAT/LC08/C02/T1_L2`, `LC09/C02/T1_L2`

| Property | Value |
|---|---|
| Resolution | 30 m optical, 100 m thermal (resampled to 30 m) |
| Revisit | 16 days each, 8 days combined |
| Archive | L8 2013→, L9 2021→ |
| Scale | ×0.0000275 − 0.2 |

**Why two reasons:** fills Sentinel-2 gaps with an independent acquisition schedule, and — more important — it's the **only 30 m thermal source**. MODIS LST is 1 km, far too coarse for field-level advisory.

**Downstream:** NDVI/NDWI as S2 backup; `ST_B10` → LST at 30 m for stress detection.

**LST conversion:** `ST_B10 × 0.00341802 + 149 − 273.15` → °C

---

### MODIS LST — `MODIS/061/MOD11A2`

| Property | Value |
|---|---|
| Resolution | 1 km |
| Cadence | 8-day composite (native) |
| Scale | ×0.02, Kelvin |

**Why:** The panel called LST out specifically — most submissions omit it, and it adds real value to moisture-stress estimation. Native 8-day cadence lines up with your compositing window with no resampling.

**Downstream:** Thermal stress indicator; cumulative thermal indices (which the panel noted work well in Rabi); crop water stress index if you go that route.

**Caveat:** 1 km over a district with fragmented fields means one pixel spans many farms. Treat it as regional context, not field-level truth. Landsat thermal is the field-level version.

---

### MODIS NDVI — `MODIS/061/MOD13Q1` ⚠️

| Property | Value |
|---|---|
| Resolution | 250 m |
| Cadence | 16-day |
| Archive | 2000 → present |

**Status: end-of-life.** The panel flagged this. Keep it for **one purpose only** — the 25-year archive is the best available source for multi-year NDVI min/max, which VCI requires. For anything forward-looking, use VIIRS.

---

### VIIRS NDVI — `NOAA/VIIRS/001/VNP13A1` ⭐

| Property | Value |
|---|---|
| Resolution | 500 m |
| Cadence | 16-day |
| Archive | 2012 → present |

**Why:** The designated MODIS successor. The panel asked for it explicitly. Including it demonstrates you understood the end-of-life warning rather than just nodding at it.

---

### CHIRPS — `UCSB-CHG/CHIRPS/DAILY`

| Property | Value |
|---|---|
| Resolution | ~5.5 km (0.05°) |
| Cadence | Daily |
| Unit | mm/day |

**Why:** Effective rainfall — the supply side of the water balance. Without it, "deficit" is meaningless.

**Downstream:** Window-sum rainfall → effective rainfall → ΔW = ETc − (P_eff + antecedent soil moisture).

**Note:** IMD gridded 0.25° is the India-specific alternative and aligns with IMD standard weeks. CHIRPS is finer and available in GEE with no download. Use CHIRPS now; add IMD if the panel prefers it.

---

### ERA5-Land — `ECMWF/ERA5_LAND/DAILY_AGGR`

| Property | Value |
|---|---|
| Resolution | ~11 km (0.1°) |
| Cadence | Daily aggregates |

**Why this is the important one for Objective 3:** it carries **every** Penman-Monteith input in a single collection. No hunting for separate temperature, wind, radiation and humidity sources.

| Band | Unit | Conversion |
|---|---|---|
| `temperature_2m` | K | −273.15 → °C |
| `dewpoint_temperature_2m` | K | −273.15 → RH |
| `u/v_component_of_wind_10m` | m/s @ 10 m | u₂ = u₁₀ × 0.748 |
| `surface_solar_radiation_downwards_sum` | J/m² | ÷10⁶ → MJ/m²/day |
| `surface_pressure` | Pa | ÷1000 → kPa |
| `total_precipitation_sum` | m | ×1000 → mm |

**This replaces the hardcoded `ET0_PER_DAY_MM = 3.0` in the current backend.** That constant is the single largest scientific gap between the deck and the code.

**Caveat:** 11 km is coarse for field advisory. It's the standard tradeoff — ET₀ varies smoothly over space, so coarse met data is more defensible than coarse vegetation data. Say this if asked rather than being caught by it.

---

### SRTM — `USGS/SRTMGL1_003`

30 m elevation → slope, aspect. Needed for drainage context, waterlogging risk in low-lying command areas, and terrain correction sanity checks.

CartoDEM (Bhoonidhi) is the Indian national DEM and a better answer for an ISRO audience — but it needs a download. SRTM now, CartoDEM later.

---

### Dynamic World — `GOOGLE/DYNAMICWORLD/V1`

10 m near-real-time land cover, 9 classes. Used as a **cropland mask** — annual mode, keep agricultural classes, drop water/trees/built/snow.

**Why it matters:** without it, advisory pixels land on rooftops, roads and water bodies. Every exported stack in this module is masked by it.

**Do not use it as crop-type labels.** Dynamic World has one "crops" class. It cannot distinguish sugarcane from paddy from ragi. That distinction is exactly what Objective 1 asks for, and it requires the ISRO ground-truth polygons.

---

## 3. What is deliberately NOT here

| Dataset | Why not |
|---|---|
| LISS-III / LISS-IV / AWiFS | Not on GEE. Bhoonidhi portal → manual download → GEE asset upload |
| EOS-04 SAR | Not on GEE. Same path. Direct-download on Bhoonidhi |
| EOS-06 8-day NDVI | Not on GEE. Bhoonidhi direct-download |
| Soil maps | ISRO offered them — ask. SoilGrids is the fallback |
| Canal command boundary | India-WRIS / CNNL / ask ISRO |
| **Ground truth polygons** | **ISRO provides at hackathon only. Hard blocker for O1.** |
| NISAR | Not launched into usable archive. Future scope framing only |

---

## 4. Compositing

**Window:** parameterised. `WINDOW_DAYS = 8` (default) or `7` for IMD standard weeks.

The panel suggested IMD weeks for India-specificity. But ISRO's own EOS-06 NDVI product is natively 8-day. That tension is unresolved — build both, ask on Monday, and present it as a considered question rather than an oversight.

**Method:**

| Data | Reducer | Why |
|---|---|---|
| Optical | median | Robust to residual cloud that survived the SCL mask |
| SAR | mean in linear power | Preserves radiometry; dB mean is mathematically wrong |
| Rainfall | sum | It's a flux, not a state |
| Weather | mean | It's a state |

**Band naming:** `<VAR>_t<NN>` — `NDVI_t01`, `VH_VV_t07`, `LST_day_t12`

Rabi at 8-day → T = 23 windows. With 11 features per window (6 optical + 4 SAR + 1 thermal) that's ~253 bands per season stack.

---

## 5. Export strategy

**Do not export 10 m over the full AOI.** 11,266 km² at 10 m is ~113 billion cells per band. It will not complete.

| Approach | When |
|---|---|
| 20 m, tiled 3×3, → Drive | Full-AOI rasters for visualisation and QA |
| `Export.image.toAsset` | Iteration — stays server-side, re-samplable without re-download |
| `sampleRegions` → CSV | **The one you actually want for training** |

For classification, your model consumes `[N_pixels, T, F]` — a table, not an image. Exporting rasters and re-flattening them locally wastes days of transfer on data you immediately discard. Once ISRO hands over polygons:

```python
samples = stack.sampleRegions(
    collection=ground_truth_fc,
    properties=['crop_label'],
    scale=10,
    tileScale=4,
)
ee.batch.Export.table.toDrive(samples, fileFormat='CSV').start()
```

Until then, sample on a regular grid to validate the extraction path end to end.

---

## 6. Run order

```bash
python run_ingest.py orbits      # FIRST — pin S1_RELATIVE_ORBIT
python run_ingest.py inventory   # what exists, where the gaps are
python run_ingest.py static      # terrain + cropland mask
python run_ingest.py meteo       # CHIRPS + ERA5 → unblocks O3
python run_ingest.py stack       # the season composite stack
```

Steps 1–2 take minutes and change how you plan everything else. Step 4 is what makes Objective 3 real rather than hardcoded.

---

## 7. Known gaps

- **`S1_RELATIVE_ORBIT = None`** until you run `orbits`. Leaving it None means all orbits in the pass are mixed.
- **`GEE_PROJECT`** must match your Cloud project.
- **`AOI_FILE`** expects `data/boundaries/mandya_mysuru_aoi.geojson` from Task 1.
- **Kharif windows will be empty of optical.** That's expected and it's the point — record it in `gap_statistics.csv` and use it as evidence for why fusion is necessary.
- **10 m analysis grid** deferred until the canal command boundary settles the real extent.
