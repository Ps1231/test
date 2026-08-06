"""
Task 1 — AOI & Boundaries
Builds the dissolved Mandya+Mysuru AOI and taluk layer from KGIS shapefiles.

Input:  District.shp, Taluk.shp  (KGIS, native CRS EPSG:32643)
Output: mandya_mysuru_aoi.geojson, mandya_mysuru_aoi_utm43n.geojson,
        taluk_boundaries.geojson, crs_standard.txt, aoi_preview.png
"""
import geopandas as gpd
from shapely.ops import unary_union
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

MASTER_CRS = 'EPSG:32643'   # UTM 43N — area/distance calcs
WEB_CRS    = 'EPSG:4326'    # WGS84 — GEE / web display / GeoJSON spec

MANDYA_CODE = '22'
MYSURU_CODE = '26'

# ---------------------------------------------------------------- load
district = gpd.read_file('/home/dell/tasks/task1/District/District.shp')
taluk    = gpd.read_file('/home/dell/tasks/task1/Taluk/Taluk.shp')

assert district.crs.to_string() == MASTER_CRS
assert taluk.crs.to_string() == MASTER_CRS

# ---------------------------------------------------------- filter + dissolve
aoi_districts = district[district['KGISDistri'].isin([MANDYA_CODE, MYSURU_CODE])].copy()
print("Districts selected:", aoi_districts['KGISDist_1'].tolist())

aoi_geom = unary_union(aoi_districts.geometry)
aoi = gpd.GeoDataFrame(
    {
        'aoi_name': ['Mandya_Mysuru_AOI'],
        'districts': ['Mandya, Mysuru'],
        'area_sq_km': [aoi_geom.area / 1e6],
    },
    geometry=[aoi_geom],
    crs=MASTER_CRS,
)
print("AOI area (sq km):", aoi['area_sq_km'].iloc[0])

# ---------------------------------------------------------------- taluks
aoi_taluks = taluk[taluk['KGISDistri'].isin([MANDYA_CODE, MYSURU_CODE])].copy()
aoi_taluks = aoi_taluks.rename(columns={
    'KGISTalukC': 'taluk_code',
    'KGISTalukN': 'taluk_name',
    'KGISDistri': 'district_code',
})
district_name_map = {MANDYA_CODE: 'Mandya', MYSURU_CODE: 'Mysuru'}
aoi_taluks['district_name'] = aoi_taluks['district_code'].map(district_name_map)
print("Taluks included:", len(aoi_taluks))

# ---------------------------------------------------------------- reproject
aoi_4326    = aoi.to_crs(WEB_CRS)
taluks_4326 = aoi_taluks.to_crs(WEB_CRS)
bounds = aoi_4326.total_bounds
print("AOI bounds (lon/lat):", bounds)

# ---------------------------------------------------------------- export
aoi_4326.to_file('mandya_mysuru_aoi.geojson', driver='GeoJSON')
aoi.to_file('mandya_mysuru_aoi_utm43n.geojson', driver='GeoJSON')
taluks_4326[['taluk_code', 'taluk_name', 'district_code', 'district_name', 'geometry']].to_file(
    'taluk_boundaries.geojson', driver='GeoJSON'
)

crs_doc = f"""CRS STANDARD — Mandya-Mysuru AOI
================================================
Master (area/distance calcs): EPSG:32643  (UTM Zone 43N)
Display / GEE / web:          EPSG:4326   (WGS84 lat/lon)

Source data (KGIS District & Taluk shapefiles) native CRS: EPSG:32643 — no reprojection needed for area work.
GeoJSON deliverables are exported in EPSG:4326 (GeoJSON spec requires WGS84).

AOI extent (EPSG:4326): minx={bounds[0]:.4f}, miny={bounds[1]:.4f}, maxx={bounds[2]:.4f}, maxy={bounds[3]:.4f}
AOI area: {aoi['area_sq_km'].iloc[0]:.2f} sq km
Districts: Mandya (KGIS code {MANDYA_CODE}), Mysuru (KGIS code {MYSURU_CODE})
Taluks included: {len(aoi_taluks)}
"""
with open('crs_standard.txt', 'w') as f:
    f.write(crs_doc)

# ---------------------------------------------------------------- preview map
fig, ax = plt.subplots(figsize=(9, 9))
taluks_4326.plot(ax=ax, column='district_name', edgecolor='white', linewidth=0.8,
                  legend=True, cmap='Set2', alpha=0.85)
aoi_4326.boundary.plot(ax=ax, color='black', linewidth=2)
for _, row in taluks_4326.iterrows():
    c = row.geometry.centroid
    ax.annotate(row['taluk_name'], (c.x, c.y), fontsize=7, ha='center')
ax.set_title('Mandya + Mysuru AOI — 16 Taluks (dissolved district boundary in black)')
ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude')
plt.tight_layout()
plt.savefig('aoi_preview.png', dpi=130)

print("Done.")
