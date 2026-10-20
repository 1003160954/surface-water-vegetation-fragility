# ============================================================
# Step 2 FIXED. Re-export JRC water metrics for ALL grid points
# Fix: unmask water bands and distance bands before sampleRegions
# ============================================================

import os
import math
import pandas as pd
import ee

PROJECT_ID = "ee-abao9474"

try:
    ee.Initialize(project=PROJECT_ID)
    print("Earth Engine initialized with existing credentials.")
except Exception:
    ee.Authenticate()
    ee.Initialize(project=PROJECT_ID)
    print("Earth Engine authenticated and initialized.")

print("Using project:", PROJECT_ID)

PATCHED_DIR = "/content/drive/MyDrive/China_VegClimate_ResponseOperator_Dataset_patched"
AQUA_DIR = os.path.join(PATCHED_DIR, "_aquatic_buffering")
GRID_FP = os.path.join(AQUA_DIR, "china_grid_points_for_water_sampling.csv")

grid = pd.read_csv(GRID_FP)

print("Grid points:", len(grid))
display(grid.head())

# ----------------------------
# JRC Global Surface Water
# ----------------------------
gsw = ee.Image("JRC/GSW1_4/GlobalSurfaceWater")

# Critical fix: unmask to 0, so non-water pixels are still sampled.
occurrence = gsw.select("occurrence").unmask(0).rename("water_occurrence")
seasonality = gsw.select("seasonality").unmask(0).rename("water_seasonality")
transition = gsw.select("transition").unmask(0).rename("water_transition")

# Define target water masks as 1/0 images, not selfMask images.
perm_water = occurrence.gte(80).And(seasonality.gte(10)).rename("is_perm_water")
seasonal_water = seasonality.gte(1).And(seasonality.lt(10)).rename("is_seasonal_water")
any_water = occurrence.gte(10).rename("is_any_water")

MAX_PIXELS = 1024
PIXEL_SIZE_M = 30
MAX_DIST_M = MAX_PIXELS * PIXEL_SIZE_M

def distance_to_water(mask_img, out_name):
    # Input is 1 for target water, 0 elsewhere.
    # fastDistanceTransform returns squared pixel distance by default.
    dist = (
        mask_img
        .unmask(0)
        .fastDistanceTransform(MAX_PIXELS)
        .sqrt()
        .multiply(PIXEL_SIZE_M)
        .rename(out_name)
    )
    return dist.unmask(MAX_DIST_M).min(MAX_DIST_M)

dist_perm_m = distance_to_water(perm_water, "dist_perm_water_m")
dist_seasonal_m = distance_to_water(seasonal_water, "dist_seasonal_water_m")
dist_any_m = distance_to_water(any_water, "dist_any_water_m")

water_img = (
    occurrence
    .addBands(seasonality)
    .addBands(transition)
    .addBands(perm_water.unmask(0).rename("is_perm_water"))
    .addBands(seasonal_water.unmask(0).rename("is_seasonal_water"))
    .addBands(any_water.unmask(0).rename("is_any_water"))
    .addBands(dist_perm_m)
    .addBands(dist_seasonal_m)
    .addBands(dist_any_m)
    .toFloat()
)

def df_to_fc(df):
    feats = []

    for _, r in df.iterrows():
        geom = ee.Geometry.Point([float(r["lon"]), float(r["lat"])])
        props = {
            "grid_id": str(r["grid_id"]),
            "lon": float(r["lon"]),
            "lat": float(r["lat"]),
            "region_id": str(r["region_id"]) if "region_id" in df.columns else ""
        }
        feats.append(ee.Feature(geom, props))

    return ee.FeatureCollection(feats)

CHUNK_SIZE = 5000
DRIVE_FOLDER = "China_Aquatic_Buffering_WaterMetrics_v2"

n_chunks = math.ceil(len(grid) / CHUNK_SIZE)
print("Total chunks:", n_chunks)

for i in range(n_chunks):
    sub = grid.iloc[i * CHUNK_SIZE : (i + 1) * CHUNK_SIZE].copy()
    fc = df_to_fc(sub)

    sampled = water_img.sampleRegions(
        collection=fc,
        properties=["grid_id", "lon", "lat", "region_id"],
        scale=300,
        geometries=False,
        tileScale=4
    )

    desc = f"water_metrics_v2_chunk_{i:03d}"

    task = ee.batch.Export.table.toDrive(
        collection=sampled,
        description=desc,
        folder=DRIVE_FOLDER,
        fileNamePrefix=desc,
        fileFormat="CSV"
    )

    task.start()
    print("Started:", desc, "rows:", len(sub))

print("\nAll v2 tasks started.")
print("After export finishes, use Drive folder:")
print(DRIVE_FOLDER)