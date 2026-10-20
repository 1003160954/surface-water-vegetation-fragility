# ============================================================
# Reproduce main results for aquatic buffering and vegetation fragility
# Main framework:
#   1. Heat–dryness stress: HDSI
#   2. Vegetation negative response: EVI_neg / LAI_neg / joint_failure
#   3. Aquatic context:
#        near permanent water
#        near seasonal water
#        far or weak water
#   4. Distance-band fragility:
#        Risk(high HDSI) - Risk(low HDSI)
#   5. Regional heterogeneity
#   6. Bootstrap 95% CI
# ============================================================

import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

try:
    import geopandas as gpd
except Exception:
    gpd = None

# ============================================================
# 0. Paths
# ============================================================

PATCHED_DIR = "/content/drive/MyDrive/China_VegClimate_ResponseOperator_Dataset_patched"
AQUA_DIR = os.path.join(PATCHED_DIR, "_aquatic_buffering")

PANEL_FP = os.path.join(AQUA_DIR, "china_aquatic_buffering_panel_v2.parquet")
CHINA_BOUNDARY_PATH = "/content/drive/MyDrive/China_boundary/china_boundary_full.gpkg"

OUT_DIR = os.path.join(AQUA_DIR, "_main_results_reproduced")
FIG_DIR = os.path.join(OUT_DIR, "figs")
TAB_DIR = os.path.join(OUT_DIR, "tables")

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TAB_DIR, exist_ok=True)

print("Panel:", PANEL_FP)
print("Exists:", os.path.exists(PANEL_FP))


# ============================================================
# 1. Nature-like plotting style
# ============================================================

mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.5,
    "axes.titlesize": 10,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.dpi": 300,
})

COL = {
    "blue": "#3B6EA8",
    "light_blue": "#9BBBD8",
    "orange": "#E69F00",
    "red": "#B2182B",
    "green": "#4E9A51",
    "teal": "#009E73",
    "purple": "#7E62A3",
    "gray": "#666666",
    "light_gray": "#D9D9D9",
    "dark": "#222222",
}

CONTEXT_COLORS = {
    "near_permanent_lt1km": COL["blue"],
    "near_seasonal_lt1km": COL["orange"],
    "far_or_weak_water": COL["gray"],
}

HDSI_COLORS = {
    "low": COL["blue"],
    "mid": COL["orange"],
    "high": COL["red"],
}

def clean_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.18, linewidth=0.4)

def panel_label(ax, label, x=-0.08, y=1.04):
    ax.text(
        x, y, label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="top",
        ha="left"
    )

def savefig(fig, name):
    png = os.path.join(FIG_DIR, name + ".png")
    pdf = os.path.join(FIG_DIR, name + ".pdf")
    fig.savefig(png, bbox_inches="tight", dpi=300)
    fig.savefig(pdf, bbox_inches="tight")
    print("Saved:", png)
    print("Saved:", pdf)

def read_china_boundary():
    if gpd is None:
        return None
    if not os.path.exists(CHINA_BOUNDARY_PATH):
        print("[WARN] China boundary not found.")
        return None
    try:
        china = gpd.read_file(CHINA_BOUNDARY_PATH, layer="china_full").to_crs("EPSG:4326")
    except Exception:
        china = gpd.read_file(CHINA_BOUNDARY_PATH).to_crs("EPSG:4326")
    china = china.dissolve().reset_index(drop=True)
    return china


# ============================================================
# 2. Load data and build core variables
# ============================================================

df = pd.read_parquet(PANEL_FP)
print("Raw panel:", df.shape)

# Convert important variables
num_cols = [
    "lon", "lat", "year",
    "spring_EVI_z", "spring_LAI_z",
    "late_summer_EVI_z", "late_summer_LAI_z",
    "summer_VPD_z", "summer_LST_Day_z",
    "summer_SM_rootproxy_z", "summer_Precip_z", "summer_Tair_z",
    "elevation", "slope", "LC_Type1",
    "dist_perm_water_km", "dist_seasonal_water_km", "dist_any_water_km",
    "water_occurrence", "water_seasonality", "water_transition"
]

for c in num_cols:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

# Build HDSI if missing
if "HDSI" not in df.columns:
    df["HDSI"] = np.nanmean(
        np.vstack([
            df["summer_VPD_z"].values,
            df["summer_LST_Day_z"].values,
            -df["summer_SM_rootproxy_z"].values,
            -df["summer_Precip_z"].values
        ]),
        axis=0
    )
else:
    df["HDSI"] = pd.to_numeric(df["HDSI"], errors="coerce")

# Build negative response variables if missing
if "EVI_neg_z05" not in df.columns:
    df["EVI_neg_z05"] = (df["late_summer_EVI_z"] < -0.5).astype(int)
else:
    df["EVI_neg_z05"] = pd.to_numeric(df["EVI_neg_z05"], errors="coerce").astype("Int64")

if "LAI_neg_z05" not in df.columns:
    df["LAI_neg_z05"] = (df["late_summer_LAI_z"] < -0.5).astype(int)
else:
    df["LAI_neg_z05"] = pd.to_numeric(df["LAI_neg_z05"], errors="coerce").astype("Int64")

df["joint_failure"] = (
    (df["EVI_neg_z05"].fillna(0).astype(int) == 1) |
    (df["LAI_neg_z05"].fillna(0).astype(int) == 1)
).astype(int)

# Main valid sample
d = df[
    df["valid_model_both_y"] &
    df["HDSI"].notna() &
    df["dist_perm_water_km"].notna() &
    df["dist_seasonal_water_km"].notna() &
    df["EVI_neg_z05"].notna() &
    df["LAI_neg_z05"].notna()
].copy()

d["EVI_neg_z05"] = d["EVI_neg_z05"].astype(int)
d["LAI_neg_z05"] = d["LAI_neg_z05"].astype(int)
d["joint_failure"] = d["joint_failure"].astype(int)

# HDSI groups
q25 = d["HDSI"].quantile(0.25)
q75 = d["HDSI"].quantile(0.75)

d["HDSI_group"] = "mid"
d.loc[d["HDSI"] <= q25, "HDSI_group"] = "low"
d.loc[d["HDSI"] >= q75, "HDSI_group"] = "high"

d["HDSI_group"] = pd.Categorical(
    d["HDSI_group"],
    categories=["low", "mid", "high"],
    ordered=True
)

# Water context classification
# Priority:
#   1. near permanent water < 1 km
#   2. near seasonal water < 1 km but not near permanent
#   3. far or weak water
d["water_context"] = "far_or_weak_water"
d.loc[d["dist_perm_water_km"] < 1.0, "water_context"] = "near_permanent_lt1km"
d.loc[
    (d["dist_perm_water_km"] >= 1.0) &
    (d["dist_seasonal_water_km"] < 1.0),
    "water_context"
] = "near_seasonal_lt1km"

d["water_context"] = pd.Categorical(
    d["water_context"],
    categories=["near_permanent_lt1km", "near_seasonal_lt1km", "far_or_weak_water"],
    ordered=True
)

# Permanent-water distance bands
dist_bins = [0, 0.5, 1, 2, 5, 10, 20]
dist_labels = ["0-0.5", "0.5-1", "1-2", "2-5", "5-10", ">10"]

d["dist_band_perm"] = pd.cut(
    d["dist_perm_water_km"],
    bins=dist_bins,
    labels=dist_labels,
    include_lowest=True,
    right=False
)

d["dist_band_perm"] = pd.Categorical(
    d["dist_band_perm"],
    categories=dist_labels,
    ordered=True
)

# Macro regions
def assign_macro_region(row):
    lon = row["lon"]
    lat = row["lat"]
    elev = row["elevation"] if "elevation" in row else np.nan

    if pd.notna(elev) and elev >= 3000 and 75 <= lon <= 105 and 25 <= lat <= 42:
        return "Tibetan_Plateau"
    if lon >= 118 and lat >= 40:
        return "Northeast"
    if lon < 105 and lat >= 35:
        return "Northwest_Arid"
    if 105 <= lon < 122 and 34 <= lat < 40:
        return "North_China"
    if 110 <= lon < 123 and 27 <= lat < 34:
        return "East_China"
    if lon >= 105 and lat < 27:
        return "South_China"
    if lon < 110 and lat < 35:
        return "Southwest"
    return "Other"

if "macro_region" not in d.columns:
    d["macro_region"] = d.apply(assign_macro_region, axis=1)
else:
    d["macro_region"] = d["macro_region"].astype(str)

print("Valid sample:", d.shape)
print("Unique grids:", d["grid_id"].nunique())
print("Years:", sorted(d["year"].dropna().astype(int).unique()))
print("HDSI Q25/Q75:", q25, q75)
print("EVI neg rate:", d["EVI_neg_z05"].mean())
print("LAI neg rate:", d["LAI_neg_z05"].mean())
print("Joint failure rate:", d["joint_failure"].mean())


# ============================================================
# 3. Bootstrap helpers
# ============================================================

def bootstrap_mean_ci(values, n_boot=500, seed=42):
    values = np.asarray(values)
    values = values[np.isfinite(values)]
    n = len(values)
    if n == 0:
        return np.nan, np.nan, np.nan
    mean = np.mean(values)
    if n == 1:
        return mean, np.nan, np.nan
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        sample = rng.choice(values, size=n, replace=True)
        boots.append(np.mean(sample))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return mean, lo, hi

def bootstrap_diff_ci(a, b, n_boot=500, seed=42):
    """
    Difference: mean(a) - mean(b)
    """
    a = np.asarray(a)
    b = np.asarray(b)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) == 0 or len(b) == 0:
        return np.nan, np.nan, np.nan
    diff = np.mean(a) - np.mean(b)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        aa = rng.choice(a, size=len(a), replace=True)
        bb = rng.choice(b, size=len(b), replace=True)
        boots.append(np.mean(aa) - np.mean(bb))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return diff, lo, hi


# ============================================================
# 4. Data overview table
# ============================================================

overview_rows = [
    {"metric": "rows", "value": len(d)},
    {"metric": "unique_grids", "value": d["grid_id"].nunique()},
    {"metric": "years", "value": ",".join(map(str, sorted(d["year"].dropna().astype(int).unique())))},
    {"metric": "HDSI_Q25", "value": q25},
    {"metric": "HDSI_Q75", "value": q75},
    {"metric": "EVI_neg_rate", "value": d["EVI_neg_z05"].mean()},
    {"metric": "LAI_neg_rate", "value": d["LAI_neg_z05"].mean()},
    {"metric": "joint_failure_rate", "value": d["joint_failure"].mean()},
    {"metric": "dist_perm_mean_km", "value": d["dist_perm_water_km"].mean()},
    {"metric": "dist_perm_median_km", "value": d["dist_perm_water_km"].median()},
    {"metric": "dist_perm_p90_km", "value": d["dist_perm_water_km"].quantile(0.90)},
    {"metric": "dist_perm_max_km", "value": d["dist_perm_water_km"].max()},
]

overview = pd.DataFrame(overview_rows)
overview.to_csv(os.path.join(TAB_DIR, "main_overview.csv"), index=False)
display(overview)


# ============================================================
# 5. Main result 1:
#    Water context × HDSI group × response risk + CI
# ============================================================

response_cols = {
    "EVI": "EVI_neg_z05",
    "LAI": "LAI_neg_z05",
    "Joint": "joint_failure"
}

rows = []

for context in d["water_context"].cat.categories:
    for hgroup in d["HDSI_group"].cat.categories:
        sub = d[
            (d["water_context"] == context) &
            (d["HDSI_group"] == hgroup)
        ].copy()

        for resp_name, col in response_cols.items():
            mean, lo, hi = bootstrap_mean_ci(sub[col].values, n_boot=500, seed=42)
            rows.append({
                "water_context": context,
                "HDSI_group": hgroup,
                "response": resp_name,
                "risk": mean,
                "ci_low": lo,
                "ci_high": hi,
                "rows": len(sub),
                "grids": sub["grid_id"].nunique()
            })

water_context_risk = pd.DataFrame(rows)
water_context_risk.to_csv(os.path.join(TAB_DIR, "water_context_HDSI_risk_with_CI.csv"), index=False)
display(water_context_risk.head(20))


# ============================================================
# 6. Main result 2:
#    High-HDSI water context contrasts + CI
# ============================================================

contrast_rows = []

high = d[d["HDSI_group"] == "high"].copy()

for resp_name, col in response_cols.items():
    near_perm = high[high["water_context"] == "near_permanent_lt1km"][col].values
    near_season = high[high["water_context"] == "near_seasonal_lt1km"][col].values
    far_weak = high[high["water_context"] == "far_or_weak_water"][col].values

    diff_far, lo_far, hi_far = bootstrap_diff_ci(far_weak, near_perm, n_boot=500, seed=42)
    diff_season, lo_season, hi_season = bootstrap_diff_ci(near_season, near_perm, n_boot=500, seed=43)

    contrast_rows.append({
        "response": resp_name,
        "contrast": "far_or_weak - near_permanent",
        "risk_difference": diff_far,
        "ci_low": lo_far,
        "ci_high": hi_far,
        "n_group1": len(far_weak),
        "n_group2": len(near_perm)
    })

    contrast_rows.append({
        "response": resp_name,
        "contrast": "near_seasonal - near_permanent",
        "risk_difference": diff_season,
        "ci_low": lo_season,
        "ci_high": hi_season,
        "n_group1": len(near_season),
        "n_group2": len(near_perm)
    })

high_context_contrasts = pd.DataFrame(contrast_rows)
high_context_contrasts.to_csv(os.path.join(TAB_DIR, "high_HDSI_water_context_contrasts_with_CI.csv"), index=False)
display(high_context_contrasts)


# ============================================================
# 7. Main result 3:
#    Distance band × HDSI group risk + CI
# ============================================================

rows = []

for band in dist_labels:
    for hgroup in ["low", "mid", "high"]:
        sub = d[
            (d["dist_band_perm"].astype(str) == band) &
            (d["HDSI_group"].astype(str) == hgroup)
        ].copy()

        for resp_name, col in response_cols.items():
            mean, lo, hi = bootstrap_mean_ci(sub[col].values, n_boot=500, seed=42)
            rows.append({
                "dist_band_perm": band,
                "HDSI_group": hgroup,
                "response": resp_name,
                "risk": mean,
                "ci_low": lo,
                "ci_high": hi,
                "rows": len(sub),
                "grids": sub["grid_id"].nunique()
            })

dist_band_risk = pd.DataFrame(rows)
dist_band_risk.to_csv(os.path.join(TAB_DIR, "distance_band_HDSI_risk_with_CI.csv"), index=False)
display(dist_band_risk.head(20))


# ============================================================
# 8. Main result 4:
#    Fragility increase = Risk(high HDSI) - Risk(low HDSI)
#    by permanent-water distance band
# ============================================================

frag_rows = []

for band in dist_labels:
    sub_low = d[
        (d["dist_band_perm"].astype(str) == band) &
        (d["HDSI_group"].astype(str) == "low")
    ].copy()

    sub_high = d[
        (d["dist_band_perm"].astype(str) == band) &
        (d["HDSI_group"].astype(str) == "high")
    ].copy()

    for resp_name, col in response_cols.items():
        diff, lo, hi = bootstrap_diff_ci(
            sub_high[col].values,
            sub_low[col].values,
            n_boot=500,
            seed=42
        )

        frag_rows.append({
            "dist_band_perm": band,
            "response": resp_name,
            "fragility_increase": diff,
            "ci_low": lo,
            "ci_high": hi,
            "risk_high": np.mean(sub_high[col].values) if len(sub_high) else np.nan,
            "risk_low": np.mean(sub_low[col].values) if len(sub_low) else np.nan,
            "n_high": len(sub_high),
            "n_low": len(sub_low)
        })

fragility = pd.DataFrame(frag_rows)
fragility.to_csv(os.path.join(TAB_DIR, "fragility_increase_by_distance_band_with_CI.csv"), index=False)
display(fragility)


# ============================================================
# 9. Main result 5:
#    Aquatic Fragility Buffering, AFB
#    AFB = Fragility_far_or_weak - Fragility_near_permanent
# ============================================================

afb_rows = []

for resp_name, col in response_cols.items():
    near = d[d["water_context"] == "near_permanent_lt1km"]
    far = d[d["water_context"] == "far_or_weak_water"]
    season = d[d["water_context"] == "near_seasonal_lt1km"]

    near_low = near[near["HDSI_group"] == "low"][col].values
    near_high = near[near["HDSI_group"] == "high"][col].values

    far_low = far[far["HDSI_group"] == "low"][col].values
    far_high = far[far["HDSI_group"] == "high"][col].values

    season_low = season[season["HDSI_group"] == "low"][col].values
    season_high = season[season["HDSI_group"] == "high"][col].values

    rng = np.random.default_rng(42)
    boots_far = []
    boots_season = []

    for _ in range(500):
        np_low = rng.choice(near_low, size=len(near_low), replace=True)
        np_high = rng.choice(near_high, size=len(near_high), replace=True)

        fw_low = rng.choice(far_low, size=len(far_low), replace=True)
        fw_high = rng.choice(far_high, size=len(far_high), replace=True)

        ns_low = rng.choice(season_low, size=len(season_low), replace=True)
        ns_high = rng.choice(season_high, size=len(season_high), replace=True)

        frag_near = np.mean(np_high) - np.mean(np_low)
        frag_far = np.mean(fw_high) - np.mean(fw_low)
        frag_season = np.mean(ns_high) - np.mean(ns_low)

        boots_far.append(frag_far - frag_near)
        boots_season.append(frag_season - frag_near)

    afb_far = np.mean(boots_far)
    afb_far_lo, afb_far_hi = np.percentile(boots_far, [2.5, 97.5])

    afb_season = np.mean(boots_season)
    afb_season_lo, afb_season_hi = np.percentile(boots_season, [2.5, 97.5])

    afb_rows.append({
        "response": resp_name,
        "contrast": "far_or_weak fragility - near_permanent fragility",
        "AFB": afb_far,
        "ci_low": afb_far_lo,
        "ci_high": afb_far_hi
    })

    afb_rows.append({
        "response": resp_name,
        "contrast": "near_seasonal fragility - near_permanent fragility",
        "AFB": afb_season,
        "ci_low": afb_season_lo,
        "ci_high": afb_season_hi
    })

afb = pd.DataFrame(afb_rows)
afb.to_csv(os.path.join(TAB_DIR, "aquatic_fragility_buffering_AFB_with_CI.csv"), index=False)
display(afb)


# ============================================================
# 10. Regional results:
#     Region-level risk, fragility, and AFB
# ============================================================

region_summary = (
    d.groupby("macro_region")
    .agg(
        rows=("grid_id", "count"),
        grids=("grid_id", "nunique"),
        HDSI_mean=("HDSI", "mean"),
        HDSI_q75=("HDSI", lambda x: x.quantile(0.75)),
        dist_perm_mean=("dist_perm_water_km", "mean"),
        dist_perm_median=("dist_perm_water_km", "median"),
        EVI_neg_rate=("EVI_neg_z05", "mean"),
        LAI_neg_rate=("LAI_neg_z05", "mean"),
        joint_failure_rate=("joint_failure", "mean")
    )
    .reset_index()
    .sort_values("joint_failure_rate", ascending=False)
)

region_summary.to_csv(os.path.join(TAB_DIR, "region_summary.csv"), index=False)
display(region_summary)

# Region water context high-HDSI risks
rows = []

for region in sorted(d["macro_region"].dropna().unique()):
    dr = d[d["macro_region"] == region].copy()

    for context in d["water_context"].cat.categories:
        sub = dr[
            (dr["water_context"] == context) &
            (dr["HDSI_group"] == "high")
        ].copy()

        for resp_name, col in response_cols.items():
            mean, lo, hi = bootstrap_mean_ci(sub[col].values, n_boot=400, seed=42)
            rows.append({
                "macro_region": region,
                "water_context": context,
                "HDSI_group": "high",
                "response": resp_name,
                "risk": mean,
                "ci_low": lo,
                "ci_high": hi,
                "rows": len(sub),
                "grids": sub["grid_id"].nunique()
            })

region_high_context = pd.DataFrame(rows)
region_high_context.to_csv(os.path.join(TAB_DIR, "region_high_HDSI_water_context_risk_with_CI.csv"), index=False)
display(region_high_context.head(20))

# Region AFB for joint failure
region_afb_rows = []

for region in sorted(d["macro_region"].dropna().unique()):
    dr = d[d["macro_region"] == region].copy()

    near = dr[dr["water_context"] == "near_permanent_lt1km"]
    far = dr[dr["water_context"] == "far_or_weak_water"]

    if len(near) < 100 or len(far) < 100:
        continue

    for resp_name, col in response_cols.items():
        near_low = near[near["HDSI_group"] == "low"][col].values
        near_high = near[near["HDSI_group"] == "high"][col].values
        far_low = far[far["HDSI_group"] == "low"][col].values
        far_high = far[far["HDSI_group"] == "high"][col].values

        if min(len(near_low), len(near_high), len(far_low), len(far_high)) < 30:
            continue

        rng = np.random.default_rng(42)
        boots = []

        for _ in range(400):
            nl = rng.choice(near_low, size=len(near_low), replace=True)
            nh = rng.choice(near_high, size=len(near_high), replace=True)
            fl = rng.choice(far_low, size=len(far_low), replace=True)
            fh = rng.choice(far_high, size=len(far_high), replace=True)

            frag_near = np.mean(nh) - np.mean(nl)
            frag_far = np.mean(fh) - np.mean(fl)
            boots.append(frag_far - frag_near)

        region_afb_rows.append({
            "macro_region": region,
            "response": resp_name,
            "AFB_far_minus_near": np.mean(boots),
            "ci_low": np.percentile(boots, 2.5),
            "ci_high": np.percentile(boots, 97.5),
            "n_near": len(near),
            "n_far": len(far)
        })

region_afb = pd.DataFrame(region_afb_rows)
region_afb.to_csv(os.path.join(TAB_DIR, "region_AFB_far_minus_near_with_CI.csv"), index=False)
display(region_afb)


# ============================================================
# 11. Figure 1:
#     Data overview and framework
# ============================================================

fig = plt.figure(figsize=(11, 7))
gs = fig.add_gridspec(2, 3, width_ratios=[1.1, 1, 1], height_ratios=[1, 1], wspace=0.35, hspace=0.35)

# Panel A: framework
ax = fig.add_subplot(gs[:, 0])
ax.axis("off")
panel_label(ax, "A", x=-0.04, y=1.02)

def rounded_box(ax, xy, w, h, text, fc):
    box = mpl.patches.FancyBboxPatch(
        xy, w, h,
        boxstyle="round,pad=0.02,rounding_size=0.025",
        fc=fc, ec="#333333", lw=0.8
    )
    ax.add_patch(box)
    ax.text(xy[0] + w/2, xy[1] + h/2, text, ha="center", va="center", fontsize=8.5, linespacing=1.25)

def arrow(ax, a, b):
    ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="->", lw=1.2, color="#555555"))

rounded_box(ax, (0.10, 0.78), 0.75, 0.12, "Compound heat–dryness\nHDSI", "#FEE0D2")
rounded_box(ax, (0.10, 0.58), 0.75, 0.12, "Aquatic context\npermanent / seasonal / far-weak", "#D9EAF7")
rounded_box(ax, (0.10, 0.38), 0.75, 0.12, "Vegetation response\nEVI, LAI, joint failure", "#E5F5E0")
rounded_box(ax, (0.10, 0.18), 0.75, 0.12, "Fragility analysis\nRisk(high) - Risk(low)", "#E6E6FA")

arrow(ax, (0.48, 0.78), (0.48, 0.70))
arrow(ax, (0.48, 0.58), (0.48, 0.50))
arrow(ax, (0.48, 0.38), (0.48, 0.30))

ax.text(0.08, 0.96, "Aquatic modulation of\nvegetation fragility", fontsize=12, fontweight="bold", va="top")

# Panel B: HDSI distribution
ax = fig.add_subplot(gs[0, 1])
panel_label(ax, "B")
ax.hist(d["HDSI"], bins=50, color=COL["light_blue"], edgecolor="white")
ax.axvline(q25, color=COL["blue"], linestyle="--", linewidth=1.3, label="Q25")
ax.axvline(q75, color=COL["red"], linestyle="--", linewidth=1.3, label="Q75")
ax.set_xlabel("HDSI")
ax.set_ylabel("Count")
ax.set_title("Heat–dryness stress distribution")
ax.legend(frameon=False)
clean_ax(ax)

# Panel C: permanent water distance distribution
ax = fig.add_subplot(gs[0, 2])
panel_label(ax, "C")
ax.hist(d["dist_perm_water_km"], bins=50, color=COL["blue"], alpha=0.75, edgecolor="white")
ax.set_xlim(0, 16)
ax.set_xlabel("Distance to permanent water (km)")
ax.set_ylabel("Count")
ax.set_title("Permanent-water distance")
clean_ax(ax)

# Panel D: overall negative response rates
ax = fig.add_subplot(gs[1, 1])
panel_label(ax, "D")
vals = [
    d["EVI_neg_z05"].mean(),
    d["LAI_neg_z05"].mean(),
    d["joint_failure"].mean()
]
labels = ["EVI", "LAI", "Joint"]
colors = [COL["green"], COL["purple"], COL["red"]]
ax.bar(labels, vals, color=colors)
ax.set_ylim(0, max(vals) * 1.35)
ax.set_ylabel("Negative-response rate")
ax.set_title("Overall vegetation negative response")
clean_ax(ax)

# Panel E: regional joint failure rate
ax = fig.add_subplot(gs[1, 2])
panel_label(ax, "E")
rs = region_summary.sort_values("joint_failure_rate")
ax.barh(rs["macro_region"], rs["joint_failure_rate"], color=COL["orange"])
ax.set_xlabel("Joint failure rate")
ax.set_title("Regional risk background")
clean_ax(ax)

plt.tight_layout()
savefig(fig, "Fig1_data_framework_overview")
plt.show()


# ============================================================
# 12. Figure 2:
#     High-HDSI risk by water context
# ============================================================

fig, axes = plt.subplots(1, 3, figsize=(11, 4.2), sharey=True)

context_labels = {
    "near_permanent_lt1km": "Near permanent\n<1 km",
    "near_seasonal_lt1km": "Near seasonal\n<1 km",
    "far_or_weak_water": "Far / weak\nwater"
}

for ax, resp_name, label in zip(axes, ["EVI", "LAI", "Joint"], ["A", "B", "C"]):
    sub = water_context_risk[
        (water_context_risk["HDSI_group"] == "high") &
        (water_context_risk["response"] == resp_name)
    ].copy()

    sub["x"] = np.arange(len(sub))
    colors = [CONTEXT_COLORS[c] for c in sub["water_context"]]

    ax.bar(sub["x"], sub["risk"], color=colors, width=0.65)
    yerr = np.vstack([
        sub["risk"] - sub["ci_low"],
        sub["ci_high"] - sub["risk"]
    ])
    ax.errorbar(sub["x"], sub["risk"], yerr=yerr, fmt="none", color="#222222", linewidth=1, capsize=3)

    ax.set_xticks(sub["x"])
    ax.set_xticklabels([context_labels[c] for c in sub["water_context"]])
    ax.set_title(f"{resp_name} risk under high HDSI")
    ax.set_ylabel("Negative-response risk")
    ax.set_ylim(0, max(0.55, sub["ci_high"].max() * 1.25))
    panel_label(ax, label)
    clean_ax(ax)

plt.tight_layout()
savefig(fig, "Fig2_high_HDSI_water_context_risk")
plt.show()


# ============================================================
# 13. Figure 3:
#     Distance-band fragility increase
# ============================================================

fig, ax = plt.subplots(figsize=(9, 5.2))

x = np.arange(len(dist_labels))
offsets = {"EVI": -0.22, "LAI": 0.0, "Joint": 0.22}
colors = {"EVI": COL["green"], "LAI": COL["purple"], "Joint": COL["red"]}

for resp_name in ["EVI", "LAI", "Joint"]:
    sub = fragility[fragility["response"] == resp_name].copy()
    sub["dist_band_perm"] = pd.Categorical(sub["dist_band_perm"], categories=dist_labels, ordered=True)
    sub = sub.sort_values("dist_band_perm")

    xx = x + offsets[resp_name]
    ax.errorbar(
        xx,
        sub["fragility_increase"],
        yerr=np.vstack([
            sub["fragility_increase"] - sub["ci_low"],
            sub["ci_high"] - sub["fragility_increase"]
        ]),
        fmt="o",
        markersize=6,
        linewidth=1.2,
        capsize=3,
        color=colors[resp_name],
        label=resp_name
    )
    ax.plot(xx, sub["fragility_increase"], color=colors[resp_name], linewidth=1.2, alpha=0.75)

ax.axhline(0, color="#333333", linewidth=0.8)
ax.set_xticks(x)
ax.set_xticklabels(dist_labels)
ax.set_xlabel("Distance to permanent water (km)")
ax.set_ylabel("Fragility increase\nRisk(high HDSI) - Risk(low HDSI)")
ax.set_title("Vegetation fragility increases with distance from permanent water")
ax.legend(frameon=False, ncol=3)
clean_ax(ax)

plt.tight_layout()
savefig(fig, "Fig3_distance_band_fragility_increase")
plt.show()


# ============================================================
# 14. Figure 4:
#     Aquatic Fragility Buffering index
# ============================================================

fig, ax = plt.subplots(figsize=(8.5, 4.6))

plot_afb = afb.copy()
plot_afb["label"] = plot_afb["response"] + "\n" + plot_afb["contrast"].str.replace(" fragility - near_permanent fragility", "", regex=False)
plot_afb = plot_afb.sort_values(["response", "contrast"])

x = np.arange(len(plot_afb))
colors = [COL["red"] if "far_or_weak" in c else COL["orange"] for c in plot_afb["contrast"]]

ax.bar(x, plot_afb["AFB"], color=colors, width=0.65)
ax.errorbar(
    x,
    plot_afb["AFB"],
    yerr=np.vstack([
        plot_afb["AFB"] - plot_afb["ci_low"],
        plot_afb["ci_high"] - plot_afb["AFB"]
    ]),
    fmt="none",
    color="#222222",
    linewidth=1,
    capsize=3
)

ax.axhline(0, color="#333333", linewidth=0.8)
ax.set_xticks(x)
ax.set_xticklabels(plot_afb["label"], rotation=30, ha="right")
ax.set_ylabel("Aquatic Fragility Buffering, AFB")
ax.set_title("How much more heat–dryness sensitivity occurs away from permanent water?")
clean_ax(ax)

plt.tight_layout()
savefig(fig, "Fig4_AFB_index")
plt.show()


# ============================================================
# 15. Figure 5:
#     Regional heterogeneity
# ============================================================

# Panel A: regional joint failure rate
# Panel B: regional AFB_joint
# Panel C: high-HDSI water-context risk in key regions

fig = plt.figure(figsize=(11, 7))
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.15], height_ratios=[1, 1], wspace=0.35, hspace=0.35)

ax = fig.add_subplot(gs[0, 0])
panel_label(ax, "A")
rs = region_summary.sort_values("joint_failure_rate")
ax.barh(rs["macro_region"], rs["joint_failure_rate"], color=COL["orange"])
ax.set_xlabel("Joint failure rate")
ax.set_title("Regional background risk")
clean_ax(ax)

ax = fig.add_subplot(gs[1, 0])
panel_label(ax, "B")
reg_joint = region_afb[region_afb["response"] == "Joint"].copy()
reg_joint = reg_joint.sort_values("AFB_far_minus_near")
y = np.arange(len(reg_joint))
ax.errorbar(
    reg_joint["AFB_far_minus_near"],
    y,
    xerr=np.vstack([
        reg_joint["AFB_far_minus_near"] - reg_joint["ci_low"],
        reg_joint["ci_high"] - reg_joint["AFB_far_minus_near"]
    ]),
    fmt="o",
    color=COL["red"],
    capsize=3
)
ax.axvline(0, color="#333333", linewidth=0.8)
ax.set_yticks(y)
ax.set_yticklabels(reg_joint["macro_region"])
ax.set_xlabel("Regional AFB for joint failure")
ax.set_title("Regional aquatic fragility buffering")
clean_ax(ax)

ax = fig.add_subplot(gs[:, 1])
panel_label(ax, "C")
# high-HDSI joint risk by water context per region
rh = region_high_context[
    (region_high_context["response"] == "Joint") &
    (region_high_context["HDSI_group"] == "high")
].copy()

# keep regions with enough rows in all contexts
pivot = rh.pivot_table(
    index="macro_region",
    columns="water_context",
    values="risk"
)

# order by far/weak risk
pivot = pivot.sort_values("far_or_weak_water", ascending=True)

y = np.arange(len(pivot))
for context, offset in zip(["near_permanent_lt1km", "near_seasonal_lt1km", "far_or_weak_water"], [-0.2, 0.0, 0.2]):
    vals = pivot[context].values if context in pivot.columns else np.full(len(pivot), np.nan)
    ax.scatter(vals, y + offset, s=45, color=CONTEXT_COLORS[context], label=context_labels[context].replace("\n", " "))

ax.set_yticks(y)
ax.set_yticklabels(pivot.index)
ax.set_xlabel("High-HDSI joint failure risk")
ax.set_title("Water-context risk varies by region")
ax.legend(frameon=False, loc="lower right")
clean_ax(ax)

plt.tight_layout()
savefig(fig, "Fig5_regional_heterogeneity")
plt.show()


# ============================================================
# 16. Figure 6:
#     Optional map of observed risk and water distance
# ============================================================

china = read_china_boundary()

grid = (
    d.groupby(["grid_id", "lon", "lat"], as_index=False)
    .agg(
        HDSI_mean=("HDSI", "mean"),
        dist_perm_water_km=("dist_perm_water_km", "mean"),
        EVI_neg_rate=("EVI_neg_z05", "mean"),
        LAI_neg_rate=("LAI_neg_z05", "mean"),
        joint_failure_rate=("joint_failure", "mean")
    )
)

def plot_points(ax, data, col, cmap, title, vmin=None, vmax=None):
    if china is not None:
        china.boundary.plot(ax=ax, linewidth=0.5, color="#333333")

    if gpd is not None:
        gdf = gpd.GeoDataFrame(
            data.copy(),
            geometry=gpd.points_from_xy(data["lon"], data["lat"]),
            crs="EPSG:4326"
        )
        gdf.plot(
            ax=ax,
            column=col,
            markersize=1.0,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            legend=True,
            legend_kwds={"shrink": 0.62}
        )
    else:
        sc = ax.scatter(
            data["lon"], data["lat"],
            c=data[col],
            s=1,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax
        )
        plt.colorbar(sc, ax=ax, shrink=0.62)

    ax.set_xlim(70, 138)
    ax.set_ylim(15, 56)
    ax.set_title(title)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.grid(alpha=0.15, lw=0.3)

cmap_risk = mpl.colors.LinearSegmentedColormap.from_list(
    "risk",
    ["#f7fbff", "#c6dbef", "#6baed6", "#fdae61", "#d7191c"]
)
cmap_dist = mpl.colors.LinearSegmentedColormap.from_list(
    "dist",
    ["#08306b", "#2171b5", "#6baed6", "#f7f7f7"]
)

fig, axes = plt.subplots(2, 2, figsize=(10.5, 8))

plot_points(axes[0, 0], grid, "HDSI_mean", cmap_risk, "Mean HDSI")
panel_label(axes[0, 0], "A")

plot_points(axes[0, 1], grid, "dist_perm_water_km", cmap_dist, "Distance to permanent water", vmin=0, vmax=8)
panel_label(axes[0, 1], "B")

plot_points(axes[1, 0], grid, "joint_failure_rate", cmap_risk, "Joint failure frequency", vmin=0, vmax=0.6)
panel_label(axes[1, 0], "C")

plot_points(axes[1, 1], grid, "LAI_neg_rate", cmap_risk, "LAI negative-response frequency", vmin=0, vmax=0.6)
panel_label(axes[1, 1], "D")

plt.tight_layout()
savefig(fig, "Fig6_spatial_patterns")
plt.show()


# ============================================================
# 17. Compact printed summary
# ============================================================

print("\n==================== MAIN RESULT TABLES SAVED ====================")
for fn in sorted(os.listdir(TAB_DIR)):
    print(os.path.join(TAB_DIR, fn))

print("\n==================== MAIN FIGURES SAVED ====================")
for fn in sorted(os.listdir(FIG_DIR)):
    print(os.path.join(FIG_DIR, fn))

print("\nKey expected tables:")
print("1) water_context_HDSI_risk_with_CI.csv")
print("2) high_HDSI_water_context_contrasts_with_CI.csv")
print("3) fragility_increase_by_distance_band_with_CI.csv")
print("4) aquatic_fragility_buffering_AFB_with_CI.csv")
print("5) region_AFB_far_minus_near_with_CI.csv")

print("\nKey expected figures:")
print("Fig1_data_framework_overview.png")
print("Fig2_high_HDSI_water_context_risk.png")
print("Fig3_distance_band_fragility_increase.png")
print("Fig4_AFB_index.png")
print("Fig5_regional_heterogeneity.png")
print("Fig6_spatial_patterns.png")