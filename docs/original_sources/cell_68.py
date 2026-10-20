# ============================================================
# Core analysis and visualization:
# Surface-water permanence, distance decay, and vegetation loss
# under compound heat–dryness events
#
# Main outputs:
#   Fig1: HDSI risk gradient + national map
#   Fig2: Distance-decay of fragility increase
#   Fig3: Water stability and risk amplification
#   Fig4: Lightweight mechanism diagnostics
#   Fig5: Seasonal propagation
#   Fig6: Regional boundary conditions
#
# ============================================================

import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

try:
    from scipy.stats import spearmanr, pearsonr
except Exception:
    spearmanr = None
    pearsonr = None

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

# Optional China boundary, change if needed
CHINA_BOUNDARY_FP = "/content/drive/MyDrive/China_boundary/china_boundary_full.gpkg"

OUT_DIR = os.path.join(AQUA_DIR, "_final_core_story_results")
TAB_DIR = os.path.join(OUT_DIR, "tables")
FIG_DIR = os.path.join(OUT_DIR, "figs")

os.makedirs(TAB_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

print("Panel:", PANEL_FP)
print("Exists:", os.path.exists(PANEL_FP))


# ============================================================
# 1. Plot style
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
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

COL = {
    "blue": "#3B6EA8",
    "orange": "#E69F00",
    "red": "#B2182B",
    "green": "#4E9A51",
    "purple": "#7E62A3",
    "gray": "#666666",
    "dark": "#222222",
    "lightgray": "#D9D9D9",
    "teal": "#009E73",
    "brown": "#8C510A"
}

RESP_COLORS = {
    "EVI": COL["green"],
    "LAI": COL["purple"],
    "Joint": COL["red"],
}

AQUA_COLORS = {
    "near_permanent": COL["blue"],
    "near_seasonal": COL["orange"],
    "other_weak": COL["gray"],
}

STRESS_COLORS = {
    "HDSI": COL["dark"],
    "SM_deficit": COL["blue"],
    "Precip_deficit": COL["teal"],
    "LST_Day": COL["red"],
    "VPD": COL["orange"],
    "Tair": COL["purple"],
}

AQUA_LABELS = {
    "near_permanent": "Near permanent",
    "near_seasonal": "Near seasonal",
    "other_weak": "Other / weak",
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

def save_table(df, name):
    fp = os.path.join(TAB_DIR, name)
    df.to_csv(fp, index=False)
    print("Saved:", fp)
    return fp


# ============================================================
# 2. Load and prepare data
# ============================================================

df = pd.read_parquet(PANEL_FP)
print("Raw panel:", df.shape)

num_cols = [
    "lon", "lat", "year",
    "spring_EVI_z", "spring_LAI_z",
    "late_summer_EVI_z", "late_summer_LAI_z",
    "summer_VPD_z", "summer_LST_Day_z",
    "summer_SM_rootproxy_z", "summer_Precip_z", "summer_Tair_z",
    "elevation", "slope", "LC_Type1",
    "dist_perm_water_km", "dist_seasonal_water_km",
    "dist_any_water_km",
    "water_occurrence", "water_seasonality", "water_transition"
]

for c in num_cols:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

# HDSI
if "HDSI" not in df.columns:
    df["HDSI"] = np.nanmean(
        np.vstack([
            df["summer_VPD_z"].values,
            df["summer_LST_Day_z"].values,
            -df["summer_SM_rootproxy_z"].values,
            -df["summer_Precip_z"].values,
        ]),
        axis=0
    )
else:
    df["HDSI"] = pd.to_numeric(df["HDSI"], errors="coerce")

# Stress components
df["SM_deficit"] = -df["summer_SM_rootproxy_z"]
df["Precip_deficit"] = -df["summer_Precip_z"]
df["LST_Day"] = df["summer_LST_Day_z"]
df["VPD"] = df["summer_VPD_z"]
df["Tair"] = df["summer_Tair_z"]

# Vegetation loss
df["EVI_loss"] = (df["late_summer_EVI_z"] < -0.5).astype(int)
df["LAI_loss"] = (df["late_summer_LAI_z"] < -0.5).astype(int)
df["Joint_loss"] = ((df["EVI_loss"] == 1) | (df["LAI_loss"] == 1)).astype(int)

# Seasonal propagation
df["delta_EVI"] = df["late_summer_EVI_z"] - df["spring_EVI_z"]
df["delta_LAI"] = df["late_summer_LAI_z"] - df["spring_LAI_z"]

# Macro region fallback
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

if "macro_region" not in df.columns:
    df["macro_region"] = df.apply(assign_macro_region, axis=1)
else:
    df["macro_region"] = df["macro_region"].astype(str)

# Valid sample
valid_mask = df["valid_model_both_y"] if "valid_model_both_y" in df.columns else True

required = [
    "grid_id", "lon", "lat", "year", "macro_region",
    "HDSI",
    "SM_deficit", "Precip_deficit", "LST_Day", "VPD", "Tair",
    "dist_perm_water_km", "dist_seasonal_water_km",
    "spring_EVI_z", "spring_LAI_z",
    "late_summer_EVI_z", "late_summer_LAI_z",
    "EVI_loss", "LAI_loss", "Joint_loss",
    "delta_EVI", "delta_LAI"
]

d = df[valid_mask].copy()
d = d.dropna(subset=required).copy()

# Aquatic context
d["aquatic_context"] = "other_weak"
d.loc[d["dist_perm_water_km"] < 1.0, "aquatic_context"] = "near_permanent"
d.loc[
    (d["dist_perm_water_km"] >= 1.0) &
    (d["dist_seasonal_water_km"] < 1.0),
    "aquatic_context"
] = "near_seasonal"

aquatic_order = ["near_permanent", "near_seasonal", "other_weak"]
d["aquatic_context"] = pd.Categorical(d["aquatic_context"], categories=aquatic_order, ordered=True)

# HDSI groups
q25 = d["HDSI"].quantile(0.25)
q75 = d["HDSI"].quantile(0.75)
q90 = d["HDSI"].quantile(0.90)

d["HDSI_group"] = "mid"
d.loc[d["HDSI"] <= q25, "HDSI_group"] = "low"
d.loc[d["HDSI"] >= q75, "HDSI_group"] = "high"
d["HDSI_group"] = pd.Categorical(d["HDSI_group"], categories=["low", "mid", "high"], ordered=True)

# Distance bands
dist_bins = [0, 0.25, 0.5, 1, 2, 3, 5, 7.5, 10, 15, np.inf]
dist_labels = ["0-0.25", "0.25-0.5", "0.5-1", "1-2", "2-3", "3-5", "5-7.5", "7.5-10", "10-15", ">15"]

d["dist_perm_band"] = pd.cut(
    d["dist_perm_water_km"],
    bins=dist_bins,
    labels=dist_labels,
    include_lowest=True,
    right=False
)

dist_mid = {
    "0-0.25": 0.125,
    "0.25-0.5": 0.375,
    "0.5-1": 0.75,
    "1-2": 1.5,
    "2-3": 2.5,
    "3-5": 4.0,
    "5-7.5": 6.25,
    "7.5-10": 8.75,
    "10-15": 12.5,
    ">15": 17.5,
}

loss_vars = {
    "EVI": "EVI_loss",
    "LAI": "LAI_loss",
    "Joint": "Joint_loss",
}

stress_vars = {
    "HDSI": "HDSI",
    "SM_deficit": "SM_deficit",
    "Precip_deficit": "Precip_deficit",
    "LST_Day": "LST_Day",
    "VPD": "VPD",
    "Tair": "Tair",
}

overview = pd.DataFrame([
    {"metric": "rows", "value": len(d)},
    {"metric": "unique_grids", "value": d["grid_id"].nunique()},
    {"metric": "years", "value": ",".join(map(str, sorted(d["year"].astype(int).unique())))},
    {"metric": "HDSI_Q25", "value": q25},
    {"metric": "HDSI_Q75", "value": q75},
    {"metric": "HDSI_Q90", "value": q90},
    {"metric": "EVI_loss_rate", "value": d["EVI_loss"].mean()},
    {"metric": "LAI_loss_rate", "value": d["LAI_loss"].mean()},
    {"metric": "Joint_loss_rate", "value": d["Joint_loss"].mean()},
    {"metric": "near_permanent_rows", "value": int((d["aquatic_context"] == "near_permanent").sum())},
    {"metric": "near_seasonal_rows", "value": int((d["aquatic_context"] == "near_seasonal").sum())},
    {"metric": "other_weak_rows", "value": int((d["aquatic_context"] == "other_weak").sum())},
    {"metric": "dist_perm_mean", "value": d["dist_perm_water_km"].mean()},
    {"metric": "dist_perm_p90", "value": d["dist_perm_water_km"].quantile(0.90)},
    {"metric": "dist_perm_max", "value": d["dist_perm_water_km"].max()},
])

save_table(overview, "00_data_overview.csv")
display(overview)


# ============================================================
# 3. Cluster bootstrap helpers
#    Bootstrap by grid_id to reduce pseudo-replication concern
# ============================================================

def cluster_boot_mean(data, value_col, cluster_col="grid_id", n_boot=300, seed=42):
    dd = data[[cluster_col, value_col]].dropna()
    if len(dd) == 0:
        return np.nan, np.nan, np.nan

    point = dd[value_col].mean()
    clusters = dd[cluster_col].drop_duplicates().values

    if len(clusters) <= 1:
        return point, np.nan, np.nan

    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)

    grouped = {k: v[value_col].values for k, v in dd.groupby(cluster_col)}

    for i in range(n_boot):
        sampled = rng.choice(clusters, size=len(clusters), replace=True)
        vals = np.concatenate([grouped[c] for c in sampled if c in grouped])
        boots[i] = vals.mean() if len(vals) else np.nan

    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return point, lo, hi

def cluster_boot_diff(data_a, data_b, value_col, cluster_col="grid_id", n_boot=300, seed=42):
    # mean(a) - mean(b)
    da = data_a[[cluster_col, value_col]].dropna()
    db = data_b[[cluster_col, value_col]].dropna()

    if len(da) == 0 or len(db) == 0:
        return np.nan, np.nan, np.nan

    point = da[value_col].mean() - db[value_col].mean()

    ca = da[cluster_col].drop_duplicates().values
    cb = db[cluster_col].drop_duplicates().values

    if len(ca) <= 1 or len(cb) <= 1:
        return point, np.nan, np.nan

    ga = {k: v[value_col].values for k, v in da.groupby(cluster_col)}
    gb = {k: v[value_col].values for k, v in db.groupby(cluster_col)}

    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)

    for i in range(n_boot):
        sa = rng.choice(ca, size=len(ca), replace=True)
        sb = rng.choice(cb, size=len(cb), replace=True)
        va = np.concatenate([ga[c] for c in sa if c in ga])
        vb = np.concatenate([gb[c] for c in sb if c in gb])
        boots[i] = va.mean() - vb.mean()

    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return point, lo, hi

def cluster_boot_fragility(data, value_col, low_mask, high_mask, cluster_col="grid_id", n_boot=300, seed=42):
    low = data[low_mask]
    high = data[high_mask]
    return cluster_boot_diff(high, low, value_col, cluster_col=cluster_col, n_boot=n_boot, seed=seed)

def corr_stats(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]

    out = {"pearson_r": np.nan, "pearson_p": np.nan, "spearman_r": np.nan, "spearman_p": np.nan}

    if len(x) >= 3:
        if pearsonr is not None:
            r, p = pearsonr(x, y)
            out["pearson_r"] = r
            out["pearson_p"] = p
        if spearmanr is not None:
            r, p = spearmanr(x, y)
            out["spearman_r"] = r
            out["spearman_p"] = p
    return out


# ============================================================
# 4. Result 1: HDSI identifies vegetation loss risk
# ============================================================

r1_rows = []

for hgroup in ["low", "mid", "high"]:
    sub = d[d["HDSI_group"] == hgroup]

    for resp, col in loss_vars.items():
        m, lo, hi = cluster_boot_mean(sub, col, n_boot=300, seed=42)

        r1_rows.append({
            "HDSI_group": hgroup,
            "response": resp,
            "risk": m,
            "ci_low": lo,
            "ci_high": hi,
            "n": len(sub),
            "grids": sub["grid_id"].nunique()
        })

r1 = pd.DataFrame(r1_rows)
save_table(r1, "01_HDSI_group_loss_risk.csv")


# ============================================================
# 5. Result 2: Distance decay of fragility increase
# ============================================================

dist_rows = []

for band in dist_labels:
    db = d[d["dist_perm_band"].astype(str) == band].copy()
    if len(db) == 0:
        continue

    for resp, col in loss_vars.items():
        low_mask = db["HDSI_group"] == "low"
        high_mask = db["HDSI_group"] == "high"

        risk_low, low_lo, low_hi = cluster_boot_mean(db[low_mask], col, n_boot=300, seed=42)
        risk_high, high_lo, high_hi = cluster_boot_mean(db[high_mask], col, n_boot=300, seed=42)
        frag, frag_lo, frag_hi = cluster_boot_fragility(db, col, low_mask, high_mask, n_boot=300, seed=42)

        dist_rows.append({
            "dist_perm_band": band,
            "dist_mid_km": dist_mid.get(band, np.nan),
            "response": resp,
            "risk_low_HDSI": risk_low,
            "risk_low_ci_low": low_lo,
            "risk_low_ci_high": low_hi,
            "risk_high_HDSI": risk_high,
            "risk_high_ci_low": high_lo,
            "risk_high_ci_high": high_hi,
            "fragility_increase": frag,
            "fragility_ci_low": frag_lo,
            "fragility_ci_high": frag_hi,
            "n_total": len(db),
            "n_low": int(low_mask.sum()),
            "n_high": int(high_mask.sum()),
            "grids": db["grid_id"].nunique()
        })

dist_res = pd.DataFrame(dist_rows)

# Correlation of distance mid vs fragility
corr_rows = []
for resp in ["EVI", "LAI", "Joint"]:
    sub = dist_res[(dist_res["response"] == resp) & (dist_res["dist_perm_band"] != ">15")].dropna(subset=["dist_mid_km", "fragility_increase"])
    cs = corr_stats(sub["dist_mid_km"], sub["fragility_increase"])
    corr_rows.append({
        "response": resp,
        **cs
    })

dist_corr = pd.DataFrame(corr_rows)

save_table(dist_res, "02_distance_decay_fragility.csv")
save_table(dist_corr, "02_distance_decay_correlation.csv")


# ============================================================
# 6. Result 3: Water stability and AFB
# ============================================================

stability_rows = []
afb_rows = []

for aq in aquatic_order:
    daq = d[d["aquatic_context"] == aq]

    for resp, col in loss_vars.items():
        low = daq[daq["HDSI_group"] == "low"]
        high = daq[daq["HDSI_group"] == "high"]

        risk_low, low_lo, low_hi = cluster_boot_mean(low, col, n_boot=300, seed=42)
        risk_high, high_lo, high_hi = cluster_boot_mean(high, col, n_boot=300, seed=42)
        frag, frag_lo, frag_hi = cluster_boot_diff(high, low, col, n_boot=300, seed=42)

        stability_rows.append({
            "aquatic_context": aq,
            "response": resp,
            "risk_low_HDSI": risk_low,
            "risk_low_ci_low": low_lo,
            "risk_low_ci_high": low_hi,
            "risk_high_HDSI": risk_high,
            "risk_high_ci_low": high_lo,
            "risk_high_ci_high": high_hi,
            "fragility_increase": frag,
            "fragility_ci_low": frag_lo,
            "fragility_ci_high": frag_hi,
            "n_low": len(low),
            "n_high": len(high),
            "grids": daq["grid_id"].nunique()
        })

stability_res = pd.DataFrame(stability_rows)

# AFB = fragility(comp) - fragility(near_permanent)
for resp, col in loss_vars.items():
    near = d[d["aquatic_context"] == "near_permanent"]
    near_low = near[near["HDSI_group"] == "low"]
    near_high = near[near["HDSI_group"] == "high"]

    for comp in ["near_seasonal", "other_weak"]:
        dc = d[d["aquatic_context"] == comp]
        comp_low = dc[dc["HDSI_group"] == "low"]
        comp_high = dc[dc["HDSI_group"] == "high"]

        # bootstrap fragility difference
        rng = np.random.default_rng(42)
        boot_vals = []

        # prepare cluster dictionaries
        def prep(dd):
            dd = dd[["grid_id", col]].dropna()
            clusters = dd["grid_id"].drop_duplicates().values
            grouped = {k: v[col].values for k, v in dd.groupby("grid_id")}
            return clusters, grouped

        cl_nl, g_nl = prep(near_low)
        cl_nh, g_nh = prep(near_high)
        cl_cl, g_cl = prep(comp_low)
        cl_ch, g_ch = prep(comp_high)

        if min(len(cl_nl), len(cl_nh), len(cl_cl), len(cl_ch)) <= 1:
            afb, lo, hi = np.nan, np.nan, np.nan
        else:
            for i in range(300):
                snl = rng.choice(cl_nl, size=len(cl_nl), replace=True)
                snh = rng.choice(cl_nh, size=len(cl_nh), replace=True)
                scl = rng.choice(cl_cl, size=len(cl_cl), replace=True)
                sch = rng.choice(cl_ch, size=len(cl_ch), replace=True)

                v_nl = np.concatenate([g_nl[c] for c in snl])
                v_nh = np.concatenate([g_nh[c] for c in snh])
                v_cl = np.concatenate([g_cl[c] for c in scl])
                v_ch = np.concatenate([g_ch[c] for c in sch])

                frag_near = v_nh.mean() - v_nl.mean()
                frag_comp = v_ch.mean() - v_cl.mean()
                boot_vals.append(frag_comp - frag_near)

            boot_vals = np.asarray(boot_vals)
            afb = np.nanmean(boot_vals)
            lo, hi = np.nanpercentile(boot_vals, [2.5, 97.5])

        afb_rows.append({
            "response": resp,
            "contrast": f"{comp} fragility - near_permanent fragility",
            "AFB": afb,
            "ci_low": lo,
            "ci_high": hi
        })

afb_res = pd.DataFrame(afb_rows)

save_table(stability_res, "03_water_stability_risk_and_fragility.csv")
save_table(afb_res, "03_aquatic_fragility_buffering_AFB.csv")


# ============================================================
# 7. Result 4: Lightweight mechanism
#    4a. Stress-channel fragility
# ============================================================

stress_frag_rows = []

for stress_name, stress_col in stress_vars.items():
    s25 = d[stress_col].quantile(0.25)
    s75 = d[stress_col].quantile(0.75)

    for aq in aquatic_order:
        daq = d[d["aquatic_context"] == aq]

        low_mask = daq[stress_col] <= s25
        high_mask = daq[stress_col] >= s75

        for resp, col in loss_vars.items():
            risk_low, _, _ = cluster_boot_mean(daq[low_mask], col, n_boot=300, seed=42)
            risk_high, _, _ = cluster_boot_mean(daq[high_mask], col, n_boot=300, seed=42)
            frag, lo, hi = cluster_boot_diff(daq[high_mask], daq[low_mask], col, n_boot=300, seed=42)

            stress_frag_rows.append({
                "stress": stress_name,
                "stress_q25": s25,
                "stress_q75": s75,
                "aquatic_context": aq,
                "response": resp,
                "risk_low_stress": risk_low,
                "risk_high_stress": risk_high,
                "fragility_increase": frag,
                "ci_low": lo,
                "ci_high": hi,
                "n_low": int(low_mask.sum()),
                "n_high": int(high_mask.sum())
            })

stress_frag = pd.DataFrame(stress_frag_rows)

# Stress-channel AFB
stress_afb_rows = []

for stress_name in stress_vars.keys():
    for resp in ["EVI", "LAI", "Joint"]:
        sub = stress_frag[(stress_frag["stress"] == stress_name) & (stress_frag["response"] == resp)]
        near = sub[sub["aquatic_context"] == "near_permanent"]

        if len(near) == 0:
            continue

        near_frag = near["fragility_increase"].iloc[0]

        for comp in ["near_seasonal", "other_weak"]:
            cc = sub[sub["aquatic_context"] == comp]
            if len(cc) == 0:
                continue

            stress_afb_rows.append({
                "stress": stress_name,
                "response": resp,
                "contrast": f"{comp} - near_permanent",
                "AFB_simple": cc["fragility_increase"].iloc[0] - near_frag,
                "fragility_near_permanent": near_frag,
                "fragility_comparison": cc["fragility_increase"].iloc[0]
            })

stress_afb = pd.DataFrame(stress_afb_rows)

save_table(stress_frag, "04_stress_channel_fragility.csv")
save_table(stress_afb, "04_stress_channel_AFB_simple.csv")


# ============================================================
# 8. Result 4b: Seasonal propagation
# ============================================================

season_rows = []
delta_vars = {
    "delta_EVI": "delta_EVI",
    "delta_LAI": "delta_LAI"
}

for aq in aquatic_order:
    daq = d[d["aquatic_context"] == aq]

    for name, col in delta_vars.items():
        low = daq[daq["HDSI_group"] == "low"]
        high = daq[daq["HDSI_group"] == "high"]

        mean_low, low_lo, low_hi = cluster_boot_mean(low, col, n_boot=300, seed=42)
        mean_high, high_lo, high_hi = cluster_boot_mean(high, col, n_boot=300, seed=42)
        deterioration, det_lo, det_hi = cluster_boot_diff(high, low, col, n_boot=300, seed=42)

        season_rows.append({
            "aquatic_context": aq,
            "delta_variable": name,
            "mean_low_HDSI": mean_low,
            "mean_low_ci_low": low_lo,
            "mean_low_ci_high": low_hi,
            "mean_high_HDSI": mean_high,
            "mean_high_ci_low": high_lo,
            "mean_high_ci_high": high_hi,
            "high_minus_low_delta": deterioration,
            "ci_low": det_lo,
            "ci_high": det_hi,
            "n_low": len(low),
            "n_high": len(high)
        })

season_res = pd.DataFrame(season_rows)

# Compare deterioration: comp - near_permanent
season_diff_rows = []

for name, col in delta_vars.items():
    near = season_res[(season_res["aquatic_context"] == "near_permanent") & (season_res["delta_variable"] == name)]
    near_val = near["high_minus_low_delta"].iloc[0]

    for comp in ["near_seasonal", "other_weak"]:
        cc = season_res[(season_res["aquatic_context"] == comp) & (season_res["delta_variable"] == name)]
        if len(cc) == 0:
            continue
        season_diff_rows.append({
            "delta_variable": name,
            "contrast": f"{comp} deterioration - near_permanent deterioration",
            "difference_simple": cc["high_minus_low_delta"].iloc[0] - near_val,
            "near_permanent_deterioration": near_val,
            "comparison_deterioration": cc["high_minus_low_delta"].iloc[0]
        })

season_diff = pd.DataFrame(season_diff_rows)

save_table(season_res, "05_seasonal_propagation_delta.csv")
save_table(season_diff, "05_seasonal_propagation_difference_simple.csv")


# ============================================================
# 9. Result 4c: High-HDSI stress components
# ============================================================

component_rows = []
component_diff_rows = []

components = {
    "SM_deficit": "Water deficit",
    "Precip_deficit": "Precipitation deficit",
    "LST_Day": "Surface heat",
    "VPD": "Atmospheric dryness",
    "Tair": "Air temperature"
}

dh = d[d["HDSI_group"] == "high"].copy()

for comp, label in components.items():
    for aq in aquatic_order:
        sub = dh[dh["aquatic_context"] == aq]
        m, lo, hi = cluster_boot_mean(sub, comp, n_boot=300, seed=42)

        component_rows.append({
            "component": comp,
            "component_label": label,
            "aquatic_context": aq,
            "mean_high_HDSI": m,
            "ci_low": lo,
            "ci_high": hi,
            "n": len(sub),
            "grids": sub["grid_id"].nunique()
        })

    near = dh[dh["aquatic_context"] == "near_permanent"]
    for comp_aq in ["near_seasonal", "other_weak"]:
        cc = dh[dh["aquatic_context"] == comp_aq]
        diff, lo, hi = cluster_boot_diff(cc, near, comp, n_boot=300, seed=42)

        component_diff_rows.append({
            "component": comp,
            "component_label": label,
            "contrast": f"{comp_aq} - near_permanent",
            "difference": diff,
            "ci_low": lo,
            "ci_high": hi,
            "n_comparison": len(cc),
            "n_near_permanent": len(near)
        })

component_state = pd.DataFrame(component_rows)
component_diff = pd.DataFrame(component_diff_rows)

save_table(component_state, "06_high_HDSI_component_state.csv")
save_table(component_diff, "06_high_HDSI_component_difference.csv")


# ============================================================
# 10. Result 5: Regional boundary conditions
# ============================================================

regional_rows = []
regional_afb_rows = []

for region in sorted(d["macro_region"].dropna().unique()):
    dr = d[d["macro_region"] == region].copy()
    if len(dr) < 1000:
        continue

    for aq in aquatic_order:
        daq = dr[dr["aquatic_context"] == aq]

        for resp, col in loss_vars.items():
            low = daq[daq["HDSI_group"] == "low"]
            high = daq[daq["HDSI_group"] == "high"]

            if len(low) < 50 or len(high) < 50:
                frag, lo, hi = np.nan, np.nan, np.nan
                risk_high = np.nan
            else:
                risk_high, _, _ = cluster_boot_mean(high, col, n_boot=200, seed=42)
                frag, lo, hi = cluster_boot_diff(high, low, col, n_boot=200, seed=42)

            regional_rows.append({
                "macro_region": region,
                "aquatic_context": aq,
                "response": resp,
                "risk_high_HDSI": risk_high,
                "fragility_increase": frag,
                "ci_low": lo,
                "ci_high": hi,
                "n_low": len(low),
                "n_high": len(high),
                "grids": daq["grid_id"].nunique()
            })

regional_res = pd.DataFrame(regional_rows)

for region in sorted(regional_res["macro_region"].dropna().unique()):
    for resp in ["EVI", "LAI", "Joint"]:
        sub = regional_res[(regional_res["macro_region"] == region) & (regional_res["response"] == resp)]
        near = sub[sub["aquatic_context"] == "near_permanent"]

        if len(near) == 0:
            continue

        near_frag = near["fragility_increase"].iloc[0]

        for comp in ["near_seasonal", "other_weak"]:
            cc = sub[sub["aquatic_context"] == comp]
            if len(cc) == 0:
                continue

            regional_afb_rows.append({
                "macro_region": region,
                "response": resp,
                "contrast": f"{comp} - near_permanent",
                "AFB_simple": cc["fragility_increase"].iloc[0] - near_frag,
                "near_permanent_fragility": near_frag,
                "comparison_fragility": cc["fragility_increase"].iloc[0],
                "n_comp_high": cc["n_high"].iloc[0],
                "n_near_high": near["n_high"].iloc[0]
            })

regional_afb = pd.DataFrame(regional_afb_rows)

save_table(regional_res, "07_regional_fragility_by_context.csv")
save_table(regional_afb, "07_regional_AFB_simple.csv")


# ============================================================
# 11. National map summaries
# ============================================================

grid_summary = (
    d.groupby("grid_id", as_index=False)
    .agg(
        lon=("lon", "first"),
        lat=("lat", "first"),
        macro_region=("macro_region", "first"),
        aquatic_context=("aquatic_context", "first"),
        dist_perm_water_km=("dist_perm_water_km", "first"),
        HDSI_mean=("HDSI", "mean"),
        HDSI_q75=("HDSI", lambda x: x.quantile(0.75)),
        high_HDSI_freq=("HDSI_group", lambda x: np.mean(x == "high")),
        EVI_loss_rate=("EVI_loss", "mean"),
        LAI_loss_rate=("LAI_loss", "mean"),
        Joint_loss_rate=("Joint_loss", "mean")
    )
)

save_table(grid_summary, "08_grid_level_map_summary.csv")


# ============================================================
# 12. Optional boundary
# ============================================================

china_boundary = None
if gpd is not None and os.path.exists(CHINA_BOUNDARY_FP):
    try:
        china_boundary = gpd.read_file(CHINA_BOUNDARY_FP)
        china_boundary = china_boundary.to_crs("EPSG:4326")
        print("Loaded China boundary:", china_boundary.shape)
    except Exception as e:
        print("Boundary load failed:", e)
        china_boundary = None


# ============================================================
# 13. Figures
# ============================================================

# ------------------------------------------------------------
# Fig 1: HDSI risk gradient + map
# ------------------------------------------------------------

fig = plt.figure(figsize=(13.5, 7))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.05], width_ratios=[1, 1, 1])

axes_top = [fig.add_subplot(gs[0, i]) for i in range(3)]
ax_map = fig.add_subplot(gs[1, :])

for ax, resp, label in zip(axes_top, ["EVI", "LAI", "Joint"], ["A", "B", "C"]):
    sub = r1[r1["response"] == resp].copy()
    sub["HDSI_group"] = pd.Categorical(sub["HDSI_group"], categories=["low", "mid", "high"], ordered=True)
    sub = sub.sort_values("HDSI_group")

    x = np.arange(len(sub))
    ax.bar(x, sub["risk"], color=RESP_COLORS[resp], width=0.65)
    ax.errorbar(
        x, sub["risk"],
        yerr=np.vstack([sub["risk"] - sub["ci_low"], sub["ci_high"] - sub["risk"]]),
        fmt="none", color="#222222", capsize=3, linewidth=1
    )
    ax.set_xticks(x)
    ax.set_xticklabels(["Low", "Mid", "High"])
    ax.set_ylabel("Loss risk")
    ax.set_title(f"{resp} loss increases with HDSI")
    panel_label(ax, label)
    clean_ax(ax)

# Map: grid-level high-HDSI frequency
if china_boundary is not None:
    china_boundary.boundary.plot(ax=ax_map, linewidth=0.5, color="#555555")

sc = ax_map.scatter(
    grid_summary["lon"],
    grid_summary["lat"],
    c=grid_summary["high_HDSI_freq"],
    s=3,
    cmap="YlOrRd",
    alpha=0.75,
    linewidths=0
)

ax_map.set_title("D. National distribution of high-HDSI exposure frequency")
ax_map.set_xlabel("Longitude")
ax_map.set_ylabel("Latitude")
clean_ax(ax_map)
cb = plt.colorbar(sc, ax=ax_map, fraction=0.025, pad=0.01)
cb.set_label("High-HDSI frequency")

plt.tight_layout()
savefig(fig, "Fig1_HDSI_risk_gradient_and_exposure_map")
plt.show()


# ------------------------------------------------------------
# Fig 2: Distance decay
# ------------------------------------------------------------

fig, axes = plt.subplots(1, 3, figsize=(14, 4.6), sharey=False)

for ax, resp, label in zip(axes, ["EVI", "LAI", "Joint"], ["A", "B", "C"]):
    sub = dist_res[(dist_res["response"] == resp) & (dist_res["dist_perm_band"] != ">15")].copy()
    sub["dist_perm_band"] = pd.Categorical(sub["dist_perm_band"], categories=dist_labels[:-1], ordered=True)
    sub = sub.sort_values("dist_perm_band")

    x = sub["dist_mid_km"].values
    y = sub["fragility_increase"].values

    ax.errorbar(
        x, y,
        yerr=np.vstack([y - sub["fragility_ci_low"].values, sub["fragility_ci_high"].values - y]),
        fmt="o",
        color=RESP_COLORS[resp],
        ecolor="#222222",
        elinewidth=1,
        capsize=3,
        markersize=5,
        label="Distance-band estimate"
    )

    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() >= 3:
        coef = np.polyfit(x[ok], y[ok], deg=1)
        xs = np.linspace(np.nanmin(x[ok]), np.nanmax(x[ok]), 100)
        ys = coef[0] * xs + coef[1]
        ax.plot(xs, ys, color="#222222", lw=1.2, ls="--", label="Linear trend")

        cs = corr_stats(x[ok], y[ok])
        txt = f"Spearman r={cs['spearman_r']:.2f}"
        ax.text(0.04, 0.93, txt, transform=ax.transAxes, ha="left", va="top", fontsize=8)

    ax.set_xlabel("Distance to permanent water (km)")
    ax.set_ylabel("Fragility increase\nRisk(high HDSI) − Risk(low HDSI)")
    ax.set_title(f"{resp}: distance-decay of aquatic buffering")
    panel_label(ax, label)
    clean_ax(ax)
    ax.legend(frameon=False, loc="lower right", fontsize=7)

plt.tight_layout()
savefig(fig, "Fig2_distance_decay_fragility_with_trend")
plt.show()


# ------------------------------------------------------------
# Fig 3: Water stability
# ------------------------------------------------------------

fig, axes = plt.subplots(2, 3, figsize=(14, 8.2), sharey="row")

for j, resp in enumerate(["EVI", "LAI", "Joint"]):
    # High-HDSI risk
    ax = axes[0, j]
    sub = stability_res[stability_res["response"] == resp].copy()
    sub["aquatic_context"] = pd.Categorical(sub["aquatic_context"], categories=aquatic_order, ordered=True)
    sub = sub.sort_values("aquatic_context")

    x = np.arange(len(sub))
    colors = [AQUA_COLORS[a] for a in sub["aquatic_context"].astype(str)]

    ax.bar(x, sub["risk_high_HDSI"], color=colors, width=0.65)
    ax.errorbar(
        x, sub["risk_high_HDSI"],
        yerr=np.vstack([sub["risk_high_HDSI"] - sub["risk_high_ci_low"], sub["risk_high_ci_high"] - sub["risk_high_HDSI"]]),
        fmt="none", color="#222222", capsize=3, linewidth=1
    )
    ax.set_xticks(x)
    ax.set_xticklabels([AQUA_LABELS[a] for a in sub["aquatic_context"].astype(str)], rotation=15, ha="right")
    ax.set_title(f"{resp}: high-HDSI loss risk")
    ax.set_ylabel("Loss risk")
    panel_label(ax, ["A", "B", "C"][j])
    clean_ax(ax)

    # Fragility increase
    ax = axes[1, j]
    ax.bar(x, sub["fragility_increase"], color=colors, width=0.65)
    ax.errorbar(
        x, sub["fragility_increase"],
        yerr=np.vstack([sub["fragility_increase"] - sub["fragility_ci_low"], sub["fragility_ci_high"] - sub["fragility_increase"]]),
        fmt="none", color="#222222", capsize=3, linewidth=1
    )
    ax.axhline(0, color="#333333", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([AQUA_LABELS[a] for a in sub["aquatic_context"].astype(str)], rotation=15, ha="right")
    ax.set_title(f"{resp}: risk amplification")
    ax.set_ylabel("Risk(high) − Risk(low)")
    panel_label(ax, ["D", "E", "F"][j])
    clean_ax(ax)

plt.tight_layout()
savefig(fig, "Fig3_water_stability_high_risk_and_fragility")
plt.show()


# ------------------------------------------------------------
# Fig 4: Lightweight mechanism
# ------------------------------------------------------------

fig, axes = plt.subplots(1, 3, figsize=(14, 4.8))

# A: stress-channel AFB for other_weak - near_permanent
ax = axes[0]
plot = stress_afb[(stress_afb["contrast"] == "other_weak - near_permanent") & (stress_afb["response"] == "Joint")].copy()
stress_order = ["HDSI", "SM_deficit", "Precip_deficit", "LST_Day", "VPD", "Tair"]
plot["stress"] = pd.Categorical(plot["stress"], categories=stress_order, ordered=True)
plot = plot.sort_values("stress")
x = np.arange(len(plot))
colors = [STRESS_COLORS[s] for s in plot["stress"].astype(str)]
ax.bar(x, plot["AFB_simple"], color=colors, width=0.65)
ax.axhline(0, color="#333333", lw=0.8)
ax.set_xticks(x)
ax.set_xticklabels(plot["stress"], rotation=35, ha="right")
ax.set_ylabel("AFB simple\nFragility(other/weak) − Fragility(near permanent)")
ax.set_title("Joint loss: pathway-specific fragility")
panel_label(ax, "A")
clean_ax(ax)

# B: high-HDSI component difference
ax = axes[1]
plot = component_diff[component_diff["contrast"] == "other_weak - near_permanent"].copy()
plot["component"] = pd.Categorical(plot["component"], categories=["SM_deficit", "Precip_deficit", "LST_Day", "VPD", "Tair"], ordered=True)
plot = plot.sort_values("component")
x = np.arange(len(plot))
colors = [STRESS_COLORS.get(s, COL["gray"]) for s in plot["component"].astype(str)]
ax.bar(x, plot["difference"], color=colors, width=0.65)
ax.errorbar(
    x, plot["difference"],
    yerr=np.vstack([plot["difference"] - plot["ci_low"], plot["ci_high"] - plot["difference"]]),
    fmt="none", color="#222222", capsize=3, linewidth=1
)
ax.axhline(0, color="#333333", lw=0.8)
ax.set_xticks(x)
ax.set_xticklabels(plot["component"], rotation=35, ha="right")
ax.set_ylabel("Component difference\nOther/weak − near permanent")
ax.set_title("High-HDSI stress components")
panel_label(ax, "B")
clean_ax(ax)

# C: seasonal propagation
ax = axes[2]
plot = season_res[season_res["delta_variable"] == "delta_LAI"].copy()
plot["aquatic_context"] = pd.Categorical(plot["aquatic_context"], categories=aquatic_order, ordered=True)
plot = plot.sort_values("aquatic_context")
x = np.arange(len(plot))
colors = [AQUA_COLORS[a] for a in plot["aquatic_context"].astype(str)]
ax.bar(x, plot["high_minus_low_delta"], color=colors, width=0.65)
ax.errorbar(
    x, plot["high_minus_low_delta"],
    yerr=np.vstack([plot["high_minus_low_delta"] - plot["ci_low"], plot["ci_high"] - plot["high_minus_low_delta"]]),
    fmt="none", color="#222222", capsize=3, linewidth=1
)
ax.axhline(0, color="#333333", lw=0.8)
ax.set_xticks(x)
ax.set_xticklabels([AQUA_LABELS[a] for a in plot["aquatic_context"].astype(str)], rotation=15, ha="right")
ax.set_ylabel("ΔLAI response\nHigh HDSI − Low HDSI")
ax.set_title("Seasonal propagation of LAI deterioration")
panel_label(ax, "C")
clean_ax(ax)

plt.tight_layout()
savefig(fig, "Fig4_lightweight_mechanism_pathways")
plt.show()


# ------------------------------------------------------------
# Fig 5: Regional boundary conditions
# ------------------------------------------------------------

fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharex=True)

for ax, resp, label in zip(axes, ["EVI", "LAI", "Joint"], ["A", "B", "C"]):
    plot = regional_afb[
        (regional_afb["response"] == resp) &
        (regional_afb["contrast"] == "other_weak - near_permanent")
    ].copy()

    plot = plot.dropna(subset=["AFB_simple"])
    plot = plot.sort_values("AFB_simple")
    y = np.arange(len(plot))

    ax.scatter(plot["AFB_simple"], y, color=RESP_COLORS[resp], s=45)
    ax.axvline(0, color="#333333", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(plot["macro_region"])
    ax.set_xlabel("Regional AFB\nFragility(other/weak) − Fragility(near permanent)")
    ax.set_title(f"{resp}: regional boundary conditions")
    panel_label(ax, label)
    clean_ax(ax)

plt.tight_layout()
savefig(fig, "Fig5_regional_boundary_conditions_AFB")
plt.show()


# ------------------------------------------------------------
# Fig 6: National map of aquatic context and Joint loss
# ------------------------------------------------------------

fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.5))

# A aquatic context
ax = axes[0]
if china_boundary is not None:
    china_boundary.boundary.plot(ax=ax, linewidth=0.5, color="#555555")

for aq in aquatic_order:
    sub = grid_summary[grid_summary["aquatic_context"].astype(str) == aq]
    ax.scatter(
        sub["lon"], sub["lat"],
        s=3,
        color=AQUA_COLORS[aq],
        alpha=0.65,
        linewidths=0,
        label=AQUA_LABELS[aq]
    )

ax.set_title("A. Aquatic context")
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
ax.legend(frameon=False, loc="lower left", fontsize=7)
clean_ax(ax)

# B joint loss rate
ax = axes[1]
if china_boundary is not None:
    china_boundary.boundary.plot(ax=ax, linewidth=0.5, color="#555555")

sc = ax.scatter(
    grid_summary["lon"], grid_summary["lat"],
    c=grid_summary["Joint_loss_rate"],
    s=3,
    cmap="Reds",
    alpha=0.75,
    linewidths=0
)
ax.set_title("B. Grid-level Joint loss rate")
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
clean_ax(ax)
cb = plt.colorbar(sc, ax=ax, fraction=0.03, pad=0.01)
cb.set_label("Joint loss rate")

plt.tight_layout()
savefig(fig, "Fig6_national_maps_aquatic_context_and_joint_loss")
plt.show()


# ============================================================
# 14. Print outputs
# ============================================================

print("\n==================== TABLES SAVED ====================")
for fn in sorted(os.listdir(TAB_DIR)):
    print(os.path.join(TAB_DIR, fn))

print("\n==================== FIGURES SAVED ====================")
for fn in sorted(os.listdir(FIG_DIR)):
    print(os.path.join(FIG_DIR, fn))

print("\nKey tables to inspect:")
print("01_HDSI_group_loss_risk.csv")
print("02_distance_decay_fragility.csv")
print("02_distance_decay_correlation.csv")
print("03_water_stability_risk_and_fragility.csv")
print("03_aquatic_fragility_buffering_AFB.csv")
print("04_stress_channel_AFB_simple.csv")
print("05_seasonal_propagation_delta.csv")
print("06_high_HDSI_component_difference.csv")
print("07_regional_AFB_simple.csv")
print("08_grid_level_map_summary.csv")