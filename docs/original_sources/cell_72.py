# ============================================================
# Complete plotting code for:
#   Fig 4 (mechanism diagnosis) – legends below subplots
#   Fig 5 (regional boundary conditions)
#   Threshold sensitivity figure
# Directly shows figures; does NOT save.
# ============================================================

import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

try:
    import seaborn as sns
    HAS_SNS = True
except ImportError:
    HAS_SNS = False

# ------------------------------------------------------------
# 0. Paths and data loading
# ------------------------------------------------------------
BASE_DIR = "/content/drive/MyDrive/China_VegClimate_ResponseOperator_Dataset_patched/_aquatic_buffering"
PANEL_FP = os.path.join(BASE_DIR, "china_aquatic_buffering_panel_v2.parquet")

df = pd.read_parquet(PANEL_FP)
print("Raw panel:", df.shape)

# ------------------------------------------------------------
# 1. Global figure style (unified fonts)
# ------------------------------------------------------------
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9.5,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8.5,
    "axes.linewidth": 0.8,
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# Color palette
COL = {
    "near_permanent": "#3B6EA8",
    "near_seasonal": "#E69F00",
    "other_weak": "#666666",
    "EVI": "#4E9A51",
    "LAI": "#7E62A3",
    "Joint": "#B2182B",
    "low": "#72A9D5",
    "high": "#D95F5F",
    "dark": "#222222",
    "gray": "#666666",
}

ctx_order = ["near_permanent", "near_seasonal", "other_weak"]
ctx_label = {
    "near_permanent": "Near permanent",
    "near_seasonal": "Near seasonal",
    "other_weak": "Other / weak",
}

def clean_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.16, linewidth=0.45)

def panel_label(ax, lab, x=-0.08, y=1.04, fontsize=12):
    ax.text(
        x, y, lab,
        transform=ax.transAxes,
        fontsize=fontsize,
        fontweight="bold",
        ha="left",
        va="top"
    )

def boot_mean_ci(x, n_boot=400, seed=42):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return np.nan, np.nan, np.nan
    m = x.mean()
    if len(x) <= 2:
        return m, np.nan, np.nan
    rng = np.random.default_rng(seed)
    boots = [rng.choice(x, size=len(x), replace=True).mean() for _ in range(n_boot)]
    return m, np.percentile(boots, 2.5), np.percentile(boots, 97.5)

def boot_diff_ci(a, b, n_boot=400, seed=42):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) == 0 or len(b) == 0:
        return np.nan, np.nan, np.nan
    d0 = a.mean() - b.mean()
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        aa = rng.choice(a, size=len(a), replace=True)
        bb = rng.choice(b, size=len(b), replace=True)
        boots.append(aa.mean() - bb.mean())
    return d0, np.percentile(boots, 2.5), np.percentile(boots, 97.5)

def sample_for_plot(x, n=2500, seed=42):
    if len(x) <= n:
        return x
    return x.sample(n=n, random_state=seed)

# ------------------------------------------------------------
# 2. Prepare variables (same as original)
# ------------------------------------------------------------
num_cols = [
    "lon", "lat", "year", "elevation", "slope",
    "HDSI",
    "dist_perm_water_km", "dist_seasonal_water_km",
    "spring_EVI_z", "spring_LAI_z",
    "late_summer_EVI_z", "late_summer_LAI_z",
    "summer_SM_rootproxy_z", "summer_Precip_z",
    "summer_LST_Day_z", "summer_VPD_z", "summer_Tair_z"
]

for c in num_cols:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

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

# Aquatic context
df["aquatic_context"] = "other_weak"
df.loc[df["dist_perm_water_km"] < 1.0, "aquatic_context"] = "near_permanent"
df.loc[
    (df["dist_perm_water_km"] >= 1.0) & (df["dist_seasonal_water_km"] < 1.0),
    "aquatic_context"
] = "near_seasonal"

# HDSI groups
q25 = df["HDSI"].quantile(0.25)
q75 = df["HDSI"].quantile(0.75)
df["HDSI_group"] = "mid"
df.loc[df["HDSI"] <= q25, "HDSI_group"] = "low"
df.loc[df["HDSI"] >= q75, "HDSI_group"] = "high"

# Loss events
df["EVI_loss"] = (df["late_summer_EVI_z"] < -0.5).astype(int)
df["LAI_loss"] = (df["late_summer_LAI_z"] < -0.5).astype(int)
df["Joint_loss"] = ((df["EVI_loss"] == 1) | (df["LAI_loss"] == 1)).astype(int)

# Seasonal deterioration
df["delta_EVI"] = df["late_summer_EVI_z"] - df["spring_EVI_z"]
df["delta_LAI"] = df["late_summer_LAI_z"] - df["spring_LAI_z"]

# Stress components (larger = worse)
df["SM_deficit"] = -df["summer_SM_rootproxy_z"]
df["Precip_deficit"] = -df["summer_Precip_z"]
df["LST_Day"] = df["summer_LST_Day_z"]
df["VPD"] = df["summer_VPD_z"]
df["Tair"] = df["summer_Tair_z"]

def assign_macro_region(row):
    lon = row.get("lon", np.nan)
    lat = row.get("lat", np.nan)
    elev = row.get("elevation", np.nan)
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

required = [
    "HDSI", "HDSI_group", "aquatic_context",
    "late_summer_EVI_z", "late_summer_LAI_z",
    "spring_EVI_z", "spring_LAI_z",
    "dist_perm_water_km", "dist_seasonal_water_km",
    "SM_deficit", "Precip_deficit", "LST_Day", "VPD", "Tair",
    "EVI_loss", "LAI_loss", "Joint_loss", "delta_LAI",
    "macro_region"
]

d = df.dropna(subset=required).copy()
d = d[d["aquatic_context"].isin(ctx_order)].copy()
print("Valid sample:", d.shape)

# ------------------------------------------------------------
# 3. Figure 4: Mechanism diagnosis (legends below subplots)
# ------------------------------------------------------------
comp_order = ["SM_deficit", "Precip_deficit", "LST_Day", "VPD", "Tair"]
comp_label = {
    "SM_deficit": "Water deficit",
    "Precip_deficit": "Precipitation deficit",
    "LST_Day": "Surface heat",
    "VPD": "Atmospheric dryness",
    "Tair": "Air temperature",
}

high = d[d["HDSI_group"] == "high"].copy()

# ---- Panel A data: long format for violin plots ----
long_list = []
for comp in comp_order:
    tmp = high[["aquatic_context", comp]].dropna().copy()
    tmp = tmp.rename(columns={comp: "value"})
    tmp["component"] = comp
    tmp["component_label"] = comp_label[comp]
    long_list.append(tmp)
long_high = pd.concat(long_list, axis=0, ignore_index=True)
long_high_plot = (
    long_high.groupby(["component", "aquatic_context"], group_keys=False)
    .apply(lambda x: sample_for_plot(x, n=2200))
    .reset_index(drop=True)
)

# ---- Panel B: component differences ----
component_diff_rows = []
for comp in comp_order:
    near = high.loc[high["aquatic_context"] == "near_permanent", comp].dropna().values
    for group in ["near_seasonal", "other_weak"]:
        vals = high.loc[high["aquatic_context"] == group, comp].dropna().values
        est, lo, hi = boot_diff_ci(vals, near, n_boot=500)
        component_diff_rows.append({
            "component": comp,
            "component_label": comp_label[comp],
            "group": group,
            "contrast": f"{ctx_label[group]} − Near permanent",
            "difference": est,
            "ci_low": lo,
            "ci_high": hi,
        })
component_diff = pd.DataFrame(component_diff_rows)

# ---- Panel C: delta LAI ----
dlai_plot = (
    d[d["HDSI_group"].isin(["low", "high"])]
    [["aquatic_context", "HDSI_group", "delta_LAI"]]
    .dropna()
    .groupby(["aquatic_context", "HDSI_group"], group_keys=False)
    .apply(lambda x: sample_for_plot(x, n=2200))
    .reset_index(drop=True)
)

# ---- Panel D: risk propagation ----
risk_rows = []
for event_col, event_name in [("LAI_loss", "LAI loss"), ("Joint_loss", "Joint loss")]:
    for ctx in ctx_order:
        sub = d[d["aquatic_context"] == ctx]
        low_vals = sub.loc[sub["HDSI_group"] == "low", event_col].values
        high_vals = sub.loc[sub["HDSI_group"] == "high", event_col].values
        low_m, low_lo, low_hi = boot_mean_ci(low_vals, n_boot=500)
        high_m, high_lo, high_hi = boot_mean_ci(high_vals, n_boot=500)
        risk_rows.append({
            "event": event_name,
            "aquatic_context": ctx,
            "low_risk": low_m,
            "low_ci_low": low_lo,
            "low_ci_high": low_hi,
            "high_risk": high_m,
            "high_ci_low": high_lo,
            "high_ci_high": high_hi,
        })
risk_prop = pd.DataFrame(risk_rows)

# ---- Draw Figure 4 with legends below each subplot ----
fig = plt.figure(figsize=(17.5, 13.0), constrained_layout=False)  # taller to fit legends
gs = fig.add_gridspec(2, 2, height_ratios=[1.12, 1.2],
                      width_ratios=[1.10, 0.90], hspace=0.50, wspace=0.30)

# Panel A
axA = fig.add_subplot(gs[0, 0])
if HAS_SNS:
    sns.violinplot(
        data=long_high_plot,
        x="value", y="component_label",
        hue="aquatic_context",
        order=[comp_label[c] for c in comp_order],
        hue_order=ctx_order,
        palette=[COL[c] for c in ctx_order],
        cut=0, inner="quartile", linewidth=0.75, scale="width",
        ax=axA
    )
    # Legend below the subplot
    axA.legend(
        loc="upper center", bbox_to_anchor=(0.5, -0.18),
        title="Aquatic context", frameon=False, ncol=3
    )
axA.axvline(0, color="black", linestyle="--", linewidth=0.9)
axA.set_title("A. High‑HDSI stress‑component distributions", loc="left", pad=8)
axA.set_xlabel("Standardized anomaly; larger = more stressful")
panel_label(axA, "A", fontsize=12)
clean_ax(axA)

# Panel B
axB = fig.add_subplot(gs[0, 1])
plotB = component_diff.copy()
plotB["component_label"] = pd.Categorical(plotB["component_label"],
                                          categories=[comp_label[c] for c in comp_order],
                                          ordered=True)
plotB = plotB.sort_values(["component_label", "group"])
y_pos, y_labs = [], []
y = 0
for comp_lab in [comp_label[c] for c in comp_order]:
    sub = plotB[plotB["component_label"] == comp_lab]
    for _, row in sub.iterrows():
        color = COL[row["group"]]
        axB.errorbar(row["difference"], y,
                     xerr=[[row["difference"] - row["ci_low"]],
                           [row["ci_high"] - row["difference"]]],
                     fmt="o", color=color, ecolor=color,
                     capsize=3.5, markersize=6, linewidth=1.7)
        y_pos.append(y)
        y_labs.append(comp_lab if row["group"] == "near_seasonal" else "")
        y += 1
    y += 0.6
axB.axvline(0, color="black", linestyle="--", linewidth=0.9)
axB.set_yticks(y_pos)
axB.set_yticklabels(y_labs)
axB.set_xlabel("Difference from near‑permanent water")
axB.set_title("B. Exposure contrasts under high HDSI", loc="left", pad=8)
panel_label(axB, "B", fontsize=12)
clean_ax(axB)
# Legend below Panel B
axB.legend(
    handles=[
        Line2D([0], [0], marker="o", color=COL["near_seasonal"],
               label="Near seasonal − Near permanent"),
        Line2D([0], [0], marker="o", color=COL["other_weak"],
               label="Other/weak − Near permanent"),
    ],
    title="Contrast", frameon=False,
    loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2
)

# Panel C
axC = fig.add_subplot(gs[1, 0])
if HAS_SNS:
    sns.violinplot(
        data=dlai_plot, x="aquatic_context", y="delta_LAI",
        hue="HDSI_group", order=ctx_order, hue_order=["low", "high"],
        palette=[COL["low"], COL["high"]], cut=0, inner=None,
        linewidth=0.75, scale="width", ax=axC
    )
    sns.boxplot(
        data=dlai_plot, x="aquatic_context", y="delta_LAI",
        hue="HDSI_group", order=ctx_order, hue_order=["low", "high"],
        width=0.22, showcaps=True,
        boxprops={"facecolor": "white", "alpha": 0.85, "linewidth": 0.7},
        whiskerprops={"linewidth": 0.7}, medianprops={"color": "black", "linewidth": 1.0},
        showfliers=False, ax=axC
    )
    # Legend below Panel C
    axC.legend(
        ["Low HDSI", "High HDSI"], title="Stress level", frameon=False,
        loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2
    )
axC.axhline(0, color="black", linestyle="--", linewidth=0.9)
axC.set_xticks(range(len(ctx_order)))
axC.set_xticklabels([ctx_label[c] for c in ctx_order])
axC.set_ylabel("ΔLAI = late‑summer LAI_z − spring LAI_z")
axC.set_title("C. Seasonal LAI deterioration", loc="left", pad=8)
panel_label(axC, "C", fontsize=12)
clean_ax(axC)

# Panel D
axD = fig.add_subplot(gs[1, 1])
plotD = risk_prop.copy()
y = 0
yticks, yticklabels = [], []
for event in ["LAI loss", "Joint loss"]:
    sub_event = plotD[plotD["event"] == event]
    for ctx in ctx_order:
        row = sub_event[sub_event["aquatic_context"] == ctx].iloc[0]
        axD.plot([row["low_risk"], row["high_risk"]], [y, y],
                 color=COL[ctx], linewidth=2.4, alpha=0.82)
        axD.errorbar(row["low_risk"], y,
                     xerr=[[row["low_risk"] - row["low_ci_low"]],
                           [row["low_ci_high"] - row["low_risk"]]],
                     fmt="o", color=COL["low"], ecolor=COL["low"],
                     capsize=3, markersize=5.8)
        axD.errorbar(row["high_risk"], y,
                     xerr=[[row["high_risk"] - row["high_ci_low"]],
                           [row["high_ci_high"] - row["high_risk"]]],
                     fmt="s", color=COL["high"], ecolor=COL["high"],
                     capsize=3, markersize=5.8)
        yticks.append(y)
        yticklabels.append(f"{event} | {ctx_label[ctx]}")
        y += 1
    y += 0.72
axD.set_yticks(yticks)
axD.set_yticklabels(yticklabels)
axD.set_xlabel("Loss probability")
axD.set_title("D. HDSI‑to‑loss propagation", loc="left", pad=8)
panel_label(axD, "D", fontsize=12)
clean_ax(axD)
# Legend below Panel D
axD.legend(
    handles=[
        Line2D([0], [0], marker="o", color=COL["low"], label="Low HDSI", linestyle="None", markersize=6.5),
        Line2D([0], [0], marker="s", color=COL["high"], label="High HDSI", linestyle="None", markersize=6.5),
    ],
    title="Stress level", frameon=False,
    loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2
)

fig.suptitle("Figure 4. Mechanism diagnostics: exposure contrast and stress‑to‑loss propagation",
             fontsize=14, fontweight="bold", y=0.98)

# Adjust margins to make room for bottom legends
plt.subplots_adjust(left=0.10, right=0.95, top=0.92, bottom=0.12)
plt.show()

# ------------------------------------------------------------
# 4. Figure 5: Regional boundary conditions (no legend needed)
# ------------------------------------------------------------
def regional_afb(sub, event_col, n_boot=500, seed=42):
    near = sub[sub["aquatic_context"] == "near_permanent"]
    other = sub[sub["aquatic_context"] == "other_weak"]
    near_low = near.loc[near["HDSI_group"] == "low", event_col].dropna().values
    near_high = near.loc[near["HDSI_group"] == "high", event_col].dropna().values
    other_low = other.loc[other["HDSI_group"] == "low", event_col].dropna().values
    other_high = other.loc[other["HDSI_group"] == "high", event_col].dropna().values
    if min(len(near_low), len(near_high), len(other_low), len(other_high)) < 20:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        nl = rng.choice(near_low, size=len(near_low), replace=True)
        nh = rng.choice(near_high, size=len(near_high), replace=True)
        ol = rng.choice(other_low, size=len(other_low), replace=True)
        oh = rng.choice(other_high, size=len(other_high), replace=True)
        frag_near = nh.mean() - nl.mean()
        frag_other = oh.mean() - ol.mean()
        boots.append(frag_other - frag_near)
    boots = np.asarray(boots)
    return boots.mean(), np.percentile(boots, 2.5), np.percentile(boots, 97.5)

region_rows = []
for region, sub in d.groupby("macro_region"):
    if len(sub) < 1000:
        continue
    for event_col, response in [("EVI_loss", "EVI"), ("LAI_loss", "LAI"), ("Joint_loss", "Joint")]:
        est, lo, hi = regional_afb(sub, event_col, n_boot=500)
        if not np.isnan(est):
            region_rows.append({
                "region": region,
                "response": response,
                "AFB": est,
                "ci_low": lo,
                "ci_high": hi
            })
regional_afb_df = pd.DataFrame(region_rows)

fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.8), sharey=False, constrained_layout=False)
for ax, response, color, lab in zip(
    axes, ["EVI", "LAI", "Joint"],
    [COL["EVI"], COL["LAI"], COL["Joint"]],
    ["A", "B", "C"]
):
    sub = regional_afb_df[regional_afb_df["response"] == response].dropna(subset=["AFB"]).copy()
    sub = sub.sort_values("AFB")
    y = np.arange(len(sub))
    ax.errorbar(sub["AFB"], y,
                xerr=[sub["AFB"] - sub["ci_low"], sub["ci_high"] - sub["AFB"]],
                fmt="o", color=color, ecolor=color,
                capsize=3.5, markersize=6, linewidth=1.8)
    ax.axvline(0, color="black", linestyle="--", linewidth=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(sub["region"])
    ax.set_title(f"{response} regional AFB", loc="left", pad=8)
    ax.set_xlabel("AFB = Fragility(other/weak) − Fragility(near permanent)")
    panel_label(ax, lab, fontsize=12)
    clean_ax(ax)

fig.suptitle("Figure 5. Regional boundary conditions of aquatic buffering",
             fontsize=14, fontweight="bold", y=0.98)
plt.subplots_adjust(left=0.08, right=0.98, top=0.88, bottom=0.14, wspace=0.36)
plt.show()

# ------------------------------------------------------------
# 5. Threshold sensitivity figure (no legend needed)
# ------------------------------------------------------------
thresholds = [-0.25, -0.5, -0.75, -1.0]
sens_rows = []
for thr in thresholds:
    tmp = d.copy()
    tmp["EVI_evt"] = (tmp["late_summer_EVI_z"] < thr).astype(int)
    tmp["LAI_evt"] = (tmp["late_summer_LAI_z"] < thr).astype(int)
    tmp["Joint_evt"] = ((tmp["EVI_evt"] == 1) | (tmp["LAI_evt"] == 1)).astype(int)
    for event_col, response in [("EVI_evt", "EVI"), ("LAI_evt", "LAI"), ("Joint_evt", "Joint")]:
        near = tmp[tmp["aquatic_context"] == "near_permanent"]
        other = tmp[tmp["aquatic_context"] == "other_weak"]
        near_low = near.loc[near["HDSI_group"] == "low", event_col].values
        near_high = near.loc[near["HDSI_group"] == "high", event_col].values
        other_low = other.loc[other["HDSI_group"] == "low", event_col].values
        other_high = other.loc[other["HDSI_group"] == "high", event_col].values
        rng = np.random.default_rng(42)
        boots = []
        if min(len(near_low), len(near_high), len(other_low), len(other_high)) > 0:
            for _ in range(500):
                nl = rng.choice(near_low, len(near_low), replace=True)
                nh = rng.choice(near_high, len(near_high), replace=True)
                ol = rng.choice(other_low, len(other_low), replace=True)
                oh = rng.choice(other_high, len(other_high), replace=True)
                boots.append((oh.mean() - ol.mean()) - (nh.mean() - nl.mean()))
            boots = np.asarray(boots)
            est = boots.mean()
            lo = np.percentile(boots, 2.5)
            hi = np.percentile(boots, 97.5)
        else:
            est, lo, hi = np.nan, np.nan, np.nan
        sens_rows.append({"threshold": thr, "response": response, "AFB": est, "ci_low": lo, "ci_high": hi})
sens_df = pd.DataFrame(sens_rows)

fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0), sharey=True, constrained_layout=False)
for ax, response, color, lab in zip(
    axes, ["EVI", "LAI", "Joint"],
    [COL["EVI"], COL["LAI"], COL["Joint"]],
    ["A", "B", "C"]
):
    sub = sens_df[sens_df["response"] == response].copy()
    ax.errorbar(sub["threshold"], sub["AFB"],
                yerr=[sub["AFB"] - sub["ci_low"], sub["ci_high"] - sub["AFB"]],
                fmt="o-", color=color, capsize=3.5, linewidth=1.8, markersize=6)
    ax.axhline(0, color="black", linestyle="--", linewidth=0.9)
    ax.set_title(response, loc="left", pad=8)
    ax.set_xlabel("Loss threshold")
    panel_label(ax, lab, fontsize=12)
    clean_ax(ax)
axes[0].set_ylabel("AFB = Fragility(other/weak) − Fragility(near permanent)")
fig.suptitle("Threshold sensitivity of aquatic buffering",
             fontsize=14, fontweight="bold", y=0.98)
plt.subplots_adjust(left=0.08, right=0.98, top=0.85, bottom=0.20, wspace=0.25)
plt.show()