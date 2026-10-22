# Computation sequence

Run commands from the repository root after installing `requirements.txt`. Input access is described in [DATA_ACCESS.md](DATA_ACCESS.md). The [source inventory](data_sources.csv) records product identities and provenance gaps. Large source panels are supplied separately; set `AQUATIC_PANEL` and, for anomaly reconstruction, `RAW_CORE_PANEL` to their local locations.

## 1. Source data and seasonal summaries

The supplied seasonal panel contains 244,306 grid-year records from 55,418 grids before validity filtering. The analysis retains 236,995 records from 48,843 grids. Source code is archived in `docs/original_sources`; `src/audit_pipeline.py` records filtering and field checks. The field dictionary and input checksum are in `analysis_evidence`.

Seasonal positions use 16-day midpoints: spring DOY 60–151, summer 152–243 and late summer 213–273. Summer exposure and late-summer response overlap on DOY 213–243. EVI and LAI validity ranges are −0.2–1.2 and 0–12. At each grid and seasonal position, 2018–2023 anomalies use the available mean and sample standard deviation (ddof=1); standard deviations at or below 1e−6 are masked and standardized values clipped to [−5,5]. The reference contains at most six annual observations.

## 2. Exposure and vegetation events

HDSI = (VPD_z + LST_z − soil_moisture_z − precipitation_z) / 4. Fixed national quartiles are −0.2974858281280519 and 0.2882174368150908. Strictly lower/higher observations form the low/high groups, with 59,249 records each. EVI and LAI events require a late-summer anomaly below −0.5; combined events mean either indicator meets that criterion.

`analysis_evidence/worked_grid_year_examples.csv` shows individual calculations. These examples illustrate the arithmetic and are not the estimation sample.

## 3. Aquatic contexts and distances

JRC Global Surface Water v1.4 permanent water uses occurrence ≥80% and seasonality ≥10 months. Seasonal water uses seasonality ≥1 and <10 months. Distances are calculated in 30 m UTM grids. The 1 km classification gives permanent context priority when both water types are nearby. Contexts contain 2,549 permanent, 15,262 seasonal and 31,032 other grids. See `gee/corrected_near_contexts.py` and `experiments/E01_spatial_support`.

Continuous-distance exports cover 17,056 grids in UTM zones 43–48. Excluding 15 grids outside the 61.44 km search range leaves 17,041 grids and 80,265 records. Bands are [0,1], (1,5], (5,10], (10,20], (20,40] and (40,61.44] km, with 962, 4,769, 4,379, 4,346, 2,280 and 305 grids. These are centroid-distance associations for a restricted geographic sample.

## 4. Risk, fragility and context contrasts

Event risk is event count divided by exposure-group sample size. Fragility is high-exposure risk minus low-exposure risk. AFB is comparison-context fragility minus permanent-context fragility. Integer counts and full-precision arithmetic are in `analysis_evidence/risk_numerators_and_arithmetic.csv` and `AFB_arithmetic.csv`. Whole-grid bootstrap resampling uses 300 replicates, seed 20261020 and fixed national exposure cutpoints. Run `src/analyze_corrected_1km.py` for context analyses.

## 5. Measured-covariate adjustment and landscape comparisons

Logistic models include exposure-by-context interactions, elevation, slope, modal 2018 land cover, background annual temperature, log-transformed precipitation, year and quadratic geographic terms. Standardized predictions use a common fitted population; grid-cluster covariance supports uncertainty estimates. The national adjusted analysis uses 118,498 low/high records from 42,060 grids.

Common-landscape comparisons use land cover and temperature, precipitation and elevation tertiles, requiring at least ten observations per context/exposure cell and five grids per context. The retained 46 strata contain 39,524 records and 13,881 grids. Shared stratum weights define the comparison population. Excluding modal cropland classes 12 and 14 retains 172,463 records from 36,238 grids; the adjusted fit uses 82,743 records from 30,212 grids. Entry points are `src/harmonized_confounding.py`, `src/common_landscape.py` and `src/noncropland_sensitivity.py`.

## 6. Distance-response modelling

Run `python src/final_release/validate_distance_gradient.py` with the seasonal panel and retained exports. Outputs are written to `outputs/distance_validation`. Alternative groupings, distance curves and contrast tables are published in `final_results/data/distance`.

The adjusted distance model uses cubic B-splines of log(1 + distance), five degrees of freedom, no spline intercept and bounds 0 to log(62.44). Internal knots use the one-third and two-thirds quantiles. The model fits 39,720 records from 14,201 grids. Grid and 1°, 2° and 5° spatial clusters provide covariance estimates. Spatial blocks use floor(longitude / block size) and floor(latitude / block size), with origin 0°. Broad distance association does not identify a fine shoreline buffering radius.

## 7. Component, seasonal and regional analyses

Run `python src/final_release/restore_core_experiments.py` to compute component exposures, seasonal LAI change and regional summaries from the seasonal panel. Outputs go to `outputs/component_regional_analysis`; published tables are in `final_results/data/diagnostics`. Regional estimates include East China and South China and report sample support.

## 8. Sensitivity and spatial uncertainty

Ten stress-index definitions comprise equal weighting, PCA, four component exclusions and four single-component indices. Event thresholds are −0.25, −0.5, −0.75 and −1.0. Six year-omission analyses assess temporal influence. Grid resampling and 1°, 2° and 5° spatial blocks assess specified dependence structures. Full scenario estimates are in `final_results/data/robustness`; entry points include `src/robustness.py` and `src/analyze_corrected_1km.py`.

Run `python src/final_release/baseline_crossfit.py` to exclude each target year from its anomaly reference, requiring at least three remaining observations. Paired comparisons retain 209,985 records from 37,878 grids. Event agreement is 91.3% for EVI, 89.5% for LAI and 87.2% for combined events. This experiment is distinct from omitting a year from the outcome analysis and does not supply an independent longer reference period.

## 9. Figures and verification

Run `python src/final_release/plot_final_figures.py` to regenerate quantitative plots from published tables without cloud exports. The editable workflow drawing is supplied separately. `final_results/figure_source_map.csv` connects each plot to its data. Main figures are 1–8 and supplementary figures 1–13. The seven geographic panels use their own southern-region crop in the inset; see `final_results/assets/README.md`.

Run `python tests/validate_public_evidence.py` to verify event arithmetic, fitted-model diagnostics, source syntax and README links. Manifests record file hashes.

## Interpretation and scope

These are observational associations. Adjustment and sensitivity analyses address measured covariates and specified analytical choices, but do not directly measure irrigation, groundwater or soil properties. Exact upstream collection and quality-mask gaps are documented rather than inferred. Independent long-reference, nationwide within-cell and grid-realignment analyses are not available. Computationally intensive exports remain paused; plotting and local validation do not restart them.
