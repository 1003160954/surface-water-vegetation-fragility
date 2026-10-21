# Surface water permanence and vegetation fragility

Reproducible analysis of aquatic context and vegetation risk amplification under compound heat and dryness across China, 2018–2023.

This repository provides analysis code, data-source documentation, derived numerical tables and visualizations. Large analysis exports are stored separately; manuscript and reviewer-response working drafts are not included. Quantitative figures are generated from recorded numerical outputs. AI assistance is disclosed separately where applicable.

## Research question and principal results

Is vegetation risk amplification under summer heat and dryness weaker in permanent-water landscape contexts? Fragility is the event probability under high heat–dryness exposure minus the probability under low exposure. Events represent negative vegetation-index anomalies relative to a recent local seasonal reference, rather than mortality or permanent vegetation loss.

- The main sample contains 236,995 grid-year records from 48,843 cells at 0.1° resolution. Corrected permanent, seasonal and other contexts contain 2,549, 15,262 and 31,032 cells, respectively.
- Permanent-context fragility is 0.1292 for EVI, 0.1408 for LAI and 0.1600 for the combined event. Corresponding other-context estimates are 0.1901, 0.1972 and 0.2209. Integer event counts and denominators are provided in the [risk arithmetic table](reviewer_evidence/risk_numerators_and_arithmetic.csv).
- Other-minus-permanent contrasts retain positive grid-cluster intervals for all three responses after measured-covariate adjustment and exclusion of dominant agricultural cover. Evidence for permanent versus seasonal water is clearest for EVI; LAI remains uncertain, and combined-event evidence depends on the target population and uncertainty treatment.
- At 1° and 2° spatial block sizes, all three other-minus-permanent contrasts and the seasonal-minus-permanent EVI contrast retain intervals above zero. At 5°, all unadjusted contrast intervals include zero. All scales are reported together.

These are observational buffering associations. A 1 km grid-centroid context does not identify a protection radius for vegetation within the cell or establish a causal hydrological mechanism.

### Adjusted contrasts using corrected contexts

AFB = fragility(comparison) − fragility(permanent). Positive values indicate weaker risk amplification in permanent-water contexts. Units are absolute probability differences. Estimates use a common empirical covariate distribution; intervals use grid-cluster sandwich covariance and the delta method, separately from the spatial-block analysis.

| Comparison | Response | Adjusted AFB | 95% interval |
| --- | --- | --- | --- |
| Seasonal − permanent | EVI | 0.0378 | [0.0126, 0.0630] |
| Other − permanent | EVI | 0.0664 | [0.0421, 0.0908] |
| Seasonal − permanent | LAI | 0.0079 | [-0.0155, 0.0314] |
| Other − permanent | LAI | 0.0412 | [0.0186, 0.0638] |
| Seasonal − permanent | combined | 0.0323 | [0.0054, 0.0591] |
| Other − permanent | combined | 0.0667 | [0.0407, 0.0926] |

## From data to results

The [computation guide](docs/COMPUTATION_STEPS.md) specifies inputs, transformations, commands, outputs and assumptions. The [field dictionary](reviewer_evidence/all_panel_fields_dictionary.csv) documents fields and missingness.

| Step | Operation | Evidence | Code |
| --- | --- | --- | --- |
| 1 | Identify source products | docs/data_sources.csv; reviewer_evidence/all_panel_fields_dictionary.csv | gee/; docs/original_sources/upstream_seasonal_source_only.ipynb |
| 2 | Verify the supplied analysis export | logs/E00_manifest.json; reviewer_evidence/raw_input_manifest.json | src/audit_pipeline.py |
| 3 | Apply validity rules and seasonal anomaly definitions | experiments/E00_original_reproduction/sample_flow.csv | docs/original_sources/upstream_seasonal_source_only.ipynb |
| 4 | Assign stress exposure and events | reviewer_evidence/worked_grid_year_examples.csv | src/audit_pipeline.py |
| 5 | Correct aquatic-context assignments | experiments/E01_spatial_support/national_corrected_1km_contexts.csv | gee/corrected_near_contexts.py |
| 6 | Calculate risk, fragility and AFB | reviewer_evidence/risk_numerators_and_arithmetic.csv; reviewer_evidence/AFB_arithmetic.csv | src/audit_pipeline.py; src/analyze_corrected_1km.py |
| 7 | Adjust for measured covariates | experiments/E01_spatial_support/corrected_1km/harmonized/ | src/harmonized_confounding.py; src/analyze_corrected_1km.py |
| 8 | Assess landscape and agricultural composition | reviewer_evidence/corrected_covariate_balance.csv; experiments/E02_confounding/corrected_noncrop/ | src/common_landscape.py; src/reviewer_noncrop.py |
| 9 | Evaluate sensitivity and spatial uncertainty | experiments/E01_spatial_support/corrected_1km/robustness/ | src/robustness.py; src/analyze_corrected_1km.py |
| 10 | Generate quantitative figures | revised_figures/figure_source_map.csv; reviewer_evidence/figure_source_map.csv | src/revision_figures.py; src/plot_corrected_fragility.py; src/build_reviewer_evidence.py |

### Worked AFB calculation

For combined events in permanent contexts, low- and high-HDSI risks are 0.214610 and 0.374590, giving fragility 0.159979. Other-context risks are 0.192537 and 0.413478, giving fragility 0.220940. Thus other-minus-permanent AFB is 0.220940 − 0.159979 = 0.060961.

Full-precision values and integer numerators/denominators are in [risk_numerators_and_arithmetic.csv](reviewer_evidence/risk_numerators_and_arithmetic.csv) and [AFB_arithmetic.csv](reviewer_evidence/AFB_arithmetic.csv). Rounded figure labels should not be used to check final-digit agreement. The 27 [worked grid-year examples](reviewer_evidence/worked_grid_year_examples.csv) come from six grids selected with fixed random seeds; they illustrate the calculations and are not the inference sample.

## Data sources

| Dataset | Use in this analysis | Source and interpretation |
| --- | --- | --- |
| JRC Global Surface Water v1.4 | Water masks, centroid classification and native-pixel distance checks | [Official catalog](https://developers.google.com/earth-engine/datasets/catalog/JRC_GSW1_4_GlobalSurfaceWater); `JRC/GSW1_4/GlobalSurfaceWater`, 30 m, 1984–2021. Long-term attributes do not measure annual water availability. |
| MODIS MCD12Q1 v6.1, 2018 | Categorical land cover for adjustment and stratification | [Official catalog](https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MCD12Q1); `MODIS/061/MCD12Q1`. Categories are aggregated by spatial mode rather than averaging class codes. |
| WorldClim v1 monthly, 1960–1990 | Background annual mean temperature and annual precipitation | [Official catalog](https://developers.google.com/earth-engine/datasets/catalog/WORLDCLIM_V1_MONTHLY); `WORLDCLIM/V1/MONTHLY`. These covariates are not an independent vegetation-anomaly reference. |
| ESA WorldCover v200, 2021 | Vegetation mask for the 24-cell within-cell distance pilot | [Official catalog](https://developers.google.com/earth-engine/datasets/catalog/ESA_WorldCover_v200); `ESA/WorldCover/v200`. Not used to define the national main response. |
| Supplied vegetation–climate analysis export | EVI/LAI, summer climate indicators, validity flags and original context variables | [Input identity and SHA-256](reviewer_evidence/raw_input_manifest.json). This is a processed analysis table, not raw satellite imagery. Some exact upstream MODIS/ERA5 collection versions and QA details remain unverified and are identified as provenance gaps. |
| Existing climatology export | Check the sample size of the recent seasonal reference | At most six annual observations per reference position; see [reference metadata](experiments/E03_independent_baseline/status.json). It does not provide an independent 2001–2017 reference. |

Further details are in [data_sources.csv](docs/data_sources.csv) and [DATA_ACCESS.md](docs/DATA_ACCESS.md). Credit water maps to **EC JRC/Google** and follow each product’s attribution and usage terms. Source catalog links identify products; they do not establish complete recovery of the original processing history.

## Covariate adjustment and sensitivity analyses

The analysis uses local seasonal standardization; adjustment for terrain, land cover, background climate, year and geographic terms; supported landscape strata with common weights; exclusion of dominant cropland and mosaic classes; and grid-cluster and multi-scale spatial uncertainty calculations. The [controls and limitations table](reviewer_evidence/confounding_controls_and_limits.csv) explains what each approach addresses.

Within the restricted, commonly weighted landscape sample, absolute descriptive standardized mean differences are approximately 0.037 for elevation, 0.002 for slope, 0.028 for temperature and 0.033 for precipitation. Selection and weighting change the target population, so the differences cannot be attributed to weighting alone or interpreted as removal of all confounding.

Excluding 2018 modal cropland class 12 and cropland/natural-vegetation mosaic class 14 leaves 172,463 records from 36,238 grids. The adjusted model fits 82,743 low/high-HDSI records from 30,212 grids, retaining national exposure cutpoints. Seasonal-minus-permanent AFB is 0.0480 [0.0170, 0.0790] for EVI, 0.0410 [0.0080, 0.0739] for combined events and 0.0218 [−0.0071, 0.0507] for LAI. All three other-minus-permanent intervals remain above zero; see [full estimates](experiments/E02_confounding/corrected_noncrop/M2_AFB.csv). This restriction reduces dominant agricultural-composition differences but does not identify irrigation-free vegetation. Groundwater, soil properties, irrigation outside these classes and dynamic management remain unmeasured.

## Figures


### Main figures using corrected aquatic contexts

#### Fig1_workflow — Workflow and calculation framework

![Workflow and calculation framework](revised_figures/figures/Fig1_workflow.png)

[PDF](revised_figures/figures/Fig1_workflow.pdf) · [Main-figure source map](revised_figures/figure_source_map.csv)

#### Fig2_JRC_attributes — JRC occurrence, seasonality and categorical transitions

![JRC occurrence, seasonality and categorical transitions](revised_figures/figures/Fig2_JRC_attributes.png)

[PDF](revised_figures/figures/Fig2_JRC_attributes.pdf) · [Main-figure source map](revised_figures/figure_source_map.csv)

#### Fig3_corrected_contexts — Corrected national contexts and combined-event spatial distribution

![Corrected national contexts and combined-event spatial distribution](revised_figures/figures/Fig3_corrected_contexts.png)

[PDF](revised_figures/figures/Fig3_corrected_contexts.pdf) · [Main-figure source map](revised_figures/figure_source_map.csv)

#### Fig4_exposure_risks — Spatial distribution of heat–dryness exposure and event risks

![Spatial distribution of heat–dryness exposure and event risks](revised_figures/figures/Fig4_exposure_risks.png)

[PDF](revised_figures/figures/Fig4_exposure_risks.pdf) · [Main-figure source map](revised_figures/figure_source_map.csv)

#### Fig5_context_fragility — Vegetation risk amplification across the three contexts

![Vegetation risk amplification across the three contexts](revised_figures/figures/Fig5_context_fragility.png)

[PDF](revised_figures/figures/Fig5_context_fragility.pdf) · [Main-figure source map](revised_figures/figure_source_map.csv)

#### Fig6_corrected_contrasts — Unadjusted and adjusted AFB estimates with intervals

![Unadjusted and adjusted AFB estimates with intervals](revised_figures/figures/Fig6_corrected_contrasts.png)

[PDF](revised_figures/figures/Fig6_corrected_contrasts.pdf) · [Main-figure source map](revised_figures/figure_source_map.csv)

#### Fig7_landscape_contrasts — Permanent versus seasonal contexts within supported landscapes

![Permanent versus seasonal contexts within supported landscapes](revised_figures/figures/Fig7_landscape_contrasts.png)

[PDF](revised_figures/figures/Fig7_landscape_contrasts.pdf) · [Main-figure source map](revised_figures/figure_source_map.csv)

#### Fig8_spatial_uncertainty — Results for all 1°, 2° and 5° spatial blocks

![Results for all 1°, 2° and 5° spatial blocks](revised_figures/figures/Fig8_spatial_uncertainty.png)

[PDF](revised_figures/figures/Fig8_spatial_uncertainty.pdf) · [Main-figure source map](revised_figures/figure_source_map.csv)

### Additional corrected results and diagnostics

#### SFig12_noncrop_sensitivity

![SFig12_noncrop_sensitivity](reviewer_evidence/figures/SFig12_noncrop_sensitivity.png)

[PDF](reviewer_evidence/figures/SFig12_noncrop_sensitivity.pdf) · [Numerical sources](reviewer_evidence/figure_source_map.csv)

#### SFig13_corrected_balance

![SFig13_corrected_balance](reviewer_evidence/figures/SFig13_corrected_balance.png)

[PDF](reviewer_evidence/figures/SFig13_corrected_balance.pdf) · [Numerical sources](reviewer_evidence/figure_source_map.csv)

#### SFig14_corrected_HDSI

![SFig14_corrected_HDSI](reviewer_evidence/figures/SFig14_corrected_HDSI.png)

[PDF](reviewer_evidence/figures/SFig14_corrected_HDSI.pdf) · [Numerical sources](reviewer_evidence/figure_source_map.csv)

#### SFig15_corrected_thresholds

![SFig15_corrected_thresholds](reviewer_evidence/figures/SFig15_corrected_thresholds.png)

[PDF](reviewer_evidence/figures/SFig15_corrected_thresholds.pdf) · [Numerical sources](reviewer_evidence/figure_source_map.csv)

#### SFig16_corrected_year_omissions

![SFig16_corrected_year_omissions](reviewer_evidence/figures/SFig16_corrected_year_omissions.png)

[PDF](reviewer_evidence/figures/SFig16_corrected_year_omissions.pdf) · [Numerical sources](reviewer_evidence/figure_source_map.csv)

#### SFig17_sample_flow

![SFig17_sample_flow](reviewer_evidence/figures/SFig17_sample_flow.png)

[PDF](reviewer_evidence/figures/SFig17_sample_flow.pdf) · [Numerical sources](reviewer_evidence/figure_source_map.csv)

### Original-result reproduction and historical distance checks

These figures document historical grouping, distance calculations and algorithms. They explain the original numerical results; **historical kilometre curves do not establish ecological reach and do not replace the corrected results above**. SFig9–11 show supplementary results using the corrected 1 km classification. See the source maps for each figure’s scope.

<details>
<summary>Expand reproduction checks and spatial figures</summary>

**Fig1_verified_workflow**

![Fig1_verified_workflow](figures/Fig1_verified_workflow.png)

[PDF](figures/Fig1_verified_workflow.pdf)

**Fig2_JRC_metrics_discrete_transition**

![Fig2_JRC_metrics_discrete_transition](figures/Fig2_JRC_metrics_discrete_transition.png)

[PDF](figures/Fig2_JRC_metrics_discrete_transition.pdf)

**Fig3_context_and_combined_events**

![Fig3_context_and_combined_events](figures/Fig3_context_and_combined_events.png)

[PDF](figures/Fig3_context_and_combined_events.pdf)

**Fig4_HDSI_and_event_risks**

![Fig4_HDSI_and_event_risks](figures/Fig4_HDSI_and_event_risks.png)

[PDF](figures/Fig4_HDSI_and_event_risks.pdf)

**Fig5_numeric_distance_bands**

![Fig5_numeric_distance_bands](figures/Fig5_numeric_distance_bands.png)

[PDF](figures/Fig5_numeric_distance_bands.pdf)

**Fig6_context_risk_and_fragility**

![Fig6_context_risk_and_fragility](figures/Fig6_context_risk_and_fragility.png)

[PDF](figures/Fig6_context_risk_and_fragility.pdf)

**Fig7_exploratory_diagnostics**

![Fig7_exploratory_diagnostics](figures/Fig7_exploratory_diagnostics.png)

[PDF](figures/Fig7_exploratory_diagnostics.pdf)

**Fig8_all_regions_with_support**

![Fig8_all_regions_with_support](figures/Fig8_all_regions_with_support.png)

[PDF](figures/Fig8_all_regions_with_support.pdf)

**SFig10_corrected_1km_spatial_blocks**

![SFig10_corrected_1km_spatial_blocks](figures/SFig10_corrected_1km_spatial_blocks.png)

[PDF](figures/SFig10_corrected_1km_spatial_blocks.pdf)

**SFig11_corrected_1km_context_map**

![SFig11_corrected_1km_context_map](figures/SFig11_corrected_1km_context_map.png)

[PDF](figures/SFig11_corrected_1km_context_map.pdf)

**SFig1_same_sample_adjustment**

![SFig1_same_sample_adjustment](figures/SFig1_same_sample_adjustment.png)

[PDF](figures/SFig1_same_sample_adjustment.pdf)

**SFig2_spatial_block_uncertainty**

![SFig2_spatial_block_uncertainty](figures/SFig2_spatial_block_uncertainty.png)

[PDF](figures/SFig2_spatial_block_uncertainty.pdf)

**SFig3_HDSI_variants**

![SFig3_HDSI_variants](figures/SFig3_HDSI_variants.png)

[PDF](figures/SFig3_HDSI_variants.pdf)

**SFig4_event_thresholds**

![SFig4_event_thresholds](figures/SFig4_event_thresholds.png)

[PDF](figures/SFig4_event_thresholds.pdf)

**SFig5_leave_analysis_year_out**

![SFig5_leave_analysis_year_out](figures/SFig5_leave_analysis_year_out.png)

[PDF](figures/SFig5_leave_analysis_year_out.pdf)

**SFig6_centroid_distance_spline**

![SFig6_centroid_distance_spline](figures/SFig6_centroid_distance_spline.png)

[PDF](figures/SFig6_centroid_distance_spline.pdf)

**SFig7_within_cell_pilot**

![SFig7_within_cell_pilot](figures/SFig7_within_cell_pilot.png)

[PDF](figures/SFig7_within_cell_pilot.pdf)

**SFig8_harmonized_and_landscape**

![SFig8_harmonized_and_landscape](figures/SFig8_harmonized_and_landscape.png)

[PDF](figures/SFig8_harmonized_and_landscape.pdf)

**SFig9_corrected_1km_permanence**

![SFig9_corrected_1km_permanence](figures/SFig9_corrected_1km_permanence.png)

[PDF](figures/SFig9_corrected_1km_permanence.pdf)

</details>

Independent native-pixel distance checks are provided in[Supplementary distance validation](revised_figures/figures/SFig18_distance_validation.png). All available figures, numerical sources, analysis scopes and scripts are indexed in [FIGURE_GALLERY.md](docs/FIGURE_GALLERY.md). Only figures supported by available numerical outputs are included.

## Local reproduction

Use Python 3.12; recorded dependency versions are listed in [requirements.txt](requirements.txt). Large input tables are stored outside Git. Set `AQUATIC_PANEL` and, when needed, `AQUATIC_CLIMATOLOGY` to the supplied analysis exports. Access to these inputs is separate from access to this repository. Public derived tables support arithmetic checks and selected figure reconstruction; they do not by themselves support an independent reconstruction from raw imagery through every model.

```text
python -m pip install -r requirements.txt
# Set AQUATIC_PANEL and optionally AQUATIC_CLIMATOLOGY in your shell.
python src/original_exact.py
python src/analyze_corrected_1km.py
python src/reviewer_noncrop.py
python src/build_reviewer_evidence.py
python src/revision_figures.py
python src/plot_corrected_fragility.py
python tests/validate_evidence.py
python tests/validate_public_evidence.py
```

`revision_figures.py` reconstructs the earlier main-figure layout, including the historical Fig. 5 distance check. `plot_corrected_fragility.py` creates the current Fig. 5 context comparison. Both are retained with distinct sources. The [computation guide](docs/COMPUTATION_STEPS.md) documents the historical reproduction entry points and expected outputs. Reproduction does not require resubmitting Earth Engine exports; heavy export tasks remain paused.

Summer (DOY 152–243) and late summer (DOY 213–273) are assigned using 16-day step midpoints and overlap on DOY 213–243. The exposure–response association is partly concurrent. Intervals containing zero do not establish equivalence, and positive point estimates alone do not establish robustness or causality. Spatial uncertainty should be read across all reported block sizes. Detailed provenance and analysis scope are documented in [REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md).

## Repository contents

The public repository contains research code, source metadata and derived evidence. It excludes manuscript and response drafts, Word documents, large input tables, credentials and private local repository history. Source-only notebooks are archived in [docs/original_sources](docs/original_sources). The initial input-discovery manifest records its historical execution state; recovered sources are described in the reproducibility documentation. File identities are recorded in [public_manifest.csv](public_manifest.csv).
