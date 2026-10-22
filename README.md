# Surface water permanence and vegetation fragility

Reproducible analysis of vegetation risk amplification under compound heat and dryness across China, 2018–2023.

The [analysis package](final_results/README.md) contains the current numerical results, figure files, full uncertainty outputs, and executable plotting code. The research design compares mapped aquatic contexts, broad centroid-distance associations, component exposure, seasonal canopy change, and regional contrasts.

## Principal findings

- The national analysis includes 236,995 grid-years from 48,843 cells. Permanent, seasonal, and other contexts contain 2,549, 15,262, and 31,032 cells.
- Permanent-context fragility is 0.1292 (EVI), 0.1408 (LAI), and 0.1600 (combined events), compared with 0.1901, 0.1972, and 0.2209 in other contexts.
- The seasonal-minus-permanent EVI contrast is supported after measured-covariate adjustment and common-landscape restriction. LAI permanence evidence is weaker. All national 5° spatial-block contrast intervals include zero.
- Broad distance-related increases occur in the geographically restricted distance subset; they do not establish a shoreline protection radius.
- Regional LAI contrasts are positive in Southwest and South China and negative in East China. Mean high-stress seasonal LAI change is similar in permanent and other contexts.
- Excluding each evaluated year from its anomaly reference retains positive permanent-versus-other contrasts in paired analyses, but event agreement is 91.3% for EVI, 89.5% for LAI and 87.2% for combined events. Classification is therefore not baseline-invariant.

## Research overview

![Methodological workflow](final_results/figures/Fig1.png)
![Aquatic contexts and vegetation events](final_results/figures/Fig3.png)
![Distance association](final_results/figures/Fig5.png)
![Aquatic-context risks and fragility](final_results/figures/Fig6.png)
![Regional contrasts](final_results/figures/Fig8.png)

All [main and supplementary graphics](final_results/figures/) are available as PNG and PDF; the workflow has an editable PowerPoint source. The [analysis guide](final_results/README.md) describes regeneration and maps the evidence files to analyses. Figures use numerical data and deterministic plotting scripts. AI-assisted tools were used for language polishing and assistance with writing portions of the analysis code.

## From data to results

The [computation guide](docs/COMPUTATION_STEPS.md) specifies inputs, transformations, commands, outputs and assumptions. The [field dictionary](analysis_evidence/all_panel_fields_dictionary.csv) documents fields and missingness.

| Step | Operation | Evidence | Code |
| --- | --- | --- | --- |
| 1 | Identify source products | docs/data_sources.csv; analysis_evidence/all_panel_fields_dictionary.csv | gee/; docs/original_sources/upstream_seasonal_source_only.ipynb |
| 2 | Verify the supplied analysis export | logs/E00_manifest.json; analysis_evidence/raw_input_manifest.json | src/audit_pipeline.py |
| 3 | Apply validity rules and seasonal anomaly definitions | experiments/E00_original_reproduction/sample_flow.csv | docs/original_sources/upstream_seasonal_source_only.ipynb |
| 4 | Assign stress exposure and events | analysis_evidence/worked_grid_year_examples.csv | src/audit_pipeline.py |
| 5 | Assign aquatic contexts using projected water distances | experiments/E01_spatial_support/national_corrected_1km_contexts.csv | gee/corrected_near_contexts.py |
| 6 | Calculate risk, fragility and AFB | analysis_evidence/risk_numerators_and_arithmetic.csv; analysis_evidence/AFB_arithmetic.csv | src/audit_pipeline.py; src/analyze_corrected_1km.py |
| 7 | Adjust for measured covariates | experiments/E01_spatial_support/corrected_1km/harmonized/ | src/harmonized_confounding.py; src/analyze_corrected_1km.py |
| 8 | Assess landscape and agricultural composition | analysis_evidence/corrected_covariate_balance.csv; experiments/E02_confounding/corrected_noncrop/ | src/common_landscape.py; src/noncropland_sensitivity.py |
| 9 | Evaluate sensitivity and spatial uncertainty | experiments/E01_spatial_support/corrected_1km/robustness/ | src/robustness.py; src/analyze_corrected_1km.py |
| 10 | Recalculate target-year-excluded anomalies | final_results/data/baseline/ | src/final_release/baseline_crossfit.py |
| 11 | Generate quantitative figures | final_results/figure_source_map.csv | src/final_release/plot_final_figures.py |

### Worked AFB calculation

For combined events in permanent contexts, low- and high-HDSI risks are 0.214610 and 0.374590, giving fragility 0.159979. Other-context risks are 0.192537 and 0.413478, giving fragility 0.220940. Thus other-minus-permanent AFB is 0.220940 − 0.159979 = 0.060961.

Full-precision values and integer numerators/denominators are in [risk_numerators_and_arithmetic.csv](analysis_evidence/risk_numerators_and_arithmetic.csv) and [AFB_arithmetic.csv](analysis_evidence/AFB_arithmetic.csv). Rounded figure labels should not be used to check final-digit agreement. The 27 [worked grid-year examples](analysis_evidence/worked_grid_year_examples.csv) come from six grids selected with fixed random seeds; they illustrate the calculations and are not the inference sample.

## Data sources

| Dataset | Use in this analysis | Source and interpretation |
| --- | --- | --- |
| JRC Global Surface Water v1.4 | Water masks, centroid classification and native-pixel distance checks | [Official catalog](https://developers.google.com/earth-engine/datasets/catalog/JRC_GSW1_4_GlobalSurfaceWater); `JRC/GSW1_4/GlobalSurfaceWater`, 30 m, 1984–2021. Long-term attributes do not measure annual water availability. |
| MODIS MCD12Q1 v6.1, 2018 | Categorical land cover for adjustment and stratification | [Official catalog](https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MCD12Q1); `MODIS/061/MCD12Q1`. Categories are aggregated by spatial mode rather than averaging class codes. |
| WorldClim v1 monthly, 1960–1990 | Background annual mean temperature and annual precipitation | [Official catalog](https://developers.google.com/earth-engine/datasets/catalog/WORLDCLIM_V1_MONTHLY); `WORLDCLIM/V1/MONTHLY`. These covariates are not an independent vegetation-anomaly reference. |
| ESA WorldCover v200, 2021 | Vegetation mask for the 24-cell within-cell distance pilot | [Official catalog](https://developers.google.com/earth-engine/datasets/catalog/ESA_WorldCover_v200); `ESA/WorldCover/v200`. Not used to define the national main response. |
| Supplied vegetation–climate analysis export | EVI/LAI, summer climate indicators, validity flags and original context variables | [Input identity and SHA-256](analysis_evidence/raw_input_manifest.json). This is a processed analysis table, not raw satellite imagery. Some exact upstream MODIS/ERA5 collection versions and QA details remain unverified and are identified as provenance gaps. |
| Existing climatology export | Check the sample size of the recent seasonal reference | At most six annual observations per reference position; see [reference metadata](experiments/E03_independent_baseline/status.json). It does not provide an independent 2001–2017 reference. |

Further details are in [data_sources.csv](docs/data_sources.csv) and [DATA_ACCESS.md](docs/DATA_ACCESS.md). Credit water maps to **EC JRC/Google** and follow each product’s attribution and usage terms. Source catalog links identify products; they do not establish complete recovery of the original processing history.

## Covariate adjustment and sensitivity analyses

The analysis uses local seasonal standardization; adjustment for terrain, land cover, background climate, year and geographic terms; supported landscape strata with common weights; exclusion of dominant cropland and mosaic classes; and grid-cluster and multi-scale spatial uncertainty calculations. The [controls and limitations table](analysis_evidence/confounding_controls_and_limits.csv) explains what each approach addresses.

Within the restricted, commonly weighted landscape sample, absolute descriptive standardized mean differences are approximately 0.037 for elevation, 0.002 for slope, 0.028 for temperature and 0.033 for precipitation. Selection and weighting change the target population, so the differences cannot be attributed to weighting alone or interpreted as removal of all confounding.

Excluding 2018 modal cropland class 12 and cropland/natural-vegetation mosaic class 14 leaves 172,463 records from 36,238 grids. The adjusted model fits 82,743 low/high-HDSI records from 30,212 grids, retaining national exposure cutpoints. Seasonal-minus-permanent AFB is 0.0480 [0.0170, 0.0790] for EVI, 0.0410 [0.0080, 0.0739] for combined events and 0.0218 [−0.0071, 0.0507] for LAI. All three other-minus-permanent intervals remain above zero; see [full estimates](experiments/E02_confounding/corrected_noncrop/M2_AFB.csv). This restriction reduces dominant agricultural-composition differences but does not identify irrigation-free vegetation. Groundwater, soil properties, irrigation outside these classes and dynamic management remain unmeasured.


## Complete figure gallery

### Fig1

![Fig1](final_results/figures/Fig1.png)

[PDF](final_results/figures/Fig1.pdf)

### Fig2

![Fig2](final_results/figures/Fig2.png)

[PDF](final_results/figures/Fig2.pdf)

### Fig3

![Fig3](final_results/figures/Fig3.png)

[PDF](final_results/figures/Fig3.pdf)

### Fig4

![Fig4](final_results/figures/Fig4.png)

[PDF](final_results/figures/Fig4.pdf)

### Fig5

![Fig5](final_results/figures/Fig5.png)

[PDF](final_results/figures/Fig5.pdf)

### Fig6

![Fig6](final_results/figures/Fig6.png)

[PDF](final_results/figures/Fig6.pdf)

### Fig7

![Fig7](final_results/figures/Fig7.png)

[PDF](final_results/figures/Fig7.pdf)

### Fig8

![Fig8](final_results/figures/Fig8.png)

[PDF](final_results/figures/Fig8.pdf)

### Fig_S1

![Fig_S1](final_results/figures/Fig_S1.png)

[PDF](final_results/figures/Fig_S1.pdf)

### Fig_S2

![Fig_S2](final_results/figures/Fig_S2.png)

[PDF](final_results/figures/Fig_S2.pdf)

### Fig_S3

![Fig_S3](final_results/figures/Fig_S3.png)

[PDF](final_results/figures/Fig_S3.pdf)

### Fig_S4

![Fig_S4](final_results/figures/Fig_S4.png)

[PDF](final_results/figures/Fig_S4.pdf)

### Fig_S5

![Fig_S5](final_results/figures/Fig_S5.png)

[PDF](final_results/figures/Fig_S5.pdf)

### Fig_S6_adjusted_contrasts

![Fig_S6_adjusted_contrasts](final_results/figures/Fig_S6_adjusted_contrasts.png)

[PDF](final_results/figures/Fig_S6_adjusted_contrasts.pdf)

### Fig_S7_common_landscapes

![Fig_S7_common_landscapes](final_results/figures/Fig_S7_common_landscapes.png)

[PDF](final_results/figures/Fig_S7_common_landscapes.pdf)

### Fig_S8_covariate_balance

![Fig_S8_covariate_balance](final_results/figures/Fig_S8_covariate_balance.png)

[PDF](final_results/figures/Fig_S8_covariate_balance.pdf)

### Fig_S9_noncropland

![Fig_S9_noncropland](final_results/figures/Fig_S9_noncropland.png)

[PDF](final_results/figures/Fig_S9_noncropland.pdf)

### Fig_S10

![Fig_S10](final_results/figures/Fig_S10.png)

[PDF](final_results/figures/Fig_S10.pdf)

### Fig_S11

![Fig_S11](final_results/figures/Fig_S11.png)

[PDF](final_results/figures/Fig_S11.pdf)

### Fig_S12

![Fig_S12](final_results/figures/Fig_S12.png)

[PDF](final_results/figures/Fig_S12.pdf)

### Fig_S13

![Fig_S13](final_results/figures/Fig_S13.png)

[PDF](final_results/figures/Fig_S13.pdf)

## Reproduction

Install `requirements.txt`, then run `python src/final_release/plot_final_figures.py` from the repository root. This reconstructs all quantitative final plots from the published derived CSV files without submitting Earth Engine tasks. The methods flowchart is supplied as an editable PowerPoint file.

For the target-year-excluded anomaly experiment, set `RAW_CORE_PANEL` to the patched 16-day source panel and `AQUATIC_PANEL` to the seasonal analysis panel, then run `python src/final_release/baseline_crossfit.py`. The [analysis guide](final_results/README.md) identifies every evidence directory. The [figure source map](final_results/figure_source_map.csv) links figures to numerical inputs.

The [computation guide](docs/COMPUTATION_STEPS.md) documents the analysis sequence. Numbered main and supplementary figures are in `final_results/figures`. Computationally intensive cloud exports are disabled by default. The independent long-reference and national within-cell/grid-realignment experiments have not been completed. Exact upstream product and quality-mask gaps are documented in the source inventory rather than inferred from variable names.

Summer (DOY 152–243) and late summer (DOY 213–273), assigned by 16-day midpoints, overlap on DOY 213–243. The exposure–response association is partly concurrent. Positive point estimates alone do not establish causality; spatial uncertainty is reported across all tested block sizes.

## Public files

This repository contains research code, source metadata and derived evidence. Large input panels, credentials and private local files are not distributed. The [manifest](final_results/manifest.csv) records hashes of the final derived files. Access to original input panels is separate from access to this repository.
