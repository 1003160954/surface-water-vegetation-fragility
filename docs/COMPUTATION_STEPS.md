# Computation steps and numerical verification

This guide separates source provenance, deterministic arithmetic, model estimation and uncertainty. A script passing a check demonstrates the stated computation, not causal validity or complete upstream provenance.

## Step 1 Original data products

Evidence: `docs/data_sources.csv; reviewer_evidence/all_panel_fields_dictionary.csv`.

Executed result or boundary: Known JRC/land-cover/climate/pilot products; original vegetation/thermal IDs partly unresolved.

Code: `gee/; docs/original_sources/upstream_seasonal_source_only.ipynb`.

## Step 2 Supplied raw analysis export

Evidence: `logs/E00_manifest.json; reviewer_evidence/raw_input_manifest.json`.

Executed result or boundary: 244306 grid-years,55418grids; no duplicate grid-year keys.

Code: `src/audit_pipeline.py`.

## Step 3 Validity and seasonal anomaly recipe

Evidence: `experiments/E00_original_reproduction/sample_flow.csv`.

Executed result or boundary: 236995 rows,48843grids; at most6 annual observations per seasonal position.

Code: `docs/original_sources/upstream_seasonal_source_only.ipynb`.

## Step 4 Stress and event assignments

Evidence: `reviewer_evidence/worked_grid_year_examples.csv`.

Executed result or boundary: HDSI four-component mean; Q25/Q75; EVI/LAI<-.5; OR event.

Code: `src/audit_pipeline.py`.

## Step 5 Corrected aquatic contexts

Evidence: `experiments/E01_spatial_support/national_corrected_1km_contexts.csv`.

Executed result or boundary: 2549permanent,15262seasonal,31032other; 17056-grid binary check zero mismatches.

Code: `gee/corrected_near_contexts.py`.

## Step 6 Event risk and fragility

Evidence: `reviewer_evidence/risk_numerators_and_arithmetic.csv; reviewer_evidence/AFB_arithmetic.csv`.

Executed result or boundary: Integer event counts / denominators -> high-low risk -> signed context contrast.

Code: `src/audit_pipeline.py; src/analyze_corrected_1km.py`.

## Step 7 Measured-confounder adjustment

Evidence: `experiments/E01_spatial_support/corrected_1km/harmonized/`.

Executed result or boundary: M0/M1 same eligible sample; M2 common-covariate predictions,grid-cluster delta95%CI.

Code: `src/harmonized_confounding.py; src/analyze_corrected_1km.py`.

## Step 8 Landscape and agricultural sensitivity

Evidence: `reviewer_evidence/corrected_covariate_balance.csv; experiments/E02_confounding/corrected_noncrop/`.

Executed result or boundary: 46strata; excludes crop/mosaic in additional analysis; no full irrigation control claimed.

Code: `src/common_landscape.py; src/reviewer_noncrop.py`.

## Step 9 Robustness and spatial uncertainty

Evidence: `experiments/E01_spatial_support/corrected_1km/robustness/`.

Executed result or boundary: 10 indices,4 thresholds,6 analysis-year omissions,1/2/5degree blocks; all scenarios retained.

Code: `src/robustness.py; src/analyze_corrected_1km.py`.

## Step 10 Quantitative visualizations

Evidence: `manuscript_revision/figure_source_map.csv; reviewer_evidence/figure_source_map.csv`.

Executed result or boundary: Main8 figures plus corrected support plots; raw tables retained; no AI image generation.

Code: `src/revision_figures.py; src/build_reviewer_evidence.py`.

## Exact definitions

1. Main selection uses valid_model_both_y on the supplied China panel. There are no duplicate grid/year keys. Do not silently replace supplied validity flags with invented upstream screening.
2. Seasonal recipe: each grid/16-day position uses available2018–2023 means and sample SD; SD≤1e-6 is masked, step z is clipped to[-5,5], seasons average available step z. SpringDOY60–151; summer152–243; late summer213–273; autumn244–304, assigned by step midpoint. Reference counts are≤6 and exposure/response overlap.
3. HDSI=(summer_VPD_z+summer_LST_Day_z−summer_SM_rootproxy_z−summer_Precip_z)/4. The main Q25/Q75 are−0.2974858281280519 and0.2882174368150908. Strict inequalities define low/high tails, each59249 observations. The noncrop sensitivity retains those full-sample cuts. Other HDSI/year scenarios recalculate cuts per scenario, held fixed across its bootstrap draws.
4. EVI/LAI events are late-summer z<−0.5; combined=logicalOR. Alternative cutoffs are−0.25,−0.5,−0.75,−1.0. There is no joint-AND interpretation.
5. Permanent JRC mask is occurrence≥80 and seasonality≥10; seasonal is1≤seasonality<10. Native masks are fixed to30m northern UTM in the centroid’s6-degree zone. A64-pixel transform supports only distance≤1km binary flags. Permanent takes precedence over seasonal. The64vs2048-pixel check covers17056grids. It is not a national full-distance or vegetation-pixel validation.
6. For each context/outcome: risk_low=low_events/low_n; risk_high=high_events/high_n; fragility=risk_high−risk_low. AFB=fragility_comparison−fragility_permanent. See worked_grid_year_examples.csv for assignments, then risk_numerators_and_arithmetic.csv and AFB_arithmetic.csv for the entire main sample. These arithmetic tables can be checked without access to the large source export.
7. New unadjusted CIs use300 paired grid resamples (seed20261020), keeping years together and the same weights across outcomes/contexts. The point is the plug-in contrast, not the bootstrap mean. The original cell59 uses500 independent grid-year resamples and mean contrast, seed42; original regional routines use400 draws. Exact historical agreement is recorded separately.
8. M0 and M1 use236995 eligible rows; M2 fits118498low/high rows,42060grids. Logit predictors: high*context, elevation,slope,lon,lat,lon²,lat²,lon*lat,year indicators,categorical2018grid modalLC,WorldClimannualtemperature,log1p(annualprecipitation). Continuous predictors are fit-sample centered/scaled. Predictions for all context/exposure combinations use the same empirical fit-sample covariate distribution. Delta intervals use grid-cluster sandwich covariance, conditional on that distribution; they are not spatial-block intervals. Model rank/convergence and coefficients are included.
9. Landscape strata are LCmode×MATtertile×MAPtertile×elevationtertile. Every context/exposure cell needs≥10rows and each context≥5grids.46strata retain39524rows,13881grids,89.64%of eligible permanent/seasonal extremes. Common weights use pooled unique-grid counts. Four bootstrap draws have empty required cells, leaving296of300; no silent renormalization. Descriptive before/after balance uses different target populations and does not measure unobserved covariate balance.
10. Agricultural sensitivity removes2018modal classes12and14, keeping172463rows/36238grids. M2fits82743extremes/30212grids,32/32rank. Broad cover exclusion is not an irrigation-free filter. All six contrasts are published, including uncertain LAI permanence.
11. Spatial blocks are floor(lon/size),floor(lat/size), origin0,0; sizes1,2,5degrees,883,265,58observed blocks,300paired draws. Degrees vary in physical area by latitude. No tested block size is asserted to be the verified correlation range. All scenarios are retained; the5degreeCIs all include zero.
12. Six analysis-year omissions retain existing reference statistics. They do not reconstruct a leave-year-out anomaly baseline. HDSI variants include equal,PC1,four omissions,four singles; PC1 uses sample-standardized components and sign alignment with the equal index. Original splines used cubic log1p-distance df5 with grid-cluster pointwise delta intervals; those distances are not physical ecological reach evidence. Corrected national splines are not available.

## Reproduction order and outputs

The README commands are the primary corrected-analysis route. src/original_exact.py uses recovered cell59 and the supplied panel to verify historical risks and CIs. src/analyze_corrected_1km.py executes corrected unadjusted, harmonized, landscape and robustness outputs. src/reviewer_noncrop.py supplies fixed-cut agricultural exclusion. src/build_reviewer_evidence.py regenerates integer counts, examples, observed balance and corrected support plots. Plotting scripts consume recorded CSVs and the original panel where maps require grid attributes. tests/validate_evidence.py verifies the original/current numerical provenance; tests/validate_public_evidence.py checks the expanded evidence and public inventory.

No command in this local route submits new GEE tasks. Heavy export scripts remain guarded by heavy_pause.flag; do not remove that flag without explicit authorization. Original code archive includes historical writes, so do not execute the whole old notebook blindly. Exact upstream satellite-to-panel reproduction is incomplete until missing product/QA provenance and original remote inputs are verified.

## Missing evidence must remain missing

Independent2001–2017 ecological reference, full national physical distances/within-cell metrics, shifted response/climate reaggregation, corrected regional/pathway attribution, and direct irrigation/groundwater/soil controls are not supplied as completed experiments. Known sources and hashes are not substitutes for those outputs. The small worked example explains calculations; it does not verify sampling or effect size by itself.
