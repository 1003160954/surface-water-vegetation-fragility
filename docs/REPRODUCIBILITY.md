# Reproducibility and analytical scope

The analysis follows the sequence of source identification, seasonal anomaly calculation, exposure and event classification, aquatic-context assignment, risk estimation, covariate adjustment, distance modelling, and sensitivity analysis. Numerical tables retain full precision; figure labels are rounded for display.

Run `python src/final_release/plot_final_figures.py` from the repository root to generate quantitative figures from the published derived tables. Model fitting and anomaly reconstruction additionally require the input panels identified in `DATA_ACCESS.md`. Source-only notebooks preserve the available extraction code. Some upstream product versions and quality-mask details remain unverified; the source inventory identifies these gaps.

Grid-cluster and spatial-block uncertainty, alternative stress indices and event thresholds, agricultural exclusions, common-landscape comparisons and target-year-excluded references assess specified sources of sensitivity. They do not guarantee removal of all confounding. Irrigation, groundwater and soil properties are not directly measured. The distance analysis is geographically restricted and describes grid-level associations, not shoreline buffering ranges. An independent longer anomaly reference and nationwide within-cell and grid-alignment tests are not available.

Published material comprises code, source metadata, derived tables and figures. Large input panels require separate access. Cloud export scripts retain a pause safeguard; plotting published results does not submit remote tasks. AI-assisted tools were used for language polishing and assistance with writing portions of the analysis code.
