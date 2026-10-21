# Final analysis outputs

The numbered figures in `figures/` are the current figure set. Main figures are Fig1–Fig8; supplementary figures are Fig_S1–Fig_S13. Fig1_editable.pptx contains the methods diagram. Historical experiments elsewhere in this repository are retained for provenance and are not substitutes for these results.

## Reproduce figures

From the repository root, install requirements and run:

```sh
python src/final_release/plot_final_figures.py
```

This regenerates quantitative plots from full-precision derived CSVs without Earth Engine or the large source panels. Fig. 1 is an editable methods drawing; Fig. S1 is regenerated from the published coverage and composition tables. It does not perform new remote exports.

## Recalculate the held-out-year reference

Place the supplied aquatic panel at the path documented by `src/audit_pipeline.py`, and provide the patched 16-day core panel at `data/raw/china_panel_16day_core_patched.parquet` (or set RAW_CORE_PANEL). Then run:

```sh
python src/final_release/baseline_crossfit.py
```

The main reference is first reproduced; each target year is then excluded from the same grid/seasonal-position reference, with at least three remaining observations. Outputs include event agreement and three paired-sample AFB scenarios. This is not an independent longer climatology.

## Data and evidence

- `data/context`: nationwide risks and fragility.
- `data/distance`: distance bands, curves, contrasts, coverage, and numerical diagnostics. The authoritative final figures use the `dist_perm_corrected_km` rows; paired numerical-distance audit rows are retained as provenance.
- `data/adjusted`, `data/noncropland`, `data/covariate_balance.csv`: full model and landscape comparisons.
- `data/robustness`: every stress-index, threshold, year-omission and spatial-block estimate.
- `data/baseline`: target-year-excluded reference experiment and event-assignment agreement.
- `data/diagnostics`: component means, seasonal LAI changes and regional estimates.
- `data/regional_sample_support.csv`: original/corrected audit counts; final tables use classification `new`.

The main sample is 236,995 grid-years from 48,843 cells. The continuous-distance subset is 80,265 records from 17,041 cells. Held-out reference comparisons retain 209,985 records from 37,878 cells. All are observational landscape associations. Unavailable independent long-reference, within-cell and grid-realignment tests are not implied by these outputs. Upstream product metadata gaps are documented in the source inventory.

Manuscript and response documents are not published in this public repository. `manifest.csv` identifies derived files and their hashes.

## Model and diagnostics entry points

From the repository root, `python src/final_release/validate_distance_gradient.py` reconstructs distance bands, adjusted curves and spatial contrasts from the seasonal panel and retained distance exports. `python src/final_release/restore_core_experiments.py` recalculates exposure, seasonal-change and regional summaries. Both require the large seasonal input identified by `AQUATIC_PANEL`; they use existing distance and covariate exports and submit no remote tasks. Their output directories are `outputs/concern2_distance_validation` and `outputs/manuscript_revision_v4`.
