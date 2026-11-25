# Sampling and spatial sensitivity analyses

These analyses use the same main panel and national exposure thresholds. They evaluate adjusted spatial uncertainty, residual autocorrelation, minimum seasonal observation counts, held-out-reference sample retention, distance covariate support and water-seasonality processing.

Run scripts from the repository root with `AQUATIC_PANEL` pointing to the seasonal input. `count_validation.py` additionally requires `RAW_CORE_PANEL`, the patched 16-day input. The large record-level validation panel is not distributed. Numerical tables and diagnostic figures are provided here.

- `spatial_models.py`: adjusted cluster-sandwich contrasts with t intervals and residual Moran diagnostics.
- `count_validation.py`: seasonal counts, two/three-step restrictions and held-out sample retention.
- `common2000.py`: 2,000 common-landscape draws, including invalid-cell diagnostics.
- `analyze_jrc.py`: national contrasts for corrected 2021 seasonality using retained 1984–2021 occurrence.

Spatial widths are sensitivity settings, not estimated correlation ranges. Minimum-count restrictions retain the primary national patterns. Residual spatial dependence remains, and several adjusted spatial intervals include zero. The seasonality comparison changes 140 of 48,843 classes; it does not test the v1.5 occurrence correction. No independent longer vegetation reference or land-only reaggregation is claimed.

Tables retain full precision. The seasonal count distributions describe primary response/exposure variables. Land-cover codes use MCD12Q1 IGBP categories. Source panels remain subject to separate data access.
