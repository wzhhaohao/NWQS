# NWQS 0.2.0 (in development)

## Breaking changes

Planned for this release (none landed yet):

- `family = "clogit"` will be removed from `nwqs()` and `nwqs_boot()`. Conditional logistic regression support is being lifted out of the main pipeline so the GLM path can be audited cleanly. It will return later as a dedicated `nwqs_clogit()` function. To stay on `clogit`, pin to `0.1.0`.
- The default exposure transform will change from discrete `q`-bin quantiles (`q = 4`) to a continuous percentile rank (empirical CDF on the training distribution, ties handled by average rank). The `q`-bin transform remains available as an opt-in via `transform_type = "q_bin"`, `q = 4`.
- `n_permutation` default will change from `10` to `30` so OOB importance estimates are more stable on small samples.
- "Soft" parameters such as `min_shape_sd`, `ties`, and `custom_knots` will move into a new `nwqs_control()` helper rather than living on the main `nwqs()` signature.
- All user-facing defaults will be sourced from a package-level `NWQS_DEFAULTS` list (`R/zzz-defaults.R`), so the same value is not declared twice in `nwqs()` and `nwqs_boot()`.

A migration section will be added here once each change lands.

## New features

Planned (none landed yet):

- `predict()`, `vcov()`, `confint()` methods for `nwqs` and `nwqs_boot` objects, plus conditional `broom::tidy()` / `broom::glance()` registrations.
- New families `negbin` (via `MASS::glm.nb`) and `ordinal` (via `MASS::polr`).
- Applied-domain vignette demonstrating an environmental-exposure-to-health-outcome workflow.
- `tests/testthat/` suite with golden regression tests, parallel-plan contract tests, and per-family contract tests.
- GitHub Actions R-CMD-check workflow on macOS / Windows / Ubuntu × R-release / R-devel.
- pkgdown site.

## Bug fixes

Planned (none landed yet):

- When `rh = 1`, `nwqs()` did not store `spline_knots` and `spline_boundary` on the returned object, which broke `plot.nwqs()` for single-holdout fits. The two `rh = 1` / `rh > 1` assembly branches will be unified.

# NWQS 0.1.0

Initial public release accompanying the thesis. Supports `gaussian`, `binomial`, `poisson`, `quasipoisson`, and `clogit` families; repeated-holdout estimation and external bootstrap inference; permutation-based variable importance; and S3 `print` / `summary` / `plot` / `coef` methods.
