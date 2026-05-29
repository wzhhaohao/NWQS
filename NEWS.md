# NWQS 0.2.0 (in development)

## Breaking changes

Landed:

- **`family = "clogit"` has been removed** from `nwqs()` and `nwqs_boot()`. Conditional logistic regression has been lifted out of the main pipeline so the GLM path can be audited cleanly. A future minor release will reintroduce it as a dedicated `nwqs_clogit()` function. To keep an existing clogit workflow, pin to `0.1.0`.
- **`strata_col` parameter removed** from both `nwqs()` and `nwqs_boot()`. The argument has no meaning without the clogit family; calls that pass `strata_col = ...` are silently absorbed by `...` and have no effect.
- **`survival` package no longer in Imports.** Users who relied on `library(NWQS)` to load `survival` indirectly must now `library(survival)` themselves.
- **Default exposure transform changed from discrete `q`-bin (q = 4) to continuous percentile rank.** Each mixture column is now mapped by `u_i = rank(x_i, ties = "average") / n`, so `u_i ∈ (0, 1]`. The spline basis is evaluated on a 100-point grid over `[0, 1]`. This matches the applied-statistics convention and aligns NWQS with the percentile-rank mode of `gWQS`.
- **`trans_quantile()` signature changed.** The old `method = c("quantile", "percentile")` argument has been renamed to `type = c("percentile_rank", "q_bin")` and a new `ties = c("average", "min", "max", "random")` argument has been added. Calls using `method = ...` will now error with `unused argument`.
- **`nwqs()` and `nwqs_boot()` gained `transform_type` and `ties` parameters.** Defaults: `transform_type = "percentile_rank"`, `ties = "average"`. The `q` argument is retained (default 4): under `q_bin` it controls bin count; under `percentile_rank` it controls how many contrast points `extract_nwqs_effects()` and `nwqs_contrast()` evaluate.
- **Fit object gained three fields**: `transform_type`, `ties`, and `train_components_sorted` (the per-component sorted training values, used by the upcoming `predict.nwqs()` to map newdata onto the training empirical CDF). `nwqs_boot()` results also gained `transform_type` and `ties`.
- **New helpers in `R/utils.R`**: `apply_percentile_rank(newdata, train_x)` for mapping new data onto a training distribution, and `build_spline_basis_knots(transform_type, q, df_spline, custom_knots, custom_boundary)` for centralized, globally aligned spline knot construction.

Migration:

- Any 0.1.x call of the form
  ```r
  nwqs(data, mix_name, family = "clogit", strata_col = "match_id", ...)
  ```
  now errors with `'arg' should be one of "gaussian", "binomial", "poisson", "quasipoisson"`. To preserve the conditional logistic analysis, stay on 0.1.0 until `nwqs_clogit()` is released.
- To reproduce 0.1.x's discrete quartile behavior, pass
  ```r
  nwqs(data, mix_name, transform_type = "q_bin", q = 4, ...)
  ```
- Numerical results for the percentile-rank default will differ from 0.1.x because the underlying spline basis and `u_i` values differ. Final weights and the `nwqs` index coefficient are not directly comparable across the two transforms; interpret weights relative to each other within a single fit.

Planned for later in this release (none landed yet):

- `n_permutation` default will change from `10` to `30` so OOB importance estimates are more stable on small samples.
- "Soft" parameters such as `min_shape_sd`, `ties`, and `custom_knots` will move into a new `nwqs_control()` helper rather than living on the main `nwqs()` signature.
- All user-facing defaults will be sourced from a package-level `NWQS_DEFAULTS` list (`R/zzz-defaults.R`), so the same value is not declared twice in `nwqs()` and `nwqs_boot()`.

## New features

Landed:

- `tests/testthat/` suite with golden regression tests for both `q_bin` (locked to 0.1.x numerical baseline) and `percentile_rank` (new default), a clogit-removal regression test, and unit tests for the new transform helpers.
- GitHub Actions R-CMD-check workflow on macOS / Windows / Ubuntu × R-release / R-devel / R-oldrel-1.
- `apply_percentile_rank()` and `build_spline_basis_knots()` helpers (exported).

Planned for later in this release:

- `predict()`, `vcov()`, `confint()` methods for `nwqs` and `nwqs_boot` objects, plus conditional `broom::tidy()` / `broom::glance()` registrations.
- New families `negbin` (via `MASS::glm.nb`) and `ordinal` (via `MASS::polr`).
- Applied-domain vignette demonstrating an environmental-exposure-to-health-outcome workflow.
- pkgdown site.

# NWQS 0.1.0

Initial public release accompanying the thesis. Supports `gaussian`, `binomial`, `poisson`, `quasipoisson`, and `clogit` families; repeated-holdout estimation and external bootstrap inference; permutation-based variable importance; and S3 `print` / `summary` / `plot` / `coef` methods.
