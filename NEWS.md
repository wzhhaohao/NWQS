# NWQS 0.2.0

## Breaking changes

Landed:

- **`family = "clogit"` has been removed** from `nwqs()` and `nwqs_boot()`. Conditional logistic regression has been lifted out of the main pipeline so the GLM path can be audited cleanly. A future minor release will reintroduce it as a dedicated `nwqs_clogit()` function. To keep an existing clogit workflow, pin to `0.1.0`.
- **`strata_col` parameter removed** from both `nwqs()` and `nwqs_boot()`. The argument has no meaning without the clogit family; calls that pass `strata_col = ...` are silently absorbed by `...` and have no effect.
- **`survival` package no longer in Imports.** Users who relied on `library(NWQS)` to load `survival` indirectly must now `library(survival)` themselves.
- **Default exposure transform changed from discrete `q`-bin (q = 4) to continuous percentile rank.** Each mixture column is now mapped by `u_i = rank(x_i, ties = "average") / n`, so `u_i ∈ (0, 1]`. The spline basis is evaluated on a 100-point grid over `[0, 1]`. This matches the applied-statistics convention and aligns NWQS with the percentile-rank mode of `gWQS`.
- **Nonlinear simulation generators now default to a percentile-rank DGP.** `gen_nonlinear_data()`, `gen_nonlinear_bio_data()`, and `gen_nonlinear_count_data()` transform mixture columns with the same percentile-rank convention used by `nwqs()` unless a custom `transform_fun` is supplied. They also attach `true_effect_curve`, a full link-scale dose-response truth matrix over `0, 0.01, ..., 1`, so simulation benchmarks can recover weights from the complete percentile curve rather than only Q4 vs Q1.
- **`trans_quantile()` signature changed.** The old `method = c("quantile", "percentile")` argument has been renamed to `type = c("percentile_rank", "q_bin")` and a new `ties = c("average", "min", "max", "random")` argument has been added. Calls using `method = ...` will now error with `unused argument`.
- **`nwqs()` and `nwqs_boot()` gained `transform_type` and `ties` parameters.** Defaults: `transform_type = "percentile_rank"`, `ties = "average"`. The `q` argument is retained (default 4): under `q_bin` it controls bin count; under `percentile_rank` it controls how many contrast points `extract_nwqs_effects()` and `nwqs_contrast()` evaluate.
- **Fit object gained three fields**: `transform_type`, `ties`, and `train_components_sorted` (the per-component sorted fit-sample values, used by `predict.nwqs()` to map `newdata` onto the fitted empirical CDF). `nwqs_boot()` results also gained `transform_type` and `ties`.
- **New helpers in `R/utils.R`**: `apply_percentile_rank(newdata, train_x, ties)` for mapping new data onto the fitted empirical distribution, and `build_spline_basis_knots(transform_type, q, df_spline, custom_knots, custom_boundary)` for centralized, globally aligned spline knot construction.

Migration:

- Any 0.1.x call of the form
  ```r
  nwqs(data, mix_name, family = "clogit", strata_col = "match_id", ...)
  ```
  now errors during family validation because `clogit` is no longer part of the supported family set. To preserve the conditional logistic analysis, stay on 0.1.0 until `nwqs_clogit()` is released.
- To reproduce 0.1.x's discrete quartile behavior, pass
  ```r
  nwqs(data, mix_name, transform_type = "q_bin", q = 4, ...)
  ```
- To reproduce the legacy q-bin nonlinear simulation DGP, pass
  ```r
  gen_nonlinear_data(..., transform_type = "q_bin", q = 4)
  ```
- Numerical results for the percentile-rank default will differ from 0.1.x because the underlying spline basis and `u_i` values differ. Final weights and the `nwqs` index coefficient are not directly comparable across the two transforms; interpret weights relative to each other within a single fit.

Landed in Phase 4 (E): centralized defaults + `nwqs_control()`

- **`NWQS_DEFAULTS` package-level list** (`R/zzz-defaults.R`) holds every user-facing default value (`q`, `df_spline`, `transform_type`, `ties`, `train_prop`, `rh`, `n_permutation`, `n_boot`, `rh_inner`, `conf_level`, `seed`, `min_shape_sd`, `zero_weight_action`). The formals on `nwqs()` and `nwqs_boot()` read from `NWQS_DEFAULTS$<key>` rather than embedding literals, so a default cannot drift between the two signatures.
- **`nwqs_control()`** packages "advanced" knobs that previously had no home: `custom_knots`, `custom_boundary`, and `zero_weight_action`. The function validates each argument (numeric, length, allowed values) and returns an `nwqs_control` object. Pass it to `nwqs()` and `nwqs_boot()` via the new `control` argument. Existing soft parameters that already lived on the main signature (`min_shape_sd`, `ties`) are not moved into control in 0.2.0 — keeping them on the signature avoids a needless breaking change.
- **`nwqs()` and `nwqs_boot()` gain a `control = nwqs_control()` argument** that is forwarded to `build_spline_basis_knots()`, so `custom_knots` and `custom_boundary` actually override the derived spline knots. `nwqs_boot()` passes `control` through to its inner `nwqs()` call.
- Existing test `test-weights-engine.R` and `test-nwqs-internals.R` updated to `eval(formals(nwqs)$<key>)` because the formals are now language objects referring to the `NWQS_DEFAULTS` list.

Landed in Phase 4 (C): family = "negbin" via MASS::glm.nb

- **`nwqs()` and `nwqs_boot()` now accept `family = "negbin"`** for overdispersed count outcomes. The final NWQS regression on the validation split is fit with `MASS::glm.nb()`; for `rh = 1` the returned `$model_obj` is a `MASS::glm.nb` object, so `predict()`, `vcov()`, `confint()`, and `broom::tidy()` / `glance()` all work without any further wiring.
- **OOB importance for negbin uses a Poisson surrogate** inside `permutation_scorer()` (via `run_oob_permutation()`). The Poisson and negative binomial deviances differ by an additive theta-dependent constant that drops out of the per-component permutation delta, so the ranking is unchanged; this avoids paying the cost of a theta search on every bootstrap iteration of the inner scorer.
- **`gen_nbin_data(n_obs, n_vars, beta, intercept, theta, snr_db, seed)`** added to `R/simulate_data.R`. Returns counts with columns `y, V1, V2, ...`; used as the fixture in `tests/testthat/test-family-negbin.R`.
- The `is_exp_family` test used by `print()` / `summary()` / `plot()` to decide between additive (`Delta`) and exponentiated (`OR / RR`) display now includes `"negbin"` and reports its effects on the rate-ratio scale, matching Poisson.

Landed in Phase 4 (A+B): standard S3 methods

- **`predict.nwqs(object, newdata, type)`** and **`predict.nwqs_boot(object, newdata, type)`** added. `type` is one of `"response"` (default), `"link"`, or `"nwqs_index"`. The fit-sample empirical CDF stored as `train_components_sorted` is reused for percentile-rank transforms; the fit-sample quantile breaks are reused for `q_bin` transforms. The same globally aligned spline knots that were used during fitting are reused so newdata is scored on the same scale.
- **`vcov.nwqs()` and `confint.nwqs()`** added. For `rh = 1` they wrap the stored inner GLM directly (real sampling-variance inference); for `rh > 1` they emit the standard "rh > 1 reflects algorithmic variance only" warning and return the RH-derived covariance / quantile CI.
- **`vcov.nwqs_boot()` and `confint.nwqs_boot()`** added. They return the bootstrap-iteration covariance and the percentile bootstrap CIs, respectively — both real sampling-variance inference.
- **Conditional `broom::tidy()` / `broom::glance()` registration** via `.onLoad`. NWQS does not add `broom` to Imports; the methods register only if `generics` is available at load time. Returned objects are plain `data.frame`s in the canonical broom column layout (`term / estimate / std.error / statistic / p.value`, optional `conf.low / conf.high`).
- **Fit object enrichment**: `nwqs()` now stores `formula` (the regression formula) and, for `rh = 1`, `model_obj` (the inner GLM). `nwqs_boot()` stores `mean_coefs` (boot-averaged regression coefs), `rh_coefs_boot` / `rh_weights_boot` (per-iteration matrices), `train_components_sorted`, `mix_name`, `covariates`, `outcome`, and `formula` so that all of the new S3 methods can run without re-fitting.
- **Version / metadata update**: `DESCRIPTION` now targets `0.2.0`. `Suggests:` gained `broom`, `generics`, `knitr`, and `rmarkdown`; `VignetteBuilder: knitr` was added; and the package description was updated to match the current family set.

Landed in Phase 3:

- **`n_permutation` default raised from `10` to `30`** across `nwqs()`, `nwqs_boot()`, `permutation_scorer()`, and `run_oob_permutation()`. OOB importance estimates are more stable on small samples; the additional compute is bounded (permutation is the inner loop of an already-OOB-bounded scorer).
- **`permutation_scorer()` now reports in-bag rank deficiency instead of silently filling NA coefficients with zero.** When `glm.fit` returns any NA coefficient on a bootstrap sample, the iteration emits a warning ("rank deficiency; iteration skipped") and returns `NULL`. The outer RH loop ignores `NULL` iterations, matching the existing failure-tolerant aggregation logic. Side effect: the golden q_bin and percentile_rank snapshot values changed slightly (<5%) because previously silently-zeroed iterations are now skipped. The updated snapshots are the new ground truth.
- **`min_shape_sd` parameter exposed on `nwqs()`** (default `1e-8`, matching the prior hard-coded literal). When a component's training-set partial linear-predictor standard deviation drops below this threshold, the per-component shape normalization is bypassed and, under `quiet = FALSE`, a `message()` names the component and RH iteration.
- **`add_noise_by_snr()` documented as link-scale SNR** with the explicit formula `SNR = Var(Xβ) / Var(ε)`. A new test (`tests/testthat/test-snr.R`) pins `Var(η) / Var(noise) ≈ 10^(snr_db / 10)` to within 5% across multiple `snr_db` levels.
- **`.calc_loss(y_true, mu_pred, fam_name)` extracted** as an internal function so the Poisson `mu → 0` clipping and the binomial `p → 0/1` clipping have direct unit-test coverage.

Landed in Phase 6: percentile-rank display unification + effect-curve API

- **Display uniformly switches from Q-labels to P-labels in percentile_rank mode.** `nwqs()` and `nwqs_boot()` no longer emit `Q*_vs_Q*` strings when `transform_type = "percentile_rank"`: `extract_nwqs_effects()`, `nwqs_contrast()`, `print.nwqs()`, `print.nwqs_boot()`, `summary.nwqs_boot()`, `plot.nwqs()`, and `plot_nwqs_contrast_box()` all label contrasts as `P0`, `P25`, `P50`, `P75`, `P95`, `P100` etc. `q_bin` fits keep every legacy `Q*_vs_Q1` label.
- **`extract_nwqs_effects()` gains `contrast_points`, `ref`, and `label_style`.** Default behavior (no args) preserves the legacy numeric grid; users can now request arbitrary contrasts via `contrast_points = c(0.25, 0.75), ref = 0.5`. Validation rejects points outside `[0, 1]` for percentile_rank. Now also accepts `nwqs_boot` objects (S3 dispatch) and returns bootstrap-percentile CIs from `rh_shapes_boot`/`rh_weights_boot`/`rh_coefs_boot`.
- **`nwqs_contrast()` gains `target`, `ref`, and `label_style`.** Legacy `q_target`/`q_ref` continue to work for q_bin backward compatibility; supplying both `target` and `q_target` raises a clear warning and prefers the new arg. percentile_rank default is `target = 0.75, ref = 0.5` (paper-style P75 vs P50). The printed header now says "Percentile-rank Contrast" or "Quantile Contrast" based on `transform_type`. The invisible return list gains `target`, `ref`, `target_label`, `ref_label`, and `transform_type`.
- **`print.nwqs()` gains `contrast_points`, `ref`, and `label_style`.** In percentile_rank mode the default table columns are now `P25 vs P50`, `P75 vs P50`, `P95 vs P50` (paper-style IQR + extreme upper). q_bin keeps the legacy `Q2 vs Q1` ... `Q{q} vs Q1` columns.
- **`summary.nwqs_boot()` no longer hard-codes `Q*_vs_Q1` for its stability table.** It now reads the largest target dynamically from `boot_table$Target` and labels the per-component effect SD with the actual contrast string (e.g., `P100_vs_P0_Effect_SD`).
- **`plot_nwqs_contrast_box()` and `plot.nwqs()` x-axis labels are transform-aware.** percentile_rank fits show "Joint exposure percentile rank"; q_bin fits keep "Exposure Quantile Index" / "Quantile Index". Boxplot facet ordering is derived from the numeric percentile/quantile parsed from each `Target` string instead of a hard-coded `Q2:Qn` factor.
- **New `extract_nwqs_effect_curve(model, grid, ref, include_components, label_style)`.** Returns a tidy `data.frame` with `term`, `x`, `ref`, `estimate`, `lower`, `upper`, `transform_type`, and `inference_type` of the joint and per-component partial-effect curve relative to the chosen reference percentile. For `nwqs` inputs the band is the RH empirical quantile (algorithmic variance only, `inference_type = "repeated_holdout"`, NA when `rh = 1`); for `nwqs_boot` inputs the band is the bootstrap percentile CI (`inference_type = "bootstrap"`).
- **New `plot_nwqs_effect_curve(model, grid, ref, include_components, label_style, base_size)`.** Wraps the extractor into a ggplot. Median-centered by default (`ref = 0.5`); set `include_components = TRUE` for component overlays. The y-axis label reports the actual reference, e.g., "Partial effect change relative to P50".
- **`nwqs_boot()` stores `rh_shapes_boot`**, a per-replicate matrix of shape coefficients. This is required for bootstrap-percentile CIs on the effect curve and for `extract_nwqs_effects.nwqs_boot()` to recompute contrasts at user-specified points.
- **`NWQS_DEFAULTS` gains** `contrast_target_pr = 0.75`, `contrast_ref_pr = 0.5`, `contrast_ref_q_bin = 0`, `contrast_points_print_pr = c(0.25, 0.75, 0.95)`, `effect_curve_grid = seq(0, 1, by = 0.01)`, `label_style = "auto"`. Every new default is the single source of truth referenced by the relevant function `formals`.

Recommended usage:

- Paper main-table contrast: `nwqs_contrast(fit, target = 0.75, ref = 0.25)` (IQR).
- Paper main-figure curve: `plot_nwqs_effect_curve(fit, ref = 0.5)` (median-centered).
- The extreme `P100 vs P0` contrast is allowed but not recommended as a paper main result.

Landed in Phase 5: applied documentation + pkgdown

- **Applied vignette**: `vignettes/nwqs-applied.Rmd` walks through an environmental-mixture example end to end: simulate a metals cohort, fit `nwqs()`, obtain valid inference with `nwqs_boot()`, score `newdata` with `predict()`, and show a count-outcome example with `family = "negbin"`.
- **pkgdown configuration**: `_pkgdown.yml` now organizes reference topics, articles, and site navigation around the 0.2.0 API. A matching GitHub Actions workflow (`.github/workflows/pkgdown.yaml`) builds and deploys the site to GitHub Pages.

## New features

Landed:

- `tests/testthat/` suite with golden regression tests for both `q_bin` (locked to 0.1.x numerical baseline) and `percentile_rank` (new default), a clogit-removal regression test, and unit tests for the new transform helpers.
- GitHub Actions R-CMD-check workflow on macOS / Windows / Ubuntu × R-release / R-devel / R-oldrel-1.
- `apply_percentile_rank()` and `build_spline_basis_knots()` helpers (exported).

Deferred beyond 0.2.0:

- `family = "ordinal"` via `MASS::polr` remains intentionally deferred to a separate plan.

# NWQS 0.1.0

Initial public release accompanying the thesis. Supports `gaussian`, `binomial`, `poisson`, `quasipoisson`, and `clogit` families; repeated-holdout estimation and external bootstrap inference; permutation-based variable importance; and S3 `print` / `summary` / `plot` / `coef` methods.
