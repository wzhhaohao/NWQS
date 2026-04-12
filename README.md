# NWQS: Non-Linear Weighted Quantile Sum Regression

<!-- badges: start -->
<!-- badges: end -->

[English](README.md) | [中文](README_CN.md)

**NWQS** extends classical Weighted Quantile Sum (WQS) regression by incorporating natural cubic splines to capture non-linear dose-response relationships (e.g., threshold, U-shaped, S-shaped curves) in environmental mixture analyses.

## Key Features

- **Non-linear dose-response modeling** via natural cubic spline basis expansion, accommodating threshold, U-shaped, inverted-U, S-shaped, and other complex exposure-response patterns.
- **Permutation-based variable importance** using out-of-bag (OOB) loss changes to derive relative component weights without regularization bias.
- **Repeated holdout architecture** for robust weight and shape estimation across random data splits.
- **External bootstrap inference** (`nwqs_boot()`) providing valid percentile confidence intervals that reflect true sampling variability.
- **Multiple GLM families**: Gaussian, binomial, Poisson, quasi-Poisson, and conditional logistic regression (`clogit`).
- **Publication-quality visualizations** including dose-response curves, component weight bar charts, and bootstrap contrast boxplots.
- **Monte Carlo simulation toolkit** for benchmarking and validating model performance under controlled scenarios.
- **Automatic parallel computing** via the `future` framework with intelligent load balancing.

## Installation

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("wzhhaohao/NWQS")
```

## Quick Start

### 1. Simulate Data

```r
library(NWQS)

n_vars <- 7
mu <- rep(0, n_vars)
sigma <- generate_sigma(n_vars, mode = "medium", seed = 42)

dat <- gen_nonlinear_data(
  n_obs       = 500,
  mu_preds    = mu,
  sigma_preds = sigma,
  beta_wqs    = 1,
  beta_preds  = c(0.4, 0.3, 0.2, 0.1, 0, 0, 0),
  snr_db      = 10,
  shape       = c("linear_like", "threshold", "u_shape",
                   "s_shape", "linear_like", "linear_like", "linear_like"),
  seed        = 123
)
```

### 2. Fit NWQS Model

```r
mix_name   <- paste0("Component", 1:n_vars)
covariates <- c("x_cont", "x_bin", "x_cat")

fit <- nwqs(
  data           = dat,
  mix_name       = mix_name,
  covariates     = covariates,
  outcome        = "y",
  family         = "gaussian",
  q              = 4,
  df_spline      = 3,
  rh             = 10,
  n_permutation  = 50,
  plan_strategy  = "sequential",
  seed           = 1234
)

print(fit)
plot(fit)
```

### 3. Bootstrap Inference

```r
boot_res <- nwqs_boot(
  data           = dat,
  mix_name       = mix_name,
  covariates     = covariates,
  outcome        = "y",
  family         = "gaussian",
  n_boot         = 100,
  rh_inner       = 1,
  n_permutation  = 50,
  plan_strategy  = "multisession",
  seed           = 42
)

print(boot_res)
plot(boot_res)
```

## Supported Families

| Family | Outcome Type | Effect Scale |
|---|---|---|
| `gaussian` | Continuous | Linear (Delta Y) |
| `binomial` | Binary (0/1) | Log-Odds / Odds Ratio |
| `poisson` | Count | Log-Rate / Rate Ratio |
| `quasipoisson` | Overdispersed Count | Log-Rate / Rate Ratio |
| `clogit` | Matched Case-Control | Log-Odds / Odds Ratio |

## Core Functions

| Function | Description |
|---|---|
| `nwqs()` | Fit a Non-Linear WQS model with repeated holdout |
| `nwqs_boot()` | Bootstrap confidence intervals for valid inference |
| `nwqs_contrast()` | Joint exposure quantile contrast significance test |
| `extract_nwqs_effects()` | Extract component-specific quantile contrast effects |
| `plot.nwqs()` | Dose-response curves and weight diagnostics |
| `plot.nwqs_boot()` | Bootstrap contrast boxplots |

## Data Simulation

| Function | Description |
|---|---|
| `gen_nonlinear_data()` | Continuous outcome with non-linear spline effects |
| `gen_nonlinear_bio_data()` | Binary outcome with automatic prevalence calibration |
| `gen_nonlinear_count_data()` | Poisson count outcome |
| `generate_linear_data()` | Linear mixture effects (baseline comparison) |
| `generate_sigma()` | Correlation matrices (low/medium/high/mixed) |

## Method Overview

The NWQS algorithm follows a three-stage architecture:

1. **Weight and Shape Discovery** (Training Set): OOB permutation importance derives component weights; spline basis coefficients capture non-linear shapes.
2. **Shape Normalization**: Shapes are standardized to unit variance to decouple shape from effect magnitude.
3. **1-DoF Effect Estimation** (Validation Set): Normalized shapes and weights are projected into a single `nwqs` index, and a standard GLM estimates the overall mixture effect.

When `rh > 1`, the standard errors from `nwqs()` reflect only data-splitting variance. Use `nwqs_boot()` for valid statistical inference.

## License

GPL-3 License. See [LICENSE](LICENSE) for details.

## Authors

**Wang Zhehao** ([wangzhehao_bill@foxmail.com](mailto:wangzhehao_bill@foxmail.com)) — Department of Public Health and Preventive Medicine, School of Medicine, Jinan University, Guangzhou, 510632, China
**Chen Shirui** ([chenshr7@bjmu.edu.cn](mailto:chenshr7@bjmu.edu.cn)) — Peking University Center for Public Health and Epidemic Preparedness & Response, Peking University, Beijing 100191, China
**Lin Ziqiang** ([linziqiang0314@gmail.com](mailto:linziqiang0314@gmail.com)) — Department of Public Health and Preventive Medicine, School of Medicine, Jinan University, Guangzhou, 510632, China (Corresponding author)
