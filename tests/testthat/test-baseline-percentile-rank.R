# Golden baseline tests for the new v0.2.0 default transform: continuous
# percentile rank (empirical CDF / n, ties = "average").
#
# These complement tests/testthat/test-baseline-golden.R, which locks in
# the legacy q_bin / q = 4 numerical baseline. Together they prove that
# future refactors do not silently drift either branch of the
# transform_type switch.

make_pr_data <- function(n = 200, p_mix = 5, family = "gaussian",
                         seed = 42) {
  set.seed(seed)
  rho <- 0.5
  Sigma <- matrix(rho, p_mix, p_mix)
  diag(Sigma) <- 1
  L <- chol(Sigma)
  X <- matrix(rnorm(n * p_mix), n, p_mix) %*% L
  colnames(X) <- paste0("X", seq_len(p_mix))

  eta <- 0.6 * X[, 1] - 0.4 * X[, 2] + 0.2 * X[, 3]^2

  y <- switch(
    family,
    gaussian     = eta + rnorm(n, sd = 0.5),
    binomial     = rbinom(n, size = 1, prob = plogis(eta - 0.5)),
    poisson      = rpois(n, lambda = exp(0.5 + 0.5 * eta)),
    quasipoisson = rpois(n, lambda = exp(0.5 + 0.5 * eta)),
    stop("unsupported family: ", family)
  )

  data.frame(y = y, X)
}

pr_mix_names <- paste0("X", 1:5)

expect_pr_contract <- function(fit, family) {
  expect_s3_class(fit, "nwqs")
  expect_equal(fit$family, family)
  expect_equal(fit$transform_type, "percentile_rank")
  expect_true(all(is.finite(fit$final_weights)))
  expect_true(all(fit$final_weights >= 0))
  expect_equal(sum(fit$final_weights), 1, tolerance = 1e-10)
  expect_setequal(names(fit$final_weights), pr_mix_names)
  expect_equal(fit$spline_boundary, c(0, 1))
  expect_named(fit$train_components_sorted, pr_mix_names)
}

run_pr_baseline <- function(family) {
  d <- make_pr_data(family = family, seed = 42)
  nwqs(
    data           = d,
    mix_name       = pr_mix_names,
    outcome        = "y",
    family         = family,
    rh             = 5,
    n_permutation  = 5,
    seed           = 1234,
    quiet          = TRUE,
    transform_type = "percentile_rank",
    q              = 4,
    ties           = "average"
  )
}

# ----- Per-family snapshots ----------------------------------------------

test_that("percentile_rank baseline: gaussian", {
  fit <- run_pr_baseline("gaussian")
  expect_pr_contract(fit, "gaussian")
  expect_snapshot_value(round(fit$final_weights, 6), style = "json2")
  expect_snapshot_value(round(fit$mean_coefs, 6), style = "json2")
})

test_that("percentile_rank baseline: binomial", {
  fit <- run_pr_baseline("binomial")
  expect_pr_contract(fit, "binomial")
  expect_snapshot_value(round(fit$final_weights, 6), style = "json2")
  expect_snapshot_value(round(fit$mean_coefs, 6), style = "json2")
})

test_that("percentile_rank baseline: poisson", {
  fit <- run_pr_baseline("poisson")
  expect_pr_contract(fit, "poisson")
  expect_snapshot_value(round(fit$final_weights, 6), style = "json2")
  expect_snapshot_value(round(fit$mean_coefs, 6), style = "json2")
})

test_that("percentile_rank baseline: quasipoisson", {
  fit <- run_pr_baseline("quasipoisson")
  expect_pr_contract(fit, "quasipoisson")
  expect_snapshot_value(round(fit$final_weights, 6), style = "json2")
  expect_snapshot_value(round(fit$mean_coefs, 6), style = "json2")
})

# ----- nwqs_boot path under percentile_rank ------------------------------

test_that("percentile_rank baseline: gaussian nwqs_boot returns valid CI structure", {
  d <- make_pr_data(family = "gaussian", seed = 42)
  boot_fit <- nwqs_boot(
    data           = d,
    mix_name       = pr_mix_names,
    outcome        = "y",
    family         = "gaussian",
    n_boot         = 10,
    rh_inner       = 1,
    n_permutation  = 5,
    seed           = 1234,
    quiet          = TRUE,
    transform_type = "percentile_rank",
    q              = 4
  )
  expect_s3_class(boot_fit, "nwqs_boot")
  expect_equal(boot_fit$transform_type, "percentile_rank")
  expect_equal(boot_fit$family, "gaussian")
  expect_true(boot_fit$n_success > 0)
  expect_true(all(c("Term", "Target", "Boot_Mean",
                    "Boot_CI_Lower", "Boot_CI_Upper") %in%
                  names(boot_fit$ci_table)))
  expect_true(all(boot_fit$ci_table$Boot_CI_Lower <=
                  boot_fit$ci_table$Boot_Mean))
  expect_true(all(boot_fit$ci_table$Boot_Mean <=
                  boot_fit$ci_table$Boot_CI_Upper))
  expect_equal(sum(boot_fit$final_weights), 1, tolerance = 1e-10)
})
