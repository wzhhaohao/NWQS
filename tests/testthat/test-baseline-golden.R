# Baseline golden regression tests for the v0.2.0 refactor.
#
# These tests capture the numerical behavior of nwqs() and nwqs_boot() on
# every family that will remain supported after Phase 1 (which removes
# clogit). They are the principal safety net protecting the GLM path from
# regressions during the refactor.
#
# After clogit is removed, after the default transform changes to
# percentile rank, after nwqs_control() is introduced, etc., these tests
# may legitimately need their snapshots updated. Each such update MUST
# be paired with a NEWS.md entry explaining the change.
#
# On the first run of devtools::test(), expect_snapshot_value() captures
# the current values under tests/testthat/_snaps/. Subsequent runs
# compare against the snapshot.

make_baseline_data <- function(n = 200, p_mix = 5, family = "gaussian",
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
    stop("unsupported family in make_baseline_data: ", family)
  )

  data.frame(y = y, X)
}

mix_names <- paste0("X", 1:5)

# ----- Contract checks shared across families -----------------------------

expect_nwqs_contract <- function(fit, family) {
  expect_s3_class(fit, "nwqs")
  expect_equal(fit$family, family)
  expect_true(all(is.finite(fit$final_weights)))
  expect_true(all(fit$final_weights >= 0))
  expect_equal(sum(fit$final_weights), 1, tolerance = 1e-10)
  expect_equal(length(fit$final_weights), length(mix_names))
  expect_setequal(names(fit$final_weights), mix_names)
  expect_true(all(is.finite(fit$mean_coefs)))
  expect_true("nwqs" %in% names(fit$mean_coefs))
}

run_nwqs_baseline <- function(family) {
  d <- make_baseline_data(family = family, seed = 42)
  nwqs(
    data    = d,
    mix_name = mix_names,
    outcome = "y",
    family  = family,
    rh      = 5,
    n_permutation = 5,
    seed    = 1234,
    quiet   = TRUE
  )
}

# ----- Per-family golden tests --------------------------------------------

test_that("baseline: gaussian nwqs has stable weights and coefs", {
  fit <- run_nwqs_baseline("gaussian")
  expect_nwqs_contract(fit, "gaussian")
  expect_snapshot_value(round(fit$final_weights, 6), style = "json2")
  expect_snapshot_value(round(fit$mean_coefs, 6), style = "json2")
})

test_that("baseline: binomial nwqs has stable weights and coefs", {
  fit <- run_nwqs_baseline("binomial")
  expect_nwqs_contract(fit, "binomial")
  expect_snapshot_value(round(fit$final_weights, 6), style = "json2")
  expect_snapshot_value(round(fit$mean_coefs, 6), style = "json2")
})

test_that("baseline: poisson nwqs has stable weights and coefs", {
  fit <- run_nwqs_baseline("poisson")
  expect_nwqs_contract(fit, "poisson")
  expect_snapshot_value(round(fit$final_weights, 6), style = "json2")
  expect_snapshot_value(round(fit$mean_coefs, 6), style = "json2")
})

test_that("baseline: quasipoisson nwqs has stable weights and coefs", {
  fit <- run_nwqs_baseline("quasipoisson")
  expect_nwqs_contract(fit, "quasipoisson")
  expect_snapshot_value(round(fit$final_weights, 6), style = "json2")
  expect_snapshot_value(round(fit$mean_coefs, 6), style = "json2")
})

# ----- Bootstrap path -----------------------------------------------------

test_that("baseline: gaussian nwqs_boot returns valid CI structure", {
  d <- make_baseline_data(family = "gaussian", seed = 42)
  boot_fit <- nwqs_boot(
    data     = d,
    mix_name = mix_names,
    outcome  = "y",
    family   = "gaussian",
    n_boot   = 10,
    rh_inner = 1,
    n_permutation = 5,
    seed     = 1234,
    quiet    = TRUE
  )
  expect_s3_class(boot_fit, "nwqs_boot")
  expect_equal(boot_fit$family, "gaussian")
  expect_true(boot_fit$n_success > 0)

  expect_true(all(c("Term", "Target", "Boot_Mean",
                    "Boot_CI_Lower", "Boot_CI_Upper") %in%
                  names(boot_fit$ci_table)))
  expect_true(nrow(boot_fit$ci_table) > 0)
  expect_true(all(is.finite(boot_fit$ci_table$Boot_Mean)))
  expect_true(all(boot_fit$ci_table$Boot_CI_Lower <=
                  boot_fit$ci_table$Boot_Mean))
  expect_true(all(boot_fit$ci_table$Boot_Mean <=
                  boot_fit$ci_table$Boot_CI_Upper))

  expect_true(all(is.finite(boot_fit$final_weights)))
  expect_equal(sum(boot_fit$final_weights), 1, tolerance = 1e-10)
  expect_setequal(names(boot_fit$final_weights), mix_names)
})
