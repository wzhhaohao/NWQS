# Contract tests for vcov / confint methods on nwqs and nwqs_boot.
#
# Invariants pinned here:
#   - rh == 1: vcov.nwqs and confint.nwqs wrap the inner GLM directly so a
#     practitioner gets sampling-variance inference with no extra wiring.
#   - rh > 1: vcov.nwqs and confint.nwqs emit the "rh > 1 is not valid
#     sampling-variance inference; use nwqs_boot()" warning AND still
#     return the algorithmic-variance covariance / CI so existing
#     downstream tooling does not crash.
#   - nwqs_boot: vcov is the per-iteration bootstrap covariance,
#     confint is read from the ci_table (already a percentile bootstrap CI).

make_vc_data <- function(n = 200, family = "gaussian", seed = 31) {
  set.seed(seed)
  X <- matrix(rnorm(n * 3), n, 3)
  colnames(X) <- paste0("X", 1:3)
  eta <- 0.6 * X[, 1] - 0.3 * X[, 2]
  y <- switch(
    family,
    gaussian     = eta + rnorm(n, sd = 0.5),
    binomial     = rbinom(n, 1, plogis(eta)),
    poisson      = rpois(n, exp(0.5 + 0.5 * eta)),
    quasipoisson = rpois(n, exp(0.5 + 0.5 * eta))
  )
  data.frame(y = y, X)
}

# ----- vcov.nwqs ---------------------------------------------------------

test_that("vcov.nwqs with rh = 1 wraps the inner GLM vcov() matrix", {
  d <- make_vc_data()
  fit <- nwqs(
    data = d, mix_name = c("X1", "X2", "X3"), outcome = "y",
    family = "gaussian", rh = 1, n_permutation = 3, seed = 1,
    quiet = TRUE, transform_type = "q_bin", q = 4
  )
  v <- vcov(fit)
  expect_true(is.matrix(v))
  expect_equal(rownames(v), colnames(v))
  expect_true("nwqs" %in% rownames(v))
  expect_true(all(diag(v) >= 0))
})

test_that("vcov.nwqs with rh > 1 warns about algorithmic variance and returns the RH covariance", {
  d <- make_vc_data()
  fit <- nwqs(
    data = d, mix_name = c("X1", "X2", "X3"), outcome = "y",
    family = "gaussian", rh = 5, n_permutation = 3, seed = 1,
    quiet = TRUE, transform_type = "q_bin", q = 4
  )
  expect_warning(v <- vcov(fit), regexp = "rh > 1")
  expect_true(is.matrix(v))
  expect_true("nwqs" %in% rownames(v))
  expect_true(all(diag(v) >= 0))
})

# ----- vcov.nwqs_boot ----------------------------------------------------

test_that("vcov.nwqs_boot is the bootstrap-iteration covariance of regression coefs", {
  d <- make_vc_data()
  boot_fit <- nwqs_boot(
    data = d, mix_name = c("X1", "X2", "X3"), outcome = "y",
    family = "gaussian", n_boot = 8, rh_inner = 1, n_permutation = 3,
    seed = 1, quiet = TRUE, transform_type = "q_bin", q = 4
  )
  v <- vcov(boot_fit)
  expect_true(is.matrix(v))
  expect_equal(rownames(v), colnames(v))
  expect_true("nwqs" %in% rownames(v))
  expect_true(all(diag(v) >= 0))
})

# ----- confint.nwqs ------------------------------------------------------

test_that("confint.nwqs with rh = 1 wraps the inner GLM confint() matrix", {
  d <- make_vc_data()
  fit <- nwqs(
    data = d, mix_name = c("X1", "X2", "X3"), outcome = "y",
    family = "gaussian", rh = 1, n_permutation = 3, seed = 1,
    quiet = TRUE, transform_type = "q_bin", q = 4
  )
  ci <- suppressMessages(confint(fit))
  expect_true(is.matrix(ci))
  expect_true("nwqs" %in% rownames(ci))
  expect_true(ci["nwqs", 1] <= ci["nwqs", 2])
})

test_that("confint.nwqs with rh > 1 warns about algorithmic variance and returns RH-derived CI", {
  d <- make_vc_data()
  fit <- nwqs(
    data = d, mix_name = c("X1", "X2", "X3"), outcome = "y",
    family = "gaussian", rh = 5, n_permutation = 3, seed = 1,
    quiet = TRUE, transform_type = "q_bin", q = 4
  )
  expect_warning(ci <- confint(fit), regexp = "rh > 1")
  expect_true(is.matrix(ci))
  expect_true("nwqs" %in% rownames(ci))
  expect_true(ci["nwqs", 1] <= ci["nwqs", 2])
})

# ----- confint.nwqs_boot -------------------------------------------------

test_that("confint.nwqs_boot returns the bootstrap percentile CI from ci_table", {
  d <- make_vc_data()
  boot_fit <- nwqs_boot(
    data = d, mix_name = c("X1", "X2", "X3"), outcome = "y",
    family = "gaussian", n_boot = 8, rh_inner = 1, n_permutation = 3,
    seed = 1, quiet = TRUE, transform_type = "q_bin", q = 4
  )
  ci <- confint(boot_fit)
  expect_true(is.matrix(ci))
  expect_true("nwqs" %in% rownames(ci))
  expect_true(ci["nwqs", 1] <= ci["nwqs", 2])
})
