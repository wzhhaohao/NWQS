make_pr_fit_curve <- function(family = "gaussian", q = 4, rh = 5) {
  set.seed(2026)
  n <- 80
  mix <- data.frame(
    Component1 = rnorm(n),
    Component2 = rnorm(n),
    Component3 = rnorm(n)
  )
  beta <- c(0.6, 0.3, 0.1)
  eta  <- as.matrix(mix) %*% beta + rnorm(n, sd = 0.5)
  y <- switch(family,
    gaussian = as.numeric(eta),
    binomial = rbinom(n, 1, plogis(eta)),
    poisson  = rpois(n, exp(eta))
  )
  dat <- cbind(mix, y = y)
  nwqs(
    data = dat, mix_name = paste0("Component", 1:3),
    outcome = "y", family = family,
    transform_type = "percentile_rank", q = q,
    rh = rh, n_permutation = 5, seed = 1234, quiet = TRUE
  )
}

make_boot_fit <- function(family = "gaussian", q = 4, n_boot = 8) {
  set.seed(2026)
  n <- 80
  mix <- data.frame(
    Component1 = rnorm(n),
    Component2 = rnorm(n),
    Component3 = rnorm(n)
  )
  beta <- c(0.6, 0.3, 0.1)
  eta  <- as.matrix(mix) %*% beta + rnorm(n, sd = 0.5)
  y <- if (family == "gaussian") as.numeric(eta) else rbinom(n, 1, plogis(eta))
  dat <- cbind(mix, y = y)
  nwqs_boot(
    data = dat, mix_name = paste0("Component", 1:3),
    outcome = "y", family = family,
    transform_type = "percentile_rank", q = q,
    n_boot = n_boot, rh_inner = 1, n_permutation = 5,
    seed = 1234, quiet = TRUE
  )
}

test_that("nwqs_boot stores per-bootstrap shape coefficients (rh_shapes_boot)", {
  fit <- make_boot_fit(n_boot = 6)
  expect_true("rh_shapes_boot" %in% names(fit))
  expect_true(is.matrix(fit$rh_shapes_boot))
  expect_equal(nrow(fit$rh_shapes_boot), fit$n_success)
  expect_equal(ncol(fit$rh_shapes_boot),
               length(fit$final_weights) * fit$df_spline)
  expected_cols <- paste0(
    rep(names(fit$final_weights), each = fit$df_spline),
    "_B",
    rep(seq_len(fit$df_spline), times = length(fit$final_weights))
  )
  expect_setequal(expected_cols, colnames(fit$rh_shapes_boot))
})
