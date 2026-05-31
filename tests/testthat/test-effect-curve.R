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
  expect_small_boot_warning(nwqs_boot(
    data = dat, mix_name = paste0("Component", 1:3),
    outcome = "y", family = family,
    transform_type = "percentile_rank", q = q,
    n_boot = n_boot, rh_inner = 1, n_permutation = 5,
    seed = 1234, quiet = TRUE
  ))
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

test_that("extract_nwqs_effect_curve.nwqs returns tidy frame with required columns", {
  fit <- make_pr_fit_curve(q = 4, rh = 5)
  curve <- extract_nwqs_effect_curve(fit, ref = 0.5)
  required <- c("term", "x", "ref", "estimate", "lower", "upper",
                "transform_type", "inference_type")
  expect_true(all(required %in% names(curve)))
  expect_equal(unique(curve$transform_type), "percentile_rank")
  expect_equal(unique(curve$inference_type), "repeated_holdout")
  expect_true("Overall" %in% unique(curve$term))
  expect_true(all(paste0("Component", 1:3) %in% unique(curve$term)))
})

test_that("extract_nwqs_effect_curve.nwqs gives estimate=0 at x=ref for every term", {
  fit <- make_pr_fit_curve(q = 4, rh = 5)
  curve <- extract_nwqs_effect_curve(fit, grid = c(0, 0.25, 0.5, 0.75, 1), ref = 0.5)
  at_ref <- curve[curve$x == 0.5, ]
  expect_equal(at_ref$estimate, rep(0, nrow(at_ref)), tolerance = 1e-10)
})

test_that("extract_nwqs_effect_curve.nwqs has NA CI when rh = 1", {
  set.seed(7)
  n <- 80
  mix <- data.frame(C1 = rnorm(n), C2 = rnorm(n))
  y   <- mix$C1 + rnorm(n, sd = 0.5)
  fit <- nwqs(data = cbind(mix, y = y), mix_name = c("C1", "C2"),
              outcome = "y", family = "gaussian",
              transform_type = "percentile_rank", q = 4,
              rh = 1, n_permutation = 5, seed = 1234, quiet = TRUE)
  curve <- extract_nwqs_effect_curve(fit, ref = 0.5)
  expect_true(all(is.na(curve$lower)))
  expect_true(all(is.na(curve$upper)))
})

test_that("extract_nwqs_effect_curve.nwqs_boot returns bootstrap CI", {
  fit <- make_boot_fit(n_boot = 8)
  curve <- extract_nwqs_effect_curve(fit, ref = 0.5)
  expect_equal(unique(curve$inference_type), "bootstrap")
  overall <- curve[curve$term == "Overall", ]
  expect_true(all(overall$lower <= overall$estimate))
  expect_true(all(overall$estimate <= overall$upper))
})

test_that("extract_nwqs_effect_curve include_components = FALSE returns only Overall", {
  fit <- make_pr_fit_curve(q = 4, rh = 5)
  curve <- extract_nwqs_effect_curve(fit, ref = 0.5, include_components = FALSE)
  expect_equal(unique(curve$term), "Overall")
})

test_that("extract_nwqs_effect_curve rejects grid outside [0,1] in percentile_rank", {
  fit <- make_pr_fit_curve(q = 4, rh = 5)
  expect_error(
    extract_nwqs_effect_curve(fit, grid = c(0, 0.5, 1.1)),
    "\\[0, 1\\]"
  )
})

test_that("plot_nwqs_effect_curve returns a ggplot with median-centered defaults", {
  fit <- make_pr_fit_curve(q = 4, rh = 5)
  p <- plot_nwqs_effect_curve(fit)
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$x, "Joint exposure percentile rank")
  expect_match(p$labels$y, "Partial effect change relative to P50")
})

test_that("plot_nwqs_effect_curve(ref=0.25) shows P25 in y-axis label", {
  fit <- make_pr_fit_curve(q = 4, rh = 5)
  p <- plot_nwqs_effect_curve(fit, ref = 0.25)
  expect_match(p$labels$y, "P25")
})

test_that("plot_nwqs_effect_curve(include_components = TRUE) keeps component layers", {
  fit <- make_pr_fit_curve(q = 4, rh = 5)
  p <- plot_nwqs_effect_curve(fit, include_components = TRUE)
  expect_true("term" %in% names(p$data))
  expect_true(length(unique(p$data$term)) > 1)
})

test_that("plot_nwqs_effect_curve dispatches on nwqs_boot", {
  fit <- make_boot_fit(n_boot = 8)
  p <- plot_nwqs_effect_curve(fit)
  expect_s3_class(p, "ggplot")
  expect_equal(unique(p$data$inference_type), "bootstrap")
})
