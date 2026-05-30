# Display-layer contracts that protect transform-aware grids and
# exponentiated display for count families.

make_display_fit <- function(family = "gaussian",
                             transform_type = "percentile_rank",
                             nwqs_coef = 1) {
  q <- 4
  df_spline <- 3
  basis <- build_spline_basis_knots(transform_type, q = q, df_spline = df_spline)
  mean_shapes <- stats::setNames(
    c(0.2, -0.1, 0.05),
    paste0("X1_B", seq_len(df_spline))
  )
  fit <- list(
    family = family,
    transform_type = transform_type,
    q = q,
    df_spline = df_spline,
    spline_knots = basis$knots,
    spline_boundary = basis$boundary,
    mean_shapes = mean_shapes,
    final_weights = c(X1 = 1),
    mean_coefs = c("(Intercept)" = -5, nwqs = nwqs_coef),
    rh = 1,
    rh_weights = matrix(c(1), nrow = 1, dimnames = list(NULL, "X1")),
    rh_shapes = matrix(mean_shapes, nrow = 1, dimnames = list(NULL, names(mean_shapes))),
    rh_coefs = matrix(c(-5, nwqs_coef), nrow = 1,
                      dimnames = list(NULL, c("(Intercept)", "nwqs"))),
    fit = list(coefficients = data.frame(
      Estimate = c(-5, nwqs_coef),
      `Std. Error` = c(0.1, 0.1),
      `z value` = c(-50, 10),
      `Pr(>|z|)` = c(0, 0),
      row.names = c("(Intercept)", "nwqs"),
      check.names = FALSE
    )),
    call = quote(nwqs())
  )
  class(fit) <- c("nwqs", "list")
  fit
}

test_that("plot.nwqs uses a [0, 1] x grid for percentile_rank fits", {
  fit <- make_display_fit(transform_type = "percentile_rank")
  p <- plot(fit, type = "curves", y_scale = "partial")
  expect_equal(range(p$data$x), c(0, 1), tolerance = 1e-12)
})

test_that("plot.nwqs exponentiates predicted curves for negbin fits", {
  fit <- make_display_fit(family = "negbin", nwqs_coef = 0)
  p <- plot(fit, type = "curves", y_scale = "predicted")
  expect_true(all(p$data$y > 0))
  expect_equal(unique(p$data$y), exp(-5), tolerance = 1e-12)
})

test_that("nwqs_contrast reports negbin contrasts on the rate-ratio scale", {
  fit <- make_display_fit(family = "negbin")
  fit$rh <- 2
  fit$rh_coefs <- rbind(fit$rh_coefs, fit$rh_coefs)
  fit$rh_weights <- rbind(fit$rh_weights, fit$rh_weights)
  fit$rh_shapes <- rbind(fit$rh_shapes, fit$rh_shapes)

  expect_output(
    nwqs_contrast(fit, q_target = 3, q_ref = 0),
    "Rate Ratio"
  )
})

test_that("print.nwqs (percentile_rank) uses P-labelled contrast columns by default", {
  fit <- make_display_fit(transform_type = "percentile_rank")
  fit$fit$deviance <- 1.0
  fit$fit$aic <- 10
  fit$b <- 30
  out <- capture.output(print(fit))
  combined <- paste(out, collapse = "\n")
  expect_true(grepl("P25 vs P50", combined))
  expect_true(grepl("P75 vs P50", combined))
  expect_true(grepl("P95 vs P50", combined))
  expect_false(grepl("Q[0-9]+ vs Q[0-9]+", combined))
})

test_that("print.nwqs (q_bin) keeps Q-labelled columns (backward compat)", {
  fit <- make_display_fit(transform_type = "q_bin")
  fit$fit$deviance <- 1.0
  fit$fit$aic <- 10
  fit$b <- 30
  out <- capture.output(print(fit))
  combined <- paste(out, collapse = "\n")
  expect_true(grepl("Q2 vs Q1", combined))
  expect_true(grepl("Q3 vs Q1", combined))
  expect_true(grepl("Q4 vs Q1", combined))
  expect_false(grepl("P[0-9]+ vs P[0-9]+", combined))
})

test_that("print.nwqs honors user contrast_points/ref overrides", {
  fit <- make_display_fit(transform_type = "percentile_rank")
  fit$fit$deviance <- 1.0
  fit$fit$aic <- 10
  fit$b <- 30
  out <- capture.output(print(fit, contrast_points = c(0.1, 0.9), ref = 0.5))
  combined <- paste(out, collapse = "\n")
  expect_true(grepl("P10 vs P50", combined))
  expect_true(grepl("P90 vs P50", combined))
})
