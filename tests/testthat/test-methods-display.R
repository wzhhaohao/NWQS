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

make_boot_fit_display <- function(family = "gaussian", n_boot = 8) {
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
    transform_type = "percentile_rank", q = 4,
    n_boot = n_boot, rh_inner = 1, n_permutation = 5,
    seed = 1234, quiet = TRUE
  ))
}

test_that("print.nwqs_boot percentile_rank fit prints P-label columns, no Q", {
  fit <- make_boot_fit_display(n_boot = 8)
  out <- capture.output(print(fit))
  combined <- paste(out, collapse = "\n")
  expect_true(grepl("P[0-9]+_vs_P[0-9]+", combined))
  expect_false(grepl("Q[0-9]+_vs_Q[0-9]+", combined))
})

test_that("summary.nwqs_boot stability table reads largest target dynamically", {
  fit <- make_boot_fit_display(n_boot = 8)
  out <- capture.output(summary(fit))
  combined <- paste(out, collapse = "\n")
  expect_false(grepl("Q[0-9]+_Effect_SD", combined))
  # Default percentile_rank contrast is {P25,P75,P95} vs P50; largest is P95_vs_P50
  expect_true(grepl("P95_vs_P50_Effect_SD", combined))
  expect_match(combined,
               "Note: P95_vs_P50_Effect_SD = SD of P95_vs_P50 effect")
})

test_that("plot_nwqs_contrast_box uses P-labelled facets in percentile_rank", {
  fit <- make_boot_fit_display(n_boot = 8)
  p <- plot_nwqs_contrast_box(fit)
  jitter_data <- p$layers[[1]]$data
  expect_true(any(grepl("^P[0-9]+$", levels(jitter_data$Quantile))))
  expect_false(any(grepl("^Q[0-9]+$", levels(jitter_data$Quantile))))
  expect_match(p$labels$x, "Percentile|percentile")
})

test_that("plot.nwqs percentile_rank x-axis label is percentile rank", {
  fit <- make_display_fit(transform_type = "percentile_rank")
  p <- plot(fit, type = "curves")
  expect_match(p$labels$x, "Percentile rank|percentile rank|Percentile Rank")
})

# ----- S2: rh>1 must not present algorithmic variance as Wald inference ----

.s2_fit <- function(rh) {
  set.seed(7)
  n  <- 150
  df <- data.frame(y = rnorm(n), M1 = rnorm(n), M2 = rnorm(n), M3 = rnorm(n))
  df$y <- 0.6 * rank(df$M1) / n + rnorm(n, sd = 0.4)
  suppressWarnings(nwqs(
    df, mix_name = c("M1", "M2", "M3"), outcome = "y", family = "gaussian",
    rh = rh, n_permutation = 2, n_shuffle = 2,
    seed = 1, plan_strategy = "sequential", quiet = TRUE
  ))
}

test_that("rh>1 print/summary neutralize Wald significance (S2)", {
  fit  <- .s2_fit(4)
  outp <- paste(capture.output(print(fit)), collapse = "\n")
  outs <- paste(capture.output(summary(fit)), collapse = "\n")
  expect_false(grepl("Signif\\. codes", outs))
  expect_false(grepl("Overall Significance", outp))
  expect_true(grepl("nwqs_boot", outs))
  expect_true(grepl("RH_SD", outs))
})

test_that("rh==1 keeps the standard GLM significance table (S2)", {
  outs1 <- paste(capture.output(summary(.s2_fit(1))), collapse = "\n")
  expect_true(grepl("Signif\\. codes", outs1))
})
