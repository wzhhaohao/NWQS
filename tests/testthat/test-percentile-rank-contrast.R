test_that(".contrast_point_label formats percentile_rank points as P{round*100}", {
  expect_equal(NWQS:::.contrast_point_label(0,   "percentile_rank", "P"), "P0")
  expect_equal(NWQS:::.contrast_point_label(0.5, "percentile_rank", "P"), "P50")
  expect_equal(NWQS:::.contrast_point_label(1/3, "percentile_rank", "P"), "P33")
  expect_equal(NWQS:::.contrast_point_label(2/3, "percentile_rank", "P"), "P67")
  expect_equal(NWQS:::.contrast_point_label(1,   "percentile_rank", "P"), "P100")
})

test_that(".contrast_point_label preserves Q labels for q_bin", {
  expect_equal(NWQS:::.contrast_point_label(0, "q_bin", "Q"), "Q1")
  expect_equal(NWQS:::.contrast_point_label(3, "q_bin", "Q"), "Q4")
})

test_that(".contrast_point_label numeric style returns trimmed numeric", {
  expect_equal(NWQS:::.contrast_point_label(0.25, "percentile_rank", "numeric"), "0.25")
  expect_equal(NWQS:::.contrast_point_label(3,    "q_bin",            "numeric"), "3")
})

test_that(".contrast_pair_label joins with _vs_", {
  expect_equal(
    NWQS:::.contrast_pair_label(0.75, 0.25, "percentile_rank", "P"),
    "P75_vs_P25"
  )
  expect_equal(
    NWQS:::.contrast_pair_label(3, 0, "q_bin", "Q"),
    "Q4_vs_Q1"
  )
})

test_that(".label_style_default picks P for percentile_rank, Q for q_bin", {
  expect_equal(NWQS:::.label_style_default("percentile_rank"), "P")
  expect_equal(NWQS:::.label_style_default("q_bin"),           "Q")
})

test_that(".validate_pr_points rejects out-of-[0,1] in percentile_rank", {
  expect_error(NWQS:::.validate_pr_points(c(0.5, 1.2), "percentile_rank"), "\\[0, 1\\]")
  expect_silent(NWQS:::.validate_pr_points(c(0, 0.5, 1), "percentile_rank"))
  expect_silent(NWQS:::.validate_pr_points(c(0, 3),     "q_bin"))
})

make_pr_fit <- function(family = "gaussian", q = 4, rh = 5) {
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

make_qbin_fit <- function(q = 4, rh = 5) {
  set.seed(99)
  n <- 80
  mix <- data.frame(C1 = rnorm(n), C2 = rnorm(n))
  y   <- mix$C1 + rnorm(n, sd = 0.5)
  dat <- cbind(mix, y = y)
  nwqs(
    data = dat, mix_name = c("C1", "C2"), outcome = "y",
    family = "gaussian",
    transform_type = "q_bin", q = q, rh = rh,
    n_permutation = 5, seed = 1234, quiet = TRUE
  )
}

test_that("extract_nwqs_effects (percentile_rank) uses P labels by default", {
  fit <- make_pr_fit(q = 4, rh = 5)
  eff <- extract_nwqs_effects(fit)
  expect_true(all(grepl("^P[0-9]+_vs_P[0-9]+$", unique(eff$Target))))
  expect_false(any(grepl("Q", eff$Target)))
  expect_equal(sort(unique(eff$Target)),
               sort(c("P33_vs_P0", "P67_vs_P0", "P100_vs_P0")))
})

test_that("extract_nwqs_effects respects user contrast_points/ref", {
  fit <- make_pr_fit(q = 4, rh = 5)
  eff <- extract_nwqs_effects(fit,
                              contrast_points = c(0.25, 0.75),
                              ref = 0.5)
  expect_equal(sort(unique(eff$Target)), c("P25_vs_P50", "P75_vs_P50"))
})

test_that("extract_nwqs_effects rejects points outside [0,1] for percentile_rank", {
  fit <- make_pr_fit(q = 4, rh = 5)
  expect_error(
    extract_nwqs_effects(fit, contrast_points = c(0.5, 1.2)),
    "\\[0, 1\\]"
  )
})

test_that("extract_nwqs_effects q_bin keeps Q labels (backward compat)", {
  fit <- make_qbin_fit(q = 4, rh = 5)
  eff <- extract_nwqs_effects(fit)
  expect_equal(sort(unique(eff$Target)),
               c("Q2_vs_Q1", "Q3_vs_Q1", "Q4_vs_Q1"))
})

test_that("nwqs_contrast prints 'Percentile-rank Contrast' in percentile_rank", {
  fit <- make_pr_fit(q = 4, rh = 5)
  out <- capture.output(nwqs_contrast(fit, target = 0.75, ref = 0.5))
  combined <- paste(out, collapse = "\n")
  expect_true(grepl("Percentile-rank Contrast", combined))
  expect_true(grepl("Target P75", combined))
  expect_true(grepl("Ref P50",    combined))
  expect_false(grepl("Q[0-9]+ vs Q[0-9]+", combined))
})

test_that("nwqs_contrast default in percentile_rank = target 0.75 vs ref 0.5", {
  fit <- make_pr_fit(q = 4, rh = 5)
  out <- capture.output(nwqs_contrast(fit))
  combined <- paste(out, collapse = "\n")
  expect_true(grepl("Target P75",  combined))
  expect_true(grepl("Ref P50",     combined))
})

test_that("nwqs_contrast q_bin keeps the legacy Q-label printout", {
  fit <- make_qbin_fit(q = 4, rh = 5)
  out <- capture.output(nwqs_contrast(fit))
  combined <- paste(out, collapse = "\n")
  expect_true(grepl("Quantile Contrast", combined))
  expect_true(grepl("Target Q4",         combined))
  expect_true(grepl("Ref Q1",            combined))
})

test_that("nwqs_contrast warns when both target and q_target are supplied", {
  fit <- make_pr_fit(q = 4, rh = 5)
  tmp <- tempfile()
  conn <- file(tmp, open = "wt")
  sink(conn)
  on.exit({ sink(); close(conn); unlink(tmp) }, add = TRUE)
  expect_warning(
    nwqs_contrast(fit, target = 0.75, q_target = 3),
    "Both"
  )
})

test_that("nwqs_contrast rejects target outside [0,1] in percentile_rank", {
  fit <- make_pr_fit(q = 4, rh = 5)
  expect_error(nwqs_contrast(fit, target = 1.5), "\\[0, 1\\]")
})

make_boot_fit_pr <- function(n_boot = 8) {
  set.seed(2026)
  n <- 80
  mix <- data.frame(
    Component1 = rnorm(n),
    Component2 = rnorm(n),
    Component3 = rnorm(n)
  )
  beta <- c(0.6, 0.3, 0.1)
  eta  <- as.matrix(mix) %*% beta + rnorm(n, sd = 0.5)
  dat <- cbind(mix, y = as.numeric(eta))
  expect_small_boot_warning(nwqs_boot(
    data = dat, mix_name = paste0("Component", 1:3),
    outcome = "y", family = "gaussian",
    transform_type = "percentile_rank", q = 4,
    n_boot = n_boot, rh_inner = 1, n_permutation = 5,
    seed = 1234, quiet = TRUE
  ))
}

test_that("extract_nwqs_effects on nwqs_boot returns bootstrap CI columns", {
  fit <- make_boot_fit_pr(n_boot = 8)
  eff <- extract_nwqs_effects(fit, contrast_points = c(0.25, 0.75), ref = 0.5)
  expect_true(all(c("Target", "Term", "Estimate",
                    "Boot_CI_Lower", "Boot_CI_Upper") %in% names(eff)))
  expect_equal(sort(unique(eff$Target)), c("P25_vs_P50", "P75_vs_P50"))
  expect_true(all(eff$Boot_CI_Lower <= eff$Estimate))
  expect_true(all(eff$Estimate      <= eff$Boot_CI_Upper))
})
