# Tests for the new transform layer introduced in v0.2.0:
#   - trans_quantile() with type/ties signature
#   - apply_percentile_rank() helper for predict
#   - build_spline_basis_knots() helper for global knot alignment
#
# These tests pin the mathematical contract:
#   percentile_rank: u_i = rank(x_i, ties = "average") / n  in (0, 1]
#   q_bin (legacy): same as 0.1.x discrete quartile binning
#
# Each test must hold for both the standalone helpers and the embedded
# call from nwqs() so refactors do not silently drift.

# ----- trans_quantile: signature and math --------------------------------

test_that("trans_quantile() default is continuous percentile rank", {
  d <- data.frame(a = c(10, 30, 20, 40, 25))
  out <- trans_quantile(d)
  expect_equal(out$a, rank(d$a, ties.method = "average") / nrow(d))
  expect_true(all(out$a > 0 & out$a <= 1))
})

test_that("trans_quantile() percentile_rank respects ties = 'average'", {
  d <- data.frame(a = c(10, 30, 10, 40, 25))
  out <- trans_quantile(d, type = "percentile_rank", ties = "average")
  expected <- rank(d$a, ties.method = "average") / nrow(d)
  expect_equal(out$a, expected)
})

test_that("trans_quantile() ties = 'min' and 'max' shift tied ranks", {
  d <- data.frame(a = c(10, 10, 20, 30))
  out_min <- trans_quantile(d, type = "percentile_rank", ties = "min")
  out_max <- trans_quantile(d, type = "percentile_rank", ties = "max")
  expect_equal(out_min$a, c(1, 1, 3, 4) / 4)
  expect_equal(out_max$a, c(2, 2, 3, 4) / 4)
})

test_that("trans_quantile() q_bin reproduces 0.1.x discrete quartile output", {
  set.seed(1)
  d <- data.frame(a = rnorm(40))
  out <- trans_quantile(d, type = "q_bin", q = 4)
  expect_true(all(out$a %in% 0:3))
  expect_equal(length(unique(out$a)), 4)
})

test_that("trans_quantile() q_bin without q errors clearly", {
  expect_error(
    trans_quantile(data.frame(a = 1:10), type = "q_bin"),
    regexp = "q"
  )
})

test_that("trans_quantile() rejects the legacy method = argument", {
  expect_error(
    trans_quantile(data.frame(a = 1:5), method = "percentile"),
    regexp = "unused argument"
  )
})

test_that("trans_quantile() refuses an all-NA column", {
  expect_error(
    trans_quantile(data.frame(a = rep(NA_real_, 5))),
    regexp = "NA"
  )
})

# ----- apply_percentile_rank ---------------------------------------------

test_that("apply_percentile_rank() maps newdata into the training ECDF in (0, 1]", {
  train_x <- c(10, 20, 30, 40, 50)
  new_x   <- c(15, 35, 50)
  out <- apply_percentile_rank(new_x, train_x)
  expect_true(all(out >= 0 & out <= 1))
  # 15 is between 10 and 20: training values <= 15 are {10} → 1/5 = 0.2
  # 35 is between 30 and 40: training values <= 35 are {10,20,30} → 3/5 = 0.6
  # 50 equals max training: 5/5 = 1.0
  expect_equal(out, c(0.2, 0.6, 1.0))
})

test_that("apply_percentile_rank() clips newdata outside training range to (0, 1)", {
  train_x <- c(10, 20, 30)
  expect_equal(apply_percentile_rank(c(0, 100), train_x), c(0, 1))
})

# ----- build_spline_basis_knots ------------------------------------------

test_that("build_spline_basis_knots() percentile_rank uses [0, 1] grid", {
  res <- build_spline_basis_knots(
    transform_type = "percentile_rank", q = 4, df_spline = 3
  )
  expect_named(res, c("knots", "boundary"))
  expect_length(res$boundary, 2)
  expect_equal(res$boundary[1], 0)
  expect_equal(res$boundary[2], 1)
})

test_that("build_spline_basis_knots() q_bin uses 0:(q-1) grid (0.1.x compatibility)", {
  res <- build_spline_basis_knots(
    transform_type = "q_bin", q = 4, df_spline = 3
  )
  expect_named(res, c("knots", "boundary"))
  expect_equal(res$boundary[1], 0)
  expect_equal(res$boundary[2], 3)
})

test_that("build_spline_basis_knots() custom_knots / custom_boundary override defaults", {
  res <- build_spline_basis_knots(
    transform_type = "percentile_rank", q = 4, df_spline = 3,
    custom_knots = 0.5, custom_boundary = c(-1, 2)
  )
  expect_equal(as.numeric(res$knots), 0.5)
  expect_equal(res$boundary, c(-1, 2))
})

# ----- nwqs() new signature ----------------------------------------------

test_that("nwqs() default transform is percentile_rank and is recorded on the fit object", {
  set.seed(11)
  n <- 200
  d <- data.frame(
    y  = rnorm(n),
    X1 = rnorm(n),
    X2 = rnorm(n),
    X3 = rnorm(n)
  )
  fit <- nwqs(
    data    = d,
    mix_name = c("X1", "X2", "X3"),
    outcome = "y",
    family  = "gaussian",
    rh      = 3,
    n_permutation = 3,
    seed    = 1,
    quiet   = TRUE
  )
  expect_equal(fit$transform_type, "percentile_rank")
  expect_true("train_components_sorted" %in% names(fit))
  expect_named(fit$train_components_sorted, c("X1", "X2", "X3"))
  expect_true(all(vapply(fit$train_components_sorted, is.numeric, logical(1))))
})

test_that("nwqs() with transform_type='q_bin' records q on the fit object", {
  set.seed(11)
  n <- 200
  d <- data.frame(
    y  = rnorm(n),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )
  fit <- nwqs(
    data           = d,
    mix_name       = c("X1", "X2"),
    outcome        = "y",
    family         = "gaussian",
    rh             = 3,
    n_permutation  = 3,
    seed           = 1,
    quiet          = TRUE,
    transform_type = "q_bin",
    q              = 4
  )
  expect_equal(fit$transform_type, "q_bin")
  expect_equal(fit$q, 4)
})
