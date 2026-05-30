# Contract tests for the v0.2.0 default-management layer:
#   - NWQS_DEFAULTS centralizes every user-facing default value in one
#     list at R/zzz-defaults.R, so changing a default only requires
#     touching one place.
#   - nwqs_control() validates and packages "advanced" parameters
#     (custom_knots, custom_boundary, zero_weight_action) that live
#     outside the main nwqs() signature so the signature does not grow
#     unbounded as new knobs get added.
#
# Backward-compat note: min_shape_sd and ties stay on the main nwqs()
# signature. They were promoted to user-facing defaults in Phase 3.3 and
# Phase 2 respectively; moving them into control now would be a needless
# breaking change. The control layer is only for "soft" parameters that
# do not yet exist on the signature.

# ----- NWQS_DEFAULTS shape ----------------------------------------------

test_that("NWQS_DEFAULTS exists and contains every documented default", {
  expect_true(is.list(NWQS_DEFAULTS))
  expected_keys <- c(
    "q", "df_spline", "transform_type", "ties",
    "train_prop", "rh", "n_permutation",
    "n_boot", "rh_inner", "conf_level",
    "seed", "min_shape_sd", "zero_weight_action"
  )
  expect_true(all(expected_keys %in% names(NWQS_DEFAULTS)))
})

test_that("NWQS_DEFAULTS values match the published 0.2.0 defaults", {
  expect_equal(NWQS_DEFAULTS$q, 4)
  expect_equal(NWQS_DEFAULTS$df_spline, 3)
  expect_equal(NWQS_DEFAULTS$transform_type, "percentile_rank")
  expect_equal(NWQS_DEFAULTS$ties, "average")
  expect_equal(NWQS_DEFAULTS$train_prop, 0.6)
  expect_equal(NWQS_DEFAULTS$rh, 10)
  expect_equal(NWQS_DEFAULTS$n_permutation, 30)
  expect_equal(NWQS_DEFAULTS$n_boot, 100)
  expect_equal(NWQS_DEFAULTS$rh_inner, 1)
  expect_equal(NWQS_DEFAULTS$conf_level, 0.95)
  expect_equal(NWQS_DEFAULTS$seed, 1234)
  expect_equal(NWQS_DEFAULTS$min_shape_sd, 1e-8)
  expect_equal(NWQS_DEFAULTS$zero_weight_action, "na")
})

test_that("nwqs() reads its defaults from NWQS_DEFAULTS", {
  # Each formal default should evaluate to the corresponding NWQS_DEFAULTS
  # entry. Using eval() because the formal is a language object referring
  # to NWQS_DEFAULTS$key, not the literal value.
  expect_equal(eval(formals(nwqs)$q), NWQS_DEFAULTS$q)
  expect_equal(eval(formals(nwqs)$df_spline), NWQS_DEFAULTS$df_spline)
  expect_equal(eval(formals(nwqs)$train_prop), NWQS_DEFAULTS$train_prop)
  expect_equal(eval(formals(nwqs)$rh), NWQS_DEFAULTS$rh)
  expect_equal(eval(formals(nwqs)$n_permutation), NWQS_DEFAULTS$n_permutation)
  expect_equal(eval(formals(nwqs)$seed), NWQS_DEFAULTS$seed)
  expect_equal(eval(formals(nwqs)$min_shape_sd), NWQS_DEFAULTS$min_shape_sd)
})

test_that("nwqs_boot() reads its defaults from NWQS_DEFAULTS", {
  expect_equal(eval(formals(nwqs_boot)$n_boot), NWQS_DEFAULTS$n_boot)
  expect_equal(eval(formals(nwqs_boot)$rh_inner), NWQS_DEFAULTS$rh_inner)
  expect_equal(eval(formals(nwqs_boot)$n_permutation), NWQS_DEFAULTS$n_permutation)
  expect_equal(eval(formals(nwqs_boot)$conf_level), NWQS_DEFAULTS$conf_level)
})

# ----- nwqs_control() contract -------------------------------------------

test_that("nwqs_control() returns a list-classed nwqs_control object", {
  ctrl <- nwqs_control()
  expect_s3_class(ctrl, "nwqs_control")
  expect_true(is.list(ctrl))
})

test_that("nwqs_control() defaults: custom_knots NULL, custom_boundary NULL, zero_weight_action 'na'", {
  ctrl <- nwqs_control()
  expect_null(ctrl$custom_knots)
  expect_null(ctrl$custom_boundary)
  expect_equal(ctrl$zero_weight_action, "na")
})

test_that("nwqs_control() accepts custom_knots and custom_boundary", {
  ctrl <- nwqs_control(custom_knots = c(0.25, 0.75),
                       custom_boundary = c(0, 1))
  expect_equal(ctrl$custom_knots, c(0.25, 0.75))
  expect_equal(ctrl$custom_boundary, c(0, 1))
})

test_that("nwqs_control() rejects a custom_boundary that is not length 2", {
  expect_error(
    nwqs_control(custom_boundary = c(0)),
    regexp = "length"
  )
  expect_error(
    nwqs_control(custom_boundary = c(0, 0.5, 1)),
    regexp = "length"
  )
})

test_that("nwqs_control() rejects a non-numeric custom_knots vector", {
  expect_error(
    nwqs_control(custom_knots = c("a", "b")),
    regexp = "numeric"
  )
})

test_that("nwqs_control() validates zero_weight_action against the allowed set", {
  expect_error(
    nwqs_control(zero_weight_action = "invalid"),
    regexp = "zero_weight_action"
  )
})

# ----- nwqs() forwards custom_knots / custom_boundary to the spline basis

test_that("nwqs() honors control$custom_knots when building the spline basis", {
  set.seed(11)
  n <- 200
  d <- data.frame(
    y  = rnorm(n),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )
  # 2 internal knots match df_spline = 3 (so the natural-spline basis has
  # df = #knots + 1 = 3 columns, matching the rest of the pipeline).
  fit <- nwqs(
    data = d, mix_name = c("X1", "X2"), outcome = "y",
    family = "gaussian", rh = 2, n_permutation = 3, seed = 1,
    quiet = TRUE, transform_type = "percentile_rank", q = 4,
    control = nwqs_control(custom_knots = c(0.3, 0.7),
                           custom_boundary = c(0, 1))
  )
  expect_equal(as.numeric(fit$spline_knots), c(0.3, 0.7))
  expect_equal(fit$spline_boundary, c(0, 1))
})
