# Simulation generator contracts for the 0.2.0 percentile-rank DGP.
# These tests keep the data-generating mechanism aligned with the default
# nwqs() transform and preserve q-bin as an explicit legacy path.

make_sim_inputs <- function(p = 3) {
  list(
    n_obs = 80,
    mu_preds = rep(0, p),
    sigma_preds = diag(p),
    beta_wqs = 0.4,
    beta_preds = seq(0.5, 0.2, length.out = p),
    seed = 2026
  )
}

call_generator <- function(kind, ...) {
  args <- c(make_sim_inputs(), list(...))
  switch(
    kind,
    gaussian = do.call(gen_nonlinear_data, c(args, list(snr_db = Inf))),
    binomial = do.call(gen_nonlinear_bio_data, c(args, list(
      intercept = -0.5, target_prop = 0.4, snr_db = Inf
    ))),
    poisson = do.call(gen_nonlinear_count_data, c(args, list(
      intercept = -1.2, snr_db = Inf
    )))
  )
}

test_that("nonlinear generators default to percentile-rank DGP with full effect curves", {
  for (kind in c("gaussian", "binomial", "poisson")) {
    dat <- call_generator(kind)
    curve <- attr(dat, "true_effect_curve")

    expect_equal(attr(dat, "transform_type"), "percentile_rank", info = kind)
    expect_equal(attr(dat, "eval_points"), seq(0, 1, length.out = 4), info = kind)
    expect_equal(attr(dat, "effect_grid"), seq(0, 1, by = 0.01), info = kind)
    expect_equal(dim(curve), c(4, 101), info = kind)
    expect_equal(colnames(curve)[1], "P000", info = kind)
    expect_equal(colnames(curve)[101], "P100", info = kind)
    expect_equal(unname(curve[, "P000"]), rep(0, 4), tolerance = 1e-12, info = kind)
    expect_equal(
      unname(curve["Overall", ]),
      unname(colSums(curve[paste0("Component", 1:3), , drop = FALSE])),
      tolerance = 1e-12,
      info = kind
    )
    expect_equal(
      colnames(attr(dat, "true_effect_mat")),
      c("Q2_vs_Q1", "Q3_vs_Q1", "Q4_vs_Q1"),
      info = kind
    )
  }
})

test_that("q_bin generator path preserves the legacy evaluation grid", {
  dat <- call_generator("gaussian", transform_type = "q_bin", q = 4)

  expect_equal(attr(dat, "transform_type"), "q_bin")
  expect_equal(attr(dat, "eval_points"), 0:3)
  expect_equal(attr(dat, "spline_boundary"), c(0, 3))
  expect_equal(
    colnames(attr(dat, "true_effect_mat")),
    c("Q2_vs_Q1", "Q3_vs_Q1", "Q4_vs_Q1")
  )
})

test_that("custom transform_fun is still used before nonlinear effects are generated", {
  called <- FALSE
  custom_transform <- function(x) {
    called <<- TRUE
    trans_quantile(x, type = "q_bin", q = 4)
  }

  dat <- call_generator(
    "gaussian",
    transform_type = "percentile_rank",
    transform_fun = custom_transform
  )

  expect_true(called)
  expect_equal(attr(dat, "transform_type"), "percentile_rank")
})

test_that("generator eval_points and effect_grid validation is explicit", {
  expect_error(
    call_generator("gaussian", eval_points = c(0, 0.5)),
    regexp = "eval_points"
  )
  expect_error(
    call_generator("gaussian", eval_points = c(0, 0.5, 0.5, 1)),
    regexp = "eval_points"
  )
  expect_error(
    call_generator("gaussian", eval_points = c(0, 0.5, 0.75, 1.2)),
    regexp = "eval_points"
  )
  expect_error(
    call_generator("gaussian", effect_grid = c(0, 0.5, NA, 1)),
    regexp = "effect_grid"
  )
  expect_error(
    call_generator("gaussian", effect_grid = c(0, 0.5, 0.4, 1)),
    regexp = "effect_grid"
  )
  expect_error(
    call_generator("gaussian", effect_grid = c(-0.1, 0, 1)),
    regexp = "effect_grid"
  )
})

test_that("calc_true_importance() defaults to curve-range weights for true effect curves", {
  curve <- rbind(
    Overall = c(0, 1, 2),
    Component1 = c(0, 1, 0),
    Component2 = c(0, 0, 2),
    Component3 = c(0, 0.5, 0.5)
  )
  colnames(curve) <- c("P000", "P050", "P100")

  expected <- c(Component1 = 1, Component2 = 2, Component3 = 0.5)
  expected <- expected / sum(expected)

  expect_equal(
    calc_true_importance(curve, paste0("Component", 1:3), method = "curve_range"),
    expected
  )
  expect_equal(
    calc_true_importance(curve, paste0("Component", 1:3)),
    expected
  )
})

test_that("percentile-rank generator and fit recover curve-range weights in a moderate simulation", {
  skip_on_cran()
  set.seed(20260530)
  p <- 6
  mix <- paste0("Component", seq_len(p))
  sigma <- generate_sigma(p, mode = "mixed", seed = 20260530)

  dat <- gen_nonlinear_data(
    n_obs = 1000,
    mu_preds = rep(0, p),
    sigma_preds = sigma,
    beta_wqs = 0.8,
    beta_preds = c(0.30, 0.24, 0.18, 0.13, 0.10, 0.05),
    snr_db = Inf,
    seed = 20260530,
    shape = rep("linear_like", p)
  )

  fit <- nwqs(
    data = dat,
    mix_name = mix,
    covariates = c("x_cont", "x_bin", "x_cat"),
    outcome = "y",
    family = "gaussian",
    transform_type = "percentile_rank",
    q = 4,
    rh = 5,
    n_permutation = 15,
    train_prop = 0.6,
    seed = 20260530,
    quiet = TRUE
  )

  true_w <- calc_true_importance(attr(dat, "true_effect_curve"), mix)
  err <- calc_weight_error(fit$final_weights, true_w)
  expect_lt(err$SAE, 0.25)
})
