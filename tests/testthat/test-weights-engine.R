# Unit tests for the permutation_scorer engine and its internal loss
# function .calc_loss(). These tests:
#   1. pin the outer n_permutation (nwqs/nwqs_boot) and inner n_shuffle
#      (permutation_scorer) defaults to 30, and that n_shuffle drives the loop
#   2. lock the Poisson mu = 0 + binomial p = 0/1 boundary defenses in
#      .calc_loss so future refactors cannot silently remove them
#   3. assert in-bag rank deficiency now produces a warning and returns
#      NULL (0.2.0 change from the legacy silent NA -> 0 fallback)

# ----- Default counts: outer n_permutation = 30, inner n_shuffle = 30 -----

test_that("default outer n_permutation = 30 and inner n_shuffle = 30", {
  expect_equal(eval(formals(nwqs)$n_permutation), 30)
  expect_equal(eval(formals(nwqs_boot)$n_permutation), 30)
  expect_equal(eval(formals(permutation_scorer)$n_shuffle), 30)
  # nwqs()/nwqs_boot() $n_shuffle defaults are asserted in test-nwqs-internals.R
  # (Task 3, where those formals are added).
})

test_that("permutation_scorer: n_shuffle controls the inner shuffle count", {
  set.seed(123)
  n  <- 80
  xm <- cbind(
    "(Intercept)" = 1,
    A_B1 = rnorm(n), A_B2 = rnorm(n), A_B3 = rnorm(n),
    B_B1 = rnorm(n), B_B2 = rnorm(n), B_B3 = rnorm(n)
  )
  yv <- 0.7 * xm[, "A_B1"] + 0.3 * xm[, "B_B2"] + rnorm(n, sd = 0.5)
  sv <- c("A_B1", "A_B2", "A_B3", "B_B1", "B_B2", "B_B3")

  set.seed(7)
  r1 <- permutation_scorer(xm, yv, mix_name = c("A", "B"), spline_vars = sv,
                           family = gaussian(), n_shuffle = 1)
  set.seed(7)
  r2 <- permutation_scorer(xm, yv, mix_name = c("A", "B"), spline_vars = sv,
                           family = gaussian(), n_shuffle = 60)

  # Same in-bag draw (same seed) but a different number of shuffles must
  # produce a different averaged permutation-importance estimate.
  expect_false(isTRUE(all.equal(unname(r1$weights), unname(r2$weights))))
})

# ----- .calc_loss boundary contracts --------------------------------------

test_that(".calc_loss() stays finite when binomial mu_pred is at 0 or 1", {
  res <- NWQS:::.calc_loss(
    y_true   = c(0, 1, 0, 1),
    mu_pred  = c(0, 1, 0, 1),
    fam_name = "binomial"
  )
  expect_true(is.finite(res))
})

test_that(".calc_loss() stays finite when Poisson mu_pred is at 0 with mixed-y data", {
  res <- NWQS:::.calc_loss(
    y_true   = c(0, 1, 5, 10),
    mu_pred  = c(0, 0, 0, 0),
    fam_name = "poisson"
  )
  expect_true(is.finite(res))
})

test_that(".calc_loss() returns 0 for a perfect Gaussian prediction", {
  res <- NWQS:::.calc_loss(
    y_true   = c(1, 2, 3, 4),
    mu_pred  = c(1, 2, 3, 4),
    fam_name = "gaussian"
  )
  expect_equal(res, 0)
})

test_that(".calc_loss() quasipoisson matches poisson on the same input", {
  yp <- c(0, 1, 3, 7)
  mu <- c(0.5, 1.2, 3.1, 5.5)
  expect_equal(
    NWQS:::.calc_loss(yp, mu, "poisson"),
    NWQS:::.calc_loss(yp, mu, "quasipoisson")
  )
})

test_that("zero_weight_action = 'uniform' turns zero importance into equal weights", {
  importance <- c(X1 = 0, X2 = 0, X3 = 0)
  weights <- NWQS:::.importance_to_weights(
    importance_scores = importance,
    zero_weight_action = "uniform"
  )
  expect_equal(weights, c(X1 = 1 / 3, X2 = 1 / 3, X3 = 1 / 3))
})

test_that("zero_weight_action = 'na' keeps zero-importance weights as NA", {
  importance <- c(X1 = 0, X2 = 0)
  weights <- NWQS:::.importance_to_weights(
    importance_scores = importance,
    zero_weight_action = "na"
  )
  expect_true(all(is.na(weights)))
  expect_named(weights, names(importance))
})

# ----- In-bag rank deficiency -> warning + NULL --------------------------

test_that("permutation_scorer warns and returns NULL on rank-deficient in-bag fit", {
  set.seed(1)
  n <- 80
  # Two perfectly collinear spline columns force glm.fit to drop one and
  # set its coefficient to NA. Under 0.2.0 this should be a warning and a
  # NULL return rather than the silent NA -> 0 of 0.1.x.
  base <- rnorm(n)
  x <- cbind(
    "(Intercept)" = rep(1, n),
    "X1_B1"       = base,
    "X1_B2"       = base
  )
  y <- rnorm(n)

  expect_warning(
    res <- permutation_scorer(
      x            = x,
      y            = y,
      mix_name     = "X1",
      spline_vars  = c("X1_B1", "X1_B2"),
      family       = gaussian(),
      n_shuffle    = 2
    ),
    regexp = "rank deficiency"
  )
  expect_null(res)
})
