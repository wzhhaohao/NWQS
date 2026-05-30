# Unit tests for the permutation_scorer engine and its internal loss
# function .calc_loss(). These tests:
#   1. pin n_permutation default to 30 across nwqs/nwqs_boot/permutation_scorer
#   2. lock the Poisson mu = 0 + binomial p = 0/1 boundary defenses in
#      .calc_loss so future refactors cannot silently remove them
#   3. assert in-bag rank deficiency now produces a warning and returns
#      NULL (0.2.0 change from the legacy silent NA -> 0 fallback)

# ----- Default n_permutation = 30 ----------------------------------------

test_that("n_permutation default is 30 across nwqs(), nwqs_boot(), permutation_scorer()", {
  expect_equal(eval(formals(nwqs)$n_permutation), 30)
  expect_equal(eval(formals(nwqs_boot)$n_permutation), 30)
  expect_equal(formals(permutation_scorer)$n_permutation, 30)
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
      n_permutation = 2
    ),
    regexp = "rank deficiency"
  )
  expect_null(res)
})
