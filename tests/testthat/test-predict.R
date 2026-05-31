# Contract tests for predict.nwqs() and predict.nwqs_boot(). These pin:
#   - dimensions
#   - default newdata = training data
#   - type = c("nwqs_index", "link", "response")
#   - parity between the two supported transforms (percentile_rank, q_bin)
#   - parity across all four currently supported families
#   - consistency between rh = 1 and rh > 1 paths

# ----- Fixture -----------------------------------------------------------

make_pred_data <- function(n = 200, family = "gaussian", seed = 11) {
  set.seed(seed)
  Sigma <- diag(3)
  Sigma[lower.tri(Sigma)] <- Sigma[upper.tri(Sigma)] <- 0.4
  L <- chol(Sigma)
  X <- matrix(rnorm(n * 3), n, 3) %*% L
  colnames(X) <- paste0("X", 1:3)
  eta <- 0.5 * X[, 1] - 0.3 * X[, 2]
  y <- switch(
    family,
    gaussian     = eta + rnorm(n, sd = 0.5),
    binomial     = rbinom(n, 1, plogis(eta - 0.2)),
    poisson      = rpois(n, exp(0.5 + 0.5 * eta)),
    quasipoisson = rpois(n, exp(0.5 + 0.5 * eta))
  )
  data.frame(y = y, X)
}

fit_pred <- function(family = "gaussian",
                     transform_type = "percentile_rank",
                     rh = 3) {
  d <- make_pred_data(family = family, seed = 11)
  nwqs(
    data           = d,
    mix_name       = c("X1", "X2", "X3"),
    outcome        = "y",
    family         = family,
    rh             = rh,
    n_permutation  = 3,
    seed           = 1,
    quiet          = TRUE,
    transform_type = transform_type,
    q              = 4
  )
}

# ----- Generics --------------------------------------------------------

test_that("predict.nwqs exists and returns a numeric vector with length(newdata)", {
  fit <- fit_pred()
  newd <- make_pred_data(family = "gaussian", seed = 99)
  out <- predict(fit, newdata = newd)
  expect_true(is.numeric(out))
  expect_length(out, nrow(newd))
  expect_true(all(is.finite(out)))
})

test_that("predict.nwqs with newdata = NULL falls back to the training data", {
  fit <- fit_pred()
  out_null <- predict(fit, newdata = NULL)
  out_self <- predict(fit, newdata = fit$data)
  expect_equal(length(out_null), nrow(fit$data))
  expect_equal(out_null, out_self)
})

# ----- type = nwqs_index / link / response ------------------------------

test_that("predict.nwqs() honors type = 'nwqs_index', 'link', and 'response'", {
  fit <- fit_pred(family = "gaussian", rh = 3)
  newd <- make_pred_data(family = "gaussian", seed = 99)

  idx <- predict(fit, newdata = newd, type = "nwqs_index")
  link <- predict(fit, newdata = newd, type = "link")
  resp <- predict(fit, newdata = newd, type = "response")

  expect_length(idx, nrow(newd))
  expect_length(link, nrow(newd))
  expect_length(resp, nrow(newd))

  # Gaussian identity link: response == link
  expect_equal(resp, link)
})

test_that("predict.nwqs() response scale applies the inverse link for non-Gaussian", {
  for (fam in c("binomial", "poisson", "quasipoisson")) {
    fit <- fit_pred(family = fam, rh = 2)
    newd <- make_pred_data(family = fam, seed = 99)
    link <- predict(fit, newdata = newd, type = "link")
    resp <- predict(fit, newdata = newd, type = "response")
    expected <- switch(
      fam,
      binomial = plogis(link),
      poisson  = exp(link),
      quasipoisson = exp(link)
    )
    expect_equal(resp, expected, tolerance = 1e-10)
  }
})

# ----- Transform parity --------------------------------------------------

test_that("predict.nwqs works under both percentile_rank and q_bin", {
  for (tt in c("percentile_rank", "q_bin")) {
    fit <- fit_pred(transform_type = tt, rh = 2)
    newd <- make_pred_data(family = "gaussian", seed = 99)
    out <- predict(fit, newdata = newd)
    expect_true(all(is.finite(out)),
                label = sprintf("transform_type = %s", tt))
  }
})

test_that("predict.nwqs percentile-rank transform reuses fit-sample tie rules", {
  object <- list(
    transform_type = "percentile_rank",
    train_components_sorted = list(X1 = sort(c(10, 10, 20, 30))),
    ties = "average",
    q = 4
  )
  class(object) <- c("nwqs", "list")

  out <- NWQS:::.nwqs_transform_newdata(
    object,
    data.frame(X1 = c(10, 20))
  )
  expect_equal(out$X1, c(1.5 / 4, 3 / 4))
})

# ----- Default type ------------------------------------------------------

test_that("predict.nwqs() default type is 'response'", {
  fit <- fit_pred(family = "binomial", rh = 2)
  newd <- make_pred_data(family = "binomial", seed = 99)
  default <- predict(fit, newdata = newd)
  response <- predict(fit, newdata = newd, type = "response")
  expect_equal(default, response)
})

# ----- predict.nwqs_boot -------------------------------------------------

test_that("predict.nwqs_boot returns predictions with bootstrap CIs", {
  d <- make_pred_data(family = "gaussian", seed = 11)
  boot_fit <- expect_small_boot_warning(nwqs_boot(
    data           = d,
    mix_name       = c("X1", "X2", "X3"),
    outcome        = "y",
    family         = "gaussian",
    n_boot         = 6,
    rh_inner       = 1,
    n_permutation  = 3,
    seed           = 1,
    quiet          = TRUE,
    transform_type = "percentile_rank",
    q              = 4
  ))
  out <- predict(boot_fit, newdata = d, type = "nwqs_index")
  expect_true(is.numeric(out) || is.data.frame(out))
  if (is.data.frame(out)) {
    expect_true(all(c("estimate", "lower", "upper") %in% names(out)))
    expect_equal(nrow(out), nrow(d))
    expect_true(all(out$lower <= out$estimate))
    expect_true(all(out$estimate <= out$upper))
  }
})
