# Contract tests for the new family = "negbin" path that dispatches to
# MASS::glm.nb. The expectations mirror the Poisson path, plus a few
# negbin-specific facts (theta on the fit, exp-scale interpretation).

skip_if_not_installed("MASS")

make_nb_data <- function(n = 300, theta = 2, seed = 51) {
  set.seed(seed)
  X <- matrix(rnorm(n * 3), n, 3)
  colnames(X) <- paste0("X", 1:3)
  eta <- 0.4 + 0.5 * X[, 1] - 0.3 * X[, 2]
  mu <- exp(eta)
  y <- MASS::rnegbin(n, mu = mu, theta = theta)
  data.frame(y = y, X)
}

# ----- Fit succeeds ------------------------------------------------------

test_that("nwqs() accepts family = 'negbin' and produces a well-formed fit", {
  d <- make_nb_data()
  fit <- nwqs(
    data = d, mix_name = c("X1", "X2", "X3"), outcome = "y",
    family = "negbin", rh = 3, n_permutation = 3, seed = 1,
    quiet = TRUE, transform_type = "percentile_rank", q = 4
  )
  expect_s3_class(fit, "nwqs")
  expect_equal(fit$family, "negbin")
  expect_true(all(fit$final_weights >= 0))
  expect_equal(sum(fit$final_weights), 1, tolerance = 1e-10)
  expect_true("nwqs" %in% names(fit$mean_coefs))
})

# ----- rh = 1: inner model_obj is a MASS::glm.nb fit --------------------

test_that("nwqs() with rh = 1 and family = 'negbin' stores a MASS::glm.nb model on $model_obj", {
  d <- make_nb_data()
  fit <- nwqs(
    data = d, mix_name = c("X1", "X2", "X3"), outcome = "y",
    family = "negbin", rh = 1, n_permutation = 3, seed = 1,
    quiet = TRUE, transform_type = "q_bin", q = 4
  )
  expect_true(inherits(fit$model_obj, "negbin"))
  expect_true(is.numeric(fit$model_obj$theta))
})

# ----- predict / vcov / confint flow through unchanged -------------------

test_that("predict.nwqs(type = 'response') returns positive predictions for negbin", {
  d <- make_nb_data()
  fit <- nwqs(
    data = d, mix_name = c("X1", "X2", "X3"), outcome = "y",
    family = "negbin", rh = 3, n_permutation = 3, seed = 1,
    quiet = TRUE, transform_type = "percentile_rank", q = 4
  )
  pred <- predict(fit, newdata = d, type = "response")
  expect_true(all(pred > 0))
  expect_length(pred, nrow(d))
})

test_that("vcov.nwqs and confint.nwqs work on a negbin rh = 1 fit", {
  d <- make_nb_data()
  fit <- nwqs(
    data = d, mix_name = c("X1", "X2", "X3"), outcome = "y",
    family = "negbin", rh = 1, n_permutation = 3, seed = 1,
    quiet = TRUE, transform_type = "q_bin", q = 4
  )
  v <- vcov(fit)
  expect_true(is.matrix(v))
  expect_true("nwqs" %in% rownames(v))
  ci <- suppressMessages(confint(fit))
  expect_true("nwqs" %in% rownames(ci))
})

# ----- nwqs_boot path ----------------------------------------------------

test_that("nwqs_boot() with family = 'negbin' returns a valid CI table", {
  d <- make_nb_data()
  boot_fit <- expect_small_boot_warning(nwqs_boot(
    data = d, mix_name = c("X1", "X2", "X3"), outcome = "y",
    family = "negbin", n_boot = 8, rh_inner = 1, n_permutation = 3,
    seed = 1, quiet = TRUE, transform_type = "q_bin", q = 4
  ))
  expect_s3_class(boot_fit, "nwqs_boot")
  expect_equal(boot_fit$family, "negbin")
  expect_true(boot_fit$n_success > 0)
  expect_true(nrow(boot_fit$ci_table) > 0)
  expect_true(all(boot_fit$ci_table$Boot_CI_Lower <=
                  boot_fit$ci_table$Boot_Mean))
  expect_true(all(boot_fit$ci_table$Boot_Mean <=
                  boot_fit$ci_table$Boot_CI_Upper))
})

# ----- Generator: gen_nbin_data() ---------------------------------------

test_that("gen_nbin_data() returns counts with a theta-controlled overdispersion", {
  skip_if_not_installed("MASS")
  out <- gen_nbin_data(n_obs = 500, n_vars = 3, theta = 2,
                       beta = c(0.5, -0.3, 0.0), intercept = 0.5,
                       snr_db = Inf, seed = 11)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("y", paste0("V", 1:3)) %in% names(out)))
  expect_true(all(out$y >= 0))
  expect_true(all(out$y == floor(out$y)))
})
