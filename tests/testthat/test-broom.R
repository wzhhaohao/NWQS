# Tests for the conditional broom hooks. We test the underlying NWQS:::
# functions directly (no broom dependency) and, when broom is installed,
# also verify that the S3 dispatch goes through the broom generics
# without crashing.

# ----- Fixture -----------------------------------------------------------

make_broom_data <- function(n = 150, seed = 41) {
  set.seed(seed)
  X <- matrix(rnorm(n * 3), n, 3)
  colnames(X) <- paste0("X", 1:3)
  eta <- 0.5 * X[, 1] - 0.3 * X[, 2]
  y <- eta + rnorm(n, sd = 0.5)
  data.frame(y = y, X)
}

fit_broom <- function(rh = 3) {
  d <- make_broom_data()
  nwqs(
    data           = d,
    mix_name       = c("X1", "X2", "X3"),
    outcome        = "y",
    family         = "gaussian",
    rh             = rh,
    n_permutation  = 3,
    seed           = 1,
    quiet          = TRUE,
    transform_type = "q_bin",
    q              = 4
  )
}

# ----- tidy / glance for nwqs -------------------------------------------

test_that("tidy_nwqs returns a data.frame with the broom canonical column names", {
  fit <- fit_broom()
  out <- NWQS:::tidy_nwqs(fit)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("term", "estimate", "std.error", "statistic", "p.value") %in%
                  names(out)))
  expect_true("nwqs" %in% out$term)
  expect_true(all(is.finite(out$estimate)))
})

test_that("tidy_nwqs(conf.int = TRUE) appends conf.low and conf.high columns", {
  fit <- fit_broom()
  out <- NWQS:::tidy_nwqs(fit, conf.int = TRUE)
  expect_true(all(c("conf.low", "conf.high") %in% names(out)))
  expect_true(all(out$conf.low <= out$estimate))
  expect_true(all(out$estimate <= out$conf.high))
})

test_that("glance_nwqs returns a one-row data.frame with n / family / rh", {
  fit <- fit_broom()
  out <- NWQS:::glance_nwqs(fit)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1)
  expect_true(all(c("n", "family", "transform_type", "rh") %in% names(out)))
  expect_equal(out$family, "gaussian")
})

# ----- tidy / glance for nwqs_boot --------------------------------------

test_that("tidy_nwqs_boot returns estimate / std.error / conf.low / conf.high", {
  d <- make_broom_data()
  boot_fit <- expect_small_boot_warning(nwqs_boot(
    data = d, mix_name = c("X1", "X2", "X3"), outcome = "y",
    family = "gaussian", n_boot = 8, rh_inner = 1, n_permutation = 3,
    seed = 1, quiet = TRUE, transform_type = "q_bin", q = 4
  ))
  out <- NWQS:::tidy_nwqs_boot(boot_fit)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("term", "estimate", "std.error",
                    "conf.low", "conf.high") %in% names(out)))
  expect_true("nwqs" %in% out$term)
})

test_that("glance_nwqs_boot includes n_boot and n_success", {
  d <- make_broom_data()
  boot_fit <- expect_small_boot_warning(nwqs_boot(
    data = d, mix_name = c("X1", "X2", "X3"), outcome = "y",
    family = "gaussian", n_boot = 8, rh_inner = 1, n_permutation = 3,
    seed = 1, quiet = TRUE, transform_type = "q_bin", q = 4
  ))
  out <- NWQS:::glance_nwqs_boot(boot_fit)
  expect_equal(nrow(out), 1)
  expect_true(all(c("n", "family", "transform_type", "n_boot", "n_success",
                    "conf_level") %in% names(out)))
  expect_equal(out$n_boot, 8)
})

# ----- broom S3 dispatch when broom is installed -------------------------

test_that("broom::tidy dispatches to tidy_nwqs when broom is available", {
  skip_if_not_installed("broom")
  fit <- fit_broom()
  out <- broom::tidy(fit)
  expect_s3_class(out, "data.frame")
  expect_true("nwqs" %in% out$term)
})

test_that("broom::glance dispatches to glance_nwqs when broom is available", {
  skip_if_not_installed("broom")
  fit <- fit_broom()
  out <- broom::glance(fit)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1)
})
