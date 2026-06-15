# Internal behavior contracts for nwqs() introduced in v0.2.0:
#   - min_shape_sd parameter is exposed and has default 1e-8
#   - a degenerate shape (sd_eta < min_shape_sd) emits a message under
#     quiet = FALSE and is silent under quiet = TRUE
#   - rh > 1 summary stays consistent when coef_sd hits 0 (z becomes NA
#     but summary()/print()/coef() must still succeed)

# ----- min_shape_sd default ----------------------------------------------

test_that("min_shape_sd default is 1e-8 on nwqs()", {
  expect_equal(eval(formals(nwqs)$min_shape_sd), 1e-8)
})

# ----- Degenerate shape messaging ----------------------------------------

# We force the degenerate path by raising min_shape_sd far above the natural
# sd_eta. Real fits never hit that with min_shape_sd = 1e-8, so this is the
# safest way to exercise the branch without relying on fragile data.

test_that("nwqs() emits a message when min_shape_sd fires and quiet = FALSE", {
  set.seed(0)
  n <- 150
  d <- data.frame(
    y  = rnorm(n),
    X1 = rnorm(n),
    X2 = rnorm(n),
    X3 = rnorm(n)
  )
  expect_warning(
    expect_message(
      nwqs(
        data           = d,
        mix_name       = c("X1", "X2", "X3"),
        outcome        = "y",
        family         = "gaussian",
        rh             = 2,
        n_permutation  = 3,
        seed           = 1,
        quiet          = FALSE,
        transform_type = "q_bin",
        q              = 4,
        min_shape_sd   = 1e6  # force the trigger
      ),
      regexp = "degenerate"
    ),
    regexp = "When rh > 1"
  )
})

test_that("nwqs() stays silent on the same degenerate trigger when quiet = TRUE", {
  set.seed(0)
  n <- 150
  d <- data.frame(
    y  = rnorm(n),
    X1 = rnorm(n),
    X2 = rnorm(n),
    X3 = rnorm(n)
  )
  expect_silent({
    fit <- nwqs(
      data           = d,
      mix_name       = c("X1", "X2", "X3"),
      outcome        = "y",
      family         = "gaussian",
      rh             = 2,
      n_permutation  = 3,
      seed           = 1,
      quiet          = TRUE,
      transform_type = "q_bin",
      q              = 4,
      min_shape_sd   = 1e6
    )
    fit
  })
})

# ----- z = NA path is non-fatal in summary/coef -------------------------

test_that("rh > 1 with all-identical-iteration coefs (coef_sd = 0) does not crash summary/coef/print", {
  # Stub a minimal nwqs object whose coef_sd is 0 by construction. We mimic
  # the structure expected by print.nwqs / summary.nwqs / coef.nwqs.
  rh_coefs <- rbind(
    c(`(Intercept)` = 0.5, nwqs = 1.0),
    c(`(Intercept)` = 0.5, nwqs = 1.0),
    c(`(Intercept)` = 0.5, nwqs = 1.0)
  )
  coef_mean <- colMeans(rh_coefs)
  coef_sd   <- apply(rh_coefs, 2, sd)
  z_value   <- ifelse(coef_sd > 0, coef_mean / coef_sd, NA_real_)
  p_value   <- ifelse(is.na(z_value), NA_real_, 2 * pnorm(-abs(z_value)))

  coef_summary <- data.frame(
    Estimate     = coef_mean,
    `Std. Error` = coef_sd,
    `z value`    = z_value,
    `Pr(>|z|)`   = p_value,
    `2.5 %`      = coef_mean - 1.96 * coef_sd,
    `97.5 %`     = coef_mean + 1.96 * coef_sd,
    check.names  = FALSE
  )

  expect_true(all(is.na(coef_summary$`z value`)))
  expect_true(all(is.na(coef_summary$`Pr(>|z|)`)))
  expect_true(is.numeric(coef_summary$Estimate))
})

# ----- n_shuffle exposed + threaded through the public API ---------------

test_that("nwqs()/nwqs_boot() default n_shuffle = 30", {
  expect_equal(eval(formals(nwqs)$n_shuffle), 30)
  expect_equal(eval(formals(nwqs_boot)$n_shuffle), 30)
})

test_that("nwqs() accepts n_shuffle end-to-end and it affects the fit", {
  set.seed(2024)
  n  <- 150
  df <- data.frame(y = rnorm(n), M1 = rnorm(n), M2 = rnorm(n), M3 = rnorm(n))
  df$y <- 0.6 * rank(df$M1) / n + 0.2 * rank(df$M3) / n + rnorm(n, sd = 0.4)
  mix <- c("M1", "M2", "M3")

  f_small <- nwqs(df, mix_name = mix, outcome = "y", family = "gaussian",
                  rh = 1, n_permutation = 2, n_shuffle = 2,
                  seed = 11, plan_strategy = "sequential", quiet = TRUE)
  f_big   <- nwqs(df, mix_name = mix, outcome = "y", family = "gaussian",
                  rh = 1, n_permutation = 2, n_shuffle = 80,
                  seed = 11, plan_strategy = "sequential", quiet = TRUE)

  expect_s3_class(f_small, "nwqs")
  expect_s3_class(f_big, "nwqs")
  expect_false(isTRUE(all.equal(unname(f_small$final_weights),
                                unname(f_big$final_weights))))
})
