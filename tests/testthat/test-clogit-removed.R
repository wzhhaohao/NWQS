# Verifies that the clogit family — which used to be supported in 0.1.x —
# now produces a clear error from match.arg(). This is the load-bearing
# guarantee that anyone porting a 0.1.x clogit workflow to 0.2.0 sees an
# immediate failure rather than silently falling through to the GLM path.

test_that("nwqs() rejects family = 'clogit'", {
  d <- data.frame(y = rnorm(20), X1 = rnorm(20), X2 = rnorm(20))
  expect_error(
    nwqs(d, mix_name = c("X1", "X2"), outcome = "y", family = "clogit"),
    regexp = "'arg' should be one of"
  )
})

test_that("nwqs_boot() rejects family = 'clogit'", {
  d <- data.frame(y = rnorm(20), X1 = rnorm(20), X2 = rnorm(20))
  expect_error(
    nwqs_boot(d, mix_name = c("X1", "X2"), outcome = "y", family = "clogit"),
    regexp = "'arg' should be one of"
  )
})


