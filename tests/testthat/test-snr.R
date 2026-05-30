# Mathematical contract for add_noise_by_snr(): noise is injected on the
# link / linear-predictor scale so that
#   var(signal_vec) / var(noise) ~~ 10^(snr_db / 10)
# holds regardless of which family will subsequently consume the noisy
# linear predictor. Documenting and pinning this contract is the
# precondition for downstream tests and for any future SNR claims in the
# applied-domain vignette.

test_that("add_noise_by_snr() matches the link-scale SNR formula to within 5%", {
  set.seed(123)
  n <- 50000  # large n keeps Monte Carlo error well below the 5% tolerance
  eta <- rnorm(n, mean = 0, sd = 2)

  for (snr_db in c(0, 5, 10, 20)) {
    set.seed(2026)
    eta_noisy <- add_noise_by_snr(eta, snr_db = snr_db)
    noise <- eta_noisy - eta
    observed_ratio <- var(eta) / var(noise)
    target_ratio <- 10^(snr_db / 10)
    rel_err <- abs(observed_ratio - target_ratio) / target_ratio
    expect_lt(rel_err, 0.05,
              label = sprintf("snr_db = %s, observed ratio = %.3f, target = %.3f",
                              snr_db, observed_ratio, target_ratio))
  }
})

test_that("add_noise_by_snr() returns the input unchanged at snr_db = Inf", {
  set.seed(7)
  eta <- rnorm(100, mean = 1, sd = 3)
  expect_identical(add_noise_by_snr(eta, snr_db = Inf), eta)
})

test_that("add_noise_by_snr() returns the input unchanged when the signal is constant", {
  expect_identical(add_noise_by_snr(rep(5, 100), snr_db = 10), rep(5, 100))
})

test_that("add_noise_by_snr() refuses a length-1 vector", {
  expect_error(add_noise_by_snr(1, snr_db = 10), regexp = "length")
})
