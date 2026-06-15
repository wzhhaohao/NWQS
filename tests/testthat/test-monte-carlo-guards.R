# Guardrail tests: the RH-based coverage helpers must warn that their CIs
# reflect algorithmic (data-splitting) variance only and are not valid
# inference. Valid coverage is via nwqs_boot() + check_boot_coverage().

test_that("check_coverage warns that RH-based coverage is not valid inference (S3)", {
  true_mat <- matrix(
    c(0.5, 0.3), nrow = 2,
    dimnames = list(c("Overall", "M1"), c("P75_vs_P50"))
  )
  est_df <- data.frame(
    Target             = "P75_vs_P50",
    Term               = c("Overall", "M1"),
    Estimate           = c(0.40, 0.25),
    Wald_CI_Lower      = c(0.10, 0.00),
    Wald_CI_Upper      = c(0.70, 0.50),
    Empirical_CI_Lower = c(0.05, -0.05),
    Empirical_CI_Upper = c(0.75, 0.55),
    stringsAsFactors   = FALSE
  )
  expect_warning(
    check_coverage(est_df, true_mat),
    "algorithmic|not valid|nwqs_boot"
  )
})

test_that("evaluate_sim_performance warns about RH-based coverage (S3)", {
  expect_warning(
    try(
      evaluate_sim_performance(data.frame(), data.frame(),
                               c(M1 = 1), matrix(0, 1, 1)),
      silent = TRUE
    ),
    "algorithmic|not valid|nwqs_boot"
  )
})
