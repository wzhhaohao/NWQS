#' @keywords internal
"_PACKAGE"

#' @importFrom stats rnorm rbinom rpois runif cor cov2cor reshape uniroot
NULL

utils::globalVariables(c(

  # dplyr / tidyverse NSE
  "Target", "Term", "True_Value", "Estimate", "Bias",
  "Wald_CI_Lower", "Wald_CI_Upper", "Covered_Wald",
  "Empirical_CI_Lower", "Empirical_CI_Upper", "Covered_Empirical",
  "Abs_Bias",

  # ggplot2 aes() variables
  "Component", "Weight", "Lower", "Upper",
  "y", "ymin", "ymax", "lower", "middle", "upper",
  "Model", "Deviance", "SAE", "Estimated_Weight",
  "Quantile", "Effect"
))
