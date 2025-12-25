test_that("multiplication works", {
  data <- gen_nonlinear_data(
  n_obs = 1000,
  mu_preds = rep(0, n_vars),
  sigma_preds = sigma_preds,
  beta_preds = beta_preds,
  beta_wqs = 3,
  snr_db = 10,
  df_spline = 3,
  seed = NULL
  )
})
