#' @title Add Gaussian Noise on the Linear-Predictor Scale at a Target SNR
#'
#' @description
#' Injects independent \eqn{N(0, \sigma^2)} noise onto a clean signal vector
#' so that the resulting signal-to-noise ratio (on whatever scale
#' \code{signal_vec} represents) matches the requested target.
#'
#' @details
#' The signal-to-noise ratio used by this function is the **link-scale SNR**:
#' \deqn{\mathrm{SNR} \;=\; \frac{\mathrm{Var}(X\beta)}{\mathrm{Var}(\varepsilon)},
#'       \qquad \mathrm{SNR}_{\mathrm{dB}} \;=\; 10\,\log_{10}\!\bigl(\mathrm{SNR}\bigr).}
#' Given a target \code{snr_db}, the noise standard deviation is solved from
#' \deqn{\sigma_{\varepsilon} \;=\; \sqrt{\mathrm{Var}(X\beta) \;/\; 10^{\mathrm{snr\_db}/10}}.}
#'
#' Callers are expected to pass the linear predictor \eqn{\eta = X\beta} as
#' \code{signal_vec}. The result is then \eqn{\eta + \varepsilon}; the caller
#' is responsible for applying the family's inverse link and any
#' family-specific sampling step. For Gaussian \code{family} with identity
#' link the link-scale SNR coincides with the outcome-scale SNR. For
#' non-Gaussian families it does not: the response variance also includes
#' the family's inherent sampling variability (Bernoulli, Poisson, ...).
#'
#' This contract is verified by the cross-family test
#' \code{tests/testthat/test-snr.R}, which asserts
#' \eqn{\mathrm{Var}(\eta)/\mathrm{Var}(\varepsilon) \approx 10^{\mathrm{snr\_db}/10}}
#' at the noise injection site for each generator path.
#'
#' @param signal_vec Numeric vector. The noise-free linear predictor
#'   \eqn{\eta = X\beta}. Length must be at least 2.
#' @param snr_db Numeric. Target link-scale signal-to-noise ratio in
#'   decibels. Higher values mean cleaner signal; \code{Inf} returns
#'   \code{signal_vec} unchanged.
#'
#' @return Numeric vector of the same length as \code{signal_vec}: the noisy
#'   linear predictor.
#'
#' @export
add_noise_by_snr <- function(signal_vec, snr_db) {
  stopifnot(is.numeric(signal_vec), length(signal_vec) > 1)

  power_signal <- mean((signal_vec - mean(signal_vec))^2)

  if (power_signal == 0) {
    return(signal_vec)
  }

  sigma_noise <- sqrt(power_signal / 10^(snr_db / 10))
  noise_vec <- rnorm(length(signal_vec), mean = 0, sd = sigma_noise)

  return(signal_vec + noise_vec)
}

#' @title Generate Covariance Matrix with Specific Correlation Structure
#'
#' @description
#' Generates a positive-definite symmetric correlation/covariance matrix with
#' a specified correlation pattern, designed for simulating high-dimensional
#' collinear exposure data.
#'
#' @param n_vars Integer. Number of mixture variables.
#' @param mode Character. Correlation pattern: \code{"low"}, \code{"medium"},
#'   \code{"high"}, or \code{"mixed"} (block-diagonal structure).
#' @param rho Numeric. Base correlation strength (for certain modes). Default
#'   is 0.7.
#' @param seed Integer or \code{NULL}. Random seed for reproducibility.
#'
#' @return A positive-definite symmetric correlation matrix.
#'
#' @export
generate_sigma <- function(n_vars, mode = c("medium", "low", "high", "mixed"),
                           rho = 0.7, seed = NULL) {
  mode <- match.arg(mode)
  if (!is.null(seed)) set.seed(seed)

  if (mode == "low") {
    A <- diag(n_vars) + matrix(runif(n_vars^2, -0.1, 0.1), nrow = n_vars)
    sigma <- cov2cor(t(A) %*% A)
  } else if (mode == "medium") {
    A <- matrix(runif(n_vars^2, -1, 1), ncol = n_vars)
    sigma <- cov2cor(t(A) %*% A)
  } else if (mode == "high") {
    sigma <- matrix(rho, nrow = n_vars, ncol = n_vars)
    diag(sigma) <- 1

    noise <- matrix(runif(n_vars^2, -0.05, 0.05), nrow = n_vars)
    sigma <- sigma + (noise + t(noise)) / 2

    eig <- eigen(sigma)
    val <- pmax(eig$values, 0.01)
    sigma <- cov2cor(eig$vectors %*% diag(val) %*% t(eig$vectors))
  } else if (mode == "mixed") {
    split_idx <- floor(n_vars / 2)
    s1 <- split_idx
    s2 <- n_vars - split_idx

    B1 <- matrix(0.8, nrow = s1, ncol = s1)
    diag(B1) <- 1

    A2 <- matrix(runif(s2^2, -1, 1), ncol = s2)
    B2 <- cov2cor(t(A2) %*% A2)

    sigma <- matrix(0, nrow = n_vars, ncol = n_vars)
    sigma[1:s1, 1:s1] <- B1
    sigma[(s1 + 1):n_vars, (s1 + 1):n_vars] <- B2

    noise <- matrix(runif(n_vars^2, -0.1, 0.1), nrow = n_vars)
    sigma <- sigma + (noise + t(noise)) / 2

    eig <- eigen(sigma)
    val <- pmax(eig$values, 0.01)
    sigma <- cov2cor(eig$vectors %*% diag(val) %*% t(eig$vectors))
  }

  return(sigma)
}

#' @title Generate Covariates for Epidemiological Simulation
#'
#' @description
#' Randomly generates a covariate dataset containing continuous, binary, and
#' categorical variables, and computes their true linear effect contribution.
#'
#' @param n_obs Integer. Sample size.
#' @param beta_cont Numeric. True regression coefficient for the continuous
#'   variable.
#' @param beta_bin Numeric. True regression coefficient for the binary
#'   variable.
#' @param beta_cat Numeric vector. True regression coefficients for the
#'   categorical variable levels.
#' @param prob_bin Numeric. Probability for the binary variable.
#' @param prob_cat Numeric vector. Level probabilities for the categorical
#'   variable.
#' @param Intercept Numeric. Baseline intercept.
#'
#' @return A list with elements: \code{original} (raw covariate data.frame),
#'   \code{mm} (design matrix form), and \code{eta_cov} (true linear predictor
#'   contribution).
#'
#' @export
generate_covariates <- function(n_obs = 1000,
                                beta_cont = 0.5,
                                beta_bin = -0.8,
                                beta_cat = c(0, -0.5, 0.7),
                                prob_bin = 0.5,
                                prob_cat = c(1 / 3, 1 / 3, 1 / 3),
                                Intercept = 0) {
  x_cont <- rnorm(n_obs, 0, 1)
  x_bin_raw <- rbinom(n_obs, 1, prob_bin)
  x_cat_raw <- sample(1:3, n_obs, replace = TRUE, prob = prob_cat)

  x_bin <- factor(x_bin_raw, levels = c(0, 1))
  x_cat <- factor(x_cat_raw, levels = 1:3)

  eta_cov <- beta_cont * x_cont + beta_bin * x_bin_raw + beta_cat[x_cat_raw] + Intercept

  df_raw <- data.frame(x_cont, x_bin, x_cat)

  df_result <- as.data.frame(cbind(eta_cov = eta_cov, df_raw))

  list(mm = df_result, original = df_raw, eta_cov = eta_cov)
}


#' @title Generate Linear Mixture Effect Data
#'
#' @description
#' Generates continuous outcome data with linear mixture effects for Monte
#' Carlo simulation. Exposure components are drawn from a multivariate normal
#' distribution, optionally quantile-transformed, combined with covariate
#' effects, and corrupted by noise according to the specified SNR.
#'
#' @param n_obs Integer. Sample size.
#' @param mu_preds Numeric vector. Mean vector for multivariate normal
#'   exposure generation.
#' @param sigma_preds Matrix. Covariance matrix for exposure generation.
#' @param beta_wqs Numeric. Overall mixture effect coefficient.
#' @param beta_preds Numeric vector. Component-specific effect weights.
#' @param snr_db Numeric. Signal-to-noise ratio in dB. Default is 10.
#' @param transform_fun Function or \code{NULL}. Custom transformation for
#'   exposures.
#' @param seed Integer or \code{NULL}. Random seed.
#' @param ... Additional arguments passed to \code{generate_covariates}.
#'
#' @export
generate_linear_data <- function(n_obs = 1000,
                                 mu_preds,
                                 sigma_preds,
                                 beta_wqs = 1,
                                 beta_preds,
                                 snr_db = 10,
                                 transform_fun = NULL,
                                 seed = NULL,
                                 ...) {
  if (!is.null(seed)) set.seed(seed)
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' required")
  }

  preds_raw <- MASS::mvrnorm(
    n_obs,
    mu = mu_preds,
    Sigma = sigma_preds
  )

  preds_scaled <- as.data.frame(scale(preds_raw))
  names(preds_scaled) <- paste0("Component", seq_len(ncol(preds_scaled)))

  if (!is.null(transform_fun) && is.function(transform_fun)) {
    preds_final <- transform_fun(preds_scaled)
  } else {
    preds_final <- preds_scaled
  }

  preds_final <- as.data.frame(preds_final)
  names(preds_final) <- names(preds_scaled)

  cov_list <- generate_covariates(n_obs = n_obs, ...)

  beta_preds <- beta_wqs * beta_preds

  y_clean <- as.matrix(preds_final) %*% beta_preds + cov_list$eta_cov
  y_observed <- add_noise_by_snr(as.vector(y_clean), snr_db = snr_db)

  cols_cov <- setdiff(names(cov_list$mm), "eta_cov")

  final_df <- cbind(
    y = y_observed,
    preds_scaled,
    cov_list$mm[, cols_cov, drop = FALSE]
  )

  as.data.frame(final_df)
}


#' @title Generate Non-Linear Spline Dose-Response Data (Continuous Outcome)
#'
#' @description
#' Advanced data generation function that allows specifying component-specific
#' non-linear dose-response trajectories (e.g., U-shaped, S-shaped, threshold
#' effects). Uses natural cubic splines to precisely map these shapes and
#' synthesizes a continuous outcome variable.
#'
#' @param n_obs Integer. Sample size.
#' @param mu_preds Numeric vector. Mean vector for multivariate normal
#'   exposure generation.
#' @param sigma_preds Matrix. Covariance matrix for exposure generation.
#' @param beta_wqs Numeric. Overall mixture effect coefficient.
#' @param beta_preds Numeric vector. Component-specific effect weights.
#' @param snr_db Numeric. Signal-to-noise ratio in dB. Default is 10.
#' @param transform_fun Function or \code{NULL}. Custom transformation for
#'   exposures.
#' @param q Integer. Number of quantile bins. Default is 4.
#' @param df_spline Integer. Degrees of freedom for natural cubic splines.
#'   Default is 3.
#' @param seed Integer or \code{NULL}. Random seed.
#' @param shape Character or character vector. Controls the true causal
#'   dose-response curve shape for each component. Options: \code{"linear_like"},
#'   \code{"u_shape"}, \code{"inv_u_shape"}, \code{"s_shape"},
#'   \code{"threshold"}, \code{"neg_linear"}, \code{"inv_threshold"}.
#' @param ... Additional arguments passed to \code{generate_covariates}.
#' @export
gen_nonlinear_data <- function(n_obs = 1000,
                               mu_preds,
                               sigma_preds,
                               beta_wqs = 1,
                               beta_preds,
                               snr_db = 10,
                               transform_fun = NULL,
                               q = 4,
                               df_spline = 3,
                               seed = NULL,
                               shape = "linear_like",
                               ...) {
  if (!requireNamespace("splines", quietly = TRUE)) stop("Package 'splines' required")
  if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' required")
  if (!is.null(seed)) set.seed(seed)

  preds_raw <- MASS::mvrnorm(n_obs, mu = mu_preds, Sigma = sigma_preds)
  preds_scaled <- as.data.frame(scale(preds_raw))
  n_vars <- ncol(preds_scaled)
  names(preds_scaled) <- paste0("Component", 1:n_vars)

  if (!is.null(transform_fun) && is.function(transform_fun)) {
    preds_trans <- transform_fun(preds_scaled)
  } else {
    preds_trans <- preds_scaled
  }

  mat_spline_list <- lapply(preds_trans, function(x) splines::ns(x, df = df_spline))

  if (length(beta_preds) != n_vars) stop("Length of 'beta_preds' must match n_vars.")

  cov_list <- generate_covariates(n_obs = n_obs, ...)

  eval_pts <- 0:(q - 1)

  temp_spline <- splines::ns(eval_pts, df = df_spline)
  global_knots <- attr(temp_spline, "knots")
  global_boundary <- attr(temp_spline, "Boundary.knots")

  mat_spline_list <- lapply(preds_trans, function(x) {
    splines::ns(x, df = df_spline, knots = global_knots, Boundary.knots = global_boundary)
  })

  basis_std_true <- splines::ns(eval_pts, df = df_spline, knots = global_knots, Boundary.knots = global_boundary, intercept = FALSE)

  if (length(shape) == 1) shape <- rep(shape, n_vars)

  eta_components_raw <- matrix(0, nrow = n_obs, ncol = n_vars)
  baseline_components <- numeric(n_vars)

  true_eff_mat <- matrix(0, nrow = n_vars + 1, ncol = q - 1)
  rownames(true_eff_mat) <- c("Overall", names(preds_scaled))
  colnames(true_eff_mat) <- paste0("Q", 2:q, "_vs_Q1")

  for (i in 1:n_vars) {
    current_shape <- shape[i]
    b <- beta_preds[i] * beta_wqs
    min_x <- min(preds_trans[, i])

    if (current_shape == "linear_like") {
      pattern <- c(1, 2, 3)
    } else if (current_shape == "neg_linear") {
      pattern <- c(-1, -2, -3)
    } else if (current_shape == "u_shape") {
      pattern <- c(1.5, -3.0, 1.5)
    } else if (current_shape == "inv_u_shape") {
      pattern <- c(-1.5, 3.0, -1.5)
    } else if (current_shape == "s_shape") {
      pattern <- c(1, -1, 1)
    } else if (current_shape == "threshold") {
      pattern <- c(0, 0.5, 4.0)
    } else if (current_shape == "inv_threshold") {
      pattern <- c(0, -0.5, -4.0)
    } else {
      pattern <- rep(1, df_spline)
    }

    comp_beta <- b * pattern
    eta_components_raw[, i] <- as.vector(mat_spline_list[[i]] %*% comp_beta)
    baseline_components[i] <- as.vector(predict(mat_spline_list[[i]], newx = min_x) %*% comp_beta)

    for (k in 2:q) {
      b_diff <- basis_std_true[k, ] - basis_std_true[1, ]
      true_eff_mat[i + 1, k - 1] <- sum(b_diff * comp_beta)
    }
  }

  for (k in 2:q) {
    true_eff_mat[1, k - 1] <- sum(true_eff_mat[2:(n_vars + 1), k - 1])
  }

  eta_spline_raw <- rowSums(eta_components_raw)
  baseline_effect <- sum(baseline_components)

  eta_spline_adjusted <- eta_spline_raw - baseline_effect
  y_clean <- eta_spline_adjusted + cov_list$eta_cov

  y_observed <- add_noise_by_snr(as.vector(y_clean), snr_db = snr_db)

  cols_cov <- setdiff(names(cov_list$mm), "eta_cov")
  final_df <- cbind(y = y_observed, preds_scaled, cov_list$mm[, cols_cov, drop = FALSE])

  attr(final_df, "true_effect_mat") <- true_eff_mat
  attr(final_df, "spline_knots") <- global_knots
  attr(final_df, "spline_boundary") <- global_boundary

  return(as.data.frame(final_df))
}

#' @title Generate Non-Linear Binary Outcome Data
#'
#' @description
#' Generates non-linear binary outcome data using a specified link function
#' (logit, probit, or cloglog). Supports automatic intercept calibration to
#' control rare event incidence rates. The \code{true_effect_mat} attribute
#' stores true log-odds ratios.
#'
#' @param n_obs Integer. Sample size.
#' @param mu_preds Numeric vector. Mean vector for multivariate normal
#'   exposure generation.
#' @param sigma_preds Matrix. Covariance matrix for exposure generation.
#' @param beta_wqs Numeric. Overall mixture effect coefficient.
#' @param beta_preds Numeric vector. Component-specific effect weights.
#' @param intercept Numeric. Baseline model intercept on the log-odds scale.
#' @param target_prop Numeric in (0, 1) or \code{NULL}. Target disease
#'   incidence rate. If provided, the intercept is automatically calibrated.
#' @param link Character. Link function: \code{"logit"}, \code{"probit"}, or
#'   \code{"cloglog"}.
#' @param snr_db Numeric. Signal-to-noise ratio in dB. Default is
#'   \code{Inf} (no noise).
#' @param transform_fun Function or \code{NULL}. Custom transformation for
#'   exposures.
#' @param q Integer. Number of quantile bins. Default is 4.
#' @param df_spline Integer. Degrees of freedom for natural cubic splines.
#'   Default is 3.
#' @param seed Integer or \code{NULL}. Random seed.
#' @param shape Character or character vector. Dose-response curve shape(s).
#' @param ... Additional arguments passed to \code{generate_covariates}.
#'
#' @export
gen_nonlinear_bio_data <- function(n_obs = 1000, mu_preds, sigma_preds,
                                   beta_wqs = 1, beta_preds,
                                   intercept = 0, target_prop = NULL,
                                   link = c("logit", "probit", "cloglog"),
                                   snr_db = Inf, transform_fun = NULL,
                                   q = 4, df_spline = 3, seed = NULL,
                                   shape = "linear_like", ...) {
  if (!requireNamespace("splines", quietly = TRUE)) stop("Package 'splines' required")
  if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' required")

  link <- match.arg(link)
  if (!is.null(seed)) set.seed(seed)

  preds_raw <- MASS::mvrnorm(n_obs, mu = mu_preds, Sigma = sigma_preds)
  preds_scaled <- as.data.frame(scale(preds_raw))
  n_vars <- ncol(preds_scaled)
  names(preds_scaled) <- paste0("Component", 1:ncol(preds_scaled))

  if (!is.null(transform_fun) && is.function(transform_fun)) {
    preds_trans <- transform_fun(preds_scaled)
  } else {
    preds_trans <- preds_scaled
  }

  if (length(beta_preds) != n_vars) stop("Length of 'beta_preds' must match n_vars.")

  cov_list <- generate_covariates(n_obs = n_obs, ...)

  eval_pts <- 0:(q - 1)
  temp_spline <- splines::ns(eval_pts, df = df_spline)
  global_knots <- attr(temp_spline, "knots")
  global_boundary <- attr(temp_spline, "Boundary.knots")

  mat_spline_list <- lapply(preds_trans, function(x) {
    splines::ns(x, df = df_spline, knots = global_knots, Boundary.knots = global_boundary)
  })

  basis_std_true <- splines::ns(eval_pts, df = df_spline, knots = global_knots, Boundary.knots = global_boundary, intercept = FALSE)

  if (length(shape) == 1) shape <- rep(shape, n_vars)

  eta_components_raw <- matrix(0, nrow = n_obs, ncol = n_vars)
  baseline_components <- numeric(n_vars)

  true_eff_mat <- matrix(0, nrow = n_vars + 1, ncol = q - 1)
  rownames(true_eff_mat) <- c("Overall", names(preds_scaled))
  colnames(true_eff_mat) <- paste0("Q", 2:q, "_vs_Q1")

  for (i in 1:n_vars) {
    current_shape <- shape[i]
    b <- beta_preds[i] * beta_wqs
    min_x <- min(preds_trans[, i])

    if (current_shape == "pure_linear") {
      eta_components_raw[, i] <- preds_trans[, i] * b
      baseline_components[i] <- min_x * b

      for (k in 2:q) {
        true_eff_mat[i + 1, k - 1] <- (eval_pts[k] - eval_pts[1]) * b
      }
    } else {
      if (current_shape == "linear_like") {
        pattern <- c(1, 2, 3)
      } else if (current_shape == "neg_linear") {
        pattern <- c(-1, -2, -3)
      } else if (current_shape == "u_shape") {
        pattern <- c(1.5, -3.0, 1.5)
      } else if (current_shape == "inv_u_shape") {
        pattern <- c(-1.5, 3.0, -1.5)
      } else if (current_shape == "s_shape") {
        pattern <- c(1, -1, 1)
      } else if (current_shape == "threshold") {
        pattern <- c(0, 0.5, 4.0)
      } else if (current_shape == "inv_threshold") {
        pattern <- c(0, -0.5, -4.0)
      } else {
        pattern <- rep(1, df_spline)
      }

      comp_beta <- b * pattern

      eta_components_raw[, i] <- as.vector(mat_spline_list[[i]] %*% comp_beta)
      baseline_components[i] <- as.vector(predict(mat_spline_list[[i]], newx = min_x) %*% comp_beta)

      for (k in 2:q) {
        b_diff <- basis_std_true[k, ] - basis_std_true[1, ]
        true_eff_mat[i + 1, k - 1] <- sum(b_diff * comp_beta)
      }
    }
  }

  for (k in 2:q) {
    true_eff_mat[1, k - 1] <- sum(true_eff_mat[2:(n_vars + 1), k - 1])
  }

  eta_spline_raw <- rowSums(eta_components_raw)
  baseline_effect <- sum(baseline_components)

  eta_spline_adjusted <- eta_spline_raw - baseline_effect
  eta_partial <- eta_spline_adjusted + cov_list$eta_cov

  if (!is.null(snr_db) && is.finite(snr_db)) {
    eta_noisy_partial <- add_noise_by_snr(eta_partial, snr_db)
  } else {
    eta_noisy_partial <- eta_partial
  }

  final_intercept <- intercept
  if (!is.null(target_prop)) {
    calc_mean_prob_diff <- function(b0) {
      eta_temp <- b0 + eta_noisy_partial
      if (link == "logit") {
        p <- 1 / (1 + exp(-eta_temp))
      } else if (link == "probit") {
        p <- pnorm(eta_temp)
      } else if (link == "cloglog") p <- 1 - exp(-exp(eta_temp))
      return(mean(p) - target_prop)
    }
    tryCatch(
      {
        final_intercept <- uniroot(calc_mean_prob_diff, interval = c(-50, 50))$root
      },
      error = function(e) {}
    )
  }

  eta_final <- final_intercept + eta_noisy_partial

  if (link == "logit") {
    probs <- 1 / (1 + exp(-eta_final))
  } else if (link == "probit") {
    probs <- pnorm(eta_final)
  } else if (link == "cloglog") probs <- 1 - exp(-exp(eta_final))

  y_binary <- rbinom(n_obs, size = 1, prob = probs)

  cols_cov <- setdiff(names(cov_list$mm), "eta_cov")
  final_df <- cbind(y = y_binary, preds_scaled, cov_list$mm[, cols_cov, drop = FALSE])

  attr(final_df, "true_effect_mat") <- true_eff_mat
  attr(final_df, "true_prob") <- probs
  attr(final_df, "spline_knots") <- global_knots
  attr(final_df, "spline_boundary") <- global_boundary
  return(as.data.frame(final_df))
}


#' @title Generate Non-Linear Poisson Count Data
#'
#' @description
#' Maps non-linear spline features of mixture exposures to expected event
#' rates (\eqn{\lambda}) via a log link and generates discrete count outcomes
#' from a Poisson process. The \code{true_effect_mat} attribute stores true
#' log rate ratios (log-RR).
#'
#' @param n_obs Integer. Sample size.
#' @param mu_preds Numeric vector. Mean vector for multivariate normal
#'   exposure generation.
#' @param sigma_preds Matrix. Covariance matrix for exposure generation.
#' @param beta_wqs Numeric. Overall mixture effect coefficient.
#' @param beta_preds Numeric vector. Component-specific effect weights.
#' @param intercept Numeric. Model intercept on the log scale, directly
#'   controlling the baseline Poisson rate \eqn{\lambda}.
#' @param snr_db Numeric. Signal-to-noise ratio in dB. Default is
#'   \code{Inf} (no noise).
#' @param transform_fun Function or \code{NULL}. Custom transformation for
#'   exposures.
#' @param q Integer. Number of quantile bins. Default is 4.
#' @param df_spline Integer. Degrees of freedom for natural cubic splines.
#'   Default is 3.
#' @param seed Integer or \code{NULL}. Random seed.
#' @param shape Character or character vector. Dose-response curve shape(s).
#' @param ... Additional arguments passed to \code{generate_covariates}.
#'
#' @export
gen_nonlinear_count_data <- function(n_obs = 1000, mu_preds, sigma_preds,
                                     beta_wqs = 1, beta_preds,
                                     intercept = 0, snr_db = Inf,
                                     transform_fun = NULL, q = 4,
                                     df_spline = 3, seed = NULL,
                                     shape = "linear_like", ...) {
  if (!requireNamespace("splines", quietly = TRUE)) stop("Package 'splines' required")
  if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' required")
  if (!is.null(seed)) set.seed(seed)

  preds_raw <- MASS::mvrnorm(n_obs, mu = mu_preds, Sigma = sigma_preds)
  preds_scaled <- as.data.frame(scale(preds_raw))
  n_vars <- ncol(preds_scaled)
  names(preds_scaled) <- paste0("Component", 1:n_vars)

  if (!is.null(transform_fun) && is.function(transform_fun)) {
    preds_trans <- transform_fun(preds_scaled)
  } else {
    preds_trans <- preds_scaled
  }

  if (length(beta_preds) != n_vars) stop("Length of 'beta_preds' must match n_vars.")

  cov_list <- generate_covariates(n_obs = n_obs, ...)

  eval_pts <- 0:(q - 1)
  temp_spline <- splines::ns(eval_pts, df = df_spline)
  global_knots <- attr(temp_spline, "knots")
  global_boundary <- attr(temp_spline, "Boundary.knots")

  mat_spline_list <- lapply(preds_trans, function(x) {
    splines::ns(x, df = df_spline, knots = global_knots, Boundary.knots = global_boundary)
  })

  basis_std_true <- splines::ns(eval_pts, df = df_spline, knots = global_knots, Boundary.knots = global_boundary, intercept = FALSE)

  if (length(shape) == 1) shape <- rep(shape, n_vars)

  eta_components_raw <- matrix(0, nrow = n_obs, ncol = n_vars)
  baseline_components <- numeric(n_vars)

  true_eff_mat <- matrix(0, nrow = n_vars + 1, ncol = q - 1)
  rownames(true_eff_mat) <- c("Overall", names(preds_scaled))
  colnames(true_eff_mat) <- paste0("Q", 2:q, "_vs_Q1")

  for (i in 1:n_vars) {
    current_shape <- shape[i]
    b <- beta_preds[i] * beta_wqs
    min_x <- min(preds_trans[, i])

    if (current_shape == "pure_linear") {
      eta_components_raw[, i] <- preds_trans[, i] * b
      baseline_components[i] <- min_x * b

      for (k in 2:q) {
        true_eff_mat[i + 1, k - 1] <- (eval_pts[k] - eval_pts[1]) * b
      }
    } else {
      if (current_shape == "linear_like") {
        pattern <- c(1, 2, 3)
      } else if (current_shape == "neg_linear") {
        pattern <- c(-1, -2, -3)
      } else if (current_shape == "u_shape") {
        pattern <- c(1.5, -3.0, 1.5)
      } else if (current_shape == "inv_u_shape") {
        pattern <- c(-1.5, 3.0, -1.5)
      } else if (current_shape == "s_shape") {
        pattern <- c(1, -1, 1)
      } else if (current_shape == "threshold") {
        pattern <- c(0, 0.5, 4.0)
      } else if (current_shape == "inv_threshold") {
        pattern <- c(0, -0.5, -4.0)
      } else {
        pattern <- rep(1, df_spline)
      }

      comp_beta <- b * pattern

      eta_components_raw[, i] <- as.vector(mat_spline_list[[i]] %*% comp_beta)
      baseline_components[i] <- as.vector(predict(mat_spline_list[[i]], newx = min_x) %*% comp_beta)

      for (k in 2:q) {
        b_diff <- basis_std_true[k, ] - basis_std_true[1, ]
        true_eff_mat[i + 1, k - 1] <- sum(b_diff * comp_beta)
      }
    }
  }

  for (k in 2:q) {
    true_eff_mat[1, k - 1] <- sum(true_eff_mat[2:(n_vars + 1), k - 1])
  }

  eta_spline_raw <- rowSums(eta_components_raw)
  baseline_effect <- sum(baseline_components)

  eta_spline_adjusted <- eta_spline_raw - baseline_effect
  eta_partial <- eta_spline_adjusted + cov_list$eta_cov

  if (!is.null(snr_db) && is.finite(snr_db)) {
    eta_noisy_partial <- add_noise_by_snr(eta_partial, snr_db)
  } else {
    eta_noisy_partial <- eta_partial
  }

  eta_final <- intercept + eta_noisy_partial
  lambda <- exp(eta_final)

  if (any(lambda > 10000)) warning("Extremely high lambda values.")

  y_count <- rpois(n_obs, lambda = lambda)

  cols_cov <- setdiff(names(cov_list$mm), "eta_cov")
  final_df <- cbind(y = y_count, preds_scaled, cov_list$mm[, cols_cov, drop = FALSE])

  attr(final_df, "true_effect_mat") <- true_eff_mat
  attr(final_df, "spline_knots") <- global_knots
  attr(final_df, "spline_boundary") <- global_boundary
  return(as.data.frame(final_df))
}


#' @title Generate Negative-Binomial Outcome Data with a Linear Predictor
#'
#' @description
#' Lightweight generator for overdispersed count data with a known true
#' linear predictor. Used as a fixture in \code{tests/testthat/test-family-negbin.R}
#' and as a starting point for the upcoming applied-domain vignette.
#'
#' @details
#' Samples independent standard-normal predictors \eqn{V_1, \ldots, V_p},
#' constructs the linear predictor
#' \deqn{\eta = \alpha + \sum_j \beta_j V_j,}
#' optionally injects Gaussian noise on the link scale at the requested
#' \code{snr_db} (see \code{\link{add_noise_by_snr}}), then draws
#' \eqn{y \sim \mathrm{NB}(\mu = \exp(\eta),\, \theta)} via
#' \code{MASS::rnegbin}.
#'
#' Predictors are returned as columns \code{V1}, \code{V2}, ... to keep
#' the fixture stable across releases.
#'
#' @param n_obs Integer. Sample size.
#' @param n_vars Integer. Number of predictors. Default is 3.
#' @param beta Numeric vector of length \code{n_vars}. Coefficients on
#'   the link scale.
#' @param intercept Numeric. Intercept on the link scale. Default is 0.
#' @param theta Numeric. Negative-binomial dispersion parameter (smaller =
#'   more overdispersion). Default is 2.
#' @param snr_db Numeric. Link-scale signal-to-noise ratio in dB used by
#'   \code{add_noise_by_snr()}. Default is \code{Inf} (no noise).
#' @param seed Integer or \code{NULL}. Random seed.
#'
#' @return A \code{data.frame} with columns \code{y, V1, V2, ..., V_{n_vars}}.
#'
#' @export
gen_nbin_data <- function(n_obs = 500, n_vars = 3, beta = rep(0, n_vars),
                          intercept = 0, theta = 2, snr_db = Inf,
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for gen_nbin_data().")
  }
  if (length(beta) != n_vars) {
    stop(sprintf("length(beta) must equal n_vars (got %d vs %d).",
                 length(beta), n_vars))
  }

  X <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
  colnames(X) <- paste0("V", seq_len(n_vars))

  eta <- intercept + as.numeric(X %*% beta)
  eta_obs <- if (is.finite(snr_db)) add_noise_by_snr(eta, snr_db) else eta
  mu <- exp(eta_obs)

  y <- MASS::rnegbin(n = n_obs, mu = mu, theta = theta)

  data.frame(y = y, X)
}
