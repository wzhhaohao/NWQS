#' Fast GLM-Permutation Scorer for Spline WQS
#'
#' @description
#' Core computational engine for a single internal NWQS iteration. It fits an unpenalized
#' generalized linear model (GLM) on an in-bag bootstrap sample, evaluates the Out-Of-Bag
#' (OOB) loss, and then extracts grouped permutation importance for each mixture component.
#'
#' @details
#' \strong{Workflow:}
#' \enumerate{
#'   \item Fits a standard GLM on the in-bag data to estimate spline basis coefficients.
#'   \item Evaluates baseline OOB loss on the out-of-bag sample.
#'   \item Computes grouped permutation importance by jointly shuffling all spline basis
#'         columns belonging to the same mixture component and recalculating the OOB loss.
#' }
#'
#' \strong{Weight Derivation:}
#' \deqn{Weights_i = \frac{\sqrt{\Delta Loss_i}}{\sum \sqrt{\Delta Loss}}}
#'
#' @param x Numeric matrix. Design matrix containing spline bases and covariates.
#' @param y Numeric vector. Outcome variable.
#' @param mix_name Character vector. Names of original mixture components.
#' @param spline_vars Character vector. Names of all spline basis columns in \code{x}.
#' @param family List. A GLM family object (e.g., \code{gaussian()}, \code{binomial()}).
#' @param n_permutation Integer. Number of permutations to stabilize importance scores.
#'   Defaults to 100.
#' @param ... Additional arguments (currently unused, for compatibility).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{weights}: Normalized importance-based weights for each mixture component.
#'   \item \code{shapes}: Estimated spline coefficients for all mixture components.
#' }
#'
#' @importFrom stats coef predict glm.fit
#' @export
permutation_scorer <- function(x, y, mix_name, spline_vars, family, n_permutation = 100, ...) {
  n_obs <- nrow(x)

  # 1. Internal bootstrap for In-Bag and Out-Of-Bag (OOB) splitting
  idx <- sample(seq_len(n_obs), size = n_obs, replace = TRUE)
  oob_idx <- setdiff(seq_len(n_obs), idx)

  if (length(oob_idx) == 0) {
    return(NULL)
  }

  x_train <- x[idx, , drop = FALSE]
  y_train <- y[idx]
  x_oob <- x[oob_idx, , drop = FALSE]
  y_oob <- y[oob_idx]

  # Remove explicit intercept column if present (to handle manually later)
  int_col <- match("(Intercept)", colnames(x_train))
  if (!is.na(int_col)) {
    x_train_net <- x_train[, -int_col, drop = FALSE]
    x_oob_net <- x_oob[, -int_col, drop = FALSE]
  } else {
    x_train_net <- x_train
    x_oob_net <- x_oob
  }

  fam_name <- family$family
  linkinv <- family$linkinv

  # 2. Family-specific OOB loss calculation
  calc_loss <- function(y_true, mu_pred) {
    if (fam_name == "gaussian") {
      return(mean((y_true - mu_pred)^2))
    }
    if (fam_name == "binomial") {
      mu_pred <- pmax(pmin(mu_pred, 1 - 1e-7), 1e-7)
      return(-2 * mean(y_true * log(mu_pred) + (1 - y_true) * log(1 - mu_pred)))
    }
    if (fam_name %in% c("poisson", "quasipoisson")) {
      mu_pred <- pmax(mu_pred, 1e-7)
      term1 <- ifelse(y_true == 0, 0, y_true * log(y_true / mu_pred))
      return(2 * mean(term1 - (y_true - mu_pred)))
    }
    return(mean((y_true - mu_pred)^2))
  }

  # 3. Unpenalized GLM Fit on In-Bag Data
  x_train_glm <- cbind(Intercept = 1, as.matrix(x_train_net))
  x_oob_glm <- cbind(Intercept = 1, as.matrix(x_oob_net))

  fit <- stats::glm.fit(x = x_train_glm, y = y_train, family = family)
  coef_all <- fit$coefficients
  coef_all[is.na(coef_all)] <- 0 # Force singular/unestimable coefficients to 0

  intercept_val <- unname(coef_all[1])
  coefs_no_int <- coef_all[-1]

  # 4. Baseline OOB Loss
  eta_oob <- as.numeric(x_oob_glm %*% coef_all)
  mu_oob <- linkinv(eta_oob)
  base_loss <- calc_loss(y_oob, mu_oob)

  # 5. Grouped Permutation Importance
  importance_scores <- numeric(length(mix_name))
  names(importance_scores) <- mix_name

  x_oob_shuffled <- x_oob_net
  n_oob <- length(oob_idx)

  for (var in mix_name) {
    target_cols <- grep(paste0("^", var, "_B"), colnames(x_oob_net))

    if (length(target_cols) == 0) {
      warning(paste("No spline basis columns found for mixture component:", var))
      return(NULL)
    }

    shuffled_loss_list <- numeric(n_permutation)

    for (k in seq_len(n_permutation)) {
      shuffle_idx <- sample(n_oob)
      x_oob_shuffled[, target_cols] <- x_oob_net[shuffle_idx, target_cols, drop = FALSE]

      eta_shuffled <- intercept_val + as.numeric(x_oob_shuffled %*% coefs_no_int)
      mu_shuffled <- linkinv(eta_shuffled)

      shuffled_loss_list[k] <- calc_loss(y_oob, mu_shuffled)
    }

    # Restore original un-shuffled state for the next component
    x_oob_shuffled[, target_cols] <- x_oob_net[, target_cols, drop = FALSE]

    # Calculate Delta Loss (max with 0 to prevent negative importance)
    importance_scores[var] <- max(0, mean(shuffled_loss_list) - base_loss)
  }

  # 6. Weight Normalization & Shape Extraction
  if (sum(importance_scores) <= 0) {
    weights <- rep(NA_real_, length(mix_name))
    names(weights) <- mix_name
    shape_coefs <- rep(NA_real_, length(spline_vars))
    names(shape_coefs) <- spline_vars
  } else {
    weights <- sqrt(importance_scores) / sum(sqrt(importance_scores))
    shape_coefs <- coefs_no_int[spline_vars]
    shape_coefs[is.na(shape_coefs)] <- 0
  }

  return(list(weights = weights, shapes = shape_coefs))
}
