#' @title Fast GLM Permutation Scoring Engine
#'
#' @description
#' Core computational engine for a single internal iteration within the NWQS
#' framework. Fits an unpenalized GLM on in-bag bootstrap samples, evaluates
#' prediction loss on out-of-bag (OOB) samples, and derives relative
#' component importance via random permutation.
#'
#' @details
#' Workflow:
#' \enumerate{
#'   \item Adaptive sampling and GLM fitting on in-bag data via
#'     \code{stats::glm.fit}.
#'   \item Baseline loss evaluation on OOB samples (MSE for Gaussian,
#'     deviance for binomial / Poisson / quasi-Poisson).
#'   \item Grouped permutation importance: for each mixture component, all
#'     corresponding spline basis columns are jointly shuffled and the OOB
#'     loss change is measured.
#' }
#'
#' Weight derivation formula:
#' \deqn{w_i = \frac{\sqrt{\max(0, \Delta Loss_i)}}{\sum \sqrt{\max(0, \Delta Loss)}}}
#'
#' @param x Numeric matrix. Design matrix containing spline basis columns and
#'   adjustment covariates.
#' @param y Numeric vector. Outcome variable.
#' @param mix_name Character vector. Names of original mixture components.
#' @param spline_vars Character vector. Column names in \code{x} that belong to
#'   spline basis functions.
#' @param family List. GLM family object (e.g., \code{gaussian()},
#'   \code{binomial()}, \code{poisson()}, \code{quasipoisson()}).
#' @param n_permutation Integer. Number of OOB permutations for stabilizing
#'   importance scores. Default is 10.
#' @param ... Additional compatibility parameters (accepted but ignored, so
#'   callers that pass extras like \code{strata_id} from older revisions still
#'   work).
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{weights}: Normalized relative importance weights for each
#'     mixture component.
#'   \item \code{shapes}: Spline basis coefficients estimated on in-bag data.
#'     Coefficients dropped due to collinearity are set to 0.
#' }
#'
#' @importFrom stats coef predict glm.fit as.formula
#' @export
permutation_scorer <- function(x, y, mix_name, spline_vars, family,
                               n_permutation = 10, ...) {
  n_obs <- nrow(x)
  fam_name <- family$family
  linkinv <- family$linkinv

  idx <- sample(seq_len(n_obs), size = n_obs, replace = TRUE)
  oob_idx <- setdiff(seq_len(n_obs), idx)
  if (length(oob_idx) == 0) return(NULL)

  x_train <- x[idx, , drop = FALSE]
  y_train <- y[idx]
  x_oob <- x[oob_idx, , drop = FALSE]
  y_oob <- y[oob_idx]

  int_col <- match("(Intercept)", colnames(x_train))
  if (!is.na(int_col)) {
    x_train_net <- x_train[, -int_col, drop = FALSE]
    x_oob_net <- x_oob[, -int_col, drop = FALSE]
  } else {
    x_train_net <- x_train
    x_oob_net <- x_oob
  }

  calc_loss <- function(y_true, mu_pred) {
    if (fam_name == "gaussian") return(mean((y_true - mu_pred)^2))
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

  x_train_glm <- cbind(Intercept = 1, as.matrix(x_train_net))
  x_oob_glm <- cbind(Intercept = 1, as.matrix(x_oob_net))

  fit <- stats::glm.fit(x = x_train_glm, y = y_train, family = family)
  coef_all <- fit$coefficients
  coef_all[is.na(coef_all)] <- 0

  intercept_val <- unname(coef_all[1])
  coefs_no_int <- coef_all[-1]

  eta_oob <- as.numeric(x_oob_glm %*% coef_all)
  mu_oob <- linkinv(eta_oob)
  base_loss <- calc_loss(y_oob, mu_oob)

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

      eta_shuffled <- intercept_val + as.numeric(as.matrix(x_oob_shuffled) %*% coefs_no_int)
      mu_shuffled <- linkinv(eta_shuffled)
      shuffled_loss_list[k] <- calc_loss(y_oob, mu_shuffled)
    }

    x_oob_shuffled[, target_cols] <- x_oob_net[, target_cols, drop = FALSE]
    importance_scores[var] <- max(0, mean(shuffled_loss_list) - base_loss)
  }

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
