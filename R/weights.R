#' @title Fast GLM and Conditional Logistic Regression Permutation Scoring Engine
#'
#' @description
#' Core computational engine for a single internal iteration within the NWQS
#' framework. Fits an unpenalized GLM or conditional logistic regression model
#' on in-bag bootstrap samples, evaluates prediction loss on out-of-bag (OOB)
#' samples, and derives relative component importance via random permutation.
#'
#' @details
#' Workflow:
#' \enumerate{
#'   \item Adaptive sampling and model fitting on in-bag data.
#'   \item Baseline loss evaluation on OOB samples (MSE or deviance). For
#'     clogit, the engine bypasses formula parsing and directly calls
#'     \code{survival::coxph.fit} with \code{iter.max = 0} for millisecond-level
#'     partial log-likelihood evaluation.
#'   \item Grouped permutation importance: for each mixture component, all
#'     corresponding spline basis columns are jointly shuffled and the OOB loss
#'     change is measured.
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
#' @param family List. GLM family object (or a pseudo-family list with
#'   \code{$family = "clogit"}).
#' @param n_permutation Integer. Number of OOB permutations for stabilizing
#'   importance scores. Default is 100.
#' @param strata_id Vector or \code{NULL}. Stratum/matching group IDs required
#'   for conditional logistic regression.
#' @param ... Additional compatibility parameters.
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
#' @importFrom survival Surv clogit coxph.fit coxph.control
#' @export
permutation_scorer <- function(x, y, mix_name, spline_vars, family,
                               n_permutation = 100, strata_id = NULL, ...) {
  n_obs <- nrow(x)
  fam_name <- family$family

  if (fam_name != "clogit") {
    linkinv <- family$linkinv
  }

  if (fam_name == "clogit") {
    if (is.null(strata_id)) stop("strata_id must be provided for family = 'clogit'")

    unique_strata <- unique(strata_id)
    n_strata <- length(unique_strata)

    sampled_strata <- sample(unique_strata, size = n_strata, replace = TRUE)

    idx_list <- lapply(seq_along(sampled_strata), function(i) {
      orig_idx <- which(strata_id == sampled_strata[i])
      data.frame(
        orig_row = orig_idx,
        new_strata = paste0(sampled_strata[i], "_boot_", i),
        stringsAsFactors = FALSE
      )
    })
    train_map <- do.call(rbind, idx_list)

    idx <- train_map$orig_row
    strata_train <- train_map$new_strata

    oob_strata <- setdiff(unique_strata, sampled_strata)
    if (length(oob_strata) == 0) return(NULL)

    oob_idx <- which(strata_id %in% oob_strata)
    strata_oob <- strata_id[oob_idx]
  } else {
    idx <- sample(seq_len(n_obs), size = n_obs, replace = TRUE)
    oob_idx <- setdiff(seq_len(n_obs), idx)
    if (length(oob_idx) == 0) return(NULL)
  }

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

  if (fam_name == "clogit") {
    fast_clogit_loss <- function(x_new, y_new, strata_new, coefs_init) {
      y_surv <- survival::Surv(rep(1, length(y_new)), y_new)
      strata_int <- as.integer(as.factor(strata_new))

      fit_eval <- tryCatch({
        survival::coxph.fit(
          x = as.matrix(x_new),
          y = y_surv,
          strata = strata_int,
          init = coefs_init,
          control = survival::coxph.control(iter.max = 0),
          method = "exact",
          rownames = NULL
        )
      }, error = function(e) NULL)

      if (is.null(fit_eval)) return(NA_real_)
      return(-2 * fit_eval$loglik[1])
    }
  } else {
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
  }

  if (fam_name == "clogit") {
    df_train <- data.frame(y_event = y_train, strata_id = strata_train)
    df_train <- cbind(df_train, as.data.frame(x_train_net))

    x_cols <- colnames(x_train_net)
    form_str <- paste0("y_event ~ ", paste(sprintf("`%s`", x_cols), collapse = " + "), " + strata(strata_id)")

    fit <- tryCatch({
      survival::clogit(as.formula(form_str), data = df_train)
    }, error = function(e) NULL)

    if (is.null(fit)) return(NULL)

    coef_all <- coef(fit)
    coef_all[is.na(coef_all)] <- 0
    intercept_val <- 0
    coefs_no_int <- coef_all

    base_loss <- fast_clogit_loss(x_oob_net, y_oob, strata_oob, coefs_no_int)
    if (is.na(base_loss)) return(NULL)
  } else {
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
  }

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

      if (fam_name == "clogit") {
        loss_val <- fast_clogit_loss(x_oob_shuffled, y_oob, strata_oob, coefs_no_int)
        shuffled_loss_list[k] <- if (is.na(loss_val)) base_loss else loss_val
      } else {
        eta_shuffled <- intercept_val + as.numeric(as.matrix(x_oob_shuffled) %*% coefs_no_int)
        mu_shuffled <- linkinv(eta_shuffled)
        shuffled_loss_list[k] <- calc_loss(y_oob, mu_shuffled)
      }
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
