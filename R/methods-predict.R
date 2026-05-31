# predict / vcov / confint methods for nwqs and nwqs_boot objects.
#
# Design notes:
#   - Per the framework, the "nwqs" index is a single-degree-of-freedom
#     latent score; once it has been constructed for newdata using the
#     saved spline knots and weights, predictions on the link and
#     response scales reduce to a standard linear-predictor evaluation.
#   - For rh == 1 the inner GLM object is stashed on the fit so the
#     standard predict.glm / vcov.glm / confint.glm machinery can be
#     reused directly. For rh > 1 we fall back to averaged structures
#     and emit the "rh > 1 is not valid inference" warning at every
#     point that vends standard errors.

.nwqs_transform_newdata <- function(object, newdata) {
  mix_name <- names(object$train_components_sorted)
  if (!all(mix_name %in% colnames(newdata))) {
    missing <- setdiff(mix_name, colnames(newdata))
    stop("newdata is missing mixture component(s): ",
         paste(missing, collapse = ", "))
  }
  newdata_trans <- newdata
  if (object$transform_type == "percentile_rank") {
    for (comp in mix_name) {
      newdata_trans[[comp]] <- apply_percentile_rank(
        newdata[[comp]],
        object$train_components_sorted[[comp]],
        ties = if (!is.null(object$ties)) object$ties else NWQS_DEFAULTS$ties
      )
    }
  } else {
    for (comp in mix_name) {
      train_vec <- object$train_components_sorted[[comp]]
      breaks <- unique(
        stats::quantile(train_vec, probs = seq(0, 1, by = 1 / object$q), na.rm = TRUE)
      )
      if (length(breaks) == 1) {
        breaks <- c(-Inf, breaks)
      } else {
        breaks[1] <- -Inf
        breaks[length(breaks)] <- Inf
      }
      newdata_trans[[comp]] <- as.numeric(
        cut(newdata[[comp]], breaks = breaks, labels = FALSE, include.lowest = TRUE)
      ) - 1
    }
  }
  newdata_trans
}

.nwqs_compute_index <- function(object, newdata_trans) {
  mix_name <- names(object$train_components_sorted)
  spline_basis <- wqs_nonlinear_expand(
    data = newdata_trans, mix_name = mix_name,
    df_spline = object$df_spline,
    knots = object$spline_knots, boundary = object$spline_boundary
  )
  combined <- numeric(ncol(spline_basis))
  names(combined) <- colnames(spline_basis)
  for (comp in mix_name) {
    comp_cols <- paste0(comp, "_B", seq_len(object$df_spline))
    combined[comp_cols] <- object$mean_shapes[comp_cols] * object$final_weights[comp]
  }
  as.numeric(spline_basis %*% combined)
}

.nwqs_inv_link <- function(family) {
  switch(
    family,
    gaussian     = function(x) x,
    binomial     = stats::plogis,
    poisson      = exp,
    quasipoisson = exp,
    negbin       = exp,
    stop("Unsupported family for response-scale prediction: ", family)
  )
}


#' @title Predict Method for NWQS Models
#'
#' @description
#' Computes the NWQS latent index, the linear predictor, or the response-scale
#' prediction for new data. The fit-sample empirical distribution (saved as
#' \code{train_components_sorted}) and the globally aligned spline knots are
#' reused so that newdata is interpreted on the same scale used during fitting.
#'
#' @param object An object of class \code{"nwqs"}.
#' @param newdata Optional \code{data.frame}. If \code{NULL}, the training
#'   data stored on \code{object$data} is used.
#' @param type Character. One of \code{"response"} (default),
#'   \code{"link"}, or \code{"nwqs_index"}.
#'   \itemize{
#'     \item \code{"nwqs_index"} returns the 1-DoF latent NWQS score for
#'       each row of newdata.
#'     \item \code{"link"} returns the full linear predictor
#'       \eqn{X\beta} including the intercept, the \code{nwqs} term, and any
#'       covariates that were in the original fit.
#'     \item \code{"response"} returns the inverse-link transform of the
#'       linear predictor (identity for Gaussian, \code{plogis} for binomial,
#'       \code{exp} for Poisson / quasi-Poisson / negbin).
#'   }
#' @param ... Currently unused.
#'
#' @return A numeric vector with one element per row of newdata.
#'
#' @export
predict.nwqs <- function(object, newdata = NULL,
                         type = c("response", "link", "nwqs_index"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) newdata <- object$data

  newdata_trans <- .nwqs_transform_newdata(object, newdata)
  nwqs_index <- .nwqs_compute_index(object, newdata_trans)

  if (type == "nwqs_index") return(nwqs_index)

  newdata_with_nwqs <- newdata
  newdata_with_nwqs$nwqs <- nwqs_index

  if (is.null(object$formula)) {
    stop("Cannot compute link/response prediction: object$formula is missing.")
  }

  mm <- stats::model.matrix(
    stats::delete.response(stats::terms(object$formula)),
    data = newdata_with_nwqs
  )

  beta <- object$mean_coefs[colnames(mm)]
  if (any(is.na(beta))) {
    missing_coefs <- colnames(mm)[is.na(beta)]
    stop(sprintf(
      "Missing coefficient(s) for %s on the fit object.",
      paste(missing_coefs, collapse = ", ")
    ))
  }

  eta <- as.numeric(mm %*% beta)
  if (type == "link") return(eta)

  inv_link <- .nwqs_inv_link(object$family)
  inv_link(eta)
}


#' @title Predict Method for Bootstrap NWQS Models
#'
#' @description
#' Computes the NWQS latent index, the linear predictor, or the
#' response-scale prediction for new data using the bootstrap-averaged
#' weights, shapes, and regression coefficients.
#'
#' @details
#' This is a point-estimate predictor. The bootstrap-iteration-level
#' predictive distribution (which would give per-row confidence intervals
#' on the prediction) is not yet implemented; pass \code{keep_fits = TRUE}
#' to \code{nwqs_boot()} and call \code{predict()} on each entry of
#' \code{boot_fits} to roll your own if you need them.
#'
#' @param object An object of class \code{"nwqs_boot"}.
#' @param newdata Optional \code{data.frame}. If \code{NULL}, the training
#'   data stored on \code{object$data} is used.
#' @param type Character. One of \code{"response"} (default),
#'   \code{"link"}, or \code{"nwqs_index"}.
#' @param ... Currently unused.
#'
#' @return A numeric vector with one element per row of newdata.
#'
#' @export
predict.nwqs_boot <- function(object, newdata = NULL,
                              type = c("response", "link", "nwqs_index"), ...) {
  type <- match.arg(type)
  stub <- list(
    transform_type          = object$transform_type,
    train_components_sorted = object$train_components_sorted,
    q                       = object$q,
    df_spline               = object$df_spline,
    spline_knots            = object$spline_knots,
    spline_boundary         = object$spline_boundary,
    mean_shapes             = object$mean_shapes,
    final_weights           = object$final_weights,
    mean_coefs              = object$mean_coefs,
    formula                 = object$formula,
    family                  = object$family,
    ties                    = if (!is.null(object$ties)) object$ties else NWQS_DEFAULTS$ties,
    data                    = object$data
  )
  class(stub) <- c("nwqs", "list")
  predict.nwqs(stub, newdata = newdata, type = type)
}


.rh_inference_warning <- function() {
  warning(
    "Standard errors / CIs from a fit with rh > 1 reflect algorithmic ",
    "(data-splitting) variance only, NOT sampling variance. Use ",
    "nwqs_boot() for valid sampling-variance inference."
  )
}


#' @title Covariance Matrix for NWQS Fits
#'
#' @description
#' For \code{rh = 1}, returns the covariance matrix of the inner GLM, which
#' is the standard sampling-variance covariance. For \code{rh > 1},
#' returns the empirical covariance of the regression coefficients across
#' RH iterations and emits a warning that this is algorithmic variance
#' only.
#'
#' @param object An object of class \code{"nwqs"}.
#' @param ... Currently unused.
#'
#' @return A symmetric numeric matrix with one row/column per regression
#'   coefficient (\code{"(Intercept)"}, \code{"nwqs"}, covariates).
#'
#' @export
vcov.nwqs <- function(object, ...) {
  if (isTRUE(object$rh == 1) && !is.null(object$model_obj)) {
    return(stats::vcov(object$model_obj))
  }
  .rh_inference_warning()
  v <- stats::cov(object$rh_coefs, use = "complete.obs")
  v
}


#' @title Covariance Matrix for Bootstrap NWQS Fits
#'
#' @description
#' Returns the empirical covariance of the per-iteration regression
#' coefficients across the bootstrap replicates. This IS a valid
#' sampling-variance covariance.
#'
#' @param object An object of class \code{"nwqs_boot"}.
#' @param ... Currently unused.
#'
#' @return A symmetric numeric matrix.
#'
#' @export
vcov.nwqs_boot <- function(object, ...) {
  if (is.null(object$rh_coefs_boot)) {
    stop("object$rh_coefs_boot is missing; refit with NWQS >= 0.2.0.")
  }
  stats::cov(object$rh_coefs_boot, use = "complete.obs")
}


#' @title Confidence Intervals for NWQS Fits
#'
#' @description
#' For \code{rh = 1}, returns \code{stats::confint()} on the inner GLM.
#' For \code{rh > 1}, returns the empirical \code{(1 - level)/2} and
#' \code{1 - (1 - level)/2} quantiles of the per-iteration regression
#' coefficients and emits the rh > 1 warning.
#'
#' @param object An object of class \code{"nwqs"}.
#' @param parm Optional. Which parameters to report; defaults to all.
#' @param level Numeric. Confidence level. Default is 0.95.
#' @param ... Forwarded to \code{stats::confint()} when \code{rh = 1}.
#'
#' @return A two-column matrix with rownames matching the regression
#'   coefficients and columns \code{"2.5 \%"} / \code{"97.5 \%"} at the
#'   default level.
#'
#' @export
confint.nwqs <- function(object, parm, level = 0.95, ...) {
  if (isTRUE(object$rh == 1) && !is.null(object$model_obj)) {
    if (missing(parm)) {
      return(stats::confint(object$model_obj, level = level, ...))
    }
    return(stats::confint(object$model_obj, parm = parm, level = level, ...))
  }
  .rh_inference_warning()
  alpha <- 1 - level
  probs <- c(alpha / 2, 1 - alpha / 2)
  ci <- apply(object$rh_coefs, 2, stats::quantile, probs = probs, na.rm = TRUE)
  ci <- t(ci)
  colnames(ci) <- sprintf("%g %%", probs * 100)
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}


#' @title Confidence Intervals for Bootstrap NWQS Fits
#'
#' @description
#' Returns the empirical \code{(1 - level)/2} and
#' \code{1 - (1 - level)/2} quantiles of the per-iteration regression
#' coefficients across bootstrap replicates. These are valid percentile
#' bootstrap CIs.
#'
#' @param object An object of class \code{"nwqs_boot"}.
#' @param parm Optional. Which parameters to report; defaults to all.
#' @param level Numeric. Confidence level. Default uses the fit's stored
#'   \code{conf_level}.
#' @param ... Currently unused.
#'
#' @return A two-column matrix with rownames matching the regression
#'   coefficients.
#'
#' @export
confint.nwqs_boot <- function(object, parm, level = NULL, ...) {
  if (is.null(object$rh_coefs_boot)) {
    stop("object$rh_coefs_boot is missing; refit with NWQS >= 0.2.0.")
  }
  if (is.null(level)) {
    level <- if (!is.null(object$conf_level)) object$conf_level else 0.95
  }
  alpha <- 1 - level
  probs <- c(alpha / 2, 1 - alpha / 2)
  ci <- apply(object$rh_coefs_boot, 2, stats::quantile, probs = probs, na.rm = TRUE)
  ci <- t(ci)
  colnames(ci) <- sprintf("%g %%", probs * 100)
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}
