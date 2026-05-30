#' @title Control Options for NWQS Fits
#'
#' @description
#' Bundles "advanced" knobs that do not appear on the main \code{nwqs()}
#' signature so the signature stays short as new options accrue. Pass an
#' \code{nwqs_control} object to \code{nwqs()} via its \code{control}
#' argument.
#'
#' @details
#' The parameters managed here are intended for sensitivity analyses
#' rather than everyday use:
#' \itemize{
#'   \item \code{custom_knots}: override the internal knot vector that
#'     \code{build_spline_basis_knots()} would otherwise derive from the
#'     transform's evaluation grid. Useful if a paper specifies particular
#'     percentile / bin positions.
#'   \item \code{custom_boundary}: override the boundary knots. Must be a
#'     length-2 numeric vector.
#'   \item \code{zero_weight_action}: how \code{permutation_scorer()}
#'     should handle the degenerate case where every component's OOB loss
#'     delta is non-positive. Default \code{"na"} returns NA weights for
#'     that iteration (matching 0.1.x). \code{"uniform"} fills with
#'     \code{1 / n_mix} so the iteration still contributes; use with
#'     caution since this hides a fitting failure.
#' }
#'
#' Parameters that already live on the main \code{nwqs()} signature
#' (\code{min_shape_sd}, \code{ties}, \code{transform_type}, \code{q})
#' stay on the signature in 0.2.0 for backward compatibility; only "new"
#' soft parameters are routed through \code{nwqs_control}.
#'
#' @param custom_knots Numeric vector or \code{NULL}. If non-NULL,
#'   overrides the internal spline knots.
#' @param custom_boundary Numeric vector of length 2 or \code{NULL}. If
#'   non-NULL, overrides the boundary knots.
#' @param zero_weight_action Character. One of \code{"na"} (default) or
#'   \code{"uniform"}.
#'
#' @return An object of class \code{c("nwqs_control", "list")}.
#'
#' @export
nwqs_control <- function(custom_knots = NULL,
                         custom_boundary = NULL,
                         zero_weight_action = NWQS_DEFAULTS$zero_weight_action) {
  if (!is.null(custom_knots)) {
    if (!is.numeric(custom_knots)) {
      stop("`custom_knots` must be numeric or NULL.")
    }
  }
  if (!is.null(custom_boundary)) {
    if (!is.numeric(custom_boundary) || length(custom_boundary) != 2) {
      stop("`custom_boundary` must be a numeric vector of length 2 or NULL.")
    }
  }
  if (!zero_weight_action %in% c("na", "uniform")) {
    stop("`zero_weight_action` must be one of 'na' or 'uniform'.")
  }

  structure(
    list(
      custom_knots       = custom_knots,
      custom_boundary    = custom_boundary,
      zero_weight_action = zero_weight_action
    ),
    class = c("nwqs_control", "list")
  )
}
