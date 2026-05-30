#' @title Extract NWQS Partial-Effect Curve Relative to a Reference
#'
#' @description
#' Computes the joint (and optionally per-component) partial-effect change for
#' an \code{nwqs} or \code{nwqs_boot} model along a 0--1 grid, relative to a
#' chosen reference percentile. For \code{nwqs} objects the uncertainty bands
#' are empirical quantiles across repeated holdout iterations (algorithmic
#' variance only); for \code{nwqs_boot} objects they are bootstrap percentile
#' intervals (valid sampling-variance inference).
#'
#' @details
#' For a grid point \eqn{x}, the curve for component \eqn{c} is
#' \deqn{\hat\eta_c(x) = \hat\beta_{nwqs}\,\hat w_c\,(B(x) - B(x_{\mathrm{ref}}))^\top \hat\theta_c}
#' where \eqn{B(\cdot)} is the natural cubic spline basis evaluated on the
#' model's globally aligned knots. The overall curve is the sum over
#' components. At \eqn{x = x_{\mathrm{ref}}} the estimate is identically zero.
#'
#' @param model An object of class \code{"nwqs"} or \code{"nwqs_boot"}.
#' @param grid Numeric vector. Evaluation points in \code{[0, 1]} for
#'   percentile_rank fits. Default \code{seq(0, 1, by = 0.01)} (101 points).
#' @param ref Numeric scalar. Reference point. Default \code{0.5} (median)
#'   for percentile_rank, \code{0} for q_bin.
#' @param include_components Logical. If \code{TRUE} (default), the returned
#'   data.frame includes one row per component plus an \dQuote{Overall} row.
#'   If \code{FALSE}, only the \dQuote{Overall} row is returned.
#' @param label_style One of \code{"auto"}, \code{"P"}, \code{"Q"},
#'   \code{"numeric"}. Reserved for future per-row labelling extensions; the
#'   curve grid is always returned as numeric \code{x}.
#'
#' @return A \code{data.frame} with columns \code{term}, \code{x}, \code{ref},
#'   \code{estimate}, \code{lower}, \code{upper}, \code{transform_type},
#'   \code{inference_type}. When \code{rh == 1} for an \code{nwqs} input,
#'   \code{lower}/\code{upper} are \code{NA}.
#'
#' @importFrom splines ns
#' @importFrom stats quantile
#' @export
extract_nwqs_effect_curve <- function(model,
                                      grid               = NWQS_DEFAULTS$effect_curve_grid,
                                      ref                = NULL,
                                      include_components = TRUE,
                                      label_style        = c("auto", "P", "Q", "numeric")) {
  UseMethod("extract_nwqs_effect_curve")
}

#' @export
extract_nwqs_effect_curve.nwqs <- function(model,
                                           grid               = NWQS_DEFAULTS$effect_curve_grid,
                                           ref                = NULL,
                                           include_components = TRUE,
                                           label_style        = c("auto", "P", "Q", "numeric")) {
  label_style <- match.arg(label_style)
  .validate_pr_points(grid, model$transform_type)
  ref_pt <- .resolve_ref(model, ref)
  comps  <- names(model$final_weights)
  df_spline <- model$df_spline

  eval_seq <- c(ref_pt, grid)
  basis <- splines::ns(eval_seq, df = df_spline,
                       knots = model$spline_knots, Boundary.knots = model$spline_boundary,
                       intercept = FALSE)
  ref_row <- basis[1L, ]
  b_diff_mat <- sweep(basis[-1L, , drop = FALSE], 2, ref_row, FUN = "-")

  rh <- model$rh
  curves_by_term <- vector("list", length(comps) + 1L)
  names(curves_by_term) <- c("Overall", comps)
  for (term in names(curves_by_term)) {
    curves_by_term[[term]] <- matrix(NA_real_, nrow = rh, ncol = length(grid))
  }

  for (i in seq_len(rh)) {
    beta_i <- model$rh_coefs[i, "nwqs"]
    if (!is.finite(beta_i)) next
    overall_val <- numeric(length(grid))
    for (comp in comps) {
      theta_cols <- paste0(comp, "_B", seq_len(df_spline))
      theta_i <- model$rh_shapes[i, theta_cols]
      w_i     <- model$rh_weights[i, comp]
      contrib <- beta_i * w_i * as.vector(b_diff_mat %*% theta_i)
      curves_by_term[[comp]][i, ] <- contrib
      overall_val <- overall_val + contrib
    }
    curves_by_term[["Overall"]][i, ] <- overall_val
  }

  terms_to_keep <- if (include_components) c("Overall", comps) else "Overall"

  out_rows <- lapply(terms_to_keep, function(term) {
    mat <- curves_by_term[[term]]
    est <- colMeans(mat, na.rm = TRUE)
    if (rh > 1) {
      lower <- apply(mat, 2, quantile, 0.025, na.rm = TRUE)
      upper <- apply(mat, 2, quantile, 0.975, na.rm = TRUE)
    } else {
      lower <- rep(NA_real_, length(grid))
      upper <- rep(NA_real_, length(grid))
    }
    data.frame(
      term            = term,
      x               = grid,
      ref             = ref_pt,
      estimate        = est,
      lower           = lower,
      upper           = upper,
      transform_type  = model$transform_type,
      inference_type  = "repeated_holdout",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out_rows)
  rownames(out) <- NULL
  out
}

#' @title Plot the NWQS Partial-Effect Curve Relative to a Reference Percentile
#'
#' @description
#' Renders the joint partial-effect curve from
#' \code{extract_nwqs_effect_curve()} as a ggplot. By default only the overall
#' curve is drawn (median-centered with \code{ref = 0.5}). For \code{nwqs_boot}
#' inputs the ribbon is a valid bootstrap percentile band; for \code{nwqs}
#' inputs with \code{rh > 1} the ribbon is the empirical RH band (algorithmic
#' variance only -- not valid sampling inference).
#'
#' @param model An object of class \code{"nwqs"} or \code{"nwqs_boot"}.
#' @param grid Numeric vector. Default \code{seq(0, 1, by = 0.01)}.
#' @param ref Numeric scalar. Default \code{0.5}.
#' @param include_components Logical. If \code{FALSE} (default), only the
#'   Overall curve is plotted. If \code{TRUE}, component curves are overlaid.
#' @param label_style One of \code{"auto"}, \code{"P"}, \code{"Q"},
#'   \code{"numeric"} (mirrors the extractor signature).
#' @param base_size Integer. Base font size for the ggplot theme. Default 12.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plot_nwqs_effect_curve <- function(model,
                                   grid               = NWQS_DEFAULTS$effect_curve_grid,
                                   ref                = NULL,
                                   include_components = FALSE,
                                   label_style        = c("auto", "P", "Q", "numeric"),
                                   base_size          = 12) {
  label_style <- match.arg(label_style)
  ref_pt    <- .resolve_ref(model, ref)
  ref_label <- .contrast_point_label(
    ref_pt,
    transform_type = model$transform_type,
    label_style    = .label_style_default(model$transform_type)
  )

  df <- extract_nwqs_effect_curve(
    model,
    grid               = grid,
    ref                = ref_pt,
    include_components = include_components,
    label_style        = label_style
  )

  x_label <- if (identical(model$transform_type, "percentile_rank")) {
    "Joint exposure percentile rank"
  } else {
    "Quantile index"
  }
  y_label <- sprintf("Partial effect change relative to %s", ref_label)

  has_ci <- any(!is.na(df$lower))

  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = .data$x, y = .data$estimate,
                                    color = .data$term, fill = .data$term,
                                    group = .data$term)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "#95A5A6", linewidth = 0.5) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::labs(title = "NWQS Partial-Effect Curve",
                  x = x_label, y = y_label)

  if (has_ci) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower,
                                               ymax = .data$upper),
                                  alpha = 0.15, color = NA)
  }
  p <- p + ggplot2::geom_line(linewidth = 1.1)

  if (!include_components) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  p
}

#' @export
extract_nwqs_effect_curve.nwqs_boot <- function(model,
                                                grid               = NWQS_DEFAULTS$effect_curve_grid,
                                                ref                = NULL,
                                                include_components = TRUE,
                                                label_style        = c("auto", "P", "Q", "numeric")) {
  label_style <- match.arg(label_style)
  .validate_pr_points(grid, model$transform_type)
  if (is.null(model$rh_shapes_boot)) {
    stop("nwqs_boot object lacks `rh_shapes_boot`; refit with current package version.",
         call. = FALSE)
  }
  ref_pt <- .resolve_ref(model, ref)
  comps  <- names(model$final_weights)
  df_spline <- model$df_spline
  n_boot <- nrow(model$rh_shapes_boot)

  eval_seq <- c(ref_pt, grid)
  basis <- splines::ns(eval_seq, df = df_spline,
                       knots = model$spline_knots, Boundary.knots = model$spline_boundary,
                       intercept = FALSE)
  ref_row <- basis[1L, ]
  b_diff_mat <- sweep(basis[-1L, , drop = FALSE], 2, ref_row, FUN = "-")

  curves_by_term <- vector("list", length(comps) + 1L)
  names(curves_by_term) <- c("Overall", comps)
  for (term in names(curves_by_term)) {
    curves_by_term[[term]] <- matrix(NA_real_, nrow = n_boot, ncol = length(grid))
  }

  for (i in seq_len(n_boot)) {
    beta_i <- model$rh_coefs_boot[i, "nwqs"]
    if (!is.finite(beta_i)) next
    overall_val <- numeric(length(grid))
    for (comp in comps) {
      theta_cols <- paste0(comp, "_B", seq_len(df_spline))
      theta_i <- model$rh_shapes_boot[i, theta_cols]
      w_i     <- model$rh_weights_boot[i, comp]
      contrib <- beta_i * w_i * as.vector(b_diff_mat %*% theta_i)
      curves_by_term[[comp]][i, ] <- contrib
      overall_val <- overall_val + contrib
    }
    curves_by_term[["Overall"]][i, ] <- overall_val
  }

  terms_to_keep <- if (include_components) c("Overall", comps) else "Overall"

  out_rows <- lapply(terms_to_keep, function(term) {
    mat <- curves_by_term[[term]]
    est   <- colMeans(mat, na.rm = TRUE)
    lower <- apply(mat, 2, quantile, 0.025, na.rm = TRUE)
    upper <- apply(mat, 2, quantile, 0.975, na.rm = TRUE)
    data.frame(
      term            = term,
      x               = grid,
      ref             = ref_pt,
      estimate        = est,
      lower           = lower,
      upper           = upper,
      transform_type  = model$transform_type,
      inference_type  = "bootstrap",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out_rows)
  rownames(out) <- NULL
  out
}
