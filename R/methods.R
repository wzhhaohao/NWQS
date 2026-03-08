# ==============================================================================
# Internal Helpers
# ==============================================================================

#' @keywords internal
.format_pval <- function(p) {
  sapply(p, function(pv) {
    if (is.na(pv)) {
      return("NA")
    }
    if (pv < 0.001) {
      return("<0.001")
    }
    if (pv < 0.01) {
      return(sprintf("%.3f", pv))
    }
    if (pv < 0.05) {
      return(sprintf("%.3f", pv))
    }
    return(sprintf("%.3f", pv))
  })
}

#' @keywords internal
.print_coef_table <- function(coef_table, digits = 4) {
  if (!any(grepl("Pr\\(", colnames(coef_table)))) {
    z_stat <- coef_table$Estimate / coef_table$`Std. Error`
    p_values <- 2 * stats::pnorm(-abs(z_stat))
    coef_table$`Pr(>|z|)` <- p_values
  }

  p_col <- grep("Pr\\(", colnames(coef_table))[1]
  p_raw <- coef_table[[p_col]]

  out <- as.data.frame(coef_table)

  num_cols <- which(sapply(out, is.numeric))
  for (j in num_cols) out[[j]] <- round(out[[j]], digits)

  out[[p_col]] <- .format_pval(p_raw)

  out$` ` <- ifelse(p_raw < 0.001, "***",
    ifelse(p_raw < 0.01, "**",
      ifelse(p_raw < 0.05, "*",
        ifelse(p_raw < 0.1, ".", " ")
      )
    )
  )

  print(out)
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  invisible(out)
}


# ==============================================================================
# S3 Methods: nwqs
# ==============================================================================

#' Plot Diagnostics for Non-linear Weighted Quantile Sum (NWQS) Models
#'
#' @description
#' Generates publication-ready diagnostic plots for \code{nwqs} objects. This function
#' can produce component weight distributions and non-linear dose-response trajectories
#' (partial effects or predicted values) for the overall mixture and individual components.
#'
#' @details
#' The dose-response curves (when \code{type = "curves"} or \code{"both"}) map the estimated
#' non-linear spline functions. When \code{y_scale = "partial"}, the curves represent the
#' isolated effect of increasing exposure while holding others constant (zero-centered
#' log-OR, log-RR, or \eqn{\Delta Y}). When \code{y_scale = "predicted"}, the curves show
#' the absolute predicted scale (e.g., probabilities for binomial models).
#'
#' When \code{plot_ci = TRUE} and \code{rh > 1}, 95\% empirical confidence ribbons
#' are drawn from the Repeated Holdout iterations. \strong{Important}: These ribbons
#' reflect data-splitting (algorithmic) variance only and are NOT valid for statistical
#' inference. A warning caption is automatically appended to the plot. Use
#' \code{\link{nwqs_boot}} with \code{keep_fits = TRUE} for valid bootstrap confidence
#' intervals.
#'
#' @param x An object of class \code{"nwqs"} resulting from a call to \code{\link{nwqs}}.
#' @param type Character. Plot type: \code{"both"} (default, combines weights and curves
#'   via \pkg{patchwork}), \code{"curves"} (dose-response trajectories), or
#'   \code{"weights"} (bar plot of component weights).
#' @param y_scale Character. Y-axis scale for curves: \code{"partial"} (default, relative
#'   effect changes) or \code{"predicted"} (absolute predicted values).
#' @param components Character vector. Specific mixture components to display. If
#'   \code{NULL} (default), all components are shown.
#' @param overlay Logical. If \code{TRUE} (default), all component curves are overlaid
#'   on a single plot. If \code{FALSE}, curves are faceted by component.
#' @param plot_ci Logical. Whether to plot 95\% empirical confidence intervals from
#'   Repeated Holdout iterations. Default is \code{FALSE}. Requires \code{rh > 1}.
#' @param base_size Integer. Base font size for the \pkg{ggplot2} theme. Default is 12.
#' @param palette Character. Color palette: \code{"default"}, \code{"palette2"},
#'   \code{"palette3"}, or \code{"palette4"}.
#' @param colorblind_friendly Logical. If \code{TRUE}, forces a colorblind-safe palette.
#'   Default is \code{FALSE}.
#' @param top_n Integer or NULL. Number of top components (by weight) to display.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return A \code{ggplot} object (for \code{"weights"} or \code{"curves"}) or a
#'   \code{patchwork} composite (for \code{"both"}).
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_errorbar geom_line geom_ribbon
#'   geom_hline scale_fill_manual scale_color_manual scale_y_continuous
#'   scale_x_continuous coord_flip facet_wrap theme_minimal labs theme expansion
#'   element_blank element_line element_text unit margin
#' @importFrom splines ns
#' @importFrom stats plogis quantile
#'
#' @seealso \code{\link{plot.nwqs_boot}} for bootstrap-based plots with valid CIs.
#' @export
#' @method plot nwqs
plot.nwqs <- function(x, type = c("both", "curves", "weights"),
                      y_scale = c("partial", "predicted"),
                      components = NULL, overlay = TRUE, plot_ci = FALSE,
                      base_size = 12, palette = "default",
                      colorblind_friendly = FALSE, top_n = NULL, ...) {
  type <- match.arg(type)
  y_scale <- match.arg(y_scale)
  if (!inherits(x, "nwqs")) stop("Object must be of class 'nwqs'")

  ci_available <- plot_ci && !is.null(x$rh_weights) && x$rh > 1
  rh_ci_caption <- "Note: Shaded ribbons / error bars reflect data-splitting variance only (NOT valid for inference). Use nwqs_boot() for valid CIs."
  is_combined <- (type == "both")

  .get_palette <- function(n_colors, palette = "default") {
    default <- c(
      "#4A90C8", "#D92828", "#6EC44A", "#8B6FB8", "#00B4D8",
      "#006B3C", "#A8D8EA", "#F4B6B6", "#5BA3D0", "#E03030", "#7AD450", "#9B7FC0"
    )
    palette2 <- c(
      "#9bbf8a", "#82afda", "#f79059", "#e7dbd3", "#c2bdde",
      "#8dcec8", "#add3e2", "#3480b8", "#ffbe7a", "#fa8878", "#c82423", "#6b5b95"
    )
    palette3 <- c(
      "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
      "#E6AB02", "#A6761D", "#666666", "#4DAF4A", "#377EB8", "#FFFF33", "#984EA3"
    )
    palette4 <- c(
      "#004C90", "#E60000", "#009E73", "#E69F00", "#56B4E9",
      "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000", "#1F78B4", "#B2DF8A"
    )
    palettes <- list(default = default, palette2 = palette2, palette3 = palette3, palette4 = palette4)
    sel <- if (tolower(palette) %in% names(palettes)) {
      palettes[[tolower(palette)]]
    } else {
      warning(paste0("Palette '", palette, "' not found. Using 'default'."))
      default
    }
    if (n_colors > length(sel)) {
      rep(sel, ceiling(n_colors / length(sel)))[seq_len(n_colors)]
    } else {
      sel[seq_len(n_colors)]
    }
  }

  # -- Weights panel -----------------------------------------------------------
  if (type %in% c("weights", "both")) {
    w_names <- names(x$final_weights)
    if (!is.null(components)) w_names <- intersect(w_names, components)

    w_df <- data.frame(Component = w_names, Weight = x$final_weights[w_names])

    if (ci_available) {
      rh_w <- x$rh_weights[, w_names, drop = FALSE]
      w_df$Lower <- apply(rh_w, 2, quantile, 0.025, names = FALSE, na.rm = TRUE)
      w_df$Upper <- apply(rh_w, 2, quantile, 0.975, names = FALSE, na.rm = TRUE)
    } else {
      w_df$Lower <- w_df$Weight
      w_df$Upper <- w_df$Weight
    }

    w_df$Component <- factor(w_df$Component, levels = w_df$Component[order(w_df$Weight)])
    w_colors <- .get_palette(length(w_names), palette)

    p_w <- ggplot2::ggplot(w_df, ggplot2::aes(x = Component, y = Weight, fill = Component)) +
      ggplot2::geom_col(alpha = 0.9, width = 0.7) +
      ggplot2::scale_fill_manual(values = w_colors, guide = "none") +
      ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.02))) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::labs(title = "Component Weights", x = "", y = "Weight") +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(color = "#2C3E50", linewidth = 0.6),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = base_size + 1),
        axis.text.y = ggplot2::element_text(face = "bold")
      )

    if (ci_available) {
      p_w <- p_w + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = Lower, ymax = Upper),
        width = 0.3, color = "#2C3E50", linewidth = 0.6
      )
    }

    if (ci_available && !is_combined) {
      p_w <- p_w + ggplot2::labs(caption = rh_ci_caption) +
        ggplot2::theme(plot.caption = ggplot2::element_text(
          face = "italic", color = "red3", size = base_size - 3, hjust = 0
        ))
    }
  }

  # -- Curves panel ------------------------------------------------------------
  if (type %in% c("curves", "both")) {
    q_level <- if (!is.null(x$q)) x$q else 4
    x_seq <- seq(0, q_level - 1, length.out = 100)
    shape_names <- names(x$mean_shapes)

    pattern <- "^(.+)_B(\\d+)$"
    parsed_names <- data.frame(
      full_name = shape_names,
      component = sub(pattern, "\\1", shape_names),
      basis_idx = as.numeric(sub(pattern, "\\2", shape_names)),
      stringsAsFactors = FALSE
    )
    df_spline <- max(parsed_names$basis_idx)
    basis_mat <- splines::ns(x_seq,
      df = df_spline, knots = x$spline_knots,
      Boundary.knots = x$spline_boundary, intercept = FALSE
    )

    unique_comps <- unique(parsed_names$component)
    if (!is.null(components)) unique_comps <- intersect(unique_comps, components)

    if (!is.null(top_n) && is.numeric(top_n) && top_n > 0) {
      top_components <- names(sort(x$final_weights, decreasing = TRUE))[seq_len(min(top_n, length(x$final_weights)))]
      unique_comps <- intersect(unique_comps, top_components)
      if (length(unique_comps) == 0) stop("No components match top_n filter.")
    }

    curve_colors <- .get_palette(length(unique_comps), palette)
    plot_data_list <- list()

    for (comp in unique_comps) {
      comp_cols <- parsed_names$full_name[parsed_names$component == comp]

      if (ci_available) {
        beta_mat <- x$rh_shapes[, comp_cols, drop = FALSE]
        scaling_factor <- x$rh_coefs[, "wqs_score"] * x$rh_weights[, comp]
        beta_mat <- sweep(beta_mat, 1, scaling_factor, "*")
        y_pred_mat <- as.matrix(basis_mat) %*% t(beta_mat)
        if (y_scale == "predicted") {
          y_pred_mat <- sweep(y_pred_mat, 2, x$rh_coefs[, "(Intercept)"], "+")
        }
      } else {
        beta_mat <- matrix(x$mean_shapes[comp_cols], nrow = 1)
        scaling_factor <- x$mean_coefs["wqs_score"] * x$final_weights[comp]
        beta_mat <- beta_mat * scaling_factor
        y_pred_mat <- as.matrix(basis_mat) %*% t(beta_mat)
        if (y_scale == "predicted") {
          y_pred_mat <- y_pred_mat + x$mean_coefs["(Intercept)"]
        }
      }

      if (y_scale == "predicted") {
        if (x$family == "binomial") {
          y_pred_mat <- stats::plogis(y_pred_mat)
        } else if (x$family %in% c("poisson", "quasipoisson")) y_pred_mat <- exp(y_pred_mat)
      }

      y_stats <- t(apply(y_pred_mat, 1, function(v) {
        c(
          mean = mean(v, na.rm = TRUE),
          lower = quantile(v, 0.025, names = FALSE, na.rm = TRUE),
          upper = quantile(v, 0.975, names = FALSE, na.rm = TRUE)
        )
      }))

      plot_data_list[[comp]] <- data.frame(
        x = x_seq, y = y_stats[, "mean"], ymin = y_stats[, "lower"], ymax = y_stats[, "upper"],
        Component = comp
      )
    }

    final_df <- do.call(rbind, plot_data_list)
    final_df$Component <- factor(final_df$Component, levels = unique_comps)

    y_label <- if (y_scale == "predicted") {
      switch(x$family,
        binomial = "Predicted Probability",
        poisson = "Predicted Expected Count",
        quasipoisson = "Predicted Expected Count",
        "Predicted Value"
      )
    } else {
      switch(x$family,
        binomial = "Partial Effect (Log-OR)",
        poisson = "Partial Effect (Log-RR)",
        quasipoisson = "Partial Effect (Log-RR)",
        "Partial Effect (\u0394Y)"
      )
    }

    p_c <- ggplot2::ggplot(final_df, ggplot2::aes(x = x, y = y, color = Component, fill = Component)) +
      ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::labs(title = "Dose-Response Curves", x = "Quantile Index", y = y_label) +
      ggplot2::scale_color_manual(values = curve_colors) +
      ggplot2::scale_fill_manual(values = curve_colors) +
      ggplot2::scale_x_continuous(
        limits = c(0, q_level - 1),
        expand = ggplot2::expansion(mult = c(0, 0.02)), breaks = seq(0, q_level - 1, by = 1)
      ) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        axis.line        = ggplot2::element_line(color = "#2C3E50", linewidth = 0.6),
        legend.position  = "right",
        legend.title     = ggplot2::element_blank(),
        legend.key.size  = ggplot2::unit(0.8, "lines"),
        plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold", size = base_size + 1)
      )

    if (y_scale == "partial") {
      p_c <- p_c + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#95A5A6", linewidth = 0.5)
    }

    if (overlay) {
      if (ci_available) p_c <- p_c + ggplot2::geom_ribbon(ggplot2::aes(ymin = ymin, ymax = ymax), alpha = 0.15, color = NA)
      p_c <- p_c + ggplot2::geom_line(linewidth = 1.1)
    } else {
      if (ci_available) p_c <- p_c + ggplot2::geom_ribbon(ggplot2::aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA)
      p_c <- p_c + ggplot2::geom_line(linewidth = 1.1) +
        ggplot2::facet_wrap(~Component, scales = "free_y", ncol = 2) +
        ggplot2::theme(legend.position = "none")
    }

    if (ci_available && !is_combined) {
      p_c <- p_c + ggplot2::labs(caption = rh_ci_caption) +
        ggplot2::theme(plot.caption = ggplot2::element_text(
          face = "italic", color = "red3", size = base_size - 3, hjust = 0
        ))
    }
  }

  if (type == "weights") {
    return(p_w)
  }
  if (type == "curves") {
    return(p_c)
  }

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Install 'patchwork': install.packages('patchwork')")
  }

  p_w + p_c +
    patchwork::plot_layout(widths = c(1, 1.8)) +
    patchwork::plot_annotation(
      title = "NWQS Model Diagnostics",
      caption = if (ci_available) rh_ci_caption else NULL,
      theme = ggplot2::theme(
        plot.title   = ggplot2::element_text(size = base_size + 3, face = "bold", hjust = 0.5),
        plot.caption = ggplot2::element_text(face = "italic", color = "red3", size = base_size - 2, hjust = 0),
        plot.margin  = ggplot2::margin(10, 10, 10, 10)
      )
    )
}


#' Print Method for NWQS Objects
#'
#' @description
#' Prints a concise, publication-style summary of the Non-linear Weighted Quantile Sum (NWQS)
#' model fit. It dynamically formats and outputs the overall mixture effect, component-specific
#' relative weights, and absolute effect contrasts (e.g., Q4 vs Q1) along with their 95\%
#' empirical confidence intervals.
#'
#' @param x An object of class \code{"nwqs"}.
#' @param digits Integer. Number of significant digits to print.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return Invisibly returns the original \code{"nwqs"} object.
#'
#' @export
#' @method print nwqs
print.nwqs <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n")
  print(x$call)

  cat("\n--- Non-linear WQS (NWQS) Results ---\n")
  cat(sprintf("Family: %s | Bootstrap: %d | Repeated Holdout: %d\n", x$family, x$b, x$rh))

  coef_table <- x$fit$coefficients
  if (!any(grepl("Pr\\(", colnames(coef_table)))) {
    z_stat <- coef_table$Estimate / coef_table$`Std. Error`
    p_values <- 2 * stats::pnorm(-abs(z_stat))
    names(p_values) <- rownames(coef_table)
    coef_table$`Pr(>|z|)` <- p_values
  } else {
    p_values <- coef_table[, grep("Pr\\(", colnames(coef_table))[1]]
    names(p_values) <- rownames(coef_table)
  }

  p_wqs <- p_values["wqs_score"]
  p_wqs_str <- ifelse(p_wqs < 0.001, "<0.001", sprintf("%.3f", p_wqs))

  q_level <- if (!is.null(x$q)) x$q else 4
  df_spline <- x$df_spline
  full_basis <- splines::ns(0:(q_level - 1),
    df = df_spline, knots = x$spline_knots,
    Boundary.knots = x$spline_boundary, intercept = FALSE
  )

  comps <- names(x$final_weights)
  rh <- x$rh

  eff_mat_str <- matrix("", nrow = length(comps) + 1, ncol = q_level - 1)
  rownames(eff_mat_str) <- c("Overall", comps)
  colnames(eff_mat_str) <- paste0("Q", 2:q_level, " vs Q1")

  is_exp_family <- x$family %in% c("binomial", "poisson", "quasipoisson")
  eff_name <- ifelse(is_exp_family, ifelse(x$family == "binomial", "OR", "RR"), "Delta")

  for (q_tgt in 2:q_level) {
    b_diff <- full_basis[q_tgt, ] - full_basis[1, ]
    iter_effs <- matrix(0, nrow = rh, ncol = length(comps) + 1)
    colnames(iter_effs) <- c("Overall", comps)

    for (i in seq_len(rh)) {
      beta_i <- if (rh == 1) x$mean_coefs["wqs_score"] else x$rh_coefs[i, "wqs_score"]
      comp_effs_i <- numeric(length(comps))
      names(comp_effs_i) <- comps

      for (comp in comps) {
        theta_cols <- paste0(comp, "_B", 1:df_spline)
        theta_i <- if (rh == 1) x$mean_shapes[theta_cols] else x$rh_shapes[i, theta_cols]
        w_i <- if (rh == 1) x$final_weights[comp] else x$rh_weights[i, comp]
        comp_effs_i[comp] <- beta_i * w_i * sum(b_diff * theta_i)
      }
      iter_effs[i, "Overall"] <- sum(comp_effs_i)
      iter_effs[i, comps] <- comp_effs_i
    }

    mean_eff <- colMeans(iter_effs, na.rm = TRUE)
    if (rh > 1) {
      lci <- apply(iter_effs, 2, quantile, probs = 0.025, na.rm = TRUE)
      uci <- apply(iter_effs, 2, quantile, probs = 0.975, na.rm = TRUE)
    } else {
      lci <- uci <- mean_eff
    }

    if (is_exp_family) {
      mean_eff <- exp(mean_eff)
      lci <- exp(lci)
      uci <- exp(uci)
    }

    eff_mat_str[, q_tgt - 1] <- if (rh > 1) {
      sprintf("%.3f [%.3f, %.3f]", mean_eff, lci, uci)
    } else {
      sprintf("%.3f", mean_eff)
    }
  }

  comps_sorted <- names(sort(x$final_weights, decreasing = TRUE))
  row_order <- c("Overall", comps_sorted)
  eff_mat_str <- eff_mat_str[row_order, , drop = FALSE]

  weight_col <- c("-", sprintf("%.3f", x$final_weights[comps_sorted]))
  print_df <- data.frame(Weight = weight_col, eff_mat_str, check.names = FALSE)

  if (rh > 1) {
    cat(
      "\n[!] CI Validity Warning:\n",
      "    The 95% empirical CIs below are derived from data-splitting variance\n",
      "    across Repeated Holdout (RH) iterations. They reflect ALGORITHMIC\n",
      "    SPLITTING VARIANCE only, NOT true sampling uncertainty.\n",
      "    These CIs are systematically too narrow and must NOT be used for\n",
      "    statistical inference or hypothesis testing.\n",
      "    --> Use nwqs_boot() to obtain valid percentile bootstrap CIs.\n",
      sep = ""
    )
  }

  cat(sprintf("\n>>> Joint & Component Effects (%s with 95%% Empirical CI):\n", eff_name))
  cat(sprintf("    Overall Significance (wqs_score): P = %s\n\n", p_wqs_str))
  print(print_df, right = FALSE)

  cat("\n>>> Model Coefficients (Averaged across RH iterations):\n")
  .print_coef_table(coef_table, digits = 4)

  aic_val <- if (is.na(x$fit$aic)) "NA" else sprintf("%.2f", x$fit$aic)
  cat(sprintf("\nMean AIC: %s | Mean Residual Deviance: %.2f\n", aic_val, x$fit$deviance))

  invisible(x)
}


#' Summary Method for NWQS Objects
#'
#' @description
#' Provides a detailed statistical summary of the NWQS model, including the global
#' regression coefficients, empirical standard errors, Z-statistics, P-values,
#' and deviance metrics.
#'
#' @param object An object of class \code{"nwqs"}.
#' @param digits Integer. Number of significant digits to print.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return Invisibly returns the original \code{"nwqs"} object.
#'
#' @export
#' @method summary nwqs
summary.nwqs <- function(object, digits = max(3L, getOption("digits") - 3L), ...) {
  coef_table <- object$fit$coefficients
  if (!any(grepl("Pr\\(", colnames(coef_table)))) {
    z_stat <- coef_table$Estimate / coef_table$`Std. Error`
    p_values <- 2 * stats::pnorm(-abs(z_stat))
    coef_table$`Pr(>|z|)` <- p_values
  }

  cat("\n=======================================================\n")
  cat("     Non-linear Weighted Quantile Sum (NWQS) Summary\n")
  cat("=======================================================\n")
  cat("Call:\n")
  print(object$call)
  cat(sprintf("\nModel Family: %s\n", object$family))

  cat("\nCoefficients:\n")
  .print_coef_table(coef_table, digits = 4)

  if (object$family != "quasipoisson") {
    cat(sprintf(
      "\nMean Null Deviance: %.2f on %d degrees of freedom",
      object$fit$null.deviance, object$fit$df.null
    ))
  }
  cat(sprintf(
    "\nMean Res. Deviance: %.2f on %d degrees of freedom\n",
    object$fit$deviance, object$fit$df.residual
  ))

  if (object$family != "quasipoisson") {
    cat(sprintf("Mean AIC: %.2f\n", object$fit$aic))
  } else {
    cat("Mean AIC: NA (Quasi-likelihood model)\n")
  }

  cat("=======================================================\n")
  invisible(object)
}


#' Extract Coefficients from NWQS Objects
#'
#' @param object An object of class \code{"nwqs"}.
#' @param ... Additional arguments.
#' @return A named numeric vector of averaged global coefficients.
#' @export
#' @method coef nwqs
coef.nwqs <- function(object, ...) {
  object$mean_coefs
}


# ==============================================================================
# S3 Methods: nwqs_boot
# ==============================================================================

#' Print Method for nwqs_boot Objects
#'
#' @description
#' Prints a concise summary of the \code{nwqs_boot} result, showing the point
#' estimate alongside outer bootstrap percentile confidence intervals for all
#' terms and targets.
#'
#' @param x An object of class \code{"nwqs_boot"}.
#' @param digits Integer. Significant digits. Defaults to 3.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns \code{x}.
#' @export
#' @method print nwqs_boot
print.nwqs_boot <- function(x, digits = 3, ...) {
  n_success <- x$ci_table$N_Success[1]
  n_total <- length(unique(x$boot_table$Boot_ID))
  conf_pct <- if (!is.null(x$conf_level)) sprintf("%.0f%%", x$conf_level * 100) else "95%"

  cat("\n--- NWQS Bootstrap Results ---\n")
  cat(sprintf(
    "Outer Bootstrap: %d total | %d successful | Conf. Level: %s\n",
    n_total, n_success, conf_pct
  ))

  cat("\n>>> Point Estimate & Bootstrap Percentile CI\n")
  cat("    (Valid sampling-variance CIs from outer bootstrap resampling)\n\n")
  print(x$formatted_table, digits = digits, row.names = FALSE)

  pf <- x$point_fit
  cat(sprintf("\n>>> Point Estimate Model (inner rh=%d fit on original data):\n", pf$rh))
  cat(sprintf(
    "    Family: %s | AIC: %s | Residual Deviance: %.2f\n\n",
    pf$family,
    ifelse(is.na(pf$fit$aic), "NA", sprintf("%.2f", pf$fit$aic)),
    pf$fit$deviance
  ))

  .print_coef_table(pf$fit$coefficients, digits = digits)
  invisible(x)
}


#' Summary Method for nwqs_boot Objects
#'
#' @description
#' Provides a detailed summary of the \code{nwqs_boot} result including bootstrap
#' CI table, component weight stability metrics, and point-estimate model summary.
#'
#' @param object An object of class \code{"nwqs_boot"}.
#' @param digits Integer. Significant digits. Defaults to 4.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns \code{object}.
#' @export
#' @method summary nwqs_boot
summary.nwqs_boot <- function(object, digits = 4, ...) {
  n_success <- object$ci_table$N_Success[1]
  n_total <- length(unique(object$boot_table$Boot_ID))
  conf_pct <- if (!is.null(object$conf_level)) sprintf("%.0f%%", object$conf_level * 100) else "95%"

  cat("\n=======================================================\n")
  cat("   NWQS Bootstrap Inference Summary (nwqs_boot)\n")
  cat("=======================================================\n")

  pf <- object$point_fit
  cat(sprintf("Family        : %s\n", pf$family))
  cat(sprintf("Outer Boots   : %d total (%d successful)\n", n_total, n_success))
  cat(sprintf("Inner RH      : rh = %d (per bootstrap replicate)\n", pf$rh))
  cat(sprintf("Conf. Level   : %s\n", conf_pct))
  cat(sprintf("Components    : %d\n", length(pf$final_weights)))

  cat("\n--- Bootstrap Percentile Confidence Intervals ---\n")
  cat("    (Outer bootstrap: valid for sampling-variance inference)\n\n")
  print(object$formatted_table, digits = digits, row.names = FALSE)

  cat("\n--- Component Weight Stability (across bootstrap replicates) ---\n")
  point_w <- pf$final_weights
  point_w_df <- data.frame(
    Component = names(point_w),
    Point_Weight = round(point_w, digits), row.names = NULL
  )

  max_target <- paste0("Q", object$q, "_vs_Q1")

  top_effs <- object$boot_table[
    object$boot_table$Target == max_target &
      object$boot_table$Term != "Overall",
  ]

  if (nrow(top_effs) > 0) {
    eff_sd <- aggregate(Estimate ~ Term, data = top_effs, FUN = sd)
    col_label <- paste0("Q", object$q, "_Effect_SD")
    names(eff_sd)[2] <- col_label
    eff_sd[[col_label]] <- round(eff_sd[[col_label]], digits)
    stability_df <- merge(point_w_df, eff_sd, by.x = "Component", by.y = "Term", all.x = TRUE)
    stability_df <- stability_df[order(-stability_df$Point_Weight), ]
    print(stability_df, row.names = FALSE)
    cat(sprintf(
      "    Note: %s = SD of Q%d vs Q1 effect across bootstrap replicates\n",
      col_label, object$q
    ))
    cat("    (Set keep_fits=TRUE in nwqs_boot() for full weight distributions)\n")
  } else {
    print(point_w_df[order(-point_w_df$Point_Weight), ], row.names = FALSE)
  }

  cat("\n--- Point Estimate: Model Coefficients ---\n")
  .print_coef_table(pf$fit$coefficients, digits = digits)

  cat(sprintf(
    "\nPoint AIC: %s | Point Residual Deviance: %.4f\n",
    ifelse(is.na(pf$fit$aic), "NA (quasi)", sprintf("%.4f", pf$fit$aic)),
    pf$fit$deviance
  ))
  cat("=======================================================\n")
  invisible(object)
}


#' Plot Diagnostics for Bootstrap NWQS Models
#'
#' @description
#' Generates diagnostic plots for \code{nwqs_boot} objects: component weights,
#' non-linear dose-response curves with bootstrap percentile confidence intervals,
#' or a combined panel of both.
#'
#' @details
#' When \code{type = "curves"} or \code{"both"}, the function extracts the fitted
#' shapes, weights, and \code{wqs_score} coefficients from each bootstrap replicate
#' (stored when \code{keep_fits = TRUE} in \code{\link{nwqs_boot}}) and projects them
#' onto the spline basis to compute per-component dose-response curves. The 95\%
#' confidence ribbons are derived from the 2.5th and 97.5th percentiles of the
#' bootstrap curve distribution — these are valid sampling-variance confidence intervals.
#'
#' If \code{keep_fits = FALSE}, only the point estimate curve (from the original data fit)
#' is shown without confidence intervals.
#'
#' @param x An object of class \code{"nwqs_boot"} resulting from \code{\link{nwqs_boot}}.
#' @param type Character. Plot type: \code{"both"} (default, weights + curves),
#'   \code{"curves"} (dose-response only), or \code{"weights"} (bar chart only).
#' @param y_scale Character. Y-axis scale: \code{"partial"} (default) or \code{"predicted"}.
#' @param components Character vector. Specific components to display. \code{NULL} for all.
#' @param overlay Logical. If \code{TRUE} (default), curves are overlaid; if \code{FALSE},
#'   faceted by component.
#' @param base_size Integer. Base font size. Default is 12.
#' @param palette Character. Color palette name.
#' @param top_n Integer or NULL. Number of top components (by weight) to show.
#' @param ... Additional arguments.
#'
#' @return A \code{ggplot} object or \code{patchwork} composite.
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_line geom_ribbon geom_hline
#'   scale_fill_manual scale_color_manual scale_y_continuous scale_x_continuous
#'   coord_flip facet_wrap theme_minimal labs theme expansion element_blank
#'   element_line element_text unit
#' @importFrom splines ns
#' @importFrom stats plogis quantile
#'
#' @seealso \code{\link{nwqs_boot}}, \code{\link{plot.nwqs}}
#' @export
#' @method plot nwqs_boot
plot.nwqs_boot <- function(x, type = c("both", "curves", "weights"),
                           y_scale = c("partial", "predicted"),
                           components = NULL, overlay = TRUE,
                           base_size = 12, palette = "default",
                           top_n = NULL, ...) {
  type <- match.arg(type)
  y_scale <- match.arg(y_scale)

  .get_palette <- function(n, pal = "default") {
    cols <- list(
      default = c(
        "#4A90C8", "#D92828", "#6EC44A", "#8B6FB8", "#00B4D8",
        "#006B3C", "#A8D8EA", "#F4B6B6", "#5BA3D0", "#E03030", "#7AD450", "#9B7FC0"
      ),
      palette2 = c(
        "#9bbf8a", "#82afda", "#f79059", "#e7dbd3", "#c2bdde",
        "#8dcec8", "#add3e2", "#3480b8", "#ffbe7a", "#fa8878", "#c82423", "#6b5b95"
      )
    )
    p <- if (pal %in% names(cols)) cols[[pal]] else cols[["default"]]
    rep(p, ceiling(n / length(p)))[seq_len(n)]
  }

  pf <- x$point_fit
  is_combined <- (type == "both")

  # -- Weight bar chart --------------------------------------------------------
  if (type %in% c("weights", "both")) {
    w <- pf$final_weights
    w_df <- data.frame(Component = factor(names(w), levels = names(sort(w))), Weight = as.numeric(w))
    pal_w <- .get_palette(nrow(w_df), palette)

    p_w <- ggplot2::ggplot(w_df, ggplot2::aes(x = Component, y = Weight, fill = Component)) +
      ggplot2::geom_col(alpha = 0.9, width = 0.7) +
      ggplot2::scale_fill_manual(values = pal_w, guide = "none") +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::labs(
        title = "Component Weights", subtitle = "Point estimate (original data)",
        x = NULL, y = "Weight"
      ) +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(color = "#2C3E50", linewidth = 0.5),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = base_size + 1),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = base_size - 1, color = "#7F8C8D")
      )
  }

  # -- Dose-Response Curves with Bootstrap CI ----------------------------------
  if (type %in% c("curves", "both")) {
    boot_fits <- x$boot_fits
    valid_fits <- NULL
    if (!is.null(boot_fits)) {
      valid_fits <- Filter(Negate(is.null), boot_fits)
      if (length(valid_fits) == 0) valid_fits <- NULL
    }
    has_boot_ci <- !is.null(valid_fits)

    if (!has_boot_ci) {
      message("No bootstrap fits available (keep_fits=FALSE or all failed). Showing point estimate only.")
    }

    q_level <- if (!is.null(pf$q)) pf$q else 4
    df_spline <- pf$df_spline
    x_seq <- seq(0, q_level - 1, length.out = 100)
    basis_mat <- splines::ns(x_seq,
      df = df_spline, knots = pf$spline_knots,
      Boundary.knots = pf$spline_boundary, intercept = FALSE
    )

    shape_names <- names(pf$mean_shapes)
    pattern <- "^(.+)_B(\\d+)$"
    parsed_names <- data.frame(
      full_name = shape_names,
      component = sub(pattern, "\\1", shape_names),
      basis_idx = as.numeric(sub(pattern, "\\2", shape_names)),
      stringsAsFactors = FALSE
    )
    unique_comps <- unique(parsed_names$component)
    if (!is.null(components)) unique_comps <- intersect(unique_comps, components)

    if (!is.null(top_n) && is.numeric(top_n) && top_n > 0) {
      top_components <- names(sort(pf$final_weights, decreasing = TRUE))[seq_len(min(top_n, length(pf$final_weights)))]
      unique_comps <- intersect(unique_comps, top_components)
      if (length(unique_comps) == 0) stop("No components match top_n filter.")
    }

    curve_colors <- .get_palette(length(unique_comps), palette)
    plot_data_list <- list()

    for (comp in unique_comps) {
      comp_cols <- parsed_names$full_name[parsed_names$component == comp]

      if (has_boot_ci) {
        y_boot_list <- lapply(valid_fits, function(bf) {
          theta <- bf$mean_shapes[comp_cols]
          beta_wqs <- bf$mean_coefs["wqs_score"]
          w <- bf$final_weights[comp]
          y_vec <- as.vector(basis_mat %*% (theta * beta_wqs * w))
          if (y_scale == "predicted") y_vec <- y_vec + bf$mean_coefs["(Intercept)"]
          return(y_vec)
        })
        y_pred_mat <- do.call(cbind, y_boot_list)
      } else {
        theta <- pf$mean_shapes[comp_cols]
        beta_wqs <- pf$mean_coefs["wqs_score"]
        w <- pf$final_weights[comp]
        y_vec <- as.vector(basis_mat %*% (theta * beta_wqs * w))
        if (y_scale == "predicted") y_vec <- y_vec + pf$mean_coefs["(Intercept)"]
        y_pred_mat <- matrix(y_vec, ncol = 1)
      }

      if (y_scale == "predicted") {
        if (pf$family == "binomial") {
          y_pred_mat <- stats::plogis(y_pred_mat)
        } else if (pf$family %in% c("poisson", "quasipoisson")) y_pred_mat <- exp(y_pred_mat)
      }

      y_stats <- t(apply(y_pred_mat, 1, function(v) {
        c(
          mean = mean(v, na.rm = TRUE),
          lower = quantile(v, 0.025, names = FALSE, na.rm = TRUE),
          upper = quantile(v, 0.975, names = FALSE, na.rm = TRUE)
        )
      }))

      plot_data_list[[comp]] <- data.frame(
        x = x_seq, y = y_stats[, "mean"],
        ymin = y_stats[, "lower"], ymax = y_stats[, "upper"],
        Component = comp
      )
    }

    final_df <- do.call(rbind, plot_data_list)
    final_df$Component <- factor(final_df$Component, levels = unique_comps)

    y_label <- if (y_scale == "predicted") {
      switch(pf$family,
        binomial = "Predicted Probability",
        poisson = "Predicted Expected Count",
        quasipoisson = "Predicted Expected Count",
        "Predicted Value"
      )
    } else {
      switch(pf$family,
        binomial = "Partial Effect (Log-OR)",
        poisson = "Partial Effect (Log-RR)",
        quasipoisson = "Partial Effect (Log-RR)",
        "Partial Effect (\u0394Y)"
      )
    }

    ci_subtitle <- if (has_boot_ci) {
      sprintf("Bootstrap Percentile CI (%d replicates)", length(valid_fits))
    } else {
      "Point Estimate (no CI available)"
    }

    p_c <- ggplot2::ggplot(final_df, ggplot2::aes(x = x, y = y, color = Component, fill = Component)) +
      ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::labs(
        title = "Dose-Response Curves", subtitle = ci_subtitle,
        x = "Quantile Index", y = y_label
      ) +
      ggplot2::scale_color_manual(values = curve_colors) +
      ggplot2::scale_fill_manual(values = curve_colors) +
      ggplot2::scale_x_continuous(
        limits = c(0, q_level - 1),
        expand = ggplot2::expansion(mult = c(0, 0.02)), breaks = seq(0, q_level - 1, by = 1)
      ) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        axis.line        = ggplot2::element_line(color = "#2C3E50", linewidth = 0.6),
        legend.position  = "right",
        legend.title     = ggplot2::element_blank(),
        legend.key.size  = ggplot2::unit(0.8, "lines"),
        plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold", size = base_size + 1),
        plot.subtitle    = ggplot2::element_text(hjust = 0.5, size = base_size - 1, color = "#7F8C8D")
      )

    if (y_scale == "partial") {
      p_c <- p_c + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#95A5A6", linewidth = 0.5)
    }

    if (overlay) {
      if (has_boot_ci) p_c <- p_c + ggplot2::geom_ribbon(ggplot2::aes(ymin = ymin, ymax = ymax), alpha = 0.15, color = NA)
      p_c <- p_c + ggplot2::geom_line(linewidth = 1.1)
    } else {
      if (has_boot_ci) p_c <- p_c + ggplot2::geom_ribbon(ggplot2::aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA)
      p_c <- p_c + ggplot2::geom_line(linewidth = 1.1) +
        ggplot2::facet_wrap(~Component, scales = "free_y", ncol = 2) +
        ggplot2::theme(legend.position = "none")
    }
  }

  if (type == "weights") {
    return(p_w)
  }
  if (type == "curves") {
    return(p_c)
  }

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Install 'patchwork': install.packages('patchwork')")
  }

  p_w + p_c +
    patchwork::plot_layout(widths = c(1, 1.6)) +
    patchwork::plot_annotation(
      title = "NWQS Bootstrap Diagnostics",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = base_size + 3, face = "bold", hjust = 0.5)
      )
    )
}


#' Extract Coefficients from nwqs_boot Objects
#'
#' @param object An object of class \code{"nwqs_boot"}.
#' @param ... Additional arguments.
#' @return A named numeric vector of point-estimate coefficients.
#' @export
#' @method coef nwqs_boot
coef.nwqs_boot <- function(object, ...) {
  object$point_fit$mean_coefs
}
