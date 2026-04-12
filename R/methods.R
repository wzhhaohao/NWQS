#' Format P-values
#'
#' @param p Numeric vector of P-values.
#' @return Formatted character vector.
#' @keywords internal
#' @noRd
.format_pval <- function(p) {
  sapply(p, function(pv) {
    if (is.na(pv)) return("NA")
    if (pv < 0.001) return("<0.001")
    return(sprintf("%.3f", pv))
  })
}

#' Print Model Coefficient Table
#'
#' @param coef_table Data frame or matrix with model coefficients.
#' @param digits Integer. Number of significant digits. Default is 4.
#' @return Invisibly returns the printed data frame.
#' @keywords internal
#' @noRd
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


#' @title Plot NWQS Model Diagnostics and Dose-Response Curves
#'
#' @description
#' Generates publication-quality diagnostic plots for an \code{nwqs} object,
#' including component weight distributions and non-linear dose-response
#' trajectories.
#'
#' @param x An object of class \code{"nwqs"}.
#' @param type Character. Plot type: \code{"both"} (default), \code{"curves"},
#'   or \code{"weights"}.
#' @param y_scale Character. Y-axis scale for curves: \code{"partial"}
#'   (default) or \code{"predicted"}.
#' @param components Character vector. Specific components to display. If
#'   \code{NULL}, all are shown.
#' @param overlay Logical. If \code{TRUE} (default), all curves are overlaid.
#'   If \code{FALSE}, faceted by component.
#' @param plot_ci Logical. Whether to plot 95\% empirical confidence ribbons
#'   from repeated holdout iterations. Default is \code{FALSE}.
#' @param base_size Integer. Base font size. Default is 12.
#' @param palette Character. Color palette: \code{"default"}, \code{"palette2"},
#'   \code{"palette3"}, or \code{"palette4"}.
#' @param colorblind_friendly Logical. If \code{TRUE}, uses a colorblind-safe
#'   palette. Default is \code{FALSE}.
#' @param top_n Integer or \code{NULL}. Show only the top \code{top_n}
#'   components by weight.
#' @param ... Additional arguments passed to other methods.
#'
#' @return A \code{ggplot} or \pkg{patchwork} composite object.
#'
#' @seealso \code{\link{plot.nwqs_boot}} for bootstrap-based inference plots.
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
        scaling_factor <- x$rh_coefs[, "nwqs"] * x$rh_weights[, comp]
        beta_mat <- sweep(beta_mat, 1, scaling_factor, "*")
        y_pred_mat <- as.matrix(basis_mat) %*% t(beta_mat)
        if (y_scale == "predicted") {
          y_pred_mat <- sweep(y_pred_mat, 2, x$rh_coefs[, "(Intercept)"], "+")
        }
      } else {
        beta_mat <- matrix(x$mean_shapes[comp_cols], nrow = 1)
        scaling_factor <- x$mean_coefs["nwqs"] * x$final_weights[comp]
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

  if (type == "weights") return(p_w)
  if (type == "curves") return(p_c)

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Install 'patchwork': install.packages('patchwork')")
  }

  p_w + p_c +
    patchwork::plot_layout(widths = c(1, 1)) +
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


#' @title Print NWQS Model Summary
#'
#' @description
#' Prints a concise publication-quality summary of the NWQS model including
#' overall mixture effects, component weights, and quantile contrast effects
#' with 95\% empirical confidence intervals.
#'
#' @param x An object of class \code{"nwqs"}.
#' @param digits Integer. Number of significant digits.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns the \code{"nwqs"} object.
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

  p_wqs <- p_values["nwqs"]
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
      beta_i <- if (rh == 1) x$mean_coefs["nwqs"] else x$rh_coefs[i, "nwqs"]
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
  cat(sprintf("    Overall Significance (nwqs): P = %s\n\n", p_wqs_str))
  print(print_df, right = FALSE)

  cat("\n>>> Model Coefficients (Averaged across RH iterations):\n")
  .print_coef_table(coef_table, digits = 4)

  aic_val <- if (is.na(x$fit$aic)) "NA" else sprintf("%.2f", x$fit$aic)
  cat(sprintf("\nMean AIC: %s | Mean Residual Deviance: %.2f\n", aic_val, x$fit$deviance))

  invisible(x)
}


#' @title Detailed Statistical Summary for NWQS Model
#'
#' @description
#' Provides a detailed statistical summary of the NWQS model including
#' regression coefficients, standard errors, Z-statistics, P-values, and
#' goodness-of-fit metrics (deviance, AIC).
#'
#' @param object An object of class \code{"nwqs"}.
#' @param digits Integer. Number of significant digits.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns the \code{"nwqs"} object.
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


#' @title Print Bootstrap NWQS Inference Results
#'
#' @description
#' Prints the \code{nwqs_boot} summary reflecting true sampling variance.
#' Shows all terms (overall and component-specific effects) with point
#' estimates and bootstrap percentile confidence intervals. For exponential
#' family models, coefficients are automatically converted to OR or RR.
#'
#' @param x An object of class \code{"nwqs_boot"}.
#' @param digits Integer. Number of significant digits. Default is 3.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns \code{x}.
#' @export
#' @method print nwqs_boot
print.nwqs_boot <- function(x, digits = 3, ...) {
  ci_tab <- x$ci_table

  if (is.null(ci_tab) || !is.data.frame(ci_tab) || nrow(ci_tab) == 0) {
    stop("`x$ci_table` is missing or empty.")
  }

  n_success <- x$n_success
  n_total <- x$n_boot
  conf_pct <- if (!is.null(x$conf_level)) sprintf("%.0f%%", x$conf_level * 100) else "95%"

  is_exp_family <- !is.null(x$family) && x$family %in% c("binomial", "poisson", "quasipoisson", "clogit")

  cat("\n--- NWQS Bootstrap Ensemble Results ---\n")
  cat(sprintf(
    "Bootstrap Replicates: %d total | %d successful | Conf. Level: %s\n",
    n_total, n_success, conf_pct
  ))
  cat(sprintf(
    "Model Family: %s | Quantiles (q): %d | Inner RH: %d\n",
    x$family, x$q, x$rh_inner
  ))

  cat("\n>>> Ensemble Effects & Bootstrap Percentile CI\n")
  if (is_exp_family) {
    cat("    (Displayed on exponentiated scale: OR/RR)\n\n")
  } else {
    cat("    (Displayed on the original coefficient scale)\n\n")
  }

  print(x$formatted_table, digits = digits, row.names = FALSE)

  invisible(x)
}


#' @title Detailed Summary for NWQS Bootstrap Ensemble Model
#'
#' @description
#' Provides a detailed summary of \code{nwqs_boot} results including
#' bootstrap confidence interval tables and component weight stability
#' metrics across bootstrap replicates.
#'
#' @param object An object of class \code{"nwqs_boot"}.
#' @param digits Integer. Number of significant digits. Default is 4.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns \code{object}.
#' @export
#' @method summary nwqs_boot
summary.nwqs_boot <- function(object, digits = 4, ...) {
  n_success <- object$n_success
  n_total <- object$n_boot
  conf_pct <- if (!is.null(object$conf_level)) sprintf("%.0f%%", object$conf_level * 100) else "95%"

  cat("\n=======================================================\n")
  cat("   NWQS Bootstrap Ensemble Inference Summary\n")
  cat("=======================================================\n")

  cat(sprintf("Family        : %s\n", object$family))
  cat(sprintf("Outer Boots   : %d total (%d successful)\n", n_total, n_success))
  cat(sprintf("Inner RH      : rh = %d (per bootstrap replicate)\n", object$rh_inner))
  cat(sprintf("Conf. Level   : %s\n", conf_pct))
  cat(sprintf("Components    : %d\n", length(object$final_weights)))

  cat("\n--- Ensemble Effects & Bootstrap Percentile CI ---\n")
  cat("    (Valid for sampling-variance inference)\n\n")
  print(object$formatted_table, digits = digits, row.names = FALSE)

  cat("\n--- Component Weight Stability (across bootstrap replicates) ---\n")
  point_w <- object$final_weights
  point_w_df <- data.frame(
    Component = names(point_w),
    Ensemble_Weight = round(point_w, digits), row.names = NULL
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
    stability_df <- stability_df[order(-stability_df$Ensemble_Weight), ]
    print(stability_df, row.names = FALSE)
    cat(sprintf(
      "    Note: %s = SD of Q%d vs Q1 effect across bootstrap replicates\n",
      col_label, object$q
    ))
    cat("    (Set keep_fits=TRUE in nwqs_boot() for full weight distributions)\n")
  } else {
    print(point_w_df[order(-point_w_df$Ensemble_Weight), ], row.names = FALSE)
  }

  cat("=======================================================\n")
  invisible(object)
}


#' @title Plot NWQS Bootstrap Diagnostic Charts
#'
#' @description
#' Generates inference-level diagnostic plots for \code{nwqs_boot} objects.
#' Visualizes quantile contrasts as boxplots based on bootstrap distributions,
#' combined with component weight bar charts.
#'
#' @param x An object of class \code{"nwqs_boot"}.
#' @param type Character. Plot type: \code{"both"} (default), \code{"curves"}
#'   (contrast boxplots), or \code{"weights"}.
#' @param y_scale Character. Effect scale: \code{"partial"}, \code{"contrast_or"},
#'   or \code{"predicted"}.
#' @param components Character vector. Specific components to display.
#' @param overlay Logical. Retained for backward compatibility.
#' @param show_ci Logical. Retained for backward compatibility.
#' @param base_size Integer. Base font size. Default is 12.
#' @param palette Character. Discrete color palette.
#' @param top_n Integer or \code{NULL}. Show only the top \code{top_n}
#'   components by weight.
#' @param ylim Numeric vector. Force Y-axis limits.
#' @param y_step Numeric. Y-axis tick spacing.
#' @param free_y Logical. If \code{TRUE}, each facet has free Y-axis scaling.
#' @param fill_alpha Numeric. Box fill transparency. Default is 0.16.
#' @param exponentiate Logical or \code{NULL}. Whether to exponentiate effects.
#'   If \code{NULL}, auto-detected from model family.
#' @param ... Additional arguments passed to \code{plot_nwqs_contrast_box}.
#'
#' @return A \code{ggplot} or \pkg{patchwork} composite object.
#' @export
#' @method plot nwqs_boot
plot.nwqs_boot <- function(x,
                           type = c("both", "curves", "weights"),
                           y_scale = c("partial", "contrast_or", "predicted"),
                           components = NULL,
                           overlay = TRUE,
                           show_ci = TRUE,
                           base_size = 12,
                           palette = "default",
                           top_n = NULL,
                           ylim = NULL,
                           y_step = NULL,
                           free_y = TRUE,
                           fill_alpha = 0.16,
                           exponentiate = NULL,
                           ...) {
  type <- match.arg(type)
  y_scale <- match.arg(y_scale)

  .clean_name <- function(nm) {
    nm <- gsub("_adj$", "", nm)
    nm <- gsub("^(ln|log|log10|log2|scale)_", "", nm, ignore.case = TRUE)
    nm <- gsub("URX", "", nm)
    nm
  }

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

  .resolve_selected_raw <- function(x, components = NULL, top_n = NULL) {
    if (!is.null(x$final_weights) && length(x$final_weights) > 0) {
      selected <- names(sort(x$final_weights, decreasing = TRUE))
    } else {
      selected <- character(0)
    }

    if (!is.null(components)) {
      keep <- selected %in% components | .clean_name(selected) %in% components
      selected <- selected[keep]
    }

    if (!is.null(top_n) && is.numeric(top_n) && top_n > 0) {
      selected <- selected[seq_len(min(top_n, length(selected)))]
    }

    selected
  }

  if (is.null(exponentiate)) {
    is_exp_family <- x$family %in% c("binomial", "poisson", "quasipoisson")
    exponentiate <- if (y_scale == "partial") FALSE else is_exp_family
  }

  is_combined <- (type == "both")

  selected_raw <- .resolve_selected_raw(x, components = components, top_n = top_n)
  selected_clean <- .clean_name(selected_raw)

  global_colors <- .get_palette(length(selected_clean), palette)
  names(global_colors) <- selected_clean

  if (type %in% c("weights", "both")) {
    w <- x$final_weights

    if (!is.null(components)) {
      keep <- names(w) %in% components | .clean_name(names(w)) %in% components
      w <- w[keep]
    }
    if (!is.null(top_n) && is.numeric(top_n) && top_n > 0) {
      w <- sort(w, decreasing = TRUE)[seq_len(min(top_n, length(w)))]
    } else {
      w <- sort(w, decreasing = TRUE)
    }

    if (length(w) == 0) stop("No components remain after filtering.")

    clean_names_w <- .clean_name(names(w))
    clean_levels_w <- .clean_name(names(sort(w, decreasing = TRUE)))

    w_df <- data.frame(
      Component = factor(clean_names_w, levels = rev(clean_levels_w)),
      Weight = as.numeric(w)
    )

    p_w <- ggplot2::ggplot(
      w_df,
      ggplot2::aes(x = Component, y = Weight, fill = Component)
    ) +
      ggplot2::geom_col(alpha = 0.58, width = 0.68) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.3f", Weight)),
        hjust = -0.10,
        size = base_size / 3.4,
        color = "black"
      ) +
      ggplot2::scale_fill_manual(values = global_colors, guide = "none") +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.12))
      ) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::labs(
        title = "Component Ensemble Weights",
        subtitle = "Averaged across bootstraps",
        x = NULL,
        y = "Weight"
      ) +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(color = "#2C3E50", linewidth = 0.5),
        axis.text = ggplot2::element_text(color = "black"),
        axis.title = ggplot2::element_text(color = "black"),
        plot.title = ggplot2::element_text(
          hjust = 0.5, face = "bold", size = base_size + 1
        ),
        plot.subtitle = ggplot2::element_text(
          hjust = 0.5, size = base_size - 1, color = "#7F8C8D"
        )
      )
  }

  if (type %in% c("curves", "both")) {
    p_c <- plot_nwqs_contrast_box(
      model = x,
      exponentiate = exponentiate,
      free_y = free_y,
      base_size = base_size,
      fill_alpha = fill_alpha,
      palette = palette,
      components = components,
      top_n = top_n,
      ylim = ylim,
      y_step = y_step
    )
  }

  if (type == "weights") return(p_w)
  if (type == "curves") return(p_c)

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Please install 'patchwork': install.packages('patchwork')")
  }

  p_w + p_c +
    patchwork::plot_layout(widths = c(1, 1)) +
    patchwork::plot_annotation(
      title = "NWQS Model Diagnostics",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = base_size + 3,
          face = "bold",
          hjust = 0.5
        ),
        plot.margin = ggplot2::margin(10, 10, 10, 10)
      )
    )
}


#' @title Extract Coefficients from NWQS Bootstrap Object
#'
#' @description
#' Extracts point estimate coefficients from the original full-data fit, not
#' bootstrap means.
#'
#' @param object An object of class \code{"nwqs_boot"}.
#' @param ... Additional arguments.
#'
#' @return A named numeric vector of global model coefficients.
#' @export
#' @method coef nwqs_boot
coef.nwqs_boot <- function(object, ...) {
  object$point_fit$mean_coefs
}
