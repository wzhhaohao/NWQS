#' Plot Diagnostics for Non-linear Weighted Quantile Sum (NWQS) Models
#'
#' @description 
#' Generates publication-ready diagnostic plots for `nwqs` objects. This function 
#' can produce component weight distributions and non-linear dose-response trajectories 
#' (partial effects or predicted values) for the overall mixture and individual components.
#'
#' @details 
#' The dose-response curves (when `type = "curves"` or `"both"`) map the estimated non-linear 
#' spline functions. When `y_scale = "partial"`, the curves represent the isolated effect 
#' of increasing exposure while holding others constant (zero-centered log-OR, log-RR, or \eqn{\Delta Y}). 
#' When `y_scale = "predicted"`, the curves show the absolute predicted scale (e.g., probabilities 
#' for binomial models or expected counts for Poisson models).
#' 
#' @param x An object of class \code{"nwqs"} resulting from a call to \code{\link{nwqs}}.
#' @param type Character. The type of diagnostic plot to generate. Options are 
#'   \code{"both"} (default, combines weights and curves using \pkg{patchwork}), 
#'   \code{"curves"} (dose-response trajectories), or \code{"weights"} (bar plot of component weights).
#' @param y_scale Character. The scale for the y-axis in the curve plots. Options are 
#'   \code{"partial"} (default, shows relative effect changes like \eqn{\Delta \eta}) or 
#'   \code{"predicted"} (shows absolute predicted values like probabilities or counts).
#' @param components Character vector. Specific mixture components to display in the plots. 
#'   If \code{NULL} (default), all components are shown.
#' @param overlay Logical. If \code{TRUE} (default), all component dose-response curves 
#'   are overlaid on a single plot. If \code{FALSE}, curves are faceted by component.
#' @param plot_ci Logical. Whether to plot the 95\% empirical confidence intervals 
#'   (shaded ribbons or error bars) derived from Repeated Holdout (RH) iterations. 
#'   Default is \code{FALSE}. Note: Requires \code{rh > 1} in the \code{nwqs} model.
#' @param base_size Integer. Base font size for the \pkg{ggplot2} theme. Default is 12.
#' @param palette Character. Color palette to use. Options include \code{"default"} 
#'   (medical/epidemiology optimized), \code{"palette2"} (pastel), \code{"palette3"}, 
#'   or \code{"palette4"}.
#' @param colorblind_friendly Logical. If \code{TRUE}, forces the use of a colorblind-safe 
#'   palette (overrides the \code{palette} argument). Default is \code{FALSE}.
#' @param top_n Integer. Number of top components (ranked by assigned weight) to display 
#'   in the curve plot. If \code{NULL} (default), all specified components are shown.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return A \code{ggplot} object (if \code{type} is "weights" or "curves") or a 
#'   \code{patchwork} composite object (if \code{type = "both"}).
#' 
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual scale_y_continuous coord_flip theme_minimal labs theme expansion element_blank element_line element_text geom_errorbar geom_hline geom_ribbon geom_line facet_wrap
#' @importFrom splines ns
#' @importFrom stats plogis quantile
#' 
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

  .get_palette <- function(n_colors, palette = "default") {

    default <- c(
      "#4A90C8", "#D92828", "#6EC44A", "#8B6FB8", "#00B4D8",
      "#006B3C", "#A8D8EA", "#F4B6B6", "#5BA3D0", "#E03030",
      "#7AD450", "#9B7FC0"
    )

    palette2 <- c(
      "#9bbf8a", "#82afda", "#f79059", "#e7dbd3", "#c2bdde",
      "#8dcec8", "#add3e2", "#3480b8", "#ffbe7a", "#fa8878",
      "#c82423", "#6b5b95"
    )

    palette3 <- c(
      "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
      "#E6AB02", "#A6761D", "#666666", "#4DAF4A", "#377EB8",
      "#FFFF33", "#984EA3"
    )

    palette4 <- c(
      "#004C90", "#E60000", "#009E73", "#E69F00", "#56B4E9",
      "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000",
      "#1F78B4", "#B2DF8A"
    )

    palettes <- list(
      "default" = default,
      "palette2" = palette2,
      "palette3" = palette3,
      "palette4" = palette4
    )

    palette_lower <- tolower(palette)
    if (palette_lower %in% names(palettes)) {
      selected_palette <- palettes[[palette_lower]]
    } else {
      selected_palette <- palettes[["default"]]
      warning(paste("Palette '", palette, "' not found. Using 'default' instead.", sep = ""))
    }


    # 颜色数量适配（超过 12 个开始循环）
    if (n_colors > length(selected_palette)) {
      n_cycles <- ceiling(n_colors / length(selected_palette))
      selected_palette <- rep(selected_palette, n_cycles)[1:n_colors]
    } else {
      selected_palette <- selected_palette[1:n_colors]
    }

    return(selected_palette)
  }

  # 模块 1: 权重图 (Weights)
  if (type %in% c("weights", "both")) {
    w_names <- names(x$final_weights)
    if (!is.null(components)) w_names <- intersect(w_names, components)

    w_df <- data.frame(Component = w_names, Weight = x$final_weights[w_names])

    if (!is.null(x$rh_weights) && x$rh > 1 && plot_ci) {
      rh_w <- x$rh_weights[, w_names, drop = FALSE]
      w_df$Lower <- apply(rh_w, 2, function(v) quantile(v, 0.025, names = FALSE, na.rm = TRUE))
      w_df$Upper <- apply(rh_w, 2, function(v) quantile(v, 0.975, names = FALSE, na.rm = TRUE))
    } else {
      w_df$Lower <- w_df$Weight
      w_df$Upper <- w_df$Weight
    }

    w_df$Component <- factor(w_df$Component, levels = w_df$Component[order(w_df$Weight)])
    w_colors <- .get_palette(length(w_names), palette = palette)

    p_w <- ggplot2::ggplot(w_df, ggplot2::aes(x = Component, y = Weight, fill = Component)) +
      ggplot2::geom_col(alpha = 0.9, width = 0.7) +
      ggplot2::scale_fill_manual(values = w_colors, guide = "none") +
      ggplot2::scale_y_continuous(
        limits = c(0, NA),
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
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

    if (plot_ci && x$rh > 1) {
      p_w <- p_w + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = Lower, ymax = Upper),
        width = 0.3, color = "#2C3E50", linewidth = 0.6
      )
    }
  }

  # 模块 2: 剂量反应曲线 (Curves)
  if (type %in% c("curves", "both")) {
    q_level <- if (!is.null(x$call$q)) x$call$q else 4
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

    # 只展示权重前 top_n 的污染物
    if (!is.null(top_n) && is.numeric(top_n) && top_n > 0) {
      weight_order <- sort(x$final_weights, decreasing = TRUE)
      top_components <- names(weight_order)[1:min(top_n, length(weight_order))]
      unique_comps <- intersect(unique_comps, top_components)

      if (length(unique_comps) == 0) {
        stop("No components match the specified top_n and components filter.")
      }

      message(sprintf(
        "Showing top %d components by weight: %s",
        length(unique_comps),
        paste(unique_comps, collapse = ", ")
      ))
    }

    n_components <- length(unique_comps)
    curve_colors <- .get_palette(n_components, palette = palette)

    plot_data_list <- list()
    for (comp in unique_comps) {
      comp_cols <- parsed_names$full_name[parsed_names$component == comp]

      if (!is.null(x$rh_shapes) && plot_ci && x$rh > 1) {
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
        x = x_seq, y = y_stats[, "mean"],
        ymin = y_stats[, "lower"], ymax = y_stats[, "upper"],
        Component = comp
      )
    }
    final_df <- do.call(rbind, plot_data_list)
    final_df$Component <- factor(final_df$Component, levels = unique_comps)

    if (y_scale == "predicted") {
      y_label <- switch(x$family,
        "binomial" = "Predicted Probability",
        "poisson" = "Predicted Expected Count",
        "quasipoisson" = "Predicted Expected Count",
        "gaussian" = "Predicted Value"
      )
    } else {
      y_label <- switch(x$family,
        "binomial" = "Partial Effect (Log-OR)",
        "poisson" = "Partial Effect (Log-RR)",
        "quasipoisson" = "Partial Effect (Log-RR)",
        "gaussian" = "Partial Effect (ΔY)"
      )
    }

    p_c <- ggplot2::ggplot(final_df, ggplot2::aes(x = x, y = y, color = Component, fill = Component)) +
      ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::labs(title = "Dose-Response Curves", x = "Quantile Index", y = y_label) +
      ggplot2::scale_color_manual(values = curve_colors) +
      ggplot2::scale_fill_manual(values = curve_colors) +
      ggplot2::scale_x_continuous(
        limits = c(0, q_level - 1),
        expand = ggplot2::expansion(mult = c(0, 0.02)),
        breaks = seq(0, q_level - 1, by = 1)
      ) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(color = "#2C3E50", linewidth = 0.6),
        legend.position = "right",
        legend.title = ggplot2::element_blank(),
        legend.key.size = ggplot2::unit(0.8, "lines"),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = base_size + 1)
      )

    if (y_scale == "partial") {
      p_c <- p_c + ggplot2::geom_hline(
        yintercept = 0, linetype = "dashed",
        color = "#95A5A6", linewidth = 0.5
      )
    }

    if (overlay) {
      if (plot_ci && x$rh > 1) {
        p_c <- p_c + ggplot2::geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.15, color = NA)
      }
      p_c <- p_c + ggplot2::geom_line(linewidth = 1.1)
    } else {
      if (plot_ci && x$rh > 1) {
        p_c <- p_c + ggplot2::geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA)
      }
      p_c <- p_c + ggplot2::geom_line(linewidth = 1.1) +
        ggplot2::facet_wrap(~Component, scales = "free_y", ncol = 2) +
        ggplot2::theme(legend.position = "none")
    }
  }

  # Output
  if (type == "weights") {
    return(p_w)
  }
  if (type == "curves") {
    return(p_c)
  }
  if (type == "both") {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop("Install 'patchwork' package: install.packages('patchwork')")
    }
    combined_plot <- p_w + p_c +
      patchwork::plot_layout(widths = c(1, 1.8)) +
      patchwork::plot_annotation(
        title = "NWQS Model Diagnostics",
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(size = base_size + 3, face = "bold", hjust = 0.5),
          plot.margin = ggplot2::margin(10, 10, 10, 10)
        )
      )
    return(combined_plot)
  }
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
    p_values <- 2 * pnorm(-abs(z_stat))
    names(p_values) <- rownames(coef_table)
    coef_table$`Pr(>|z|)` <- p_values
  } else {
    p_values <- coef_table[, grep("Pr\\(", colnames(coef_table))[1]]
    names(p_values) <- rownames(coef_table)
  }

  p_wqs <- p_values["wqs_score"]
  p_wqs_str <- ifelse(p_wqs < 0.001, "<0.001", sprintf("%.3f", p_wqs))

  q_level <- ifelse(is.null(x$q), 4, eval(x$q))
  df_spline <- x$df_spline
  full_basis <- splines::ns(0:(q_level - 1), df = df_spline, knots = x$spline_knots, Boundary.knots = x$spline_boundary, intercept = FALSE)

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
      lci <- mean_eff
      uci <- mean_eff
    }

    if (is_exp_family) {
      mean_eff <- exp(mean_eff)
      lci <- exp(lci)
      uci <- exp(uci)
    }

    if (rh > 1) {
      formatted_str <- sprintf("%.3f [%.3f, %.3f]", mean_eff, lci, uci)
    } else {
      formatted_str <- sprintf("%.3f", mean_eff)
    }
    eff_mat_str[, q_tgt - 1] <- formatted_str
  }

  comps_sorted <- names(sort(x$final_weights, decreasing = TRUE))
  row_order <- c("Overall", comps_sorted)
  eff_mat_str <- eff_mat_str[row_order, , drop = FALSE]

  weight_col <- c("-", sprintf("%.3f", x$final_weights[comps_sorted]))
  print_df <- data.frame(Weight = weight_col, eff_mat_str, check.names = FALSE)

  cat(sprintf("\n>>> Joint & Component Effects (%s with 95%% Empirical CI):\n", eff_name))
  cat(sprintf("    Overall Significance (wqs_score): P = %s\n\n", p_wqs_str))
  print(print_df, right = FALSE)

  cat("\n>>> Model Coefficients (Averaged across RH iterations):\n")
  stats::printCoefmat(as.matrix(coef_table), digits = 4, has.Pvalue = TRUE, P.values = TRUE, na.print = "NA")

  aic_val <- if (is.na(x$fit$aic)) "NA" else sprintf("%.2f", x$fit$aic)
  cat(sprintf("\nMean AIC: %s | Mean Residual Deviance: %.2f\n", aic_val, x$fit$deviance))

  invisible(x)
}

#' Summary Method for NWQS Objects
#'
#' @description 
#' Provides a detailed statistical summary of the NWQS model, including the global 
#' regression coefficients, empirical standard errors (derived from outer bootstrap 
#' or repeated holdout iterations), Z-statistics, P-values, and deviance metrics.
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
    p_values <- 2 * pnorm(-abs(z_stat))
    coef_table$`Pr(>|z|)` <- p_values
  }

  cat("\n=======================================================\n")
  cat("     Non-linear Weighted Quantile Sum (NWQS) Summary\n")
  cat("=======================================================\n")
  cat("Call:\n")
  print(object$call)
  cat(sprintf("\nModel Family: %s\n", object$family))

  cat("\nCoefficients (Empirical Bootstrap SE used for inference):\n")
  stats::printCoefmat(as.matrix(coef_table), digits = 4, has.Pvalue = TRUE, P.values = TRUE, na.print = "NA")

  if (object$family != "quasipoisson") {
    cat(sprintf("\nMean Null Deviance: %.2f on %d degrees of freedom", object$fit$null.deviance, object$fit$df.null))
  }
  cat(sprintf("\nMean Res. Deviance: %.2f on %d degrees of freedom\n", object$fit$deviance, object$fit$df.residual))

  if (object$family != "quasipoisson") {
    cat(sprintf("Mean AIC: %.2f\n", object$fit$aic))
  } else {
    cat("Mean AIC: NA (Quasi-likelihood model)\n")
  }
  cat("=======================================================\n")
  invisible(object)
}

#' Extract Global Coefficients from NWQS Objects
#'
#' @description 
#' Extracts the averaged global regression coefficients (intercept, overall \code{wqs_score} 
#' effect, and covariates) from the fitted NWQS model. These coefficients represent the 
#' ensemble average across all Repeated Holdout (RH) or Bootstrap iterations.
#'
#' @param object An object of class \code{"nwqs"}.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return A named numeric vector of the averaged global coefficients.
#' 
#' @export
#' @method coef nwqs
coef.nwqs <- function(object, ...) {
  return(object$mean_coefs)
}
