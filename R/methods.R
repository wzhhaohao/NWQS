# ==============================================================================
# Internal Helpers (内部辅助函数)
# ==============================================================================

#' 格式化 P 值 (内部辅助函数)
#'
#' @description 将 P 值格式化为符合医学期刊发表标准（如 <0.001）的字符串。
#' @param p 数值向量。包含需要格式化的 P 值。
#' @return 格式化后的字符向量。
#' @keywords internal
#' @noRd
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

#' 打印模型系数表 (内部辅助函数)
#'
#' @description 格式化输出 GLM 系数表，自动计算并附加 Z 检验的 P 值及显著性星号。
#' @param coef_table 数据框或矩阵。包含模型系数、标准误等。
#' @param digits 整数。保留的有效数字位数，默认为 4。
#' @return 隐式返回打印的数据框。
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


# ==============================================================================
# S3 Methods: nwqs (点估计与经验算法推断)
# ==============================================================================

#' @title 绘制非线性加权分位数和 (NWQS) 模型的诊断与剂量反应图
#'
#' @description
#' 为 \code{nwqs} 对象生成达到出版标准的诊断图。该函数能够可视化混合物各组分的权重分布，
#' 以及整体混合物或单个组分的非线性剂量反应轨迹（部分效应或预测绝对值）。
#'
#' @details
#' \strong{剂量反应曲线解读:} \cr
#' 当设置 \code{type = "curves"} 或 \code{"both"} 时，函数将映射估计的非线性样条曲线。
#' \itemize{
#'   \item 若 \code{y_scale = "partial"}：曲线代表在保持其他变量不变的情况下，暴露水平增加的独立相对效应
#'     （例如：零中心化的对数比值比 Log-OR、对数相对危险度 Log-RR 或连续型结局的变化量 \eqn{\Delta Y}）。
#'   \item 若 \code{y_scale = "predicted"}：曲线展示绝对预测尺度（例如：二分类模型的绝对预测概率）。
#' }
#'
#' \strong{关于置信带 (Confidence Ribbons) 的严正警告:} \cr
#' 当 \code{plot_ci = TRUE} 且 \code{rh > 1} 时，图中的 95\% 经验置信带来源于重复保留 (Repeated Holdout) 迭代。
#' 请务必注意：**这些置信带仅反映了数据拆分带来的算法方差 (Algorithmic Variance)，绝对不能用于正式的统计推断或假设检验。**
#' 它们通常比真实的置信带更窄。若需发表具备统计学效力的带有真实置信区间的图表，请务必使用 \code{\link{plot.nwqs_boot}}。
#'
#' @param x \code{"nwqs"} 类的对象，由 \code{\link{nwqs}} 函数拟合产生。
#' @param type Character。绘图类型：\code{"both"}（默认，使用 \pkg{patchwork} 拼接权重与曲线）、\code{"curves"}（非线性剂量反应曲线）或 \code{"weights"}（组分权重条形图）。
#' @param y_scale Character。曲线的 Y 轴尺度：\code{"partial"}（默认，相对偏效应）或 \code{"predicted"}（绝对预测值）。
#' @param components Character vector。需要显示的特定混合物组分名称。若为 \code{NULL}（默认），则显示所有组分。
#' @param overlay Logical。若为 \code{TRUE}（默认），所有组分的曲线将叠加在同一图层中。若为 \code{FALSE}，则按组分进行分面 (Faceted) 绘图。
#' @param plot_ci Logical。是否绘制来自 Repeated Holdout 迭代的 95\% 经验置信带。默认为 \code{FALSE}。需要 \code{rh > 1}。
#' @param base_size Integer。\pkg{ggplot2} 主题的基础字体大小，默认为 12。
#' @param palette Character。配色方案，可选：\code{"default"}, \code{"palette2"}, \code{"palette3"}, 或 \code{"palette4"}。
#' @param colorblind_friendly Logical。若为 \code{TRUE}，则强制使用对色盲友好的调色板，默认为 \code{FALSE}。
#' @param top_n Integer 或 \code{NULL}。按权重降序排列，仅显示前 \code{top_n} 个最重要的组分。
#' @param ... 传递给其他方法的额外参数。
#'
#' @return 返回一个 \code{ggplot} 对象（当 \code{type} 为 \code{"weights"} 或 \code{"curves"} 时），
#'   或一个 \pkg{patchwork} 复合拼接对象（当 \code{type = "both"} 时）。
#'
#' @seealso \code{\link{plot.nwqs_boot}} 以获取基于 Bootstrap 的具有真实统计学效力的诊断图。
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


#' @title 打印 NWQS 模型摘要 (基于数据拆分经验推断)
#'
#' @description
#' 打印非线性加权分位数和 (NWQS) 模型的简洁出版级摘要。
#' 自动格式化并输出整体混合物效应、特定组分的相对权重以及绝对效应对比（如最高分位数 Q4 vs 最低分位数 Q1）的 95\% 经验置信区间。
#'
#' @param x \code{"nwqs"} 类的对象。
#' @param digits Integer。需打印的有效数字位数。
#' @param ... 传递给其他方法的额外参数。
#'
#' @return 隐式返回原始的 \code{"nwqs"} 对象。
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


#' @title NWQS 模型的详细统计摘要
#'
#' @description
#' 提供 NWQS 模型的详尽统计摘要，包含全局 GLM 回归系数、经验标准误、Z 统计量、P 值以及偏差 (Deviance) 和 AIC 等拟合优度指标。
#'
#' @param object \code{"nwqs"} 类的对象。
#' @param digits Integer。保留的有效数字位数。
#' @param ... 传递给其他方法的额外参数。
#'
#' @return 隐式返回原始的 \code{"nwqs"} 对象。
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
# S3 Methods: nwqs_boot (正式的 Bootstrap 统计推断)
# ==============================================================================

#' @title 打印基于 Bootstrap 的 NWQS 模型正式推断结果
#'
#' @description
#' 打印 \code{nwqs_boot} 的摘要结果。此输出反映了模型真实参数的抽样变异 (Sampling Variance)，
#' 展示了所有项（整体效应与组分特定效应）在各分位数对比下的点估计值和外部 Bootstrap 百分位置信区间。
#' 针对流行病学研究，会自动将指数族模型的系数转换为比值比 (OR) 或相对危险度 (RR)。
#'
#' @param x \code{"nwqs_boot"} 类的对象。
#' @param digits Integer。有效数字位数，默认为 3。
#' @param ... 额外参数。
#'
#' @return 隐式返回 \code{x}。
#' @export
#' @method print nwqs_boot
print.nwqs_boot <- function(x, digits = 3, ...) {
  ci_tab <- x$ci_table

  if (is.null(ci_tab) || !is.data.frame(ci_tab) || nrow(ci_tab) == 0) {
    stop("`x$ci_table` is missing or empty.")
  }

  # 我们在新的 nwqs_boot 中已经直接导出了这些参数，无需再去猜
  n_success <- x$n_success
  n_total <- x$n_boot
  conf_pct <- if (!is.null(x$conf_level)) sprintf("%.0f%%", x$conf_level * 100) else "95%"

  # 注意这里直接从 x$family 读取，并加上了 clogit
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

  # 打印核心的汇总表格（包含平均权重和效应值）
  print(x$formatted_table, digits = digits, row.names = FALSE)

  # （注：删除了原本打印单一模型 AIC、Deviance 和底层系数表的代码，
  # 因为作为纯集成模型，这里展示的 formatted_table 就是全部的最终统计推断结果）

  invisible(x)
}



#' @title NWQS Bootstrap 集成模型的详细摘要
#'
#' @description
#' 提供 \code{nwqs_boot} 结果的详尽摘要，除了包含具备统计学效力的 Bootstrap 置信区间表外，
#' 还专门提供了跨 Bootstrap 样本的“组分权重稳定性指标 (Component Weight Stability)”，这对于评估
#' 高度共线性的混合物暴露在重抽样下的表现（如是否存在权重的不稳定跳转）具有重要科学价值。
#'
#' @param object \code{"nwqs_boot"} 类的对象。
#' @param digits Integer。有效数字位数，默认为 4。
#' @param ... 额外参数。
#'
#' @return 隐式返回 \code{object}。
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




#' @title 绘制 NWQS Bootstrap 诊断图 (权重与对比箱线图)
#'
#' @description
#' 为 \code{nwqs_boot} 对象生成推断级的诊断图。与点估计模型使用曲线图不同，Bootstrap 推断由于其集成特性，
#' 这里将剂量反应对比（如 \eqn{Q_i} vs \eqn{Q_1}）可视化为基于大量 Bootstrap 样本分布的箱线图 (Boxplots)。
#' 这种可视化能够极其直观地呈现真实样本变异下的置信度。
#'
#' @param x \code{"nwqs_boot"} 类的对象。
#' @param type Character。绘图类型：\code{"both"}（默认，权重图与对比箱图拼接）、\code{"curves"}（这里实际渲染为对比箱线图）或 \code{"weights"}（平均权重分布）。
#' @param y_scale Character。对比分布的尺度：\code{"partial"}（不进行指数转换的系数尺度）、\code{"contrast_or"}（指数转换后的 OR/RR 尺度）或 \code{"predicted"}（预测尺度）。
#' @param components Character vector。需要显示的特定混合物组分名称。若为 \code{NULL}，则显示所有组分。
#' @param overlay Logical。为兼容旧接口保留，在箱线图模式下不生效。
#' @param show_ci Logical。为兼容旧接口保留，在箱线图模式下不生效（箱线图自带分布信息）。
#' @param base_size Integer。基础字体大小，默认为 12。
#' @param palette Character。离散调色板。
#' @param top_n Integer 或 \code{NULL}。按权重降序排列，仅显示前 \code{top_n} 组分。
#' @param ylim Numeric vector。限制 Y 轴范围（如 \code{c(min, max)}），用于截断极端离群值。
#' @param y_step Numeric。Y 轴刻度步长。
#' @param free_y Logical。若为 \code{TRUE}，各个组分的 Y 轴尺度自由适配。
#' @param fill_alpha Numeric。箱线图的填充透明度，默认为 0.16。
#' @param exponentiate Logical 或 \code{NULL}。是否对 Y 轴效应量进行指数化（计算 OR/RR）。若为 \code{NULL}，则根据模型族自动推断。
#' @param ... 传递给 \code{plot_nwqs_contrast_box} 的额外参数。
#'
#' @return 返回 \code{ggplot} 对象或 \pkg{patchwork} 拼接图。
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

  # 为了兼容旧接口保留 y_scale：
  # partial -> 不指数化；其他 -> 对 exp family 默认指数化
  if (is.null(exponentiate)) {
    is_exp_family <- x$family %in% c("binomial", "poisson", "quasipoisson")
    exponentiate <- if (y_scale == "partial") FALSE else is_exp_family
  }

  # overlay / show_ci 仅为兼容旧接口保留，当前箱图版本不再使用
  is_combined <- (type == "both")

  selected_raw <- .resolve_selected_raw(x, components = components, top_n = top_n)
  selected_clean <- .clean_name(selected_raw)

  global_colors <- .get_palette(length(selected_clean), palette)
  names(global_colors) <- selected_clean

  # -- Weight bar chart --------------------------------------------------------
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

    if (length(w) == 0) stop("筛选后没有可绘制的权重")

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

  # -- Contrast box plot -------------------------------------------------------
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

  if (type == "weights") {
    return(p_w)
  }
  if (type == "curves") {
    return(p_c)
  }

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("请先安装 'patchwork'：install.packages('patchwork')")
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


#' @title 从 NWQS Bootstrap 对象中提取基础系数
#'
#' @description 注意：此函数提取的是在完整原始数据集上拟合的点估计基础系数，并非 Bootstrap 均值。
#' @param object \code{"nwqs_boot"} 类的对象。
#' @param ... 额外参数。
#' @return 包含全局模型基础系数的命名数值向量。
#' @export
#' @method coef nwqs_boot
coef.nwqs_boot <- function(object, ...) {
  object$point_fit$mean_coefs
}


