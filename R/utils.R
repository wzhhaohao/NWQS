# ==============================================================================
# utils.R — NWQS 模型核心工具函数（清理版）
# ==============================================================================

#' @title 分位数或百分位数数据转换 (Quantile or Percentile Transformation)
#'
#' @description
#' 将连续的混合物暴露变量转换为离散的分位数区间 (Quantile Bins) 或连续的百分位秩 (Percentile Ranks)。
#' 这一步是加权分位数和 (WQS) 及其扩展方法标准化不同量纲暴露特征的基石。
#'
#' @details
#' \strong{方法学考量与潜在偏倚风险:}
#' \describe{
#'   \item{\code{"quantile"}}{经典的 gWQS 分箱逻辑，将数据划分为 \code{q} 个组（如 0 到 \code{q-1}）。
#'     内部通过指定 \code{-Inf} 和 \code{Inf} 作为边界，极大地增强了对极端异常值 (Outliers) 的稳健性。
#'     *流行病学警示:* 虽然分箱能抵抗异常值，但如果真实的剂量反应关系在某个分箱内部存在急剧的非线性变化，
#'     分类操作可能引入一定程度的残余混杂 (Residual Confounding)。此外，分位点的确定高度依赖于当前样本，
#'     若样本存在选择偏倚 (Selection Bias)，则截断点可能无法代表目标人群的真实暴露水平。}
#'   \item{\code{"percentile"}}{通过 \eqn{rank(x) / (n + 1)} 将连续变量映射到 (0, 1) 区间。保留了比分类更多的
#'     等级信息，适用于样本量较小或需要精细平滑非线性曲线的场景。}
#' }
#'
#' @param data \code{data.frame}。包含需要转换的混合物变量的数据框。
#' @param method Character。转换方法：\code{"quantile"}（默认分位数）或 \code{"percentile"}（百分位秩）。
#' @param q Integer。分位数的分箱数量（仅在 \code{method = "quantile"} 时生效）。默认为 4（四分位数）。
#'
#' @return 返回一个与输入同维度和同列名的 \code{data.frame}，包含转换后的无量纲数值。
#'
#' @importFrom stats quantile
#' @export
trans_quantile <- function(data, method = c("quantile", "percentile"), q = 4) {
  data <- as.data.frame(data)
  method <- match.arg(method)

  if (method == "percentile") {
    transform_func <- function(x) {
      rank(x) / (length(x) + 1)
    }
    res_list <- lapply(data, transform_func)
  } else {
    if (!is.numeric(q) || q < 1) stop("'q' must be a positive number")

    transform_func <- function(x) {
      breaks <- unique(quantile(x, probs = seq(0, 1, by = 1 / q), na.rm = TRUE))
      if (length(breaks) == 1) {
        breaks <- c(-Inf, breaks)
      } else {
        breaks[1] <- -Inf
        breaks[length(breaks)] <- Inf
      }
      as.numeric(cut(x, breaks = breaks, labels = FALSE, include.lowest = TRUE)) - 1
    }
    res_list <- lapply(data, transform_func)
  }

  res_df <- as.data.frame(res_list)
  names(res_df) <- names(data)
  return(res_df)
}


# -------------------------------------------------------------------------

#' @title WQS 混合物组分的非线性样条扩展
#'
#' @description
#' 使用全局固定的内部节点 (Knots) 和边界节点 (Boundary Knots)，将（已完成分位数转换的）混合物变量
#' 转换为自然三次样条 (Natural Cubic Splines) 基矩阵。
#'
#' @details
#' \strong{严谨性说明 (基函数对齐):} \cr
#' 在机器学习和统计建模结合的架构中（如 NWQS 采用的 Repeated Holdout 或 Bootstrap），
#' 绝不能在不同的子样本中动态计算样条节点。此函数强制要求传入预先在全样本上确定的 \code{knots}
#' 和 \code{boundary}，这确保了训练集、验证集和重抽样集之间基函数的绝对对齐。若不进行这种全局约束，
#' 极易引发空间漂移偏倚 (Spatial Drift Bias)，导致交叉验证中提取的形状系数在验证集上失效。
#'
#' @param data \code{data.frame}。包含混合物变量（通常已经过转换）的数据框。
#' @param mix_name Character vector。需要进行样条扩展的混合物组分列名。
#' @param df_spline Integer。自然样条的自由度，默认为 3。
#' @param knots Numeric vector。内部节点位置。为了确保全局尺度对齐，必须显式提供。
#' @param boundary Numeric vector (长度为 2)。边界节点位置。必须显式提供以固定外推边界。
#'
#' @return 一个数值型矩阵，列名格式为 \code{{Component}_B{BasisIndex}}
#'   （例如，\code{Component1_B1}, \code{Component1_B2}）。
#'
#' @importFrom splines ns
#' @export
wqs_nonlinear_expand <- function(data, mix_name, df_spline = 3,
                                 knots = NULL, boundary = NULL) {
  trans_data <- data[, mix_name, drop = FALSE]

  if (is.null(knots) || is.null(boundary)) {
    stop("'knots' and 'boundary' must be provided to ensure global scale alignment.")
  }

  mat_spline_list <- lapply(trans_data, function(x) {
    splines::ns(x, df = df_spline, knots = knots, Boundary.knots = boundary)
  })

  mat_spline_full <- do.call(cbind, mat_spline_list)
  cols_per_mix <- ncol(mat_spline_list[[1]])

  colnames(mat_spline_full) <- paste0(
    rep(mix_name, each = cols_per_mix),
    "_B", rep(seq_len(cols_per_mix), times = length(mix_name))
  )

  return(mat_spline_full)
}


# -------------------------------------------------------------------------

#' @title 计算 NWQS 联合暴露分位数对比效应
#'
#' @description
#' 计算整体混合物非线性效应的全局显著性，并评估联合暴露的分位数对比
#' （例如：所有暴露组分同时处于最高分位数 Q4 相较于均处于最低分位数 Q1 时的效应变化）。
#' 对于逻辑回归和泊松回归，函数会自动将其转换为临床易解的比值比 (OR) 和相对危险度 (RR)。
#'
#' @param model \code{"nwqs"} 类的对象。
#' @param q_target Integer。目标分位数索引（基于 0 起始）。例如，3 代表 Q4。若为 \code{NULL}，将自动推断为最大分位数。
#' @param q_ref Integer。参考分位数索引（基于 0 起始）。默认为 0 (Q1)。
#'
#' @return 隐式返回一个列表，包含计算的绝对偏效应变化 \code{delta_eta} 以及置信区间上下限 \code{lower} 和 \code{upper}。
#'   （注意：当 \code{rh = 1} 时，由于无法估计变异，上下限返回 \code{NA}）。
#'
#' @importFrom splines ns
#' @importFrom stats quantile coef
#' @export
nwqs_contrast <- function(model, q_target = NULL, q_ref = 0) {
  # 【修改核心】：动态获取模型的最大分位数，不再写死 Q4(3)
  if (is.null(q_target)) {
    q_target <- if (!is.null(model$q)) model$q - 1 else 3
  }

  is_binomial <- FALSE
  if (inherits(model, "glm")) {
    is_binomial <- model$family$family == "binomial"
    rh <- 1
  } else {
    is_binomial <- model$family == "binomial"
    rh <- model$rh
  }

  cat("\n======================================================\n")
  cat("      NWQS Overall Mixture Effect Significance\n")
  cat("======================================================\n")

  if (rh == 1) {
    wqs_sum <- summary(model)$coefficients["nwqs", , drop = FALSE]
    print(round(wqs_sum, 4))
    if (wqs_sum[1, 4] < 0.05) {
      cat("\nConclusion: The overall NWQS latent risk score has a significant joint effect on the outcome (P < 0.05).\n")
    } else {
      cat("\nConclusion: The overall NWQS latent risk score does not have a significant joint effect on the outcome (P >= 0.05).\n")
    }
  } else {
    wqs_betas <- model$rh_coefs[, "nwqs"]
    wqs_betas <- wqs_betas[is.finite(wqs_betas)]

    if (length(wqs_betas) == 0) {
      stop("All RH 'nwqs' coefficients are NA/NaN/Inf; cannot compute overall mixture significance.")
    }

    mean_beta <- mean(wqs_betas, na.rm = TRUE)
    ci_lower <- quantile(wqs_betas, 0.025, na.rm = TRUE)
    ci_upper <- quantile(wqs_betas, 0.975, na.rm = TRUE)
    p_val_emp <- 2 * min(
      mean(wqs_betas > 0, na.rm = TRUE),
      mean(wqs_betas < 0, na.rm = TRUE)
    )

    cat(sprintf("Mean Beta across RH iterations :  %.4f\n", mean_beta))
    cat(sprintf("Empirical 95%% CI             : [%.4f, %.4f]\n", ci_lower, ci_upper))
    cat(sprintf("Empirical P-value            :  %.4f\n", p_val_emp))

    if (isTRUE(ci_lower > 0 || ci_upper < 0)) {
      cat("\nConclusion: The overall joint mixture effect is significant (95% CI excludes 0).\n")
    } else {
      cat("\nConclusion: The overall joint mixture effect is NOT significant (95% CI includes 0).\n")
    }
  }

  cat("\n======================================================\n")
  cat(sprintf(" Joint Exposure Quantile Contrast: Target Q%d vs. Ref Q%d\n", q_target + 1, q_ref + 1))
  cat("======================================================\n")

  if (rh == 1) {
    shapes_vec <- model$mean_shapes
    weights_vec <- model$final_weights
    beta_wqs <- coef(model)["nwqs"]
  } else {
    shapes_vec <- model$mean_shapes
    weights_vec <- model$final_weights
    beta_wqs <- model$mean_coefs["nwqs"]
  }

  df_spline <- max(as.numeric(sub("^.+_B(\\d+)$", "\\1", names(shapes_vec))))
  b_target <- splines::ns(c(q_target, q_ref),
    df = df_spline,
    knots = model$spline_knots, Boundary.knots = model$spline_boundary,
    intercept = FALSE
  )

  calc_diff <- function(b_mat, w_vec, s_vec) {
    diff_val <- 0
    for (comp in names(w_vec)) {
      comp_cols <- paste0(comp, "_B", 1:df_spline)
      theta <- matrix(s_vec[comp_cols], ncol = 1)
      w <- w_vec[comp]
      s_tgt <- b_mat[1, , drop = FALSE] %*% theta
      s_ref <- b_mat[2, , drop = FALSE] %*% theta
      diff_val <- diff_val + (s_tgt - s_ref) * w
    }
    return(as.numeric(diff_val))
  }

  if (rh > 1) {
    diff_list <- numeric(rh)
    for (i in 1:rh) {
      beta_i <- model$rh_coefs[i, "nwqs"]
      shape_i <- model$rh_shapes[i, ]
      weight_i <- model$rh_weights[i, ]
      score_diff_i <- calc_diff(b_target, weight_i, shape_i)
      diff_list[i] <- score_diff_i * beta_i
    }
    delta_eta <- mean(diff_list)
    ci_lower_eta <- quantile(diff_list, 0.025)
    ci_upper_eta <- quantile(diff_list, 0.975)
  } else {
    score_diff <- calc_diff(b_target, weights_vec, shapes_vec)
    delta_eta <- score_diff * beta_wqs
    ci_lower_eta <- NA
    ci_upper_eta <- NA
  }

  cat(sprintf("Absolute Partial Effect Change (\u0394 Eta) :  %.4f\n", delta_eta))
  if (rh > 1) {
    cat(sprintf("95%% CI (\u0394 Eta)                      : [%.4f, %.4f]\n", ci_lower_eta, ci_upper_eta))
  }

  is_rate_family <- if (inherits(model, "glm")) {
    model$family$family %in% c("poisson", "quasipoisson")
  } else {
    model$family %in% c("poisson", "quasipoisson")
  }

  if (is_binomial) {
    cat("\n----------------- Converted to Odds Ratio (OR) -----------------\n")
    cat(sprintf("Overall Joint OR :  %.4f\n", exp(delta_eta)))
    if (rh > 1) {
      cat(sprintf("95%% CI (OR)      : [%.4f, %.4f]\n", exp(ci_lower_eta), exp(ci_upper_eta)))
    }
    cat(sprintf(
      "\n[Interpretation]: When all components at Q%d vs Q%d, OR = %.4f.\n",
      q_target + 1, q_ref + 1, exp(delta_eta)
    ))
  } else if (is_rate_family) {
    cat("\n----------------- Converted to Rate Ratio (RR) -----------------\n")
    cat(sprintf("Overall Joint RR :  %.4f\n", exp(delta_eta)))
    if (rh > 1) {
      cat(sprintf("95%% CI (RR)      : [%.4f, %.4f]\n", exp(ci_lower_eta), exp(ci_upper_eta)))
    }
    cat(sprintf(
      "\n[Interpretation]: When all components at Q%d vs Q%d, RR = %.4f.\n",
      q_target + 1, q_ref + 1, exp(delta_eta)
    ))
  } else {
    cat("\n[Interpretation]: Absolute predicted increment on the response scale.\n")
  }

  invisible(list(delta_eta = delta_eta, lower = ci_lower_eta, upper = ci_upper_eta))
}


# -------------------------------------------------------------------------

#' @title 提取 NWQS 模型的详细分位数对比效应
#'
#' @description
#' 计算所有可能的分位数对比效应（例如 Q2 vs Q1，Q3 vs Q1，Q4 vs Q1）。
#' 函数不仅会输出整体混合物的效应，还会将整体效应分解至每个单独的混合物组分，
#' 并从 Repeated Holdout 迭代中提取标准误 (SE) 和经验置信区间。
#'
#' @param model_res \code{"nwqs"} 类的对象。
#' @param return_raw Logical。当前未使用，为未来扩展保留。
#'
#' @return 一个 \code{data.frame}，包含以下列：Target (对比目标), Term (变量名/整体),
#'   Estimate (估计值), SE (经验标准误), Wald_CI_Lower, Wald_CI_Upper,
#'   Empirical_CI_Lower, Empirical_CI_Upper。
#'
#' @importFrom splines ns
#' @export
extract_nwqs_effects <- function(model_res, return_raw = FALSE) {
  q_level <- if (!is.null(model_res$q)) model_res$q else 4
  df_spline <- model_res$df_spline
  comps <- colnames(model_res$rh_weights)
  rh <- model_res$rh
  model_knots <- model_res$spline_knots
  model_boundary <- model_res$spline_boundary
  eval_points_std <- 0:(q_level - 1)

  basis_std <- splines::ns(eval_points_std,
    df = df_spline, knots = model_knots, Boundary.knots = model_boundary, intercept = FALSE
  )

  res_list <- list()

  for (q_tgt in 2:q_level) {
    b_diff <- basis_std[q_tgt, ] - basis_std[1, ]
    iter_effects <- matrix(0, nrow = rh, ncol = length(comps) + 1)
    colnames(iter_effects) <- c("Overall", comps)

    iter_effects <- matrix(NA_real_, nrow = rh, ncol = length(comps) + 1)
    colnames(iter_effects) <- c("Overall", comps)

    for (i in seq_len(rh)) {
      beta_i <- model_res$rh_coefs[i, "nwqs"]
      if (!is.finite(beta_i)) next
      comp_effs_i <- numeric(length(comps))
      names(comp_effs_i) <- comps

      for (comp in comps) {
        theta_cols <- paste0(comp, "_B", 1:df_spline)
        theta_i <- model_res$rh_shapes[i, theta_cols]
        w_i <- model_res$rh_weights[i, comp]
        comp_effs_i[comp] <- beta_i * w_i * sum(b_diff * theta_i)
      }

      iter_effects[i, "Overall"] <- sum(comp_effs_i)
      iter_effects[i, comps] <- comp_effs_i
    }

    est_vec <- colMeans(iter_effects, na.rm = TRUE)
    se_vec <- apply(iter_effects, 2, sd, na.rm = TRUE)
    emp_ci <- apply(iter_effects, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

    res_list[[q_tgt - 1]] <- data.frame(
      Target = paste0("Q", q_tgt, "_vs_Q1"),
      Term = names(est_vec),
      Estimate = est_vec,
      SE = se_vec,
      Wald_CI_Lower = est_vec - 1.96 * se_vec,
      Wald_CI_Upper = est_vec + 1.96 * se_vec,
      Empirical_CI_Lower = emp_ci[1, ],
      Empirical_CI_Upper = emp_ci[2, ],
      stringsAsFactors = FALSE
    )
  }

  final_df <- do.call(rbind, res_list)
  rownames(final_df) <- NULL
  return(final_df)
}


#' @title 绘制 NWQS 组分特异性对比箱线图 (Bootstrap 分布图)
#'
#' @description
#' 专为医学与流行病学顶级期刊设计的高质量诊断图。此函数不仅展示效应的点估计，
#' 而是通过箱线图加抖动点 (Jitter) 的形式，全景展示各暴露组分在不同 Bootstrap 样本中的效应量分布。
#'
#' @details
#' \strong{可视化透明度与科学价值:} \cr
#' 环境暴露数据通常伴随极高的多重共线性 (Multicollinearity)。单纯依赖点估计与标准误容易掩盖
#' 权重分配的不稳定性。通过将 Bootstrap 集成过程中的每一次迭代可视化：
#' \itemize{
#'   \item 若箱体极宽或分布呈现双峰，高度提示该组分的效应受到共线性干扰或严重依赖特定的极端抽样 (Extreme Samples)。
#'   \item 图中的横线代表零效应（线性回归为 0，指数族回归为 1）。只有当箱体的绝大部分（如 95% 区间）
#'   远离该基线时，我们才能对该特定组分的相对贡献抱有高度的统计信心。
#' }
#' 绘图严格遵循 \pkg{ggplot2} 的无冗余设计哲学（清晰的主题、明确的坐标轴标签及对比度良好的颜色比例）。
#'
#' @param model \code{"nwqs_boot"} 类的对象。
#' @param exponentiate Logical 或 \code{NULL}。是否对 Y 轴效应量进行指数化以显示为 OR 或 RR。若为 \code{NULL}，将自动根据误差分布族判断。
#' @param free_y Logical。若为 \code{TRUE}（默认），各个组分分面的 Y 轴尺度将自由适配，方便观察效应极小的组分分布。
#' @param base_size Integer。图表基础字体大小，默认 12，适合大多数学术出版的 A4 排版。
#' @param fill_alpha Numeric。箱体填充的透明度，默认为 0.16，以确保底层抖动点清晰可见。
#' @param palette Character。离散调色板，默认为 \code{"default"}。
#' @param components Character vector。指定需要绘制的特定组分。
#' @param top_n Integer 或 \code{NULL}。按权重降序排列，仅显示前 \code{top_n} 个最重要的组分。
#' @param ylim Numeric vector (长度为 2)。强行限制 Y 轴范围（如截断极端 Bootstrap 离群值影响视觉比例时使用）。
#' @param y_step Numeric。强制指定 Y 轴的刻度间距。
#'
#' @return 一个具备学术出版质量的 \code{ggplot} 分面箱线图对象。
#'
#' @export
#' @method plot nwqs_boot
plot_nwqs_contrast_box <- function(model,
                                   exponentiate = NULL,
                                   free_y = TRUE,
                                   base_size = 12,
                                   fill_alpha = 0.16,
                                   palette = "default",
                                   components = NULL,
                                   top_n = NULL,
                                   ylim = NULL,
                                   y_step = NULL) {
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
        "#006B3C", "#A8D8EA", "#F4B6B6", "#5BA3D0", "#E03030",
        "#7AD450", "#9B7FC0"
      ),
      palette2 = c(
        "#9bbf8a", "#82afda", "#f79059", "#e7dbd3", "#c2bdde",
        "#8dcec8", "#add3e2", "#3480b8", "#ffbe7a", "#fa8878",
        "#c82423", "#6b5b95"
      )
    )
    p <- if (pal %in% names(cols)) cols[[pal]] else cols[["default"]]
    rep(p, ceiling(n / length(p)))[seq_len(n)]
  }

  calc_facet_layout <- function(n) {
    if (n <= 5) {
      return(list(nc = n))
    }
    if (n == 8) {
      return(list(nc = 4))
    }
    if (n == 9) {
      return(list(nc = 3))
    }
    if (n == 10) {
      return(list(nc = 5))
    }
    if (n == 12) {
      return(list(nc = 4))
    }
    sqrt_n <- ceiling(sqrt(n))
    list(nc = min(5, sqrt_n))
  }

  .resolve_selected_raw <- function(model, components = NULL, top_n = NULL) {
    if (!is.null(model$final_weights) && length(model$final_weights) > 0) {
      selected <- names(sort(model$final_weights, decreasing = TRUE))
    } else {
      selected <- unique(setdiff(model$boot_table$Term, "Overall"))
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

  if (!inherits(model, "nwqs_boot")) {
    stop("需传入 'nwqs_boot' 对象")
  }
  if (is.null(model$boot_table) || nrow(model$boot_table) == 0) {
    stop("model$boot_table 为空，无法绘图")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("请先安装 dplyr：install.packages('dplyr')")
  }

  selected_raw <- .resolve_selected_raw(model, components = components, top_n = top_n)
  if (length(selected_raw) == 0) {
    stop("筛选后没有可绘制的组分")
  }

  boot_df <- model$boot_table
  boot_df <- boot_df[boot_df$Term %in% c("Overall", selected_raw), , drop = FALSE]

  plot_df <- data.frame(
    Iteration = boot_df$Boot_ID,
    Component_Raw = boot_df$Term,
    Quantile = sub("_vs_.*", "", boot_df$Target),
    Effect = boot_df$Estimate,
    stringsAsFactors = FALSE
  )

  selected_clean <- .clean_name(selected_raw)

  plot_df$Component <- ifelse(
    plot_df$Component_Raw == "Overall",
    "Overall",
    .clean_name(plot_df$Component_Raw)
  )
  plot_df$Component <- factor(
    plot_df$Component,
    levels = c("Overall", selected_clean)
  )

  q_levels_str <- paste0("Q", 2:(if (!is.null(model$q)) model$q else 4))
  plot_df$Quantile <- factor(plot_df$Quantile, levels = q_levels_str)

  is_exp_family <- model$family %in% c("binomial", "poisson", "quasipoisson")
  if (is.null(exponentiate)) exponentiate <- is_exp_family

  if (exponentiate) {
    plot_df$Effect <- exp(plot_df$Effect)
    y_label <- "Odds Ratio / Risk Ratio (95% Percentile CI)"
    y_intercept <- 1
  } else {
    y_label <- "Coefficient (\u03B2 with 95% Percentile CI)"
    y_intercept <- 0
  }

  stats_df <- plot_df %>%
    dplyr::group_by(Component, Quantile) %>%
    dplyr::summarise(
      ymin = quantile(Effect, 0.025, na.rm = TRUE),
      lower = quantile(Effect, 0.25, na.rm = TRUE),
      middle = median(Effect, na.rm = TRUE),
      upper = quantile(Effect, 0.75, na.rm = TRUE),
      ymax = quantile(Effect, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  final_colors <- c(
    "Overall" = "#555555",
    stats::setNames(.get_palette(length(selected_clean), palette), selected_clean)
  )
  final_colors <- final_colors[c("Overall", selected_clean)]

  layout_config <- calc_facet_layout(length(c("Overall", selected_clean)))

  p <- ggplot2::ggplot() +
    ggplot2::geom_jitter(
      data = plot_df,
      ggplot2::aes(x = Quantile, y = Effect, fill = Component),
      shape = 21,
      color = "white",
      alpha = 0.35,
      width = 0.22,
      size = 0.85,
      stroke = 0.15
    ) +
    ggplot2::geom_boxplot(
      data = stats_df,
      ggplot2::aes(
        x = Quantile,
        ymin = ymin,
        lower = lower,
        middle = middle,
        upper = upper,
        ymax = ymax,
        fill = Component
      ),
      stat = "identity",
      color = "black",
      linewidth = 0.6,
      width = 0.55,
      alpha = fill_alpha
    ) +
    ggplot2::geom_hline(
      yintercept = y_intercept,
      linetype = "dashed",
      color = "firebrick",
      linewidth = 0.6
    ) +
    ggplot2::facet_wrap(
      ~Component,
      scales = ifelse(free_y, "free_y", "fixed"),
      ncol = layout_config$nc,
      axes = "all_x",
      axis.labels = "all_x"
    ) +
    ggplot2::scale_fill_manual(values = final_colors, drop = FALSE) +
    ggplot2::labs(
      title = "Bootstrap Specificity and Stability Analysis",
      subtitle = sprintf(
        "Boxplots: Median [IQR]; Whiskers: 2.5th to 97.5th Percentiles (95%% CI); n=%d",
        model$n_success
      ),
      x = "Exposure Quantile Index",
      y = y_label
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold",
        color = "black",
        size = base_size
      ),
      axis.line = ggplot2::element_line(
        color = "black",
        linewidth = 0.5
      ),
      axis.text.x = ggplot2::element_text(color = "black"),
      axis.ticks.x = ggplot2::element_line(color = "black"),
      axis.title.x = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5,
        size = base_size + 2
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        color = "gray30",
        size = base_size - 1
      ),
      panel.spacing = ggplot2::unit(1.2, "lines")
    )

  if (!is.null(y_step) && !is.null(ylim)) {
    y_breaks <- seq(ylim[1], ylim[2], by = y_step)
    p <- p + ggplot2::scale_y_continuous(breaks = y_breaks)
  } else if (!is.null(y_step)) {
    p <- p + ggplot2::scale_y_continuous(breaks = scales::breaks_width(y_step))
  }

  if (!is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(ylim = ylim)
  }

  p
}
