# ==============================================================================
# monte_carlo.R — Monte Carlo 模拟评估专用工具函数
# 从 utils.R 中分离，这些函数仅用于模拟研究，不属于 NWQS 包核心功能
# ==============================================================================

#' @title 计算权重分配的误差与相似度指标 (Calculate Weight Allocation Error Metrics)
#'
#' @description
#' 计算模型估计的混合物组分权重与数据生成时的真实（基准）权重之间的多种误差和相似度指标。
#' 这对于评估模型在存在多重共线性或未测量混杂时，精准识别关键毒性物质的能力至关重要。
#'
#' @details
#' \strong{指标解析:}
#' \itemize{
#'   \item \strong{SAE (绝对误差和):} $\sum |w_{est} - w_{true}|$。在成分数据 (Compositional Data) 框架下（权重和为1），SAE 是评估整体分配偏差的最直观指标。
#'   \item \strong{Cosine Similarity (余弦相似度):} 衡量高维权重向量在方向上的对齐程度，对权重的绝对数值缩放具有不变性。
#'   \item \strong{Spearman/Pearson:} 评估模型对暴露物相对重要性排序的恢复能力。
#' }
#'
#' @param w_est Numeric vector。模型估计出的混合物相对权重向量。
#' @param w_true Numeric vector。数据生成机制 (DGM) 设定的真实基准权重向量。必须与 \code{w_est} 长度一致。若两者均有命名属性，函数将自动按名称对齐。
#'
#' @return 返回一个包含以下元素的命名列表：
#'   \code{SAE}, \code{MAE}, \code{Pearson}, \code{Spearman}, \code{CosSim}。
#'
#' @export
calc_weight_error <- function(w_est, w_true) {
  if (length(w_est) != length(w_true)) stop("Lengths of estimated and true weights must match.")

  if (!is.null(names(w_est)) && !is.null(names(w_true))) {
    w_true <- w_true[names(w_est)]
  }

  error_diff <- w_est - w_true
  sae <- sum(abs(error_diff))
  mae <- mean(abs(error_diff))
  cor_pearson <- suppressWarnings(cor(w_est, w_true, method = "pearson"))
  cor_spearman <- suppressWarnings(cor(w_est, w_true, method = "spearman"))

  if (is.na(cor_pearson)) cor_pearson <- 0
  if (is.na(cor_spearman)) cor_spearman <- 0

  dot_prod <- sum(w_est * w_true)
  norm_est <- sqrt(sum(w_est^2))
  norm_true <- sqrt(sum(w_true^2))
  cos_sim <- if (norm_est > 0 && norm_true > 0) dot_prod / (norm_est * norm_true) else 0

  return(list(
    SAE = sae, MAE = mae,
    Pearson = cor_pearson, Spearman = cor_spearman, CosSim = cos_sim
  ))
}


#' @title 检验单次模拟的置信区间覆盖率 (Check Single-Simulation Coverage)
#'
#' @description
#' 将模型提取的暴露效应估计值（如特定分位数的对比效应）与基准真实效应矩阵进行比对，
#' 判定 Wald 置信区间和经验 (Empirical) 置信区间是否成功覆盖了真实值。
#'
#' @details
#' \strong{覆盖率与推断有效性:} \cr
#' 在流行病学推断中，如果一个 95% 置信区间是无偏的，那么在成千上万次 Monte Carlo 模拟中，
#' 它包含真实参数的比例应当趋近于 95%。此函数在单次模拟层面上打下布尔值标签 (\code{TRUE}/\code{FALSE})，
#' 为后续计算宏观的名义覆盖率 (Nominal Coverage Probability) 提供基础。
#'
#' @param est_df \code{data.frame}。由 \code{\link{extract_nwqs_effects}} 提取的效应估计结果表。
#' @param true_mat Matrix。真实的效应矩阵（通常从数据生成函数如 \code{\link{gen_nonlinear_data}} 的属性 \code{true_effect_mat} 中获取）。
#'
#' @return 返回一个 \code{data.frame}，包含对比目标 (\code{Target})、变量名 (\code{Term})、真实值、估计值、绝对偏差 (\code{Bias})，
#'   以及 Wald 和经验置信区间的上下限及其覆盖指示变量 (\code{Covered_Wald}, \code{Covered_Empirical})。
#'
#' @importFrom dplyr select arrange %>%
#' @export
check_coverage <- function(est_df, true_mat) {
  true_df <- as.data.frame(as.table(true_mat))
  colnames(true_df) <- c("Term", "Target", "True_Value")
  true_df$Term <- as.character(true_df$Term)
  true_df$Target <- as.character(true_df$Target)

  merged_df <- merge(est_df, true_df, by = c("Target", "Term"), all.x = TRUE)
  merged_df$Bias <- merged_df$Estimate - merged_df$True_Value

  merged_df$Covered_Wald <- (merged_df$True_Value >= merged_df$Wald_CI_Lower) &
    (merged_df$True_Value <= merged_df$Wald_CI_Upper)
  merged_df$Covered_Empirical <- (merged_df$True_Value >= merged_df$Empirical_CI_Lower) &
    (merged_df$True_Value <= merged_df$Empirical_CI_Upper)

  final_df <- merged_df %>%
    dplyr::select(
      Target, Term, True_Value, Estimate, Bias,
      Wald_CI_Lower, Wald_CI_Upper, Covered_Wald,
      Empirical_CI_Lower, Empirical_CI_Upper, Covered_Empirical
    ) %>%
    dplyr::arrange(Target, factor(Term, levels = c("Overall", setdiff(unique(Term), "Overall"))))

  return(final_df)
}


#' @title 评估 Monte Carlo 模拟的宏观统计性能 (Evaluate Macro-level Performance)
#'
#' @description
#' 跨越多次模拟迭代汇总结果，计算宏观层面的统计性能指标。涵盖区间覆盖概率、平均偏差、均方根误差 (RMSE)、
#' 统计效能 (Power) / 第一类错误率 (Type I Error)，以及变量选择的敏感性与特异性。
#'
#' @details
#' \strong{核心统计学评价维度:}
#' \itemize{
#'   \item \strong{Type I Error (假阳性率):} 当真实效应为 0 时，置信区间未覆盖 0 的模拟比例。控制在 $\alpha = 0.05$ 附近是模型稳健性的底线。
#'   \item \strong{Power (检验效能):} 当真实效应不为 0 时，模型成功拒绝零假设的比例。反映了模型在给定样本量和信噪比下发现真实关联的能力。
#'   \item \strong{敏感性/特异性 (Sensitivity/Specificity):} 评估权重分配的特征选择能力（基于设定的 \code{w_threshold}）。
#' }
#'
#' @param sim_weight_df \code{data.frame}。多次模拟权重的汇总表，每行代表一次模拟。
#' @param sim_effect_df \code{data.frame}。多次模拟效应估计的汇总表。
#' @param true_w Named numeric vector。基准真实权重向量。
#' @param true_eff_mat Matrix。数据生成时的真实效应矩阵。
#' @param w_threshold Numeric。将估计权重判定为“有效检测出信号”的硬阈值，默认为 0.01。
#'
#' @return 一个包含两部分的列表：
#' \describe{
#'   \item{\code{Weight_Metrics}}{权重分配的宏观性能（平均 SAE、敏感度、特异度）。}
#'   \item{\code{Effect_Metrics}}{按变量和对比目标汇总的效应估计性能（包含 Mean_Bias, RMSE, 覆盖概率 CP, 拒绝零假设率 Reject_H0 等）。}
#' }
#'
#' @export
evaluate_sim_performance <- function(sim_weight_df, sim_effect_df, true_w, true_eff_mat, w_threshold = 0.01) {
  true_toxics <- names(true_w)[true_w > 0]
  true_noises <- names(true_w)[true_w == 0]

  if (length(true_toxics) > 0) {
    toxic_est <- sim_weight_df[, paste0("NWQS_", true_toxics), drop = FALSE]
    sens <- mean(as.matrix(toxic_est) > w_threshold)
  } else {
    sens <- NA
  }

  if (length(true_noises) > 0) {
    noise_est <- sim_weight_df[, paste0("NWQS_", true_noises), drop = FALSE]
    spec <- mean(as.matrix(noise_est) <= w_threshold)
  } else {
    spec <- NA
  }

  mean_sae <- mean(sim_weight_df$NWQS_SAE, na.rm = TRUE)

  weight_res <- data.frame(
    Mean_SAE = mean_sae, Sensitivity = sens, Specificity = spec
  )

  true_eff_long <- as.data.frame(as.table(true_eff_mat))
  colnames(true_eff_long) <- c("Term", "Target", "True_Value")
  true_eff_long$Target <- as.character(true_eff_long$Target)
  true_eff_long$Term <- as.character(true_eff_long$Term)

  eval_df <- merge(sim_effect_df, true_eff_long, by = c("Target", "Term"))
  eval_df$Covered_Wald <- (eval_df$True_Value >= eval_df$Wald_CI_Lower) &
    (eval_df$True_Value <= eval_df$Wald_CI_Upper)
  eval_df$Covered_Empirical <- (eval_df$True_Value >= eval_df$Empirical_CI_Lower) &
    (eval_df$True_Value <= eval_df$Empirical_CI_Upper)
  eval_df$Abs_Bias <- eval_df$Estimate - eval_df$True_Value

  effect_res <- eval_df %>%
    dplyr::group_by(Target, Term) %>%
    dplyr::summarise(
      True_Value = mean(True_Value),
      Mean_Est = mean(Estimate),
      Mean_Bias = mean(Abs_Bias),
      RB_pct = ifelse(abs(mean(True_Value)) > 1e-5, mean(Abs_Bias) / abs(mean(True_Value)) * 100, NA),
      RMSE = sqrt(mean(Abs_Bias^2)),
      CP_Wald = mean(Covered_Wald, na.rm = TRUE),
      CP_Empirical = mean(Covered_Empirical, na.rm = TRUE),
      Reject_H0 = mean(Empirical_CI_Lower > 0 | Empirical_CI_Upper < 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Metric_Type = ifelse(abs(True_Value) > 1e-5, "Power", "Type I Error")
    ) %>%
    dplyr::arrange(Target, factor(Term, levels = c("Overall", names(true_w))))

  return(list(Weight_Metrics = weight_res, Effect_Metrics = effect_res))
}


#' @title 检验 Bootstrap 置信区间的覆盖率 (Check Bootstrap CI Coverage)
#'
#' @description
#' 专门用于评估经过 \code{\link{nwqs_boot}} 运行后得到的百分位 Bootstrap 置信区间是否覆盖了设定的真实因果效应值。
#' 相比于数据拆分产生的经验区间，外部 Bootstrap 区间理应提供更接近名义水平（如 95%）的真实抽样方差覆盖率。
#'
#' @param boot_res \code{"nwqs_boot"} 类的对象，包含内部的 \code{ci_table}。
#' @param true_value Numeric。用于对照的基准真实效应值。
#'
#' @return 返回原始的 \code{ci_table}，并增补了 \code{True_Value}, \code{Bias}, 及 \code{Covered_Bootstrap} 列。
#'
#' @export
check_boot_coverage <- function(boot_res, true_value) {
  if (is.null(boot_res$ci_table) || nrow(boot_res$ci_table) == 0) {
    stop("boot_res$ci_table is missing.")
  }

  out <- boot_res$ci_table
  out$True_Value <- true_value
  out$Bias <- out$Estimate - out$True_Value

  out$Covered_Bootstrap <- with(
    out,
    !is.na(Boot_CI_Lower) &
      !is.na(Boot_CI_Upper) &
      Boot_CI_Lower <= True_Value &
      Boot_CI_Upper >= True_Value
  )

  return(out)
}


#' @title 从真实效应矩阵推导相对重要性基准权重 (Derive True Importance Weights)
#'
#' @description
#' 将成分特定的绝对偏效应（从真实生成矩阵中提取）转换为归一化的、和为 1 的相对重要性权重。
#'
#' @details
#' \strong{非线性权重的合理定义:} \cr
#' 在纯线性模型中，权重通常可以直接基于真实斜率（系数）的绝对值分配（即 \code{method = "q4q1_abs"}）。
#' 但当涉及非单调剂量反应曲线（如 U 型、倒 U 型）时，最高分位数与最低分位数的对比（Q4 vs Q1）可能正好抵消为 0。
#' 为此，\code{method = "max_range"} 提取各个分位数节点处的效应极大值与极小值之差（极差），
#' 以此作为衡量该暴露组分在整个分布范围内对结局产生的最大变异贡献，从而科学地定义非线性场景下的“真实权重”。
#'
#' @param true_effect_mat Matrix。包含真实效应的矩阵。
#' @param mix_name Character vector。混合物组分的名称列表。
#' @param method Character。重要性推导方法。可选 \code{"q4q1_abs"}（仅适用单调线性场景）或 \code{"max_range"}（适用复杂非线性场景，默认）。
#'
#' @return 命名数值向量。代表各个组分的真实相对重要性基准，总和为 1。
#'
#' @export
calc_true_importance <- function(true_effect_mat, mix_name, method = "max_range") {
  
  if (method == "q4q1_abs") {
    if (!"Q4_vs_Q1" %in% colnames(true_effect_mat)) {
      stop("true_effect_mat must contain column 'Q4_vs_Q1'.")
    }
    contrib <- abs(true_effect_mat[mix_name, "Q4_vs_Q1"])
    
  } else if (method == "max_range") {
    # 提取所有对比 Q1 的列
    req_cols <- c("Q2_vs_Q1", "Q3_vs_Q1", "Q4_vs_Q1")
    if (!all(req_cols %in% colnames(true_effect_mat))) {
      stop("true_effect_mat must contain Q2_vs_Q1, Q3_vs_Q1, Q4_vs_Q1 columns.")
    }
    
    # 计算极差 (Max - Min, 记得要把参照组 Q1 的 0 包含进去)
    contrib <- sapply(mix_name, function(nm) {
      vals <- c(0, as.numeric(true_effect_mat[nm, req_cols]))
      max(vals, na.rm = TRUE) - min(vals, na.rm = TRUE)
    })
    
  } else {
    stop("Unsupported method.")
  }

  if (all(is.na(contrib)) || sum(contrib, na.rm = TRUE) == 0) {
    w_true <- rep(1 / length(mix_name), length(mix_name))
    names(w_true) <- mix_name
    return(w_true)
  }

  w_true <- contrib / sum(contrib, na.rm = TRUE)
  names(w_true) <- mix_name
  return(w_true)
}


#' @title 绘制 Monte Carlo 基准测试结果复合图 (Plot Monte Carlo Benchmark Results)
#'
#' @description
#' 为多模型基准测试 (Benchmarking) 生成达到顶级期刊发表标准的复合可视化图表面板。
#' 图表集成了雨云图 (Raincloud plots) 与分面箱线图，全景展示各候选模型的稳健性与偏差。
#'
#' @details
#' \strong{面板布局与学术解读:}
#' \itemize{
#'   \item \strong{面板 A (模型拟合误差 Deviance):} 通过雨云图展示各模型在多次模拟中的残余偏差分布。长尾分布往往暗示模型在某些模拟设置（如极端共线性或低信噪比）下发生崩溃。
#'   \item \strong{面板 B (权重提取误差 SAE):} 衡量模型降维和特征选择的准确性。分布越集中于底部（接近 0），表明模型锁定关键毒性物质的能力越稳定。
#'   \item \strong{面板 C (组分特异性权重恢复):} 以真实权重为水平参考线（虚线），展示每个模型在各个组分上的权重估计分布。能够直观揭示模型是否存在系统性的方向性偏倚（例如总是过度惩罚高方差变量）。
#' }
#' 默认采用了严谨的学术色板，并支持通过 \pkg{patchwork} 动态适应暴露组分的数量进行排版。
#'
#' @param dev_data \code{data.frame}。包含 \code{Model} 和 \code{Deviance} 列。
#' @param sae_data \code{data.frame}。包含 \code{Model} 和 \code{SAE} 列。
#' @param weight_data \code{data.frame}。包含 \code{Model}, \code{Component}, \code{Estimated_Weight}, 及 \code{True_Value} 列。
#' @param scen_name Character。模拟场景名称（用于总标题）。
#' @param family_name Character。GLM 分布族名称（用于总标题）。
#' @param custom_palette 命名字符向量。基于模型名称映射的十六进制颜色字典。若为 \code{NULL}，使用默认学术色板。
#' @param save_path Character 或 \code{NULL}。保存高分辨率图片的本地路径。若为 \code{NULL} 则仅在绘图设备中渲染。
#' @param base_size Numeric。基准字体大小，默认为 14。
#'
#' @return 隐式返回 \pkg{patchwork} 复合绘图对象。
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_point geom_hline facet_wrap scale_fill_manual scale_color_manual theme_bw labs theme element_text element_blank element_rect position_jitter position_nudge ggsave
#' @export
plot_monte_carlo_benchmark <- function(dev_data, sae_data, weight_data,
                                       scen_name = "Unknown Scenario",
                                       family_name = "Unknown Family",
                                       custom_palette = NULL, save_path = NULL, base_size = 14) {
  if (!requireNamespace("ggdist", quietly = TRUE)) stop("Please install 'ggdist' package.")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Please install 'patchwork' package.")

  if (is.null(custom_palette)) {
    custom_palette <- c(
      "NWQS" = "#4A90C8", "WQS" = "#D92828", "QGcomp" = "#6EC44A",
      "Ridge" = "#8B6FB8", "Lasso" = "#00B4D8", "ElasticNet" = "#006B3C",
      "RandomForest" = "#A8D8EA"
    )
  }

  model_levels <- names(custom_palette)

  dev_df <- as.data.frame(dev_data)
  dev_df$Model <- factor(dev_df$Model, levels = intersect(model_levels, unique(dev_df$Model)))
  sae_df <- as.data.frame(sae_data)
  sae_df$Model <- factor(sae_df$Model, levels = intersect(model_levels, unique(sae_df$Model)))
  weight_df <- as.data.frame(weight_data)
  weight_df$Model <- factor(weight_df$Model, levels = intersect(model_levels, unique(weight_df$Model)))

  n_comps <- length(unique(weight_df$Component))
  dynamic_nrow <- ceiling(n_comps / 7)
  dynamic_ncol <- ceiling(n_comps / dynamic_nrow)

  p_dev <- ggplot2::ggplot(dev_df, ggplot2::aes(x = Model, y = Deviance, fill = Model, color = Model)) +
    ggdist::stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, justification = -0.3, point_colour = NA, alpha = 0.7) +
    ggplot2::geom_boxplot(
      width = 0.2, outlier.shape = NA, alpha = 0.5, color = "black",
      position = ggplot2::position_nudge(x = -0.1)
    ) +
    ggplot2::geom_point(
      size = 1.3, alpha = 0.4,
      position = ggplot2::position_jitter(width = 0.05, height = 0)
    ) +
    ggplot2::scale_fill_manual(values = custom_palette) +
    ggplot2::scale_color_manual(values = custom_palette) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = "A. Model Fit Error (Deviance)",
      subtitle = "Lower deviance indicates better non-linear fit.",
      x = "", y = "Residual Deviance"
    )

  p_sae <- ggplot2::ggplot(sae_df, ggplot2::aes(x = Model, y = SAE, fill = Model, color = Model)) +
    ggdist::stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, justification = -0.3, point_colour = NA, alpha = 0.7) +
    ggplot2::geom_boxplot(
      width = 0.2, outlier.shape = NA, alpha = 0.5, color = "black",
      position = ggplot2::position_nudge(x = -0.1)
    ) +
    ggplot2::geom_point(
      size = 1.3, alpha = 0.4,
      position = ggplot2::position_jitter(width = 0.05, height = 0)
    ) +
    ggplot2::scale_fill_manual(values = custom_palette) +
    ggplot2::scale_color_manual(values = custom_palette) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = "B. Weight Extraction Error (SAE)",
      subtitle = "Lower SAE indicates higher accuracy.",
      x = "", y = "Sum of Absolute Errors (SAE)"
    )

  p_facet <- ggplot2::ggplot(weight_df, ggplot2::aes(x = Model, y = Estimated_Weight, fill = Model)) +
    ggplot2::geom_boxplot(alpha = 0.8, outlier.size = 0.5, color = "black", width = 0.6) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = True_Value), linetype = "dashed", color = "black", linewidth = 1) +
    ggplot2::facet_wrap(~Component, scales = "free_y", nrow = dynamic_nrow, ncol = dynamic_ncol) +
    ggplot2::scale_fill_manual(values = custom_palette) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom", legend.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#ECF0F1"),
      strip.text = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(
      title = "C. Component-Specific Weight Recovery Accuracy",
      x = "", y = "Estimated Relative Weight"
    )

  dynamic_title <- sprintf(
    "Monte Carlo Simulation Benchmark\nScenario: %s | Family: %s",
    scen_name, toupper(family_name)
  )

  dynamic_height_ratio <- 0.8 * dynamic_nrow
  final_plot <- (p_dev | p_sae) / p_facet +
    patchwork::plot_layout(heights = c(1.2, dynamic_height_ratio)) +
    patchwork::plot_annotation(
      title = dynamic_title,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = base_size + 4, face = "bold", hjust = 0.5)
      )
    )

  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path,
      plot = final_plot, width = 16,
      height = 11 + (dynamic_nrow - 1) * 3.5, dpi = 500
    )
  }

  return(final_plot)
}
