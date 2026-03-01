#' Quantile or Percentile Transformation
#'
#' @description
#' A hybrid function combining flexible ranking methods:
#' 1. "quantile": gWQS-style integer binning (handles ties and boundaries robustly).
#' 2. "percentile": Continuous percentile ranking (0 to 1).
#' 这是一个融合函数，结合了灵活的排名方法：
#' 1. "quantile": gWQS 风格的整数分箱（稳健处理重复值和边界）。
#' 2. "percentile": 连续百分位数排名（0 到 1）。
#'
#' @param data data.frame. Input data. / 输入数据。
#' @param method character. "quantile" (default) or "percentile". / 变换方法。
#' @param q integer. Number of quantiles (only used if method = "quantile"). / 分位数数量。
#' @return data.frame. Transformed data. / 变换后的数据。
#' @importFrom splines ns
#' @importFrom stats quantile
#' @export
trans_quantile <- function(data, method = c("quantile", "percentile"), q = 4) {
  data <- as.data.frame(data)
  method <- match.arg(method)

  if (method == "percentile") {
    transform_func <- function(x) {
      rank(x) / (length(x) + 1)
    }

    # Apply to all columns
    res_list <- lapply(data, transform_func)
  } else {
    if (!is.numeric(q) || q < 1) stop("'q' must be a positive number")

    transform_func <- function(x) {
      breaks <- unique(quantile(x, probs = seq(0, 1, by = 1 / q), na.rm = TRUE))

      # Handle Boundaries (gWQS style: Force -Inf / Inf for robustness)
      if (length(breaks) == 1) {
        breaks <- c(-Inf, breaks)
      } else {
        breaks[1] <- -Inf
        breaks[length(breaks)] <- Inf
      }

      # Cut (Binning): Returns integer values from 0 to q-1
      as.numeric(cut(x, breaks = breaks, labels = FALSE, include.lowest = TRUE)) - 1
    }

    # Apply to all columns
    res_list <- lapply(data, transform_func)
  }

  # Convert list back to data.frame and preserve names
  res_df <- as.data.frame(res_list)
  names(res_df) <- names(data)

  return(res_df)
}

# -------------------------------------------------------------------------
#' Nonlinear Expansion for WQS (Natural Splines) / WQS 非线性展开 (自然样条)
#'
#' @description
#' Transforms mixture variables into natural cubic spline bases to capture nonlinear effects.
#' 将混合物变量转换为自然三次样条基函数，以捕捉非线性效应。
#'
#' @details
#' By default, this function performs a quantile transformation (quartiles) before spline expansion
#' if no `transform_fun` is provided. It uses `splines::ns` for the basis expansion.
#' 默认情况下，如果未提供 `transform_fun`，该函数会在样条展开前执行四分位数转换。
#' 它使用 `splines::ns` 生成样条基底。
#'
#' @param data data.frame. The dataset containing the mixture variables.
#'   包含混合物变量的数据集。
#' @param mix_name character vector. Names of the mixture components to be expanded.
#'   需要展开的混合物组分名称。
#' @param transform_fun function. Optional custom transformation function applied before spline expansion.
#'   If NULL, applies a default quantile transformation (q=4).
#'   可选的自定义转换函数。如果为 NULL，则应用默认的四分位数转换。
#' @param df_spline integer. Degrees of freedom for the natural spline. Default is 3.
#'   自然样条的自由度。默认为 3。
#'
#' @return matrix. A matrix containing the spline basis functions for all mixture components.
#'   Column names are formatted as `{Component}_B{BasisIndex}`.
#'   返回包含所有混合物组分样条基函数的矩阵。
#' @export
wqs_nonlinear_expand <- function(data, mix_name, df_spline = 3, knots = NULL, boundary = NULL) {
  trans_data <- data[, mix_name]

  # 如果没传尺子（防呆设计），就降级回老办法
  if (is.null(knots) || is.null(boundary)) {
    stop("Error: 'knots' and 'boundary' must be provided to ensure global scale alignment.")
  }

  mat_spline_list <- lapply(trans_data, function(x) {
    splines::ns(x, df = df_spline, knots = knots, Boundary.knots = boundary)
  })
  mat_spline_full <- do.call(cbind, mat_spline_list)

  # 生成列名
  total_cols <- ncol(mat_spline_full)
  cols_per_mix <- total_cols / length(mix_name)

  colnames(mat_spline_full) <- paste0(
    rep(mix_name, each = cols_per_mix),
    "_B",
    rep(1:cols_per_mix, times = length(mix_name))
  )

  return(mat_spline_full)
}






#' Calculate NWQS Joint Exposure Quantile Contrast and Overall Significance
#'
#' @description
#' Computes the overall significance of the non-linear mixture effect and evaluates
#' the joint exposure quantile contrast (e.g., comparing all exposures at Q4 vs Q1).
#' Automatically converts to Odds Ratios (OR) for binomial models.
#'
#' @param model An object of class "nwqs" or "nwqs_result".
#' @param q_target Integer. Target quantile index (e.g., 3 represents the 4th quartile, Q4).
#' @param q_ref Integer. Reference quantile index (default is 0, representing Q1).
#' @export
#' @importFrom splines ns
#' @importFrom stats quantile coef
nwqs_contrast <- function(model, q_target = 3, q_ref = 0) {
  # 1. Identify model type (continuous vs. binomial) and RH iterations
  is_binomial <- FALSE
  if (inherits(model, "glm")) {
    is_binomial <- model$family$family == "binomial"
    rh <- 1
  } else {
    is_binomial <- model$family == "binomial"
    rh <- model$rh
  }

  # ==========================================================
  # Module 1: Extract and print overall WQS effect significance
  # ==========================================================
  cat("\n======================================================\n")
  cat("      NWQS Overall Mixture Effect Significance\n")
  cat("======================================================\n")

  if (rh == 1) {
    wqs_sum <- summary(model)$coefficients["wqs_score", , drop = FALSE]
    print(round(wqs_sum, 4))
    if (wqs_sum[1, 4] < 0.05) {
      cat("\nConclusion: The overall NWQS latent risk score has a significant joint effect on the outcome (P < 0.05).\n")
    } else {
      cat("\nConclusion: The overall NWQS latent risk score does not have a significant joint effect on the outcome (P >= 0.05).\n")
    }
  } else {
    wqs_betas <- model$rh_coefs[, "wqs_score"]
    mean_beta <- mean(wqs_betas, na.rm = TRUE)
    ci_lower <- quantile(wqs_betas, 0.025, na.rm = TRUE)
    ci_upper <- quantile(wqs_betas, 0.975, na.rm = TRUE)

    # Calculate Empirical P-value based on RH distribution
    p_val_emp <- 2 * min(mean(wqs_betas > 0), mean(wqs_betas < 0))

    cat(sprintf("Mean Beta across RH iterations :  %.4f\n", mean_beta))
    cat(sprintf("Empirical 95%% CI             : [%.4f, %.4f]\n", ci_lower, ci_upper))
    cat(sprintf("Empirical P-value            :  %.4f\n", p_val_emp))

    if (ci_lower > 0 || ci_upper < 0) {
      cat("\nConclusion: The overall joint mixture effect is significant (95% CI excludes 0).\n")
    } else {
      cat("\nConclusion: The overall joint mixture effect is NOT significant (95% CI includes 0).\n")
    }
  }

  # ==========================================================
  # Module 2: Calculate flexible joint exposure quantile contrast
  # ==========================================================
  cat("\n======================================================\n")
  cat(sprintf(" Joint Exposure Quantile Contrast: Target Q%d vs. Ref Q%d\n", q_target + 1, q_ref + 1))
  cat("======================================================\n")

  # Extract model components
  if (rh == 1) {
    shapes_vec <- model$shapes
    weights_vec <- model$final_weights
    beta_wqs <- coef(model)["wqs_score"]
  } else {
    shapes_vec <- model$mean_shapes
    weights_vec <- model$final_weights
    beta_wqs <- model$mean_coefs["wqs_score"]
  }

  # Dynamically determine the degrees of freedom for splines
  df_spline <- max(as.numeric(sub("^.+_B(\\d+)$", "\\1", names(shapes_vec))))

  # Build basis matrix for target and reference quantiles
  b_target <- splines::ns(c(q_target, q_ref), df = df_spline, intercept = FALSE)

  # Helper function to compute score difference
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

  # Compute absolute partial effect change (Delta Eta)
  if (rh > 1) {
    diff_list <- numeric(rh)
    for (i in 1:rh) {
      beta_i <- model$rh_coefs[i, "wqs_score"]
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

  # Interpret Odds Ratio (OR) if family is binomial
  if (is_binomial) {
    cat("\n----------------- Converted to Odds Ratio (OR) -----------------\n")
    cat(sprintf("Overall Joint OR :  %.4f\n", exp(delta_eta)))
    if (rh > 1) {
      cat(sprintf("95%% CI (OR)      : [%.4f, %.4f]\n", exp(ci_lower_eta), exp(ci_upper_eta)))
    }
    cat(sprintf(
      "\n[Epidemiological Interpretation]: \nWhen all mixture components are concurrently at quantile level %d (Q%d),\ncompared to all components at level %d (Q%d), the odds of the event is \nestimated to be %.4f times that of the reference group.\n",
      q_target, q_target + 1, q_ref, q_ref + 1, exp(delta_eta)
    ))
  } else {
    cat("\n[Epidemiological Interpretation]: \nThis value represents the absolute predicted increment on the continuous \nresponse scale when comparing the target joint exposure to the reference.\n")
  }

  invisible(list(delta_eta = delta_eta, lower = ci_lower_eta, upper = ci_upper_eta))
}

#' Plot Faceted Boxplot for NWQS Quantile Contrasts (Faceted by Component, with Scatter)
#' 绘制 NWQS 剂量反应轨迹：按组分分面的分组箱线图 (带散点叠加、高对比配色、独立 Y 轴)
#'
#' @param model nwqs 模型结果对象 (必须是 rh > 1 的结果)
#' @param exponentiate logical. 是否转换为 OR/RR 值 (默认为 NULL, 若为 binomial/poisson/quasipoisson 自动转换)
#' @param custom_colors 字符向量. 自定义配色板。默认提取了高对比度的绿、蓝、黄配色。
#' @param free_y logical. 各个分面的 Y 轴是否自由缩放。默认为 TRUE。
#' @export
#' @importFrom splines ns
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter geom_hline theme_bw labs scale_fill_manual element_text element_blank element_rect facet_wrap
plot_nwqs_contrast_box = function(model, exponentiate = NULL, 
                                  free_y = TRUE, 
                                  # 默认前三个颜色对标你的参考图: 绿, 蓝, 黄。后面保留安全色防报错
                                  custom_colors = c("#7DB97F", "#82B0D2", "#D92828", "#F2C05D", "#8B6FB8",
                                                    "#00B4D8", "#006B3C", "#F4B6B6", "#5BA3D0", "#E03030", 
                                                    "#7AD450", "#9B7FC0")) {
  if (model$rh < 2) {
    stop("This plot requires rh > 1 (Repeated Holdout iterations) to generate boxplots.")
  }

  is_exp_family = model$family %in% c("binomial", "poisson", "quasipoisson")
  if (is.null(exponentiate)) exponentiate = is_exp_family

  q_level = eval(model$call$q)
  if (is.null(q_level)) q_level = 4

  df_spline = max(as.numeric(sub("^.+_B(\\d+)$", "\\1", colnames(model$rh_shapes))))
  comps = colnames(model$rh_weights)
  full_basis = splines::ns(0:(q_level - 1), df = df_spline, intercept = FALSE)

  results_list = list()

  # 提取数据
  for (i in 1:model$rh) {
    beta_i = model$rh_coefs[i, "wqs_score"]

    for (q_tgt in 1:(q_level - 1)) {
      b_diff = full_basis[q_tgt + 1, ] - full_basis[1, ]
      
      current_comp_effects = numeric(length(comps))
      names(current_comp_effects) = comps

      for (comp in comps) {
        theta_comp = model$rh_shapes[i, paste0(comp, "_B", 1:df_spline)]
        w_comp = model$rh_weights[i, comp]

        comp_effect = beta_i * w_comp * sum(b_diff * theta_comp)
        current_comp_effects[comp] = comp_effect

        results_list[[length(results_list) + 1]] = data.frame(
          Iteration = i,
          Component = comp,
          Quantile = paste0("Q", q_tgt + 1), # ✨ 删掉了 "vs Q1"，直接叫 Q2, Q3, Q4
          Effect = comp_effect
        )
      }
      
      results_list[[length(results_list) + 1]] = data.frame(
        Iteration = i,
        Component = "Overall",
        Quantile = paste0("Q", q_tgt + 1),
        Effect = sum(current_comp_effects)
      )
    }
  }

  plot_df = do.call(rbind, results_list)
  plot_df$Component = factor(plot_df$Component, levels = c("Overall", comps))
  
  quantile_levels = paste0("Q", 2:q_level)
  plot_df$Quantile = factor(plot_df$Quantile, levels = quantile_levels)

  if (exponentiate) {
    plot_df$Effect = exp(plot_df$Effect)
    y_label = if (model$family == "binomial") "Odds Ratio (OR)" else "Rate Ratio (RR)"
    y_intercept = 1
  } else {
    y_label = "Absolute Effect Change (\u0394 Eta)"
    y_intercept = 0
  }

  n_facets = length(unique(plot_df$Component)) 
  dynamic_nrow = ceiling(n_facets / 7)
  dynamic_ncol = ceiling(n_facets / dynamic_nrow)
  
  n_contrasts = q_level - 1
  plot_colors = rep(custom_colors, length.out = n_contrasts)
  facet_scales = ifelse(free_y, "free_y", "fixed")


  p = ggplot2::ggplot(plot_df, ggplot2::aes(x = Quantile, y = Effect, fill = Quantile)) +
    ggplot2::geom_boxplot(
      outlier.shape = NA, 
      alpha = 0.6,          
      color = "gray20",
      linewidth = 0.5,
      width = 0.5
    ) +
    
    ggplot2::geom_jitter(
      shape = 21,          
      color = "gray30",    
      alpha = 0.7,          
      width = 0.2,        
      size = 1.2           
    ) +
    
    ggplot2::geom_hline(yintercept = y_intercept, linetype = "dashed", color = "#2C3E50", linewidth = 0.8) +
    ggplot2::facet_wrap(~ Component, scales = facet_scales, nrow = dynamic_nrow, ncol = dynamic_ncol) +
    ggplot2::scale_fill_manual(values = plot_colors) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 0, face = "bold", size = 11), # 标签变短了，无需倾斜，水平居中即可
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(), 
      strip.background = ggplot2::element_rect(fill = "#ECF0F1"),
      strip.text = ggplot2::element_text(face = "bold", size = 12),
      plot.caption = ggplot2::element_text(face = "italic", color = "gray40", size = 10, hjust = 0) # 左下角的脚注格式
    ) +
    ggplot2::labs(
      title = "Component-Specific Dose-Response Trajectories",
      subtitle = paste("Based on", model$rh, "Repeated Holdout Iterations"),
      caption = "* Note: All estimated effects (e.g., Q2, Q3) are compared to the first quartile (Q1) baseline.", # ✨ 明确的比较说明文字
      x = "Exposure Quantiles",
      y = y_label
    )

  return(p)
}

#' Calculate Weight Allocation Error Metrics
#' 计算权重分配的 SAE 和 RMSE
#'
#' @param w_est numeric vector. 模型估计出的权重 (Estimated weights).
#' @param w_true numeric vector. 真实的设定权重 (True weights).
#' @return A list containing SAE and RMSE.
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
    SAE = sae,
    MAE = mae,
    Pearson = cor_pearson,
    Spearman = cor_spearman,
    CosSim = cos_sim
  ))
}


#' 评估 100 次模拟的宏观统计学指标 (Macro-level Metrics)
#'
#' @param sim_df data.frame. 包含了 100 次模拟结果的汇总数据框
#' @param true_weights numeric vector. 真实的设定权重
#' @param true_effect numeric. 设定的真实总体效应 (例如 Q4 vs Q1 的 Delta)
evaluate_sim_performance <- function(sim_weight_df, sim_effect_df, true_w, true_eff_mat, w_threshold = 0.05) {
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

  weight_res <- data.frame(Mean_SAE = mean(sim_weight_df$NWQS_SAE, na.rm = TRUE), Sensitivity = sens, Specificity = spec)

  true_eff_long <- as.data.frame(as.table(true_eff_mat))
  colnames(true_eff_long) <- c("Term", "Target", "True_Value")
  true_eff_long$Target <- as.character(true_eff_long$Target)
  true_eff_long$Term <- as.character(true_eff_long$Term)

  sim_effect_df$CI_Lower <- sim_effect_df$Estimate - 1.96 * sim_effect_df$SE
  sim_effect_df$CI_Upper <- sim_effect_df$Estimate + 1.96 * sim_effect_df$SE

  eval_df <- merge(sim_effect_df, true_eff_long, by = c("Target", "Term"))
  eval_df$Covered <- (eval_df$True_Value >= eval_df$CI_Lower) & (eval_df$True_Value <= eval_df$CI_Upper)
  eval_df$Abs_Bias <- eval_df$Estimate - eval_df$True_Value

  effect_res <- eval_df %>%
    group_by(Target, Term) %>%
    summarise(
      True_Value = mean(True_Value),
      Mean_Est = mean(Estimate),
      Mean_Bias = mean(Abs_Bias),
      RB_pct = ifelse(abs(mean(True_Value)) > 1e-5, mean(Abs_Bias) / abs(mean(True_Value)) * 100, NA),
      RMSE = sqrt(mean(Abs_Bias^2)),
      Coverage_Prob = mean(Covered),
      Reject_H0 = mean(CI_Lower > 0 | CI_Upper < 0),
      .groups = "drop"
    ) %>%
    mutate(
      Metric_Type = ifelse(abs(True_Value) > 1e-5, "Power", "Type I Error")
    ) %>%
    arrange(Target, factor(Term, levels = c("Overall", names(true_w))))

  return(list(Weight_Metrics = weight_res, Effect_Metrics = effect_res))
}



#' Extract Effects and Empirical Bootstrap CI from NWQS Object (真·全局经验分位数版)
#' @export
extract_nwqs_effects <- function(model_res, return_raw = FALSE) {
  q_level <- ifelse(is.null(model_res$q), 4, eval(model_res$q))
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

    for (i in seq_len(rh)) {
      beta_i <- model_res$rh_coefs[i, "wqs_score"]

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




#' 检查单次模拟结果的置信区间是否覆盖真值
#'
#' @param est_df data.frame. extract_nwqs_effects() 输出的估计值数据框
#' @param true_mat matrix. 真实的效应矩阵 (true_effect_mat)
#' @return data.frame 包含估计值、真值、偏差以及两种 CI 的覆盖判断(TRUE/FALSE)
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

  library(dplyr)
  final_df <- merged_df %>%
    select(
      Target, Term, True_Value, Estimate, Bias,
      Wald_CI_Lower, Wald_CI_Upper, Covered_Wald,
      Empirical_CI_Lower, Empirical_CI_Upper, Covered_Empirical
    ) %>%
    arrange(Target, factor(Term, levels = c("Overall", setdiff(unique(Term), "Overall"))))

  return(final_df)
}


#' Evaluate Simulation Performance (Macro-level)
#' 计算 100 次模拟的宏观统计学指标 (CP, Bias, Sensitivity, Specificity)
#'
#' @param sim_weight_df data.frame. 100次模拟的权重估计集合
#' @param sim_effect_df data.frame. 100次模拟的效应量估计集合 (由 extract_nwqs_effects 组合而成)
#' @param true_w numeric vector. 真实的权重设定
#' @param true_eff_mat matrix. 真实的效应矩阵 (来自 gen_nonlinear_data 的 attribute)
#' @param w_threshold numeric. 判断权重是否被“识别”的阈值 (默认 0.01)
#' @return list 包含 Weight_Metrics 和 Effect_Metrics
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
    Mean_SAE = mean_sae,
    Sensitivity = sens,
    Specificity = spec
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
    group_by(Target, Term) %>%
    summarise(
      True_Value = mean(True_Value),
      Mean_Est = mean(Estimate),
      Mean_Bias = mean(Abs_Bias),
      RB_pct = ifelse(abs(mean(True_Value)) > 1e-5, mean(Abs_Bias) / abs(mean(True_Value)) * 100, NA),
      RMSE = sqrt(mean(Abs_Bias^2)),
      # 两种覆盖率
      CP_Wald = mean(Covered_Wald, na.rm = TRUE),
      CP_Empirical = mean(Covered_Empirical, na.rm = TRUE),
      # 统计效能 (用 Empirical CI)
      Reject_H0 = mean(Empirical_CI_Lower > 0 | Empirical_CI_Upper < 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Metric_Type = ifelse(abs(True_Value) > 1e-5, "Power", "Type I Error")
    ) %>%
    arrange(Target, factor(Term, levels = c("Overall", names(true_w))))

  return(list(Weight_Metrics = weight_res, Effect_Metrics = effect_res))
}

#' Plot Monte Carlo Benchmark Results (Deviance, SAE, and Weight Estimation)
#' @param dev_data 数据框。必须包含两列: `Model` (模型名称) 和 `Deviance` (残差拟合误差)。
#' @param sae_data 数据框。必须包含两列: `Model` (模型名称) 和 `SAE` (相对误差数值)。
#' @param weight_data 数据框。必须包含四列: `Model` (模型名称), `Component` (变量名), `Estimated_Weight` (估计权重), `True_Value` (真实权重)。
#' @param custom_palette 命名的颜色字符向量。如果为 NULL，将使用默认的高级学术色板。
#' @param save_path 字符串。保存图片的路径。如果为 NULL，则不自动保存。
#' @param base_size 基础字体大小，默认为 14。
#' @return 返回一个 patchwork 组合图形对象。
#' @export
plot_monte_carlo_benchmark <- function(dev_data, sae_data, weight_data, custom_palette = NULL, save_path = NULL, base_size = 14) {
  if (!requireNamespace("ggdist", quietly = TRUE)) stop("请安装 'ggdist' 包")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("请安装 'patchwork' 包")

  # 🔥 关键修复：这里的名字必须与引擎循环中完全一致
  if (is.null(custom_palette)) {
    custom_palette <- c(
      "NWQS"         = "#4A90C8", 
      "gWQS"         = "#D92828", 
      "QGcomp"      = "#6EC44A",
      "Ridge"        = "#8B6FB8", 
      "Lasso"        = "#00B4D8",
      "ElasticNet"   = "#006B3C", 
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
    ggplot2::geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5, color = "black", position = ggplot2::position_nudge(x = -0.1)) +
    ggplot2::geom_point(size = 1.3, alpha = 0.4, position = ggplot2::position_jitter(width = 0.05, height = 0)) +
    ggplot2::scale_fill_manual(values = custom_palette) +
    ggplot2::scale_color_manual(values = custom_palette) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, face = "bold"), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(title = "A. Model Fit Error (Deviance)", subtitle = "Lower deviance indicates better non-linear fit.", x = "", y = "Residual Deviance")

  p_sae <- ggplot2::ggplot(sae_df, ggplot2::aes(x = Model, y = SAE, fill = Model, color = Model)) +
    ggdist::stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, justification = -0.3, point_colour = NA, alpha = 0.7) +
    ggplot2::geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5, color = "black", position = ggplot2::position_nudge(x = -0.1)) +
    ggplot2::geom_point(size = 1.3, alpha = 0.4, position = ggplot2::position_jitter(width = 0.05, height = 0)) +
    ggplot2::scale_fill_manual(values = custom_palette) +
    ggplot2::scale_color_manual(values = custom_palette) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, face = "bold"), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(title = "B. Weight Extraction Error (SAE)", subtitle = "Lower SAE indicates higher accuracy.", x = "", y = "Sum of Absolute Errors (SAE)")

  p_facet <- ggplot2::ggplot(weight_df, ggplot2::aes(x = Model, y = Estimated_Weight, fill = Model)) +
    ggplot2::geom_boxplot(alpha = 0.8, outlier.size = 0.5, color = "black", width = 0.6) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = True_Value), linetype = "dashed", color = "black", linewidth = 1) +
    ggplot2::facet_wrap(~Component, scales = "free_y", nrow = dynamic_nrow, ncol = dynamic_ncol) +
    ggplot2::scale_fill_manual(values = custom_palette) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom", legend.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#ECF0F1"), strip.text = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(title = "C. Component-Specific Weight Recovery Accuracy", x = "", y = "Estimated Relative Weight")

  dynamic_height_ratio <- 0.8 * dynamic_nrow
  final_plot <- (p_dev | p_sae) / p_facet +
    patchwork::plot_layout(heights = c(1.2, dynamic_height_ratio)) +
    patchwork::plot_annotation(title = "Monte Carlo Simulation Benchmark", theme = ggplot2::theme(plot.title = ggplot2::element_text(size = base_size + 4, face = "bold", hjust = 0.5)))

  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, plot = final_plot, width = 16, height = 11 + (dynamic_nrow - 1) * 3.5, dpi = 500)
  }
  return(final_plot)
}