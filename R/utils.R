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
trans_quantile = function(data, method = c("quantile", "percentile"), q = 4) {
    data = as.data.frame(data)
    method = match.arg(method)
  
    if (method == "percentile") {
        transform_func = function(x) {
            rank(x) / (length(x) + 1)
        }
    
        # Apply to all columns
        res_list = lapply(data, transform_func)
    
    } else {
        if (!is.numeric(q) || q < 1) stop("'q' must be a positive number")

        transform_func = function(x) {
            breaks = unique(quantile(x, probs = seq(0, 1, by = 1/q), na.rm = TRUE))
      
            # Handle Boundaries (gWQS style: Force -Inf / Inf for robustness)
            if (length(breaks) == 1) {
                breaks = c(-Inf, breaks)
            } else {
                breaks[1] = -Inf
                breaks[length(breaks)] = Inf
            }
      
            # Cut (Binning): Returns integer values from 0 to q-1
            as.numeric(cut(x, breaks = breaks, labels = FALSE, include.lowest = TRUE)) - 1
        }
    
        # Apply to all columns
        res_list = lapply(data, transform_func)
    }
  
    # Convert list back to data.frame and preserve names
    res_df = as.data.frame(res_list)
    names(res_df) = names(data)
  
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
wqs_nonlinear_expand = function(data, mix_name, df_spline = 3, q = 4) {

    trans_data = data[, mix_name] 

    X = 0:(q - 1)
    temp_spline = splines::ns(X, df = df_spline)
    temp_knots = attr(temp_spline, "knots")
    temp_boundary = attr(temp_spline, "Boundary.knots")

    # 对每一列应用自然样条展开
    mat_spline_list = lapply(trans_data, function(x) splines::ns(x, knots = temp_knots, Boundary.knots = temp_boundary))
    mat_spline_full = do.call(cbind, mat_spline_list)

    # 生成列名
    total_cols = ncol(mat_spline_full)
    cols_per_mix = total_cols / length(mix_name)

    colnames(mat_spline_full) = paste0(
        rep(mix_name, each = cols_per_mix), 
        "_B", 
        rep(1:cols_per_mix, times = length(mix_name)))
    
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
    wqs_sum <- summary(model)$coefficients["wqs_score", , drop=FALSE]
    print(round(wqs_sum, 4))
    if (wqs_sum[1, 4] < 0.05) {
      cat("\nConclusion: The overall NWQS latent risk score has a significant joint effect on the outcome (P < 0.05).\n")
    } else {
      cat("\nConclusion: The overall NWQS latent risk score does not have a significant joint effect on the outcome (P >= 0.05).\n")
    }
  } else {
    wqs_betas <- model$rh_coefs[, "wqs_score"]
    mean_beta <- mean(wqs_betas, na.rm=TRUE)
    ci_lower <- quantile(wqs_betas, 0.025, na.rm=TRUE)
    ci_upper <- quantile(wqs_betas, 0.975, na.rm=TRUE)
    
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
  cat(sprintf(" Joint Exposure Quantile Contrast: Target Q%d vs. Ref Q%d\n", q_target+1, q_ref+1))
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
      s_tgt <- b_mat[1, , drop=FALSE] %*% theta
      s_ref <- b_mat[2, , drop=FALSE] %*% theta
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
    cat(sprintf("\n[Epidemiological Interpretation]: \nWhen all mixture components are concurrently at quantile level %d (Q%d),\ncompared to all components at level %d (Q%d), the odds of the event is \nestimated to be %.4f times that of the reference group.\n", 
                q_target, q_target+1, q_ref, q_ref+1, exp(delta_eta)))
  } else {
    cat("\n[Epidemiological Interpretation]: \nThis value represents the absolute predicted increment on the continuous \nresponse scale when comparing the target joint exposure to the reference.\n")
  }
  
  invisible(list(delta_eta = delta_eta, lower = ci_lower_eta, upper = ci_upper_eta))
}



# 临时放一下，一个作图的
#' Plot Grouped Boxplot for NWQS Quantile Contrasts
#' 绘制 NWQS 各个成分暴露分位数对比的分组箱线图 (常规原始尺度)
#'
#' @param model nwqs 模型结果对象 (必须是 rh > 1 的结果)
#' @param exponentiate logical. 是否转换为 OR/RR 值 (默认为 NULL, 若为 binomial/poisson/quasipoisson 自动转换)
#' @export
#' @importFrom splines ns
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_hline theme_bw labs position_dodge scale_fill_brewer
plot_nwqs_contrast_box <- function(model, exponentiate = NULL) {
  
  if (model$rh < 2) {
    stop("This plot requires rh > 1 (Repeated Holdout iterations) to generate boxplots.")
  }
  
  # 自动识别类型和参数 (支持二分类和计数型自动指数化)
  is_exp_family <- model$family %in% c("binomial", "poisson", "quasipoisson")
  if (is.null(exponentiate)) exponentiate <- is_exp_family
  
  q_level <- eval(model$call$q)
  if (is.null(q_level)) q_level <- 4
  
  df_spline <- max(as.numeric(sub("^.+_B(\\d+)$", "\\1", colnames(model$rh_shapes))))
  comps <- colnames(model$rh_weights)
  
  # 生成完整的基函数矩阵，保证节点 (knots) 划分准确
  full_basis <- splines::ns(0:(q_level-1), df = df_spline, intercept = FALSE)
  
  results_list <- list()
  
  # 遍历每一次 RH 迭代提取数据
  for (i in 1:model$rh) {
    beta_i <- model$rh_coefs[i, "wqs_score"]
    
    # 遍历不同的目标分位数 (Q2, Q3, Q4) 与 Q1 (index 0) 对比
    for (q_tgt in 1:(q_level-1)) {
      b_diff <- full_basis[q_tgt + 1, ] - full_basis[1, ]
      
      # 仅计算各个成分的独立贡献量
      for (comp in comps) {
        theta_comp <- model$rh_shapes[i, paste0(comp, "_B", 1:df_spline)]
        w_comp <- model$rh_weights[i, comp]
        
        comp_effect <- beta_i * w_comp * sum(b_diff * theta_comp)
        
        results_list[[length(results_list) + 1]] <- data.frame(
          Iteration = i,
          Component = comp,
          Quantile = paste0("Q", q_tgt + 1, " vs Q1"),
          Effect = comp_effect
        )
      }
    }
  }
  
  # 合并数据并处理因子顺序
  plot_df <- do.call(rbind, results_list)
  # 仅保留各个 Component (去掉了 Overall)
  plot_df$Component <- factor(plot_df$Component, levels = comps)
  
  # 是否转换为 OR/RR 值，并智能匹配 Y 轴标签
  if (exponentiate) {
    plot_df$Effect <- exp(plot_df$Effect)
    y_label <- if (model$family == "binomial") "Odds Ratio (OR)" else "Rate Ratio (RR)"
    y_intercept <- 1
  } else {
    y_label <- "Absolute Effect Change (\u0394 Eta)"
    y_intercept <- 0
  }
  
  # 绘图部分
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Component, y = Effect, fill = Quantile)) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(0.8), 
                          outlier.size = 0.5, 
                          alpha = 0.85, 
                          color = "gray20") +
    ggplot2::geom_hline(yintercept = y_intercept, linetype = "dashed", color = "red", linewidth = 0.8) +
    ggplot2::scale_fill_brewer(palette = "YlGnBu") +  
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold"),
      panel.grid.major.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = "Distribution of Component-specific Effects",
      subtitle = paste("Based on", model$rh, "Repeated Holdout Iterations"),
      x = "Exposure Components",
      y = y_label,
      fill = "Contrast Level:"
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
calc_weight_error = function(w_est, w_true) {
  if(length(w_est) != length(w_true)) stop("Lengths of estimated and true weights must match.")
  
  # 确保顺序对齐
  if(!is.null(names(w_est)) && !is.null(names(w_true))) {
    w_true = w_true[names(w_est)]
  }
  
  error_diff = w_est - w_true
  
  sae = sum(abs(error_diff))
  rmse = sqrt(mean(error_diff^2))
  
  return(list(SAE = sae, RMSE = rmse))
}

