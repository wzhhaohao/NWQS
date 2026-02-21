# #' Plot Diagnostics for NWQS
# #' @param x An object of class "nwqs".
# #' @param type "curves" (剂量反应曲线) 或 "weights" (权重条形图)。
# #' @param components 字符向量。指定要画的污染物名称。默认为 NULL (画全部)。
# #' @param overlay 逻辑值。TRUE表示将多条曲线画在同一张图中；FALSE则分面显示。默认为 TRUE。
# #' @param plot_ci 逻辑值。是否绘制 95% 置信区间。默认为 TRUE。
# #' @param ... 额外参数。
# #' @export
# #' @method plot nwqs
# #' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_col geom_errorbar geom_hline facet_wrap theme_minimal labs theme element_text coord_flip scale_color_viridis_d scale_fill_viridis_d
# #' @importFrom splines ns
# #' @importFrom stats quantile
# plot.nwqs <- function(x, type = c("curves", "weights"), components = NULL, 
#                       overlay = TRUE, plot_ci = TRUE, base_size = 14, ...) {
  
#   type <- match.arg(type)
#   if (!inherits(x, "nwqs")) stop("Object must be of class 'nwqs'")
  
#   # ==========================================
#   # 图形 1: 权重直条图 (Weights Barplot)
#   # ==========================================
#   if (type == "weights") {
#     w_names <- names(x$final_weights)
#     if (!is.null(components)) w_names <- intersect(w_names, components)
    
#     w_df <- data.frame(Component = w_names, Weight = x$final_weights[w_names])
    
#     if (!is.null(x$rh_weights) && x$rh > 1 && plot_ci) {
#       rh_w <- x$rh_weights[, w_names, drop = FALSE]
#       w_df$Lower <- apply(rh_w, 2, function(v) quantile(v, 0.025, names = FALSE, na.rm=TRUE))
#       w_df$Upper <- apply(rh_w, 2, function(v) quantile(v, 0.975, names = FALSE, na.rm=TRUE))
#     } else {
#       w_df$Lower <- w_df$Weight
#       w_df$Upper <- w_df$Weight
#     }
    
#     # 排序使条形图美观
#     w_df$Component <- factor(w_df$Component, levels = w_df$Component[order(w_df$Weight)])
    
#     p <- ggplot2::ggplot(w_df, ggplot2::aes(x = Component, y = Weight)) +
#       ggplot2::geom_col(fill = "#3498DB", alpha = 0.8, width = 0.6) +
#       ggplot2::coord_flip() +
#       ggplot2::theme_minimal(base_size = base_size) +
#       ggplot2::labs(title = "Component Weights (NWQS)", x = "Pollutant", y = "Estimated Weight")
    
#     if (plot_ci && x$rh > 1) {
#       p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper), width = 0.2, color = "#2C3E50")
#     }
#     return(p)
#   }
  
#   # ==========================================
#   # 图形 2: 剂量反应曲线 (Dose-Response Curves)
#   # ==========================================
#   q_level <- if (!is.null(x$call$q)) x$call$q else 4 
#   x_seq <- seq(0, q_level - 1, length.out = 100)
#   shape_names <- names(x$mean_shapes)
  
#   pattern <- "^(.+)_B(\\d+)$"
#   parsed_names <- data.frame(
#     full_name = shape_names,
#     component = sub(pattern, "\\1", shape_names),
#     basis_idx = as.numeric(sub(pattern, "\\2", shape_names)),
#     stringsAsFactors = FALSE
#   )
  
#   df_spline <- max(parsed_names$basis_idx)
#   basis_mat <- splines::ns(x_seq, df = df_spline, intercept = FALSE)
  
#   unique_comps <- unique(parsed_names$component)
#   if (!is.null(components)) unique_comps <- intersect(unique_comps, components)
  
#   plot_data_list <- list()
#   for (comp in unique_comps) {
#     comp_cols <- parsed_names$full_name[parsed_names$component == comp]
#     if (!is.null(x$rh_shapes) && plot_ci && x$rh > 1) {
#       beta_mat <- x$rh_shapes[, comp_cols, drop = FALSE]
#     } else {
#       beta_mat <- matrix(x$mean_shapes[comp_cols], nrow = 1)
#     }
#     y_pred_mat <- as.matrix(basis_mat) %*% t(beta_mat)
    
#     y_stats <- t(apply(y_pred_mat, 1, function(v) {
#       c(mean = mean(v, na.rm=TRUE), 
#         lower = quantile(v, 0.025, names = FALSE, na.rm=TRUE), 
#         upper = quantile(v, 0.975, names = FALSE, na.rm=TRUE))
#     }))
    
#     plot_data_list[[comp]] <- data.frame(
#       x = x_seq, y = y_stats[, "mean"],
#       ymin = y_stats[, "lower"], ymax = y_stats[, "upper"],
#       Component = comp
#     )
#   }
#   final_df <- do.call(rbind, plot_data_list)
  
#   # 构建基础 ggplot
#   p <- ggplot2::ggplot(final_df, ggplot2::aes(x = x, y = y, color = Component, fill = Component)) +
#     ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
#     ggplot2::theme_minimal(base_size = base_size) +
#     ggplot2::labs(title = "Non-linear Dose-Response Curves (NWQS)", 
#                   x = "Quantile Index (Exposure)", y = "Estimated Spline Effect")
  
#   # 根据是否 overlay (叠加) 来决定画法
#   if (overlay) {
#     # 叠加模式：画在同一个图里，使用科学色板
#     if (plot_ci && x$rh > 1) {
#       p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = ymin, ymax = ymax), alpha = 0.15, color = NA)
#     }
#     p <- p + ggplot2::geom_line(linewidth = 1.2) +
#              ggplot2::scale_color_viridis_d(option = "turbo") + 
#              ggplot2::scale_fill_viridis_d(option = "turbo")
#   } else {
#     # 分面模式：每个污染物一个子图
#     if (plot_ci && x$rh > 1) {
#       p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA)
#     }
#     p <- p + ggplot2::geom_line(linewidth = 1.2) +
#              ggplot2::facet_wrap(~ Component, scales = "free_y") +
#              ggplot2::theme(legend.position = "none") # 分面时自动隐藏图例
#   }
  
#   return(p)
# }
#' Plot Diagnostics for NWQS
#' @param x An object of class "nwqs".
#' @param type "both", "curves", 或 "weights".
#' @param y_scale "partial" (从0开始的偏效应) 或 "predicted" (加上截距的绝对预测值，让曲线悬浮).
#' @param components 字符向量。指定要画的污染物名称。
#' @param overlay 逻辑值。TRUE表示叠加，FALSE分面。
#' @param plot_ci 逻辑值。是否绘制置信区间。
#' @param base_size Integer. Base font size for ggplot2.
#' @param ... 额外参数。
#' @export
#' @method plot nwqs
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_col geom_errorbar geom_hline facet_wrap theme_minimal labs theme element_text coord_flip scale_color_viridis_d scale_fill_viridis_d
#' @importFrom splines ns
#' @importFrom stats quantile
plot.nwqs <- function(x, type = c("both", "curves", "weights"), 
                      y_scale = c("partial", "predicted"),
                      components = NULL, overlay = TRUE, plot_ci = TRUE, base_size = 14, ...) {
  
  type <- match.arg(type)
  y_scale <- match.arg(y_scale)
  if (!inherits(x, "nwqs")) stop("Object must be of class 'nwqs'")
  
  # ==========================================
  # 模块 1: 构建权重直条图 (p_w) (代码不变)
  # ==========================================
  if (type %in% c("weights", "both")) {
    w_names <- names(x$final_weights)
    if (!is.null(components)) w_names <- intersect(w_names, components)
    
    w_df <- data.frame(Component = w_names, Weight = x$final_weights[w_names])
    
    if (!is.null(x$rh_weights) && x$rh > 1 && plot_ci) {
      rh_w <- x$rh_weights[, w_names, drop = FALSE]
      w_df$Lower <- apply(rh_w, 2, function(v) quantile(v, 0.025, names = FALSE, na.rm=TRUE))
      w_df$Upper <- apply(rh_w, 2, function(v) quantile(v, 0.975, names = FALSE, na.rm=TRUE))
    } else {
      w_df$Lower <- w_df$Weight
      w_df$Upper <- w_df$Weight
    }
    
    w_df$Component <- factor(w_df$Component, levels = w_df$Component[order(w_df$Weight)])
    
    p_w <- ggplot2::ggplot(w_df, ggplot2::aes(x = Component, y = Weight)) +
      ggplot2::geom_col(fill = "#3498DB", alpha = 0.8, width = 0.6) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::labs(title = "Component Weights", x = "", y = "Estimated Weight")
    
    if (plot_ci && x$rh > 1) {
      p_w <- p_w + ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper), width = 0.2, color = "#2C3E50")
    }
  }
  
  # ==========================================
  # 模块 2: 构建剂量反应曲线 (p_c) - [截距悬浮修正]
  # ==========================================
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
    basis_mat <- splines::ns(x_seq, df = df_spline, intercept = FALSE)
    
    unique_comps <- unique(parsed_names$component)
    if (!is.null(components)) unique_comps <- intersect(unique_comps, components)
    
    plot_data_list <- list()
    for (comp in unique_comps) {
      comp_cols <- parsed_names$full_name[parsed_names$component == comp]
      
      # 1. 计算原始偏效应 (起点必然为0)
      if (!is.null(x$rh_shapes) && plot_ci && x$rh > 1) {
        beta_mat <- x$rh_shapes[, comp_cols, drop = FALSE]
        wqs_coefs <- x$rh_coefs[, "wqs_score"]
        comp_weights <- x$rh_weights[, comp]
        
        scaling_factor <- wqs_coefs * comp_weights
        beta_mat <- beta_mat * scaling_factor
        
        y_pred_mat <- as.matrix(basis_mat) %*% t(beta_mat)
        
        # [核心]: 如果需要悬浮，加上 RH 每次迭代各自的截距
        if (y_scale == "predicted") {
           intercepts <- x$rh_coefs[, "(Intercept)"]
           y_pred_mat <- sweep(y_pred_mat, 2, intercepts, "+")
        }
        
      } else {
        beta_mat <- matrix(x$mean_shapes[comp_cols], nrow = 1)
        scaling_factor <- x$mean_coefs["wqs_score"] * x$final_weights[comp]
        beta_mat <- beta_mat * scaling_factor
        
        y_pred_mat <- as.matrix(basis_mat) %*% t(beta_mat)
        
        # [核心]: 如果需要悬浮，加上全局平均截距
        if (y_scale == "predicted") {
           y_pred_mat <- y_pred_mat + x$mean_coefs["(Intercept)"]
        }
      }
      
      y_stats <- t(apply(y_pred_mat, 1, function(v) {
        c(mean = mean(v, na.rm=TRUE), 
          lower = quantile(v, 0.025, names = FALSE, na.rm=TRUE), 
          upper = quantile(v, 0.975, names = FALSE, na.rm=TRUE))
      }))
      
      plot_data_list[[comp]] <- data.frame(
        x = x_seq, y = y_stats[, "mean"],
        ymin = y_stats[, "lower"], ymax = y_stats[, "upper"],
        Component = comp
      )
    }
    final_df <- do.call(rbind, plot_data_list)
    
    # 根据 y_scale 动态修改 Y 轴标签
    y_label <- ifelse(y_scale == "predicted", 
                      "Predicted Effect (inc. Intercept)", 
                      "Absolute Partial Effect (from 0)")
    
    p_c <- ggplot2::ggplot(final_df, ggplot2::aes(x = x, y = y, color = Component, fill = Component)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::labs(title = "Non-linear Dose-Response", 
                    x = "Quantile Index (Exposure)", y = y_label)
    
    if (overlay) {
      if (plot_ci && x$rh > 1) {
        p_c <- p_c + ggplot2::geom_ribbon(ggplot2::aes(ymin = ymin, ymax = ymax), alpha = 0.15, color = NA)
      }
      p_c <- p_c + ggplot2::geom_line(linewidth = 1.2) +
               ggplot2::scale_color_viridis_d(option = "turbo") + 
               ggplot2::scale_fill_viridis_d(option = "turbo")
    } else {
      if (plot_ci && x$rh > 1) {
        p_c <- p_c + ggplot2::geom_ribbon(ggplot2::aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA)
      }
      p_c <- p_c + ggplot2::geom_line(linewidth = 1.2) +
               ggplot2::facet_wrap(~ Component, scales = "free_y") +
               ggplot2::theme(legend.position = "none")
    }
  }
  
  # ==========================================
  # 模块 3: 输出分发
  # ==========================================
  if (type == "weights") return(p_w)
  if (type == "curves") return(p_c)
  
  if (type == "both") {
    if (!requireNamespace("patchwork", quietly = TRUE)) stop("Install 'patchwork' package.")
    
    combined_plot <- p_w + p_c + patchwork::plot_layout(widths = c(1, 1.5)) +
      patchwork::plot_annotation(
        title = "NWQS Model Diagnostics",
        theme = ggplot2::theme(plot.title = ggplot2::element_text(size = base_size + 4, face = "bold", hjust = 0.5))
      )
    return(combined_plot)
  }
}



# 暂时也先不用，最后封包再弄
# #' Print method for NWQS objects
# #' @param x An object of class "nwqs".
# #' @param ... Further arguments.
# #' @export
# #' @method print nwqs
# print.nwqs <- function(x, ...) {
#   cat("\n=== Non-linear Weighted Quantile Sum (NWQS) Regression ===\n\n")
#   cat("Call:\n")
#   print(x$call)
  
#   cat(sprintf("\nDistribution Family: %s | Repeated Holdouts (RH): %d\n", x$family, x$rh))
  
#   cat("\n--- Pooled Coefficients ---\n")
#   print(round(x$mean_coefs, 4))
  
#   cat("\n--- Component Weights (Sorted) ---\n")
#   print(round(sort(x$final_weights, decreasing = TRUE), 4))
  
#   cat("\n")
#   invisible(x)
# }


# 这个暂时有问题，先不管
# #' Summary method for NWQS objects
# #' @param object An object of class "nwqs".
# #' @param ... Further arguments.
# #' @export
# #' @method summary nwqs
# #' @importFrom stats pnorm sd
# summary.nwqs <- function(object, ...) {
#   # 提取系数均值
#   est <- object$mean_coefs
  
#   # 利用 RH 迭代的方差计算经验标准误
#   if (object$rh > 1 && !is.null(object$rh_coefs)) {
#     se <- apply(object$rh_coefs, 2, sd, na.rm = TRUE)
#   } else {
#     warning("Cannot calculate empirical standard errors with rh = 1.")
#     se <- rep(NA, length(est))
#   }
  
#   # 计算 Z 统计量和 P 值
#   z_val <- est / se
#   p_val <- 2 * pnorm(abs(z_val), lower.tail = FALSE)
  
#   # 组合成标准回归结果表
#   coef_table <- data.frame(
#     Estimate = est,
#     `Std. Error` = se,
#     `z value` = z_val,
#     `Pr(>|z|)` = p_val,
#     check.names = FALSE
#   )
  
#   res <- list(
#     call = object$call,
#     family = object$family,
#     rh = object$rh,
#     coefficients = coef_table,
#     mean_aic = object$mean_aic,
#     df_null = object$df_null,
#     df_res = object$df_res
#   )
  
#   class(res) <- "summary.nwqs"
#   return(res)
# }

# #' Print method for summary.nwqs
# #' @export
# #' @method print summary.nwqs
# #' @importFrom stats printCoefmat
# print.summary.nwqs <- function(x, ...) {
#   cat("\nCall:\n")
#   print(x$call)
#   cat(sprintf("\nFamily: %s | RH: %d\n\n", x$family, x$rh))
  
#   cat("Coefficients:\n")
#   stats::printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE, signif.stars = TRUE)
  
#   cat(sprintf("\nMean AIC: %.2f\n", x$mean_aic))
#   invisible(x)
# }


# #' Extract Coefficients from NWQS objects
# #' @param object An object of class "nwqs".
# #' @param ... Further arguments.
# #' @export
# #' @method coef nwqs
# coef.nwqs <- function(object, ...) {
#   return(object$mean_coefs)
# }