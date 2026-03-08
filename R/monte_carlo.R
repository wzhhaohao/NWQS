# ==============================================================================
# monte_carlo.R — Monte Carlo 模拟评估专用工具函数
# 从 utils.R 中分离，这些函数仅用于模拟研究，不属于 NWQS 包核心功能
# ==============================================================================


#' Calculate Weight Allocation Error Metrics
#'
#' @description
#' Computes multiple error and similarity metrics between estimated and true
#' mixture component weights, including Sum of Absolute Errors (SAE), Mean
#' Absolute Error (MAE), Pearson and Spearman correlations, and cosine similarity.
#'
#' @param w_est Numeric vector. Estimated weights from the model.
#' @param w_true Numeric vector. True (ground-truth) weights from data generation.
#'   Must be the same length as \code{w_est}. If both vectors are named, they are
#'   aligned by name.
#'
#' @return A named list with components:
#' \describe{
#'   \item{SAE}{Sum of Absolute Errors: \eqn{\sum |w_{est} - w_{true}|}.}
#'   \item{MAE}{Mean Absolute Error.}
#'   \item{Pearson}{Pearson correlation coefficient.}
#'   \item{Spearman}{Spearman rank correlation coefficient.}
#'   \item{CosSim}{Cosine similarity between the two weight vectors.}
#' }
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


#' Check Single-Simulation Coverage of Confidence Intervals
#'
#' @description
#' Compares estimated effects from \code{\link{extract_nwqs_effects}} against a
#' ground-truth effect matrix and determines whether Wald and empirical CIs cover
#' the true values.
#'
#' @param est_df A \code{data.frame} produced by \code{\link{extract_nwqs_effects}}.
#' @param true_mat A matrix of true effects (typically the \code{true_effect_mat}
#'   attribute from data generation functions like \code{\link{gen_nonlinear_data}}).
#'
#' @return A \code{data.frame} with columns for Target, Term, True_Value, Estimate,
#'   Bias, Wald and Empirical CI bounds, and logical coverage indicators.
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


#' Evaluate Monte Carlo Simulation Performance (Macro-level Metrics)
#'
#' @description
#' Aggregates results across multiple simulation replicates to compute macro-level
#' performance metrics: coverage probability, bias, RMSE, statistical power / Type I
#' error, and weight recovery accuracy (sensitivity/specificity).
#'
#' @param sim_weight_df A \code{data.frame} with one row per simulation, containing
#'   columns named \code{NWQS_{Component}} for estimated weights and \code{NWQS_SAE}.
#' @param sim_effect_df A \code{data.frame} with columns: Sim_ID, Target, Term,
#'   Estimate, Wald_CI_Lower, Wald_CI_Upper, Empirical_CI_Lower, Empirical_CI_Upper.
#' @param true_w Named numeric vector. True weight allocation.
#' @param true_eff_mat Matrix. True effect values from data generation.
#' @param w_threshold Numeric. Threshold for classifying a weight as "detected".
#'   Defaults to 0.01.
#'
#' @return A list with:
#' \describe{
#'   \item{Weight_Metrics}{A one-row \code{data.frame} with Mean_SAE, Sensitivity,
#'     and Specificity.}
#'   \item{Effect_Metrics}{A \code{data.frame} with per-Target/Term metrics including
#'     True_Value, Mean_Est, Mean_Bias, RB_pct, RMSE, CP_Wald, CP_Empirical,
#'     Reject_H0, and Metric_Type.}
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


#' Check Bootstrap CI Coverage for nwqs_boot Object
#'
#' @description
#' Evaluates whether bootstrap percentile confidence intervals from
#' \code{\link{nwqs_boot}} cover the supplied true values.
#'
#' @param boot_res An object of class \code{"nwqs_boot"}.
#' @param true_value Numeric. The true effect value(s) to check coverage against.
#'
#' @return A \code{data.frame} augmented with True_Value, Bias, and
#'   Covered_Bootstrap columns.
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


#' Derive True Importance Weights from Ground-Truth Effect Matrix
#'
#' @description
#' Converts the absolute component-specific effects (e.g., Q4 vs Q1 contrast)
#' from the true effect matrix into normalized importance weights, serving as
#' the gold standard for weight recovery benchmarking.
#'
#' @param true_effect_mat Matrix. The true effect matrix (from data generation attributes).
#' @param mix_name Character vector. Names of mixture components (must match rownames).
#' @param method Character. Method for deriving importance. Currently only
#'   \code{"q4q1_abs"} (absolute Q4 vs Q1 contrast) is supported.
#'
#' @return A named numeric vector of true importance weights summing to 1.
#'
#' @export
calc_true_importance <- function(true_effect_mat, mix_name, method = "q4q1_abs") {
  if (method == "q4q1_abs") {
    if (!"Q4_vs_Q1" %in% colnames(true_effect_mat)) {
      stop("true_effect_mat must contain column 'Q4_vs_Q1'.")
    }
    contrib <- abs(true_effect_mat[mix_name, "Q4_vs_Q1"])
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


#' Plot Monte Carlo Benchmark Results
#'
#' @description
#' Produces a three-panel composite figure for Monte Carlo simulation benchmarking:
#' (A) Model Fit Error (Deviance) as raincloud plots, (B) Weight Extraction Error
#' (SAE) as raincloud plots, and (C) Component-Specific Weight Recovery Accuracy
#' as faceted boxplots with true-value reference lines.
#'
#' @param dev_data A \code{data.frame} with columns \code{Model} and \code{Deviance}.
#' @param sae_data A \code{data.frame} with columns \code{Model} and \code{SAE}.
#' @param weight_data A \code{data.frame} with columns \code{Model}, \code{Component},
#'   \code{Estimated_Weight}, and \code{True_Value}.
#' @param scen_name Character. Scenario label for the plot title.
#' @param family_name Character. GLM family label for the plot title.
#' @param custom_palette Named character vector of colors keyed by model name.
#'   If \code{NULL}, uses a default academic palette.
#' @param save_path Character or NULL. File path to save the plot. If \code{NULL},
#'   the plot is not saved automatically.
#' @param base_size Numeric. Base font size. Defaults to 14.
#'
#' @return A \code{patchwork} composite plot object (invisibly if saved).
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_point geom_hline facet_wrap
#'   scale_fill_manual scale_color_manual theme_bw labs theme element_text
#'   element_blank element_rect position_jitter position_nudge ggsave
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
