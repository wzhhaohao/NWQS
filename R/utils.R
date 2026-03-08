# ==============================================================================
# utils.R — NWQS 模型核心工具函数（清理版）
# 移除了所有 Monte Carlo 专用函数，这些已迁移至 monte_carlo.R
# 移除了重复定义的 evaluate_sim_performance (旧版多模型版本)
# ==============================================================================


#' Quantile or Percentile Transformation
#'
#' @description
#' Transforms continuous mixture variables into discrete quantile bins or
#' continuous percentile ranks.
#'
#' @details
#' Two methods are available:
#' \describe{
#'   \item{\code{"quantile"}}{gWQS-style integer binning into \code{q} groups
#'     (0 to \code{q-1}). Handles ties and boundary values robustly using
#'     \code{-Inf}/\code{Inf} endpoints.}
#'   \item{\code{"percentile"}}{Continuous percentile ranking scaled to (0, 1)
#'     using \code{rank(x) / (n + 1)}.}
#' }
#'
#' @param data A \code{data.frame} of mixture variables to transform.
#' @param method Character. Transformation method: \code{"quantile"} (default)
#'   or \code{"percentile"}.
#' @param q Integer. Number of quantile bins (only used when
#'   \code{method = "quantile"}). Defaults to 4 (quartiles).
#'
#' @return A \code{data.frame} with the same column names, containing
#'   transformed values.
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

#' Nonlinear Spline Expansion for WQS Mixture Components
#'
#' @description
#' Transforms (already quantile-transformed) mixture variables into natural cubic
#' spline basis matrices using globally fixed knots and boundary knots. This ensures
#' consistent basis alignment across training, validation, and bootstrap datasets.
#'
#' @param data A \code{data.frame} containing the mixture variables (already transformed).
#' @param mix_name Character vector. Column names of the mixture components to expand.
#' @param df_spline Integer. Degrees of freedom for the natural spline. Defaults to 3.
#' @param knots Numeric vector. Internal knot positions. Must be provided to ensure
#'   global scale alignment.
#' @param boundary Numeric vector of length 2. Boundary knot positions. Must be provided.
#'
#' @return A numeric matrix with columns named \code{{Component}_B{BasisIndex}}
#'   (e.g., \code{Component1_B1}, \code{Component1_B2}, ...).
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

#' Compute NWQS Joint Exposure Quantile Contrast
#'
#' @description
#' Computes the overall significance of the non-linear mixture effect and evaluates
#' joint exposure quantile contrasts (e.g., all components at Q4 vs all at Q1).
#' Automatically converts to Odds Ratios (OR) for binomial models and Rate Ratios
#' (RR) for Poisson models.
#'
#' @param model An object of class \code{"nwqs"}.
#' @param q_target Integer. Target quantile index (0-based). E.g., 3 represents Q4.
#'   Defaults to 3.
#' @param q_ref Integer. Reference quantile index (0-based). Defaults to 0 (Q1).
#'
#' @return Invisibly returns a list with \code{delta_eta}, \code{lower}, and
#'   \code{upper} (the latter two are \code{NA} when \code{rh = 1}).
#'
#' @importFrom splines ns
#' @importFrom stats quantile coef
#' @export
nwqs_contrast <- function(model, q_target = 3, q_ref = 0) {
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

  cat("\n======================================================\n")
  cat(sprintf(" Joint Exposure Quantile Contrast: Target Q%d vs. Ref Q%d\n", q_target + 1, q_ref + 1))
  cat("======================================================\n")

  if (rh == 1) {
    shapes_vec <- model$shapes
    weights_vec <- model$final_weights
    beta_wqs <- coef(model)["wqs_score"]
  } else {
    shapes_vec <- model$mean_shapes
    weights_vec <- model$final_weights
    beta_wqs <- model$mean_coefs["wqs_score"]
  }

  df_spline <- max(as.numeric(sub("^.+_B(\\d+)$", "\\1", names(shapes_vec))))
  b_target <- splines::ns(c(q_target, q_ref), df = df_spline, intercept = FALSE)

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
  } else {
    cat("\n[Interpretation]: Absolute predicted increment on the response scale.\n")
  }

  invisible(list(delta_eta = delta_eta, lower = ci_lower_eta, upper = ci_upper_eta))
}


# -------------------------------------------------------------------------

#' Plot Faceted Boxplot for NWQS Quantile Contrasts
#'
#' @description
#' Generates a faceted boxplot showing the distribution of component-specific
#' dose-response effects across Repeated Holdout iterations. Each facet represents
#' one mixture component (plus an "Overall" panel), with separate boxplots for
#' each quantile contrast (e.g., Q2, Q3, Q4 vs Q1 baseline).
#'
#' @param model An object of class \code{"nwqs"} with \code{rh > 1}.
#' @param exponentiate Logical or NULL. If \code{TRUE}, effects are exponentiated
#'   to OR/RR scale. If \code{NULL} (default), automatically set to \code{TRUE}
#'   for binomial/Poisson/quasi-Poisson families.
#' @param free_y Logical. Whether each facet has independent y-axis scaling.
#'   Defaults to \code{TRUE}.
#' @param custom_colors Character vector. Color palette for quantile groups.
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom splines ns
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter geom_hline facet_wrap
#'   scale_fill_manual theme_bw labs theme element_text element_blank element_rect
#' @export
plot_nwqs_contrast_box <- function(model, exponentiate = NULL,
                                   free_y = TRUE,
                                   custom_colors = c(
                                     "#7DB97F", "#82B0D2", "#D92828", "#F2C05D", "#8B6FB8",
                                     "#00B4D8", "#006B3C", "#F4B6B6", "#5BA3D0", "#E03030",
                                     "#7AD450", "#9B7FC0"
                                   )) {
  if (model$rh < 2) {
    stop("This plot requires rh > 1 (Repeated Holdout iterations) to generate boxplots.")
  }

  is_exp_family <- model$family %in% c("binomial", "poisson", "quasipoisson")
  if (is.null(exponentiate)) exponentiate <- is_exp_family

  q_level <- eval(model$call$q)
  if (is.null(q_level)) q_level <- 4

  df_spline <- max(as.numeric(sub("^.+_B(\\d+)$", "\\1", colnames(model$rh_shapes))))
  comps <- colnames(model$rh_weights)
  full_basis <- splines::ns(0:(q_level - 1), df = df_spline, intercept = FALSE)

  results_list <- list()

  for (i in 1:model$rh) {
    beta_i <- model$rh_coefs[i, "wqs_score"]

    for (q_tgt in 1:(q_level - 1)) {
      b_diff <- full_basis[q_tgt + 1, ] - full_basis[1, ]

      current_comp_effects <- numeric(length(comps))
      names(current_comp_effects) <- comps

      for (comp in comps) {
        theta_comp <- model$rh_shapes[i, paste0(comp, "_B", 1:df_spline)]
        w_comp <- model$rh_weights[i, comp]

        comp_effect <- beta_i * w_comp * sum(b_diff * theta_comp)
        current_comp_effects[comp] <- comp_effect

        results_list[[length(results_list) + 1]] <- data.frame(
          Iteration = i, Component = comp,
          Quantile = paste0("Q", q_tgt + 1), Effect = comp_effect
        )
      }

      results_list[[length(results_list) + 1]] <- data.frame(
        Iteration = i, Component = "Overall",
        Quantile = paste0("Q", q_tgt + 1), Effect = sum(current_comp_effects)
      )
    }
  }

  plot_df <- do.call(rbind, results_list)
  plot_df$Component <- factor(plot_df$Component, levels = c("Overall", comps))
  plot_df$Quantile <- factor(plot_df$Quantile, levels = paste0("Q", 2:q_level))

  if (exponentiate) {
    plot_df$Effect <- exp(plot_df$Effect)
    y_label <- if (model$family == "binomial") "Odds Ratio (OR)" else "Rate Ratio (RR)"
    y_intercept <- 1
  } else {
    y_label <- "Absolute Effect Change (\u0394 Eta)"
    y_intercept <- 0
  }

  n_facets <- length(unique(plot_df$Component))
  dynamic_nrow <- ceiling(n_facets / 7)
  dynamic_ncol <- ceiling(n_facets / dynamic_nrow)
  plot_colors <- rep(custom_colors, length.out = q_level - 1)
  facet_scales <- ifelse(free_y, "free_y", "fixed")

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Quantile, y = Effect, fill = Quantile)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6, color = "gray20",
                          linewidth = 0.5, width = 0.5) +
    ggplot2::geom_jitter(shape = 21, color = "gray30", alpha = 0.7, width = 0.2, size = 1.2) +
    ggplot2::geom_hline(yintercept = y_intercept, linetype = "dashed", color = "#2C3E50", linewidth = 0.8) +
    ggplot2::facet_wrap(~Component, scales = facet_scales, nrow = dynamic_nrow, ncol = dynamic_ncol) +
    ggplot2::scale_fill_manual(values = plot_colors) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 0, face = "bold", size = 11),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#ECF0F1"),
      strip.text = ggplot2::element_text(face = "bold", size = 12),
      plot.caption = ggplot2::element_text(face = "italic", color = "gray40", size = 10, hjust = 0)
    ) +
    ggplot2::labs(
      title = "Component-Specific Dose-Response Trajectories",
      subtitle = paste("Based on", model$rh, "Repeated Holdout Iterations"),
      caption = "* Note: All effects are relative to the Q1 baseline.",
      x = "Exposure Quantiles", y = y_label
    )

  return(p)
}


# -------------------------------------------------------------------------

#' Extract Effects from NWQS Object
#'
#' @description
#' Computes all quantile contrast effects (e.g., Q2 vs Q1, Q3 vs Q1, Q4 vs Q1)
#' for each mixture component and the overall mixture, along with standard errors
#' and empirical confidence intervals derived from Repeated Holdout iterations.
#'
#' @param model_res An object of class \code{"nwqs"}.
#' @param return_raw Logical. Currently unused. Reserved for future use.
#'
#' @return A \code{data.frame} with columns: Target, Term, Estimate, SE,
#'   Wald_CI_Lower, Wald_CI_Upper, Empirical_CI_Lower, Empirical_CI_Upper.
#'
#' @importFrom splines ns
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



#' Plot Component Weight Distributions for NWQS Models
#'
#' @description
#' Generates a vertical boxplot displaying the distribution of estimated mixture
#' component weights across bootstrap replicates (for \code{nwqs_boot} objects) or
#' Repeated Holdout iterations (for \code{nwqs} objects with \code{rh > 1}). Components
#' are ordered from highest to lowest by their point-estimate weight.
#'
#' @details
#' The function automatically detects the input model type and extracts weight
#' distributions from the appropriate source:
#' \itemize{
#'   \item \strong{nwqs_boot with keep_fits=TRUE}: Extracts weights from each bootstrap
#'     replicate's fitted model. Subtitle indicates "Bootstrap weight distribution".
#'   \item \strong{nwqs_boot with keep_fits=FALSE and rh_inner>1}: Falls back to the
#'     inner point fit's RH weight matrix.
#'   \item \strong{nwqs with rh>1}: Uses the RH weight matrix directly. A red caption
#'     warns that distributions reflect data-splitting variance only.
#' }
#'
#' @param model An object of class \code{"nwqs_boot"} or \code{"nwqs"}.
#'   For \code{"nwqs_boot"}, set \code{keep_fits = TRUE} during fitting to enable
#'   bootstrap-based weight distributions. For \code{"nwqs"}, \code{rh > 1} is required.
#' @param base_size Integer. Base font size for the ggplot2 theme. Default is 12.
#' @param palette Character. Color palette name. Options include \code{"default"}
#'   and \code{"palette2"}. Default is \code{"default"}.
#' @param ... Additional arguments (currently unused, reserved for future extensions).
#'
#' @return A \code{ggplot} object showing vertical boxplots of weight distributions
#'   for each mixture component, with jittered raw points and diamond-shaped mean markers.
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter stat_summary
#'   scale_fill_manual theme_bw labs theme element_blank element_line element_text
#' @export
plot_nwqs_weight_box <- function(model, base_size = 12,
                                 palette = "default", ...) {

  .get_palette <- function(n, pal="default") {
    cols <- list(
      default  = c("#4A90C8","#D92828","#6EC44A","#8B6FB8","#00B4D8",
                   "#006B3C","#A8D8EA","#F4B6B6","#5BA3D0","#E03030","#7AD450","#9B7FC0"),
      palette2 = c("#9bbf8a","#82afda","#f79059","#e7dbd3","#c2bdde",
                   "#8dcec8","#add3e2","#3480b8","#ffbe7a","#fa8878","#c82423","#6b5b95")
    )
    p <- if (pal %in% names(cols)) cols[[pal]] else cols[["default"]]
    rep(p, ceiling(n/length(p)))[seq_len(n)]
  }

  # ── 判断输入类型并提取权重矩阵 ──
  if (inherits(model, "nwqs_boot")) {
    boot_fits <- model$boot_fits
    valid_fits <- if (!is.null(boot_fits)) Filter(Negate(is.null), boot_fits) else NULL

    if (!is.null(valid_fits) && length(valid_fits) > 0) {
      w_mat <- do.call(rbind, lapply(valid_fits, function(bf) bf$final_weights))
      ci_source <- "bootstrap"
      n_iter <- length(valid_fits)
      inner_rh <- model$point_fit$rh
      if (inner_rh > 1) ci_source <- "bootstrap_with_inner_rh"
    } else {
      pf <- model$point_fit
      if (is.null(pf$rh_weights) || pf$rh <= 1)
        stop("Weight boxplot requires either keep_fits=TRUE in nwqs_boot() or rh > 1 in the inner model.")
      w_mat <- pf$rh_weights
      ci_source <- "rh_splitting"
      n_iter <- pf$rh
    }
  } else if (inherits(model, "nwqs")) {
    if (is.null(model$rh_weights) || model$rh <= 1)
      stop("Weight boxplot requires rh > 1.")
    w_mat <- model$rh_weights
    ci_source <- "rh_splitting"
    n_iter <- model$rh
  } else {
    stop("model must be of class 'nwqs' or 'nwqs_boot'.")
  }

  # ── 提取点估计权重用于排序 ──
  mix_names <- colnames(w_mat)
  n_comps   <- length(mix_names)

  if (inherits(model, "nwqs_boot")) point_w <- model$final_weights
  else point_w <- model$final_weights

  # 按点估计权重从高到低排序
  ordered_names <- names(sort(point_w, decreasing = TRUE))

  # ── 转为长格式 ──
  w_long <- data.frame(
    Component = rep(mix_names, each = nrow(w_mat)),
    Weight    = as.vector(w_mat)
  )
  w_long$Component <- factor(w_long$Component, levels = ordered_names)

  # ── caption / subtitle 逻辑 ──
  if (ci_source == "bootstrap") {
    subtitle_text <- sprintf("Bootstrap weight distribution (%d replicates)", n_iter)
    caption_text <- NULL
  } else if (ci_source == "bootstrap_with_inner_rh") {
    subtitle_text <- sprintf("Bootstrap weight distribution (%d replicates, inner rh=%d)", n_iter, model$point_fit$rh)
    caption_text <- "Note: Inner rh > 1 introduces additional splitting variance; distributions may be slightly wider than true sampling variance."
  } else {
    subtitle_text <- sprintf("Weight distribution across %d Repeated Holdout iterations", n_iter)
    caption_text <- "Note: Distributions reflect data-splitting variance only (NOT valid for inference). Use nwqs_boot(keep_fits=TRUE) for valid bootstrap distributions."
  }

  colors <- .get_palette(n_comps, palette)

  # ── 竖直箱线图 (组分在 X 轴, 权重在 Y 轴) ──
  p <- ggplot2::ggplot(w_long, ggplot2::aes(x = Component, y = Weight, fill = Component)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, color = "gray30") +
    ggplot2::geom_jitter(alpha = 0.3, width = 0.15, size = 0.8, color = "gray40") +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "#2C3E50") +
    ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::labs(
      title    = "Component Weight Distribution",
      subtitle = subtitle_text,
      x = NULL, y = "Weight"
    ) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.line     = ggplot2::element_line(color = "#2C3E50", linewidth = 0.5),
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = base_size - 1, color = "#7F8C8D"),
      axis.text.x   = ggplot2::element_text(face = "bold", angle = 45, hjust = 1)
    )

  # ── caption 警告 (仅 RH 模式) ──
  if (!is.null(caption_text))
    p <- p + ggplot2::labs(caption = caption_text) +
      ggplot2::theme(plot.caption = ggplot2::element_text(
        face = "italic", color = "red3", size = base_size - 3, hjust = 0))

  return(p)
}
