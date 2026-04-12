#' @title Quantile or Percentile Transformation
#'
#' @description
#' Transforms continuous mixture exposure variables into discrete quantile bins
#' or continuous percentile ranks. This standardization step is fundamental to
#' weighted quantile sum (WQS) and its extensions for harmonizing exposures
#' measured on different scales.
#'
#' @param data \code{data.frame}. Contains the mixture variables to transform.
#' @param method Character. Transformation method: \code{"quantile"} (default)
#'   or \code{"percentile"}.
#' @param q Integer. Number of quantile bins (only used when
#'   \code{method = "quantile"}). Default is 4 (quartiles).
#'
#' @return A \code{data.frame} with the same dimensions and column names as the
#'   input, containing the transformed dimensionless values.
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


#' @title Non-Linear Spline Expansion for WQS Mixture Components
#'
#' @description
#' Transforms (quantile-transformed) mixture variables into a natural cubic
#' spline basis matrix using globally fixed internal and boundary knots.
#'
#' @details
#' This function enforces the use of pre-computed \code{knots} and
#' \code{boundary} to ensure basis function alignment across training,
#' validation, and resampled datasets. Dynamically computing knots within
#' subsamples would cause spatial drift bias.
#'
#' @param data \code{data.frame}. Contains the mixture variables (typically
#'   already transformed).
#' @param mix_name Character vector. Column names of mixture components to
#'   expand.
#' @param df_spline Integer. Degrees of freedom for the natural spline.
#'   Default is 3.
#' @param knots Numeric vector. Internal knot positions. Must be explicitly
#'   provided for global scale alignment.
#' @param boundary Numeric vector (length 2). Boundary knot positions. Must be
#'   explicitly provided.
#'
#' @return A numeric matrix with columns named
#'   \code{{Component}_B{BasisIndex}}.
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


#' @title Compute NWQS Joint Exposure Quantile Contrast Effects
#'
#' @description
#' Evaluates the overall mixture non-linear effect significance and computes
#' joint exposure quantile contrasts (e.g., all components simultaneously at
#' their highest quantile vs. lowest). For logistic and Poisson models, results
#' are automatically converted to odds ratios (OR) and rate ratios (RR).
#'
#' @param model An object of class \code{"nwqs"}.
#' @param q_target Integer. Target quantile index (0-based). For example, 3
#'   represents Q4. If \code{NULL}, automatically set to the maximum quantile.
#' @param q_ref Integer. Reference quantile index (0-based). Default is 0 (Q1).
#'
#' @return Invisibly returns a list containing \code{delta_eta}, \code{lower},
#'   and \code{upper}. When \code{rh = 1}, confidence bounds are \code{NA}.
#'
#' @importFrom splines ns
#' @importFrom stats quantile coef
#' @export
nwqs_contrast <- function(model, q_target = NULL, q_ref = 0) {
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


#' @title Extract Detailed Quantile Contrast Effects from NWQS Model
#'
#' @description
#' Computes all possible quantile contrast effects (e.g., Q2 vs Q1, Q3 vs Q1,
#' Q4 vs Q1), decomposing the overall effect into component-specific
#' contributions. Standard errors and empirical confidence intervals are
#' derived from repeated holdout iterations.
#'
#' @param model_res An object of class \code{"nwqs"}.
#' @param return_raw Logical. Currently unused; reserved for future extensions.
#'
#' @return A \code{data.frame} with columns: Target, Term, Estimate, SE,
#'   Wald_CI_Lower, Wald_CI_Upper, Empirical_CI_Lower, Empirical_CI_Upper.
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


#' @title Plot Component-Specific Bootstrap Contrast Boxplots
#'
#' @description
#' Generates publication-quality diagnostic boxplots showing the distribution
#' of component-specific effects across bootstrap replicates. Jittered points
#' overlay the boxplots to reveal the full distribution of effect estimates.
#'
#' @param model An object of class \code{"nwqs_boot"}.
#' @param exponentiate Logical or \code{NULL}. Whether to exponentiate the
#'   Y-axis to display OR or RR. If \code{NULL}, automatically determined from
#'   the model family.
#' @param free_y Logical. If \code{TRUE} (default), each facet panel has an
#'   independently scaled Y-axis.
#' @param base_size Integer. Base font size. Default is 12.
#' @param fill_alpha Numeric. Box fill transparency. Default is 0.16.
#' @param palette Character. Discrete color palette. Default is
#'   \code{"default"}.
#' @param components Character vector. Specific components to display.
#' @param top_n Integer or \code{NULL}. Show only the top \code{top_n}
#'   components ranked by weight.
#' @param ylim Numeric vector (length 2). Force Y-axis limits.
#' @param y_step Numeric. Force Y-axis tick spacing.
#'
#' @return A publication-quality \code{ggplot} faceted boxplot object.
#'
#' @export
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
    if (n <= 5) return(list(nc = n))
    if (n == 8) return(list(nc = 4))
    if (n == 9) return(list(nc = 3))
    if (n == 10) return(list(nc = 5))
    if (n == 12) return(list(nc = 4))
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
    stop("Input must be an 'nwqs_boot' object.")
  }
  if (is.null(model$boot_table) || nrow(model$boot_table) == 0) {
    stop("model$boot_table is empty; cannot plot.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install 'dplyr': install.packages('dplyr')")
  }

  selected_raw <- .resolve_selected_raw(model, components = components, top_n = top_n)
  if (length(selected_raw) == 0) {
    stop("No components remain after filtering.")
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
