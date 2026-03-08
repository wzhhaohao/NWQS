#' Fit Non-linear Weighted Quantile Sum (NWQS) Regression
#'
#' @description
#' The main entry point for fitting Non-linear Weighted Quantile Sum (NWQS) regression models.
#' NWQS evaluates the overall joint effect and component-specific relative importance of a highly
#' correlated mixture, while flexibly accommodating non-linear dose-response relationships (e.g.,
#' threshold effects, U-shaped curves). The method separates penalized representation learning
#' from unpenalized final-stage effect estimation.
#'
#' @details
#' \strong{Algorithmic Architecture:} \cr
#' The `nwqs` function implements a robust "Repeated Holdout" (RH) framework combined with a
#' shape-decoupling architecture:
#' \enumerate{
#'   \item \strong{Outer Repeated Holdout Splitting:} For each RH iteration, the data are randomly
#'     split into a Training set (for weight/shape discovery) and a Validation set (for effect estimation).
#'   \item \strong{Weight & Shape Discovery:} On the Training set, a multi-parameter engine
#'     (specified by \code{weight_engine}) estimates the non-linear basis coefficients (shapes) and
#'     extracts relative importance (weights) via Out-Of-Bag (OOB) permutation.
#'   \item \strong{Shape Normalization:} To weaken the coupling between spline shape and effect
#'     magnitude, shapes are standardized to have a unit variance on the predictor scale using the Training set.
#'   \item \strong{1-DoF Effect Estimation:} The normalized shapes and weights are projected onto
#'     the Validation set to construct a single \code{wqs_score}. A generalized linear model (GLM)
#'     is then fitted to estimate the overall mixture effect without directly penalizing the final effect.
#' }
#'
#' By setting \code{rh > 1}, the algorithm averages the weights, shapes, and coefficients across
#' multiple random splits, producing stable empirical inferences. Note that the standard errors
#' and confidence intervals generated from RH iterations reflect data-splitting variance, not true
#' sampling variance. For robust statistical inference, use \code{nwqs_boot()}.
#'
#' @param data A \code{data.frame} containing the mixture components, covariates, and outcome variable.
#' @param mix_name Character vector. Column names of the mixture components to be evaluated.
#' @param covariates Character vector. Column names of covariates/confounders to adjust for. If
#'   \code{NULL}, an unadjusted model is fitted.
#' @param outcome Character. Column name of the dependent/outcome variable. Defaults to \code{"y"}.
#' @param weight_engine Function. The core engine used for weight and shape discovery on the training
#'   set. Defaults to \code{permutation_scorer}.
#' @param q Integer. Number of quantiles used to categorize the continuous mixture components
#'   (e.g., 4 for quartiles, 10 for deciles). Defaults to 4.
#' @param df_spline Integer. Degrees of freedom for the natural cubic splines used to model
#'   non-linearity. Defaults to 3.
#' @param transform_fun Function. A custom transformation applied to the mixture components before
#'   modeling. If \code{NULL}, the default quantile transformation (\code{trans_quantile}) is used.
#' @param split_prop Numeric between (0, 1). The proportion of the full dataset allocated to the
#'   Training set for weight/shape discovery in each RH iteration. Defaults to 0.6.
#' @param rh Integer. Number of Repeated Holdout iterations. To obtain stable distributions, set
#'   \code{rh} to at least 100. Defaults to 1.
#' @param seed Integer. Random seed for reproducible repeated holdout splitting. Defaults to 1234.
#' @param n_permutation Integer. Number of internal permutations or bootstraps passed to
#'   \code{weight_engine} to compute the variable importance metric. Defaults to 30.
#' @param family Character or function. Specifies the GLM error distribution and link function.
#'   Supported options include \code{"gaussian"}, \code{"binomial"}, \code{"poisson"}, and \code{"quasipoisson"}.
#' @param plan_strategy Character. Strategy for parallel computation powered by the \pkg{future}
#'   package. Options are \code{"sequential"}, \code{"multicore"}, or \code{"multisession"}.
#' @param n_workers Integer. Number of parallel workers to use if \code{plan_strategy != "sequential"}.
#'   If \code{NULL}, the package automatically detects and optimizes CPU cores.
#' @param quiet Logical. If \code{TRUE}, suppress the warning about RH-based standard errors
#'   when \code{rh > 1}. Useful when called internally by \code{nwqs_boot()}. Defaults to \code{FALSE}.
#' @param ... Additional arguments passed to \code{run_oob_permutation} or the selected \code{weight_engine}.
#'
#' @return An object of class \code{c("nwqs", "list")} containing the following key components:
#' \itemize{
#'   \item \code{fit}: A list containing the summarized GLM object (coefficients, Deviance, AIC).
#'   \item \code{final_weights}: A named numeric vector of the ensemble-averaged relative weights.
#'   \item \code{mean_shapes}: A named numeric vector of the ensemble-averaged normalized spline coefficients.
#'   \item \code{rh_coefs}: A matrix of the global GLM coefficients across all \code{rh} iterations.
#'   \item \code{rh_weights}: A matrix of the extracted weights across all \code{rh} iterations.
#'   \item \code{data}: The original dataset augmented with the final calculated ensemble \code{wqs_score}.
#' }
#'
#' @seealso \code{\link{nwqs_boot}}, \code{\link{plot.nwqs}}, \code{\link{summary.nwqs}}
#'
#' @importFrom stats glm coef AIC pnorm sd median as.formula
#' @importFrom splines ns
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#' @export
nwqs <- function(data, mix_name, covariates = NULL, outcome = "y",
                 weight_engine = permutation_scorer, q = 4, df_spline = 3, transform_fun = NULL,
                 split_prop = 0.6, rh = 1, seed = 1234, n_permutation = 30,
                 family = c("gaussian", "binomial", "poisson", "quasipoisson"),
                 plan_strategy = c("sequential", "multisession", "multicore"),
                 n_workers = NULL, quiet = FALSE, ...) {
  # ──────────────────────────────────────────────────────────────────────────
  # [FIX #7] quiet 参数：允许 nwqs_boot() 内部调用时静默 RH 警告
  # ──────────────────────────────────────────────────────────────────────────
  t_start <- Sys.time()
  args <- list(...)
  family <- match.arg(family)
  plan_strategy <- match.arg(plan_strategy)

  if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
    stop("Please install 'future' and 'future.apply' packages.")
  }
  if (rh < 1) stop("'rh' must be at least 1.")
  if (split_prop <= 0 || split_prop >= 1) stop("'split_prop' must be in (0, 1).")

  current_reserve <- if (!is.null(args$cpu_reserve)) args$cpu_reserve else 0.2

  # ──────────────────────────────────────────────────────────────────────────
  # [FIX #7] on.exit 保护：确保 old_plan 为 NULL 时不会报错

  # ──────────────────────────────────────────────────────────────────────────
  old_plan <- tryCatch(
    configure_parallel_plan(
      loop_number = rh, strategy = plan_strategy,
      n_workers = n_workers, reserve_cpu = current_reserve
    ),
    error = function(e) {
      warning(
        "Failed to configure parallel plan: ", conditionMessage(e),
        ". Falling back to sequential."
      )
      NULL
    }
  )
  on.exit(
    {
      if (!is.null(old_plan)) future::plan(old_plan)
    },
    add = TRUE
  )

  use_parallel <- !inherits(future::plan(), "sequential")

  if (is.null(transform_fun)) {
    transform_fun <- function(x) trans_quantile(x, q = q)
  }

  data_Q <- data
  data_Q[mix_name] <- transform_fun(data[mix_name])

  eval_points_std <- 0:(q - 1)
  temp_spline <- splines::ns(eval_points_std, df = df_spline)
  model_knots <- attr(temp_spline, "knots")
  model_boundary <- attr(temp_spline, "Boundary.knots")

  # ──────────────────────────────────────────────────────────────────────────
  # [FIX #3] 种子管理：串行模式下用 set.seed()，并行模式下完全依赖
  #          future.seed = TRUE (L'Ecuyer-CMRG)，不再手动设种子
  # ──────────────────────────────────────────────────────────────────────────
  if (!use_parallel && !is.null(seed)) set.seed(seed)

  n_obs <- nrow(data)

  if (is.null(covariates)) {
    formula_str <- paste(outcome, "~ wqs_score")
  } else {
    missing_cov <- setdiff(covariates, names(data))
    if (length(missing_cov) > 0) {
      stop(paste("Missing covariates:", paste(missing_cov, collapse = ", ")))
    }
    formula_str <- paste(outcome, "~ wqs_score +", paste(covariates, collapse = " + "))
  }
  formula_final <- as.formula(formula_str)

  # --- RH 单次迭代函数 ---
  one_rh <- function(i) {
    idx_all <- sample(seq_len(n_obs))
    n_train <- floor(n_obs * split_prop)
    train_idx <- idx_all[seq_len(n_train)]
    valid_idx <- idx_all[(n_train + 1):n_obs]

    data_train <- data_Q[train_idx, , drop = FALSE]
    data_valid <- data_Q[valid_idx, , drop = FALSE]

    boot_res <- run_oob_permutation(
      data = data_train, mix_name = mix_name, outcome = outcome,
      weight_engine = weight_engine, n_permutation = n_permutation,
      q = q, df_spline = df_spline,
      model_knots = model_knots, model_boundary = model_boundary,
      boot_strategy = "sequential",
      ...
    )
    valid_res <- Filter(Negate(is.null), boot_res)
    if (length(valid_res) == 0) {
      return(NULL)
    }

    w_matrix_iter <- do.call(rbind, lapply(valid_res, function(x) x$weights))
    s_matrix_iter <- do.call(rbind, lapply(valid_res, function(x) x$shapes))

    mean_weights_iter <- colMeans(w_matrix_iter, na.rm = TRUE)
    mean_shapes_iter <- colMeans(s_matrix_iter, na.rm = TRUE)

    if (!all(is.finite(mean_weights_iter)) || sum(mean_weights_iter, na.rm = TRUE) <= 0) {
      return(NULL)
    }

    final_weights_iter <- mean_weights_iter / sum(mean_weights_iter)

    train_trans <- wqs_nonlinear_expand(
      data_train, mix_name,
      knots = model_knots, boundary = model_boundary
    )
    valid_trans <- wqs_nonlinear_expand(
      data_valid, mix_name,
      knots = model_knots, boundary = model_boundary
    )

    # ──────────────────────────────────────────────────────────────────────
    # [FIX #5] 防御性检查：确保展开后的列名与期望的 shape 名一致
    # ──────────────────────────────────────────────────────────────────────
    expected_cols <- names(mean_shapes_iter)
    actual_train_cols <- colnames(train_trans)
    actual_valid_cols <- colnames(valid_trans)
    if (!all(expected_cols %in% actual_train_cols) ||
      !all(expected_cols %in% actual_valid_cols)) {
      warning(
        "Column mismatch in wqs_nonlinear_expand() output at RH iteration ", i,
        ". Skipping this iteration."
      )
      return(NULL)
    }

    normalized_shapes_iter <- numeric(length(mean_shapes_iter))
    names(normalized_shapes_iter) <- names(mean_shapes_iter)

    combined_coefs <- numeric(length(mean_shapes_iter))
    names(combined_coefs) <- names(mean_shapes_iter)

    for (comp in mix_name) {
      comp_cols <- paste0(comp, "_B", 1:df_spline)
      theta_raw <- mean_shapes_iter[comp_cols]

      partial_eta <- as.vector(as.matrix(train_trans[, comp_cols, drop = FALSE]) %*% theta_raw)
      sd_eta <- sd(partial_eta, na.rm = TRUE)
      if (is.na(sd_eta) || sd_eta < 1e-8) sd_eta <- 1

      theta_norm <- theta_raw / sd_eta
      normalized_shapes_iter[comp_cols] <- theta_norm
      combined_coefs[comp_cols] <- theta_norm * final_weights_iter[comp]
    }

    wqs_score <- as.matrix(valid_trans[, expected_cols, drop = FALSE]) %*% combined_coefs
    data_valid$wqs_score <- as.vector(wqs_score)

    fit <- glm(formula_final, data = data_valid, family = family)
    aic_val <- if (family == "quasipoisson") NA_real_ else AIC(fit)

    list(
      fit_obj = fit, weights = final_weights_iter, shapes = normalized_shapes_iter,
      coefs = coef(fit), aic = aic_val, null_dev = fit$null.deviance,
      res_dev = fit$deviance, df_null = fit$df.null, df_res = fit$df.residual
    )
  }

  # ──────────────────────────────────────────────────────────────────────────
  # [FIX #3] 并行时传 future.seed，串行时已经在上面 set.seed 了
  # ──────────────────────────────────────────────────────────────────────────
  if (use_parallel) {
    rh_results <- future.apply::future_lapply(
      seq_len(rh), one_rh,
      future.seed = if (!is.null(seed)) seed else TRUE
    )
  } else {
    rh_results <- lapply(seq_len(rh), one_rh)
  }

  rh_results <- Filter(Negate(is.null), rh_results)
  if (length(rh_results) == 0) stop("All iterations failed.")

  if (rh == 1) {
    final_w_global <- rh_results[[1]]$weights
    final_s_global <- rh_results[[1]]$shapes
  } else {
    weight_mat_temp <- do.call(rbind, lapply(rh_results, function(x) x$weights))
    mean_weights_temp <- colMeans(weight_mat_temp, na.rm = TRUE)
    final_w_global <- mean_weights_temp / sum(mean_weights_temp)

    shape_mat_temp <- do.call(rbind, lapply(rh_results, function(x) x$shapes))
    final_s_global <- colMeans(shape_mat_temp, na.rm = TRUE)
  }

  full_trans <- wqs_nonlinear_expand(data_Q, mix_name, knots = model_knots, boundary = model_boundary)

  # ──────────────────────────────────────────────────────────────────────────
  # [FIX #5] 防御性检查：full_trans 列名匹配
  # ──────────────────────────────────────────────────────────────────────────
  expected_cols_full <- names(final_s_global)
  if (!all(expected_cols_full %in% colnames(full_trans))) {
    stop("Column mismatch between final shapes and wqs_nonlinear_expand() output on full data.")
  }

  combined_coefs_full <- numeric(length(final_s_global))
  names(combined_coefs_full) <- names(final_s_global)

  for (comp in mix_name) {
    comp_cols <- paste0(comp, "_B", 1:df_spline)
    combined_coefs_full[comp_cols] <- final_s_global[comp_cols] * final_w_global[comp]
  }

  final_data <- data
  final_data$wqs_score <- as.vector(
    as.matrix(full_trans[, expected_cols_full, drop = FALSE]) %*% combined_coefs_full
  )

  if (rh == 1) {
    single_res <- rh_results[[1]]
    final_obj <- single_res$fit_obj
    coef_summary <- as.data.frame(summary(final_obj)$coefficients)

    fit_obj <- list(
      coefficients = coef_summary, aic = single_res$aic, deviance = single_res$res_dev,
      null.deviance = single_res$null_dev, df.residual = single_res$df_res, df.null = single_res$df_null
    )

    results <- list(
      fit = fit_obj, final_weights = single_res$weights, mean_coefs = single_res$coefs,
      mean_shapes = single_res$shapes, rh_coefs = t(as.matrix(single_res$coefs)),
      rh_weights = t(as.matrix(single_res$weights)), rh_shapes = t(as.matrix(single_res$shapes)),
      rh = 1, b = n_permutation, q = q, df_spline = df_spline, family = family,
      spline_knots = model_knots, spline_boundary = model_boundary, call = match.call(), data = final_data
    )

    class(results) <- c("nwqs", "list")
    return(results)
  }

  coef_mat <- do.call(rbind, lapply(rh_results, function(x) x$coefs))
  weight_mat <- do.call(rbind, lapply(rh_results, function(x) x$weights))
  shape_mat <- do.call(rbind, lapply(rh_results, function(x) x$shapes))

  mean_coefs <- colMeans(coef_mat, na.rm = TRUE)
  mean_weights <- colMeans(weight_mat, na.rm = TRUE)
  mean_weights <- mean_weights / sum(mean_weights)
  mean_shapes <- colMeans(shape_mat, na.rm = TRUE)

  if (family == "quasipoisson") {
    mean_aic <- NA_real_
  } else {
    mean_aic <- mean(vapply(rh_results, function(x) x$aic, numeric(1)), na.rm = TRUE)
  }

  mean_null_dev <- mean(vapply(rh_results, function(x) x$null_dev, numeric(1)), na.rm = TRUE)
  mean_res_dev <- mean(vapply(rh_results, function(x) x$res_dev, numeric(1)), na.rm = TRUE)

  df_null <- rh_results[[1]]$df_null
  df_res <- rh_results[[1]]$df_res

  coef_mean <- colMeans(coef_mat, na.rm = TRUE)
  coef_sd <- apply(coef_mat, 2, sd, na.rm = TRUE)

  # ──────────────────────────────────────────────────────────────────────────
  # [FIX #1] quiet 参数控制是否发出 RH 警告
  # ──────────────────────────────────────────────────────────────────────────
  if (!quiet) {
    warning(
      "When rh > 1, the Standard Errors and 95% CIs in `fit$coefficients` ",
      "represent ONLY the data-splitting (algorithmic) variance across holdout ",
      "iterations, NOT true sampling variance. They are inherently too narrow ",
      "and should NOT be used for statistical inference. Please use `nwqs_boot()` ",
      "to obtain valid percentile bootstrap confidence intervals."
    )
  }

  # ──────────────────────────────────────────────────────────────────────────
  # [FIX #4] coef_summary 结构统一：补充 z_value 和 p_value 列，
  #          使 rh > 1 与 rh == 1 输出的 coefficients 表结构一致
  # ──────────────────────────────────────────────────────────────────────────
  z_value <- ifelse(coef_sd > 0, coef_mean / coef_sd, NA_real_)
  p_value <- ifelse(is.na(z_value), NA_real_, 2 * pnorm(-abs(z_value)))

  coef_summary <- data.frame(
    Estimate     = coef_mean,
    `Std. Error` = coef_sd,
    `z value`    = z_value,
    `Pr(>|z|)`   = p_value,
    `2.5 %`      = coef_mean - 1.96 * coef_sd,
    `97.5 %`     = coef_mean + 1.96 * coef_sd,
    check.names  = FALSE
  )

  fit_obj <- list(
    coefficients = coef_summary, aic = mean_aic, deviance = mean_res_dev,
    null.deviance = mean_null_dev, df.residual = df_res, df.null = df_null
  )

  results <- list(
    fit = fit_obj, final_weights = mean_weights, mean_coefs = mean_coefs, mean_shapes = mean_shapes,
    rh_coefs = coef_mat, rh_weights = weight_mat, rh_shapes = shape_mat,
    rh = rh, b = n_permutation, q = q, df_spline = df_spline, family = family, call = match.call(),
    transform_fun = transform_fun, data = final_data, spline_knots = model_knots, spline_boundary = model_boundary
  )

  class(results) <- c("nwqs", "list")
  return(results)
}


#' Bootstrap CI Wrapper for NWQS
#'
#' @description
#' Performs outer bootstrap resampling for NWQS to approximate true sampling variability.
#' Instead of computing just one contrast, this function extracts ALL terms and targets
#' simultaneously and provides a comprehensive formatted table.
#'
#' @details
#' \strong{Algorithmic Note on `rh_inner`:} \cr
#' The parameter \code{rh_inner} controls the number of Repeated Holdout (RH) splits
#' inside each individual bootstrap iteration. Setting \code{rh_inner = 1} (default)
#' is highly recommended. The outer bootstrap is sufficient for deriving valid
#' confidence intervals.
#'
#' @param data data.frame. Original dataset.
#' @param mix_name Character vector. Mixture component names.
#' @param covariates Character vector or NULL. Covariates to adjust for.
#' @param outcome Character. Outcome variable name. Defaults to "y".
#' @param family Character. One of "gaussian", "binomial", "poisson", "quasipoisson".
#' @param n_boot Integer. Number of outer bootstrap replicates. Defaults to 100.
#' @param rh_inner Integer. Number of RH iterations used inside each \code{nwqs()} fit.
#'   Defaults to 1 to prevent severe computational overhead.
#' @param conf_level Numeric. Confidence level. Defaults to 0.95.
#' @param seed Integer or NULL. Random seed.
#' @param keep_fits Logical. Whether to save all bootstrap \code{nwqs} model objects
#'   in memory. Defaults to FALSE to save memory.
#' @param plan_strategy Character. Parallel strategy for outer bootstrap.
#' @param n_workers Integer or NULL. Number of workers for outer bootstrap.
#' @param ... Additional arguments passed to \code{nwqs()} (e.g., \code{q},
#'   \code{df_spline}, \code{weight_engine}).
#'
#' @return A list with class \code{c("nwqs_boot", "list")} containing:
#' \itemize{
#'   \item \code{point_fit}: NWQS fit on the original data.
#'   \item \code{point_effects}: Filtered \code{extract_nwqs_effects()} output.
#'   \item \code{ci_table}: Long-format summary table with bootstrap percentile CIs.
#'   \item \code{formatted_table}: Publication-ready wide-format table.
#'   \item \code{boot_table}: Raw estimates across all bootstrap iterations.
#'   \item \code{boot_fits}: Optional list of bootstrap model objects (if \code{keep_fits=TRUE}).
#'   \item \code{conf_level}: The confidence level used.
#'   \item \code{final_weights}: Ensemble-averaged weights from the point fit (mirrors nwqs output).
#'   \item \code{mean_shapes}: Ensemble-averaged spline shapes from the point fit.
#'   \item \code{mean_coefs}: Averaged GLM coefficients from the point fit.
#'   \item \code{family}: GLM family used.
#'   \item \code{q}: Number of quantiles used.
#'   \item \code{df_spline}: Spline degrees of freedom.
#'   \item \code{spline_knots}: Knot positions from the point fit.
#'   \item \code{spline_boundary}: Boundary knots from the point fit.
#'   \item \code{data}: Original dataset augmented with wqs_score from point fit.
#'   \item \code{call}: The matched call.
#' }
#' @seealso \code{\link{nwqs}}, \code{\link{plot.nwqs_boot}}, \code{\link{extract_nwqs_effects}}
#' @importFrom stats aggregate quantile
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#' @export
nwqs_boot <- function(data,
                      mix_name,
                      covariates = NULL,
                      outcome = "y",
                      family = c("gaussian", "binomial", "poisson", "quasipoisson"),
                      n_boot = 100,
                      rh_inner = 1,
                      conf_level = 0.95,
                      seed = 1234,
                      keep_fits = FALSE,
                      plan_strategy = c("sequential","multisession","multicore"),
                      n_workers = NULL,
                      ...) {
  start_time <- Sys.time()

  family <- match.arg(family)
  plan_strategy <- match.arg(plan_strategy)

  if (n_boot < 20) {
    warning("'n_boot' is quite small; bootstrap percentile CI may be unstable.")
  }
  if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
    stop("Please install 'future' and 'future.apply' packages.")
  }

  # ── 1. Point Estimate on Original Data ──────────────────────────────────────
  # ──────────────────────────────────────────────────────────────────────────
  # [FIX #1] 使用 quiet = TRUE 避免 point_fit 的 RH 警告（当 rh_inner > 1）
  # ──────────────────────────────────────────────────────────────────────────
  point_fit <- nwqs(
    data = data, mix_name = mix_name, covariates = covariates, outcome = outcome,
    family = family, rh = rh_inner, quiet = TRUE, ...
  )

  point_effects <- extract_nwqs_effects(point_fit)
  cols_to_keep <- c("Target", "Term", "Estimate")
  point_effects_clean <- point_effects[, names(point_effects) %in% cols_to_keep, drop = FALSE]

  # ── 2. Parallel Configuration ────────────────────────────────────────────────
  # ──────────────────────────────────────────────────────────────────────────
  # [FIX #7] on.exit 保护：tryCatch 包裹 configure_parallel_plan
  # ──────────────────────────────────────────────────────────────────────────
  old_plan <- tryCatch(
    configure_parallel_plan(
      loop_number = n_boot,
      strategy    = plan_strategy,
      n_workers   = n_workers
    ),
    error = function(e) {
      warning(
        "Failed to configure parallel plan: ", conditionMessage(e),
        ". Falling back to sequential."
      )
      NULL
    }
  )
  on.exit(
    {
      if (!is.null(old_plan)) future::plan(old_plan)
    },
    add = TRUE
  )

  # ──────────────────────────────────────────────────────────────────────────
  # [FIX #2] 种子管理：串行模式下才用 set.seed()；并行完全靠 future.seed
  # ──────────────────────────────────────────────────────────────────────────
  use_parallel <- !inherits(future::plan(), "sequential")
  if (!use_parallel && !is.null(seed)) set.seed(seed)

  n_obs <- nrow(data)
  alpha <- 1 - conf_level

  # ── 3. Single Bootstrap Iteration ───────────────────────────────────────────
  one_boot <- function(b) {
    idx_boot <- sample(seq_len(n_obs), size = n_obs, replace = TRUE)
    data_boot <- data[idx_boot, , drop = FALSE]

    # ──────────────────────────────────────────────────────────────────────
    # [FIX #2] 内部 nwqs() 不再手动传 seed，让 future.seed 管理随机性
    # [FIX #1] 内部使用 quiet = TRUE 避免每次 bootstrap 迭代都发警告
    # ──────────────────────────────────────────────────────────────────────
    fit_b <- tryCatch(
      {
        nwqs(
          data = data_boot, mix_name = mix_name, covariates = covariates, outcome = outcome,
          family = family, plan_strategy = "sequential", rh = rh_inner,
          seed = NULL, quiet = TRUE, ...
        )
      },
      error = function(e) NULL
    )

    if (is.null(fit_b)) {
      return(list(Success = FALSE, Effects = NULL, Fit = NULL))
    }

    eff_b <- tryCatch(extract_nwqs_effects(fit_b), error = function(e) NULL)

    if (is.null(eff_b)) {
      return(list(Success = FALSE, Effects = NULL, Fit = if (keep_fits) fit_b else NULL))
    }

    eff_b_clean <- eff_b[, names(eff_b) %in% cols_to_keep, drop = FALSE]
    list(Success = TRUE, Effects = eff_b_clean, Fit = if (keep_fits) fit_b else NULL)
  }

  # ── 4. Execute Bootstrap Loop ────────────────────────────────────────────────
  # ──────────────────────────────────────────────────────────────────────────
  # [FIX #2] 并行时用 future.seed = seed 确保可复现
  # ──────────────────────────────────────────────────────────────────────────
  if (use_parallel) {
    boot_results <- future.apply::future_lapply(
      seq_len(n_boot), one_boot,
      future.seed = if (!is.null(seed)) seed else TRUE
    )
  } else {
    boot_results <- lapply(seq_len(n_boot), one_boot)
  }

  boot_success <- vapply(boot_results, function(x) x$Success, logical(1))
  n_success <- sum(boot_success)

  if (n_success < max(20, ceiling(0.5 * n_boot))) {
    warning("A large proportion of bootstrap fits failed. Bootstrap CI may be unstable.")
  }
  if (n_success == 0) {
    stop("All bootstrap replicates failed.")
  }

  # ── 5. Aggregate Results ─────────────────────────────────────────────────────
  valid_effs <- lapply(boot_results[boot_success], function(x) x$Effects)

  # ──────────────────────────────────────────────────────────────────────────
  # [FIX #8] Boot_ID 使用原始的 bootstrap 迭代编号，而非成功序号
  # ──────────────────────────────────────────────────────────────────────────
  success_indices <- which(boot_success)
  for (i in seq_along(valid_effs)) {
    valid_effs[[i]]$Boot_ID <- success_indices[i]
  }
  all_effs_df <- do.call(rbind, valid_effs)

  ci_lower <- aggregate(Estimate ~ Target + Term,
    data = all_effs_df,
    FUN = function(x) quantile(x, probs = alpha / 2, na.rm = TRUE)
  )
  names(ci_lower)[names(ci_lower) == "Estimate"] <- "Boot_CI_Lower"

  ci_upper <- aggregate(Estimate ~ Target + Term,
    data = all_effs_df,
    FUN = function(x) quantile(x, probs = 1 - alpha / 2, na.rm = TRUE)
  )
  names(ci_upper)[names(ci_upper) == "Estimate"] <- "Boot_CI_Upper"

  ci_table <- merge(point_effects_clean, ci_lower, by = c("Target", "Term"), all.x = TRUE)
  ci_table <- merge(ci_table, ci_upper, by = c("Target", "Term"), all.x = TRUE)
  ci_table$N_Success <- n_success

  # ── 6. Formatted Table ───────────────────────────────────────────────────────
  ci_table$Formatted <- sprintf(
    "%.3f [%.3f, %.3f]",
    ci_table$Estimate,
    ci_table$Boot_CI_Lower,
    ci_table$Boot_CI_Upper
  )

  formatted_table <- reshape(
    ci_table[, c("Term", "Target", "Formatted")],
    idvar = "Term", timevar = "Target", direction = "wide"
  )
  names(formatted_table) <- gsub("Formatted\\.", "", names(formatted_table))

  weights <- point_fit$final_weights
  if (!is.null(weights)) {
    weight_df <- data.frame(
      Term = names(weights),
      Weight = round(weights, 3),
      stringsAsFactors = FALSE
    )
    weight_df <- rbind(data.frame(
      Term = "Overall", Weight = NA_real_,
      stringsAsFactors = FALSE
    ), weight_df)
    formatted_table <- merge(weight_df, formatted_table, by = "Term", all.y = TRUE)
    formatted_table <- formatted_table[order(
      formatted_table$Term != "Overall",
      -formatted_table$Weight
    ), ]
    rownames(formatted_table) <- NULL
  }

  boot_fits <- if (keep_fits) lapply(boot_results, function(x) x$Fit) else NULL

  # ── 7. Build Output ──────────────────────────────────────────────────────────
  out <- list(
    point_fit       = point_fit,
    point_effects   = point_effects_clean,
    ci_table        = ci_table,
    formatted_table = formatted_table,
    boot_table      = all_effs_df,
    boot_fits       = boot_fits,
    conf_level      = conf_level,

    # ── Mirror nwqs output fields for downstream compatibility ──────────────
    final_weights   = point_fit$final_weights,
    mean_shapes     = point_fit$mean_shapes,
    mean_coefs      = point_fit$mean_coefs,
    rh_weights      = point_fit$rh_weights,
    rh_shapes       = point_fit$rh_shapes,
    rh_coefs        = point_fit$rh_coefs,
    family          = point_fit$family,
    q               = point_fit$q,
    df_spline       = point_fit$df_spline,
    spline_knots    = point_fit$spline_knots,
    spline_boundary = point_fit$spline_boundary,
    data            = point_fit$data,
    call            = match.call()
  )

  class(out) <- c("nwqs_boot", "list")

  message(sprintf(
    "NWQS Bootstrap completed: %d/%d successful fits. Time taken: %.2f mins",
    n_success, n_boot,
    as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  ))

  return(out)
}
