#' Fit Non-linear Weighted Quantile Sum (NWQS) Regression
#'
#' @description
#' The main entry point for fitting Non-linear Weighted Quantile Sum (NWQS) regression models. 
#' NWQS evaluates the overall joint effect and component-specific relative importance of a 
#' highly correlated mixture, while flexibly accommodating non-linear dose-response relationships 
#' (e.g., threshold effects, U-shaped curves) without suffering from shrinkage bias.
#'
#' @details
#' \strong{Algorithmic Architecture:} \cr
#' The `nwqs` function implements a robust "Repeated Holdout" (RH) framework combined with 
#' a shape-decoupling architecture:
#' \enumerate{
#'   \item \strong{Outer Bootstrap & Splitting:} For each RH iteration, the data is bootstrapped 
#'     to capture sampling variability, then strictly split into a Training set (for weight/shape 
#'     discovery) and a Validation set (for effect estimation).
#'   \item \strong{Weight & Shape Discovery:} On the Training set, a penalized multi-parameter 
#'     engine (specified by \code{model_func}) estimates the non-linear basis coefficients (shapes) 
#'     and extracts relative importance (weights) via Out-Of-Bag (OOB) permutation.
#'   \item \strong{Shape Normalization:} To break the penalization paradox (shrinkage bias), 
#'     NWQS completely decouples the spline shape from the effect magnitude. Shapes are standardized 
#'     to have a unit variance on the predictor scale.
#'   \item \strong{1-DoF Effect Estimation:} The normalized shapes and weights are projected onto 
#'     the Validation set to construct a single \code{wqs_score}. A generalized linear model (GLM) 
#'     is then fitted to estimate the unbiased overall mixture effect and calculate nominal confidence intervals.
#' }
#' 
#' By setting \code{rh > 1}, the algorithm averages the weights, shapes, and coefficients across 
#' multiple random splits, producing highly stable empirical inferences.
#'
#' @param data A \code{data.frame} containing the mixture components, covariates, and outcome variable.
#' @param mix_name Character vector. Column names of the mixture components to be evaluated.
#' @param covariates Character vector. Column names of covariates/confounders to adjust for. 
#'   If \code{NULL}, an unadjusted model is fitted.
#' @param dependent_var Character. Column name of the dependent/outcome variable. Defaults to \code{"y"}.
#' @param model_func Function. The core engine used for weight and shape discovery on the training set. 
#'   Defaults to \code{ridge_permutation_scorer}, which utilizes internal cross-validated L2 penalization.
#' @param q Integer. Number of quantiles used to categorize the continuous mixture components 
#'   (e.g., 4 for quartiles, 10 for deciles). Defaults to 4.
#' @param df_spline Integer. Degrees of freedom for the natural cubic splines used to model 
#'   non-linearity. Defaults to 3.
#' @param split_prop Numeric between (0, 1). The proportion of the bootstrapped dataset allocated 
#'   to the Training set for weight/shape discovery. The remainder goes to the Validation set. Defaults to 0.6.
#' @param seed Integer. Random seed for reproducible data splitting and bootstrap sampling.
#' @param rh Integer. Number of Repeated Holdout (Outer Bootstrap) iterations. To obtain empirical 
#'   confidence intervals and stable distributions, set \code{rh} to at least 100. Defaults to 1.
#' @param family Character or function. Specifies the GLM error distribution and link function. 
#'   Supported options include \code{"gaussian"}, \code{"binomial"}, \code{"poisson"}, and \code{"quasipoisson"}.
#' @param transform_fun Function. A custom transformation applied to the mixture components before 
#'   modeling. If \code{NULL}, the default quantile transformation (\code{trans_quantile}) is used.
#' @param plan_strategy Character. Strategy for parallel computation powered by the \pkg{future} package. 
#'   Options are \code{"sequential"} or \code{"multicore"}.
#' @param n_workers Integer. Number of parallel workers to use if \code{plan_strategy = "multicore"}. 
#'   If \code{NULL}, the package automatically detects and optimizes CPU cores.
#' @param B Integer. Number of internal permutations or bootstraps passed to \code{model_func} 
#'   to compute the variable importance metric. Defaults to 100.
#' @param ... Additional arguments passed to \code{run_oob_permutation} or the selected \code{model_func}.
#'
#' @return An object of class \code{c("nwqs", "list")} containing the following key components:
#' \itemize{
#'   \item \code{fit}: A list containing the summarized GLM object (coefficients, Deviance, AIC).
#'   \item \code{final_weights}: A named numeric vector of the ensemble-averaged relative weights for each mixture component.
#'   \item \code{mean_shapes}: A named numeric vector of the ensemble-averaged normalized spline coefficients.
#'   \item \code{rh_coefs}: A matrix of the global GLM coefficients across all \code{rh} iterations (used for empirical CIs).
#'   \item \code{rh_weights}: A matrix of the extracted weights across all \code{rh} iterations.
#'   \item \code{data}: The original dataset augmented with the final calculated ensemble \code{wqs_score}.
#' }
#' 
#' @seealso \code{\link{plot.nwqs}}, \code{\link{summary.nwqs}}, \code{\link{print.nwqs}}
#' 
#' @importFrom stats glm coef AIC pnorm sd median as.formula
#' @importFrom splines ns
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#' @export
nwqs <- function(data,
                 mix_name,
                 covariates = NULL,
                 dependent_var = "y",
                 model_func = ridge_permutation_scorer,
                 q = 4,
                 df_spline = 3,
                 split_prop = 0.6,
                 seed = 1234,
                 rh = 1,
                 family = c("gaussian", "binomial", "poisson", "quasipoisson"),
                 transform_fun = NULL,
                 plan_strategy = c("sequential", "multicore", "multicore"),
                 n_workers = NULL,
                 B = 100,
                 ...) {
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
  old_plan <- configure_parallel_plan(loop_number = rh, strategy = plan_strategy, n_workers = n_workers, reserve_cpu = current_reserve)
  on.exit(future::plan(old_plan), add = TRUE)
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

  if (!is.null(seed)) set.seed(seed)
  n_obs <- nrow(data)

  if (is.null(covariates)) {
    formula_str <- paste(dependent_var, "~ wqs_score")
  } else {
    missing_cov <- setdiff(covariates, names(data))
    if (length(missing_cov) > 0) stop(paste("Missing covariates:", paste(missing_cov, collapse = ", ")))
    formula_str <- paste(dependent_var, "~ wqs_score +", paste(covariates, collapse = " + "))
  }
  formula_final <- as.formula(formula_str)

  # --- 3. RH 单次迭代函数 (核心重构区) ---
  one_rh <- function(i) {
    # =========================================================================
    # [统计学终极修复]: Outer Bootstrap
    # 有放回地抽取 n_obs 个样本，模拟真实的总体抽样变异
    # =========================================================================
    boot_idx <- sample(seq_len(n_obs), size = n_obs, replace = TRUE)
    data_boot <- data_Q[boot_idx, , drop = FALSE]

    # 在这个带有抽样变异的新样本集 (data_boot) 上，再切分 Train / Valid
    train_idx <- sample(seq_len(n_obs), size = floor(n_obs * split_prop))
    data_train <- data_boot[train_idx, , drop = FALSE]
    data_valid <- data_boot[-train_idx, , drop = FALSE]
    # =========================================================================

    boot_res <- run_oob_permutation(
      data = data_train, mix_name = mix_name, dependent_var = dependent_var,
      model_func = model_func, B = 100, q = q, df_spline = df_spline,
      model_knots = model_knots, model_boundary = model_boundary, ...
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

    valid_trans <- wqs_nonlinear_expand(data_valid, mix_name, knots = model_knots, boundary = model_boundary)

    # =========================================================================
    # [架构彻底重构]: Shape Normalization (解耦形状与幅度)
    # =========================================================================
    normalized_shapes_iter <- numeric(length(mean_shapes_iter))
    names(normalized_shapes_iter) <- names(mean_shapes_iter)

    combined_coefs <- numeric(length(mean_shapes_iter))
    names(combined_coefs) <- names(mean_shapes_iter)

    for (comp in mix_name) {
      comp_cols <- paste0(comp, "_B", 1:df_spline)
      theta_raw <- mean_shapes_iter[comp_cols]

      partial_eta <- as.vector(as.matrix(valid_trans[, comp_cols]) %*% theta_raw)
      sd_eta <- sd(partial_eta, na.rm = TRUE)

      if (is.na(sd_eta) || sd_eta < 1e-8) sd_eta <- 1

      theta_norm <- theta_raw / sd_eta
      normalized_shapes_iter[comp_cols] <- theta_norm
      combined_coefs[comp_cols] <- theta_norm * final_weights_iter[comp]
    }
    # =========================================================================

    wqs_score <- as.matrix(valid_trans) %*% combined_coefs
    data_valid$wqs_score <- as.vector(wqs_score)

    fit <- glm(formula_final, data = data_valid, family = family)

    list(
      fit_obj = fit,
      weights = final_weights_iter,
      shapes = normalized_shapes_iter,
      coefs = coef(fit),
      aic = AIC(fit),
      null_dev = fit$null.deviance,
      res_dev = fit$deviance,
      df_null = fit$df.null,
      df_res = fit$df.residual
    )
  }

  rh_results <- if (use_parallel) {
    future.apply::future_lapply(seq_len(rh), one_rh, future.seed = TRUE)
  } else {
    lapply(seq_len(rh), one_rh)
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

  combined_coefs_full <- numeric(length(final_s_global))
  names(combined_coefs_full) <- names(final_s_global)
  for (comp in mix_name) {
    comp_cols <- paste0(comp, "_B", 1:df_spline)
    combined_coefs_full[comp_cols] <- final_s_global[comp_cols] * final_w_global[comp]
  }

  final_data <- data
  final_data$wqs_score <- as.vector(as.matrix(full_trans) %*% combined_coefs_full)

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
      rh = 1, b = B, q = q, df_spline = df_spline, family = family,
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
    mean_aic <- NA
  } else {
    mean_aic <- mean(vapply(rh_results, function(x) x$aic, numeric(1)), na.rm = TRUE)
  }
  mean_null_dev <- mean(vapply(rh_results, function(x) x$null_dev, numeric(1)), na.rm = TRUE)
  mean_res_dev <- mean(vapply(rh_results, function(x) x$res_dev, numeric(1)), na.rm = TRUE)
  df_null <- rh_results[[1]]$df_null
  df_res <- if (!is.null(rh_results[[1]]$df_res)) rh_results[[1]]$df_res else rh_results[[1]]$df_residual

  coef_mean <- colMeans(coef_mat)
  coef_sd <- apply(coef_mat, 2, sd)
  coef_summary <- data.frame(
    Estimate = coef_mean, `Std. Error` = coef_sd,
    `2.5 %` = coef_mean - 1.96 * coef_sd, `97.5 %` = coef_mean + 1.96 * coef_sd, check.names = FALSE
  )

  fit_obj <- list(
    coefficients = coef_summary, aic = mean_aic, deviance = mean_res_dev,
    null.deviance = mean_null_dev, df.residual = df_res, df.null = df_null
  )

  results <- list(
    fit = fit_obj, final_weights = mean_weights, mean_coefs = mean_coefs, mean_shapes = mean_shapes,
    rh_coefs = coef_mat, rh_weights = weight_mat, rh_shapes = shape_mat,
    rh = rh, b = B, q = q, df_spline = df_spline, family = family, call = match.call(),
    transform_fun = transform_fun, data = final_data,
    spline_knots = model_knots, spline_boundary = model_boundary
  )
  class(results) <- c("nwqs", "list")
  return(results)
}
