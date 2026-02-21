#' Repeated Holdout Non-linear / Linear WQS Regression (Dual Mode)
#' 重复保留非线性/线性 WQS 回归 (双模式)
#'
#' @description
#' The main entry point for fitting NWQS models. It employs a "Repeated Holdout" validation framework
#' to ensure robust estimates and avoid overfitting.
#'
#' @param data data.frame. Full dataset.
#' @param mix_name character vector. Names of the mixture components.
#' @param covariates character vector. Covariate names to include in the model.
#' @param dependent_var character. Name of the outcome variable. Defaults to "y".
#' @param model_func function. Function to calculate weights (e.g., `calc_spline_wqs_weights`).
#' @param q integer. Number of quantiles for data transformation (default 4).
#' @param df_spline integer. Degrees of freedom for spline expansion (default 3).
#' @param split_prop numeric. Proportion of data used for training/weight discovery (default 0.6).
#' @param B integer. Number of bootstrap samples inside each holdout iteration (default 100).
#' @param seed integer. Random seed for reproducibility.
#' @param rh integer. Number of Repeated Holdout iterations.
#' @param family character/function. GLM family (e.g., "gaussian", "binomial").
#' @param transform_fun function. Custom transformation function. If NULL, uses quantiles.
#' @param plan_strategy character. Parallel strategy ("sequential", "multicore", "multicore").
#' @param n_workers integer. Number of workers. NULL triggers auto-optimization.
#' @param ... Additional arguments passed to `run_bootstrap` or `model_func`.
#'
#' @return An object of class `nwqs_result`. Structure depends on `rh`.
#'   Includes `final_data` which contains the original dataset with the calculated `wqs_score`.
#' @importFrom stats glm coef AIC pnorm sd median as.formula
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#' @export
nwqs = function(data,
                mix_name,
                covariates = NULL,
                dependent_var = "y",
                model_func = calc_spline_wqs_weights,
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

  t_start = Sys.time()
  args = list(...)
  family = match.arg(family)
  plan_strategy = match.arg(plan_strategy)

  if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
    stop("Please install 'future' and 'future.apply' packages.")
  }
  if (rh < 1) stop("'rh' must be at least 1.")
  if (split_prop <= 0 || split_prop >= 1) stop("'split_prop' must be in (0, 1).")

  # --- 1. 并行环境配置 (重构后) ---
  current_reserve = if (!is.null(args$cpu_reserve)) args$cpu_reserve else 0.2

  old_plan = configure_parallel_plan(
    loop_number = rh, 
    strategy = plan_strategy, 
    n_workers = n_workers,
    reserve_cpu = current_reserve
  )
  
  on.exit(future::plan(old_plan), add = TRUE)

  use_parallel = !inherits(future::plan(), "sequential")

  # --- 2. 预处理 ---
  if (is.null(transform_fun)) {
    transform_fun = function(x) trans_quantile(x, q = q)
  }

  data_Q = data
  data_Q[mix_name] = transform_fun(data[mix_name])

  if (!is.null(seed)) set.seed(seed)
  n_obs = nrow(data)

  if (is.null(covariates)) {
    formula_str = paste(dependent_var, "~ wqs_score")
  } else {
    missing_cov = setdiff(covariates, names(data))
    if (length(missing_cov) > 0) stop(paste("Missing covariates:", paste(missing_cov, collapse = ", ")))
    formula_str = paste(dependent_var, "~ wqs_score +", paste(covariates, collapse = " + "))
  }
  formula_final = as.formula(formula_str)

  # --- 3. RH 单次迭代函数 ---
  one_rh = function(i) {
    train_idx = sample(seq_len(n_obs), size = floor(n_obs * split_prop))
    data_train = data_Q[train_idx, , drop = FALSE]
    data_valid = data_Q[-train_idx, , drop = FALSE]

    boot_res = run_bootstrap(
      data = data_train,
      mix_name = mix_name,
      dependent_var = dependent_var,
      model_func = model_func,
      B = 100,
      q = q,
      df_spline = df_spline,
      ...
    )

    valid_res = Filter(Negate(is.null), boot_res)
    if (length(valid_res) == 0) return(NULL)

    w_matrix_iter = do.call(rbind, lapply(valid_res, function(x) x$weights))
    s_matrix_iter = do.call(rbind, lapply(valid_res, function(x) x$shapes))

    mean_weights_iter = colMeans(w_matrix_iter, na.rm = TRUE)
    mean_shapes_iter = colMeans(s_matrix_iter, na.rm = TRUE)

    if (!all(is.finite(mean_weights_iter)) || sum(mean_weights_iter, na.rm = TRUE) <= 0) return(NULL)
    
    final_weights_iter = mean_weights_iter / sum(mean_weights_iter)

    valid_trans = wqs_nonlinear_expand(data_valid, mix_name, df_spline = df_spline, q = q)
    weight_expanded = rep(final_weights_iter, each = df_spline)
    combined_coefs = mean_shapes_iter * weight_expanded
    
    wqs_score = as.matrix(valid_trans) %*% combined_coefs
    data_valid$wqs_score = as.vector(wqs_score)

    fit = glm(formula_final, data = data_valid, family = family)

    list(
      fit_obj = fit,
      weights = final_weights_iter,
      shapes = mean_shapes_iter, 
      coefs = coef(fit),
      aic = AIC(fit),
      null_dev = fit$null.deviance,
      res_dev = fit$deviance,
      df_null = fit$df.null,
      df_res = fit$df.residual
    )
  }

  # --- 4. 执行 RH 循环 ---
  rh_results = if (use_parallel) {
    future.apply::future_lapply(seq_len(rh), one_rh, future.seed = TRUE)
  } else {
    lapply(seq_len(rh), one_rh)
  }

  rh_results = Filter(Negate(is.null), rh_results)
  if (length(rh_results) == 0) stop("All iterations failed.")

  # =========================================================================
  # --- 5. [新增] 计算全集数据的最终 WQS Score ------------------------------
  # 提前提取出合并后的 weights 和 shapes，应用于全集数据
  # =========================================================================
  if (rh == 1) {
    final_w_global = rh_results[[1]]$weights
    final_s_global = rh_results[[1]]$shapes
  } else {
    weight_mat_temp = do.call(rbind, lapply(rh_results, function(x) x$weights))
    mean_weights_temp = colMeans(weight_mat_temp, na.rm = TRUE)
    final_w_global = mean_weights_temp / sum(mean_weights_temp)
    
    shape_mat_temp = do.call(rbind, lapply(rh_results, function(x) x$shapes))
    final_s_global = colMeans(shape_mat_temp, na.rm = TRUE)
  }

  # 对全体数据 (data_Q) 进行扩展计算
  full_trans = wqs_nonlinear_expand(data_Q, mix_name, df_spline = df_spline, q = q)
  weight_expanded_full = rep(final_w_global, each = df_spline)
  combined_coefs_full = final_s_global * weight_expanded_full
  
  # 算出来的 WQS Score 拼接回原始数据集 (data)
  final_data = data
  final_data$wqs_score = as.vector(as.matrix(full_trans) %*% combined_coefs_full)
  # =========================================================================


  # --- 6. 结果输出 (rh=1 vs rh>1) ---
  if (rh == 1) {
    single_res = rh_results[[1]]
    final_obj = single_res$fit_obj
    final_obj$final_weights = single_res$weights
    final_obj$shapes = single_res$shapes
    final_obj$final_data = final_data   # [新增] 将全集数据绑定到结果中
    final_obj$rh = 1
    final_obj$call = match.call()
    class(final_obj) = c("nwqs_result", class(final_obj))
    return(final_obj)
  }

  # Pooled Inference
  coef_mat = do.call(rbind, lapply(rh_results, function(x) x$coefs))
  weight_mat = do.call(rbind, lapply(rh_results, function(x) x$weights))
  shape_mat = do.call(rbind, lapply(rh_results, function(x) x$shapes))

  mean_coefs = colMeans(coef_mat, na.rm = TRUE)
  median_coefs = apply(coef_mat, 2, median, na.rm = TRUE)

  mean_weights = colMeans(weight_mat, na.rm = TRUE)
  mean_weights = mean_weights / sum(mean_weights)
  
  mean_shapes = colMeans(shape_mat, na.rm = TRUE)

  # 使用 vapply 确保类型安全
  mean_aic = mean(vapply(rh_results, function(x) x$aic, numeric(1)), na.rm = TRUE)
  mean_null_dev = mean(vapply(rh_results, function(x) x$null_dev, numeric(1)), na.rm = TRUE)
  mean_res_dev = mean(vapply(rh_results, function(x) x$res_dev, numeric(1)), na.rm = TRUE)

  # 处理自由度 (兼容性写法)
  df_null = rh_results[[1]]$df_null
  df_res = if(!is.null(rh_results[[1]]$df_res)) rh_results[[1]]$df_res else rh_results[[1]]$df_residual

  result = list(
    call = match.call(),
    mean_coefs = mean_coefs,
    median_coefs = median_coefs,
    final_weights = mean_weights,
    mean_aic = mean_aic,
    mean_null_dev = mean_null_dev,
    mean_res_dev = mean_res_dev,
    df_null = df_null,
    df_res = df_res,
    rh_coefs = coef_mat,
    rh_weights = weight_mat,
    rh = rh,
    mean_shapes = mean_shapes, 
    rh_shapes = shape_mat,
    family = family,
    final_data = final_data # [新增] 将全集数据绑定到结果列表中
  )

  class(result) = c("nwqs", "list")
  t_end = Sys.time()
  duration = difftime(t_end, t_start, units = "auto")
  message(sprintf("\n=== NWQS model finished in %.2f %s ===", as.numeric(duration), units(duration)))
  return(result)
}



