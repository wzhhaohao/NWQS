#' @importFrom stats glm coef AIC pnorm sd median as.formula
#' @importFrom future plan
#' @importFrom future.apply future_lapply

#' Repeated Holdout Non-linear / Linear WQS Regression (Dual Mode)
#' 重复保留非线性/线性 WQS 回归 (双模式)
#'
#' @description
#' The main entry point for fitting NWQS models. It employs a "Repeated Holdout" validation framework
#' to ensure robust estimates and avoid overfitting.
#' \cr
#' NWQS 模型拟合的主入口函数。它采用“重复保留（Repeated Holdout）”验证框架，以确保估计的稳健性并避免过拟合。
#'
#' @details
#' **Dual Mode Output (双模式输出):**
#' \itemize{
#'   \item **If rh = 1**: Behaves like a standard regression. Returns a single GLM object fitted on the validation set, with an added `final_weights` attribute.
#'   \item **If rh > 1**: Performs Pooled Inference. Returns averaged coefficients (Mean/Median), averaged weights, and pooled model performance metrics (AIC, Deviance) across all holdout iterations. This avoids "Double Dipping" (re-fitting on full data) and provides valid statistical inference.
#' }
#' \cr
#' **双模式输出：**
#' \itemize{
#'   \item **若 rh = 1**：行为类似于标准回归。返回在验证集上拟合的单个 GLM 对象，并附带 `final_weights` 属性。
#'   \item **若 rh > 1**：执行池化推断（Pooled Inference）。返回所有保留迭代中的平均系数（均值/中位数）、平均权重和汇总的模型性能指标（AIC, Deviance）。这避免了“数据二次使用（Double Dipping）”，提供了有效的统计推断。
#' }
#'
#' @param data data.frame. Full dataset.
#'   完整数据框。
#' @param mix_name character vector. Names of the mixture components.
#'   混合物组分名称。
#' @param covariates character vector. Covariate names to include in the model.
#'   模型中包含的协变量名称。
#' @param dependent_var character. Name of the outcome variable. Defaults to "y".
#'   因变量名称。
#' @param model_func function. Function to calculate weights (e.g., `calc_spline_wqs_weights`).
#'   计算权重的函数。
#' @param q integer. Number of quantiles for data transformation (default 4).
#'   分位数变换的数量（默认 4）。
#' @param df_spline integer. Degrees of freedom for spline expansion (default 3).
#'   样条展开的自由度（默认 3）。
#' @param split_prop numeric. Proportion of data used for training/weight discovery (default 0.6).
#'   用于训练/权重发现的数据比例（默认 0.6）。
#' @param B integer. Number of bootstrap samples inside each holdout iteration (default 100).
#'   每次保留迭代内部的 Bootstrap 样本数（默认 100）。
#' @param seed integer. Random seed for reproducibility.
#'   随机种子。
#' @param rh integer. Number of Repeated Holdout iterations.
#'   If 1, returns a single fit. If >1, returns pooled results.
#'   重复保留迭代次数。若为 1 返回单次拟合结果；若 >1 返回汇总结果。
#' @param family character/function. GLM family (e.g., "gaussian", "binomial").
#'   GLM 分布族。
#' @param transform_fun function. Custom transformation function. If NULL, uses quantiles.
#'   自定义变换函数。
#' @param plan_strategy character. Parallel strategy ("sequential", "multisession", "multicore").
#'   并行策略。
#' @param n_workers integer. Number of workers. NULL triggers auto-optimization.
#'   工作核心数。NULL 触发自动优化。
#' @param ... Additional arguments passed to `run_bootstrap` or `model_func`.
#'   传递给 `run_bootstrap` 或 `model_func` 的额外参数。
#'
#' @return An object of class `nwqs_result`. Structure depends on `rh`.
#'   返回 `nwqs_result` 类对象，结构取决于 `rh`。
#'
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
                family = c("gaussian", "binomial"),
                transform_fun = NULL,
                plan_strategy = c("sequential", "multisession", "multicore"),
                n_workers = NULL,
                force_inner_sequential_when_nested = TRUE,
                ...) {

  t_start = Sys.time()

  # --- 0. 环境与参数检查 ---
  family = match.arg(family)
  plan_strategy = match.arg(plan_strategy)

  if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
    stop("Please install 'future' and 'future.apply' packages.")
  }
  if (rh < 1) stop("'rh' must be at least 1.")
  if (split_prop <= 0 || split_prop >= 1) stop("'split_prop' must be in (0, 1).")

  extra_args = list(...)

  # --- 1. 并行环境配置 (重构后) ---
  # 这里将 rh 传给 loop_number 用于负载均衡计算
  old_plan = configure_parallel_plan(
    loop_number = rh, 
    strategy = plan_strategy, 
    n_workers = n_workers
  )
  
  # 注册环境还原 (函数退出时自动执行)
  on.exit(future::plan(old_plan), add = TRUE)

  use_parallel = !inherits(future::plan(), "sequential")

  # --- 2. 预处理 ---
  if (is.null(transform_fun)) {
    current_q = q
    transform_fun = function(x) trans_quantile(x, q = current_q)
  }

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
    # ... (此处逻辑保持不变，为节省篇幅省略，请保留你原有的完整逻辑) ...
    # A. 数据切分
    train_idx = sample(seq_len(n_obs), size = floor(n_obs * split_prop))
    data_train = data[train_idx, , drop = FALSE]
    data_valid = data[-train_idx, , drop = FALSE]

    # B. 训练集 bootstrap (强制串行)
    boot_res = run_bootstrap(
      data = data_train,
      mix_name = mix_name,
      dependent_var = dependent_var,
      model_func = model_func,
      B = 100,
      transform_fun = transform_fun,
      ...
    )

    valid_res = Filter(Negate(is.null), boot_res)
    if (length(valid_res) == 0) return(NULL)

    # C. 聚合权重
    w_matrix_iter = do.call(rbind, valid_res)
    mean_weights_iter = colMeans(w_matrix_iter, na.rm = TRUE)
    if (!all(is.finite(mean_weights_iter)) || sum(mean_weights_iter, na.rm = TRUE) <= 0) return(NULL)
    final_weights_iter = mean_weights_iter / sum(mean_weights_iter)

    # D. 验证集拟合
    valid_trans = wqs_nonlinear_expand(data_valid, mix_name, transform_fun = transform_fun, ...)
    wqs_score = as.matrix(valid_trans) %*% rep(final_weights_iter, each = df_spline)
    data_valid$wqs_score = as.vector(wqs_score)

    fit = glm(formula_final, data = data_valid, family = family)

    list(
      fit_obj = fit,
      weights = final_weights_iter,
      coefs = coef(fit),
      aic = AIC(fit),
      null_dev = fit$null.deviance,
      res_dev = fit$deviance,
      df_null = fit$df.null,
      df_res = fit$df.residual
    )
  }

  # --- 4. 执行 RH 循环 ---
  # use_parallel 标志位已经在 Step D 确定
  rh_results = if (use_parallel) {
    future.apply::future_lapply(seq_len(rh), one_rh, future.seed = TRUE)
  } else {
    lapply(seq_len(rh), one_rh)
  }

  rh_results = Filter(Negate(is.null), rh_results)
  if (length(rh_results) == 0) stop("All iterations failed.")

  # --- 5. 结果输出 (rh=1 vs rh>1) ---
  if (rh == 1) {
    single_res = rh_results[[1]]
    final_obj = single_res$fit_obj
    final_obj$final_weights = single_res$weights
    final_obj$rh = 1
    final_obj$call = match.call()
    class(final_obj) = c("nwqs_result", class(final_obj))
    return(final_obj)
  }

  # Pooled Inference
  coef_mat = do.call(rbind, lapply(rh_results, function(x) x$coefs))
  weight_mat = do.call(rbind, lapply(rh_results, function(x) x$weights))

  mean_coefs = colMeans(coef_mat, na.rm = TRUE)
  median_coefs = apply(coef_mat, 2, median, na.rm = TRUE)

  mean_weights = colMeans(weight_mat, na.rm = TRUE)
  mean_weights = mean_weights / sum(mean_weights)

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
    family = family
  )

  class(result) = "nwqs_result"
  t_end = Sys.time()
  duration = difftime(t_end, t_start, units = "auto")
  message(sprintf("\n=== NWQS model finished in %.2f %s ===", as.numeric(duration), units(duration)))
  return(result)
}


# TODO: Evaluate model performance against beta_preds
# Metrics needed:
# 1. MAE (Mean Absolute Error)
# 2. MSE (Mean Squared Error) = mean((est - true)^2)
# 3. Bias = mean(est - true)