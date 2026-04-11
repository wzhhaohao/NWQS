#' @title 拟合非线性加权分位数和 (NWQS) 回归模型
#'
#' @description
#' `nwqs` 是拟合非线性加权分位数和 (Non-linear Weighted Quantile Sum) 回归模型的核心函数。
#' 本方法旨在评估高共线性混合物暴露的整体联合效应 (overall joint effect) 以及各个组分的
#' 相对重要性，同时灵活地容纳非线性剂量反应关系（例如：阈值效应、U型曲线等）。
#' 算法在架构上将带惩罚的表征学习（权重与形状发现）与无惩罚的最终效应估计进行了严格分离。
#'
#' @details
#' \strong{算法架构与统计学考量:} \cr
#' 本函数实现了一个稳健的“重复保留 (Repeated Holdout, RH)”架构，并结合了形状解耦机制：
#' \enumerate{
#'   \item \strong{外部重复保留拆分 (Outer Repeated Holdout Splitting):} 在每次 RH 迭代中，数据被随机
#'     拆分为训练集（用于发现权重和非线性形状）和验证集（用于效应估计）。
#'     *医学统计学注:* 对于条件逻辑回归 (`family = "clogit"`)，算法严格在层/匹配组（Cluster/Stratum）级别
#'     进行随机拆分，这对于维持配对病例对照研究 (Matched Case-Control Studies) 的设计完整性至关重要，可避免破坏匹配结构。
#'   \item \strong{权重与形状发现:} 在训练集上，多参数引擎估计非线性基函数的系数（形状），并通过
#'     袋外 (Out-Of-Bag, OOB) 置换法提取相对重要性（权重）。
#'   \item \strong{形状归一化 (Shape Normalization):} 为了削弱样条形状与效应量之间的耦合，利用训练集将
#'     形状标准化为在预测变量尺度上具有单位方差。
#'   \item \strong{单自由度效应估计 (1-DoF Effect Estimation):} 将归一化后的形状和权重投影到验证集上，
#'     构建单一的 `nwqs` 指数。随后拟合广义线性模型 (GLM) 以估计整体混合物效应，从而避免对最终效应的直接惩罚。
#' }
#'
#' \strong{关于科学严谨性与偏倚风险的警告:} \cr
#' \itemize{
#'   \item \strong{方差膨胀与推断错误:} 当设置 \code{rh > 1} 时，算法会平均多次拆分的结果以提供稳定的经验推断。
#'     但必须严谨地认识到：此处输出的 \code{fit$coefficients} 中的标准误 (SE) 和置信区间 (CI) **仅反映了
#'     数据拆分带来的算法方差 (Algorithmic Variance)，绝不能代表真实的抽样变异 (Sampling Variance)**。
#'     直接使用该 P 值进行假设检验会导致第一类错误严重膨胀。进行可靠的推断必须使用 \code{nwqs_boot()}。
#'   \item \strong{残余混杂 (Residual Confounding):} 请务必通过 \code{covariates} 传入充分的混杂因素。
#'     如果存在未测量的混杂因素，或者模型对非线性关系的设定存在严重偏差，估计的联合效应可能产生偏倚。
#'   \item \strong{选择偏倚 (Selection Bias):} 如果原始样本中的数据缺失是“非随机缺失 (MNAR)”，或者
#'     训练集/验证集的随机拆分导致小样本下某些暴露特征的分布失衡，可能会在权重分配上引入偏倚。
#' }
#'
#' @param data \code{data.frame}。包含混合物组分、协变量、匹配变量以及结局变量的数据框。
#' @param mix_name Character vector。需要评估的混合物组分（暴露变量）的列名。
#' @param covariates Character vector。用于调整的协变量/混杂因素的列名。若为 \code{NULL}，则拟合未调整模型。
#' @param outcome Character。结局/因变量的列名，默认值为 \code{"y"}。
#' @param strata_col Character。用于条件逻辑回归的层/匹配组 ID 的列名。若 \code{family = "clogit"}，此项为必填。
#' @param weight_engine Function。用于在训练集上发现权重和形状的核心引擎，默认为 \code{permutation_scorer}。
#' @param q Integer。将连续混合物变量分类的分位数数量（如 4 代表四分位数，10 代表十分位数），默认为 4。
#' @param df_spline Integer。用于拟合非线性曲线的自然三次样条 (Natural Cubic Splines) 的自由度，默认为 3。
#' @param transform_fun Function。在建模前对混合物组分应用的自定义转换函数。若为 \code{NULL}，默认使用分位数转换 (\code{trans_quantile})。
#' @param train_prop Numeric (0, 1)。在单次 RH 迭代中，分配给训练集（用于权重/形状发现）的数据比例，默认为 0.6。
#' @param rh Integer。重复保留 (Repeated Holdout) 的迭代次数。若要获得稳定的分布，建议设置为 100 以上。默认为 1。
#' @param seed Integer。用于可重复的数据拆分的随机种子，默认为 1234。
#' @param n_permutation Integer。传递给 \code{weight_engine} 计算变量重要性的内部置换或 Bootstrap 次数，默认为 10。
#' @param family Character 包含拟合 GLM 的误差分布及链接函数。支持的选项包括 \code{"gaussian"}（线性回归）,
#'   \code{"binomial"}（逻辑回归）, \code{"poisson"}, \code{"quasipoisson"}, 以及 \code{"clogit"}（条件逻辑回归）。
#' @param plan_strategy Character。基于 \pkg{future} 包的并行计算策略。可选 \code{"sequential"}, \code{"multicore"}, 或 \code{"multisession"}。
#' @param n_workers Integer。若不采取 \code{"sequential"}，指定并行的核心数。若为 \code{NULL}，将自动检测并优化 CPU 核心分配。
#' @param quiet Logical。若为 \code{TRUE}，则静默关于 \code{rh > 1} 时标准误不代表抽样方差的警告信息。通常在 \code{nwqs_boot()} 内部调用时开启，默认为 \code{FALSE}。
#' @param ... 传递给 \code{run_oob_permutation} 或所选 \code{weight_engine} 的其他额外参数。
#'
#' @return 返回一个 \code{c("nwqs", "list")} 类的对象，包含以下核心内容：
#' \itemize{
#'   \item \code{fit}: 包含汇总 GLM 对象（系数、偏差、AIC）的列表。
#'   \item \code{final_weights}: 命名数值向量，表示多次迭代平均后的集成相对权重。
#'   \item \code{mean_shapes}: 命名数值向量，表示多次迭代平均后的标准化样条系数。
#'   \item \code{rh_coefs}: 矩阵，记录所有 \code{rh} 迭代中的全局 GLM 回归系数。
#'   \item \code{rh_weights}: 矩阵，记录所有 \code{rh} 迭代中提取的权重。
#'   \item \code{data}: 附加了最终计算出的集成 \code{nwqs} 指数得分的原始数据集。
#' }
#'
#' @seealso \code{\link{nwqs_boot}}, \code{\link{plot.nwqs}}, \code{\link{summary.nwqs}}
#'
#' @importFrom stats glm coef AIC pnorm sd median as.formula
#' @importFrom splines ns
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#' @export
nwqs <- function(data, mix_name, covariates = NULL, outcome = "y", strata_col = NULL,
                 weight_engine = permutation_scorer, q = 4, df_spline = 3, transform_fun = NULL,
                 train_prop = 0.6, rh = 100, seed = 1234, n_permutation = 100,
                 family = c("gaussian", "binomial", "poisson", "quasipoisson", "clogit"),
                 plan_strategy = c("sequential", "multisession", "multicore"),
                 n_workers = NULL, quiet = FALSE, ...) {
  family <- match.arg(family)
  plan_strategy <- match.arg(plan_strategy)
  if (length(covariates) == 0) covariates <- NULL

  # ──────────────────────────────────────────────────────────────────────────
  # [NEW] Clogit 安全检查
  # ──────────────────────────────────────────────────────────────────────────
  if (family == "clogit") {
    if (is.null(strata_col)) {
      stop("For conditional logistic regression (family = 'clogit'), 'strata_col' must be provided.")
    }
    if (!requireNamespace("survival", quietly = TRUE)) {
      stop("Please install the 'survival' package to use clogit.")
    }
    if (!(strata_col %in% names(data))) {
      stop(paste("strata_col '", strata_col, "' not found in data.", sep = ""))
    }
  }

  # ──────────────────────────────────────────────────────────────────────────
  # [FIX #7] quiet 参数：允许 nwqs_boot() 内部调用时静默 RH 警告
  # ──────────────────────────────────────────────────────────────────────────
  t_start <- Sys.time()
  args <- list(...)

  if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
    stop("Please install 'future' and 'future.apply' packages.")
  }
  if (rh < 1) stop("'rh' must be at least 1.")
  if (train_prop <= 0 || train_prop >= 1) stop("'train_prop' must be in (0, 1).")

  current_reserve <- if (!is.null(args$cpu_reserve)) args$cpu_reserve else 0.2

  # ──────────────────────────────────────────────────────────────────────────
  # [FIX #7] on.exit 保护：确保 old_plan 为 NULL 时不会报错
  # ──────────────────────────────────────────────────────────────────────────
  old_plan <- tryCatch(
    configure_parallel_plan(
      loop_number = rh, strategy = plan_strategy,
      n_workers = n_workers, reserve_cpu = current_reserve,
      verbose = !isTRUE(quiet)
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
  # [FIX #3] 种子管理
  # ──────────────────────────────────────────────────────────────────────────
  if (!use_parallel && !is.null(seed)) set.seed(seed)

  n_obs <- nrow(data)

  # ──────────────────────────────────────────────────────────────────────────
  # [NEW] Clogit 公式构建：注入 strata()
  # ──────────────────────────────────────────────────────────────────────────
  if (is.null(covariates)) {
    if (family == "clogit") {
      formula_str <- paste0(outcome, " ~ nwqs + strata(", strata_col, ")")
    } else {
      formula_str <- paste(outcome, "~ nwqs")
    }
  } else {
    missing_cov <- setdiff(covariates, names(data))
    if (length(missing_cov) > 0) {
      stop(paste("Missing covariates:", paste(missing_cov, collapse = ", ")))
    }
    if (family == "clogit") {
      formula_str <- paste0(outcome, " ~ nwqs + ", paste(covariates, collapse = " + "), " + strata(", strata_col, ")")
    } else {
      formula_str <- paste(outcome, "~ nwqs +", paste(covariates, collapse = " + "))
    }
  }
  formula_final <- as.formula(formula_str)

  # --- RH 单次迭代函数 ---
  one_rh <- function(i) {
    # ──────────────────────────────────────────────────────────────────────
    # [NEW] 如果是 Clogit，必须按组（Cluster）进行训练和验证集的拆分！
    # ──────────────────────────────────────────────────────────────────────
    if (family == "clogit") {
      unique_strata <- unique(data_Q[[strata_col]])
      n_strata <- length(unique_strata)
      shuffled_strata <- sample(unique_strata)

      n_train_strata <- max(1, floor(n_strata * train_prop))
      train_strata <- shuffled_strata[seq_len(n_train_strata)]
      valid_strata <- shuffled_strata[(n_train_strata + 1):n_strata]

      train_idx <- which(data_Q[[strata_col]] %in% train_strata)
      valid_idx <- which(data_Q[[strata_col]] %in% valid_strata)
    } else {
      idx_all <- sample(seq_len(n_obs))
      n_train <- floor(n_obs * train_prop)
      train_idx <- idx_all[seq_len(n_train)]
      valid_idx <- idx_all[(n_train + 1):n_obs]
    }

    data_train <- data_Q[train_idx, , drop = FALSE]
    data_valid <- data_Q[valid_idx, , drop = FALSE]

    vars_needed <- unique(c(mix_name, outcome, covariates))
    if (family == "clogit") vars_needed <- unique(c(vars_needed, strata_col))

    boot_res <- run_oob_permutation(
      data = data_train[, vars_needed, drop = FALSE],
      mix_name = mix_name,
      outcome = outcome,
      covariates = covariates,
      weight_engine = weight_engine,
      n_permutation = n_permutation,
      q = q,
      df_spline = df_spline,
      model_knots = model_knots,
      model_boundary = model_boundary,
      family = family,
      boot_strategy = "sequential",
      strata_col = strata_col, # 确保底层置换引擎也能识别匹配组
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
    # [FIX #5] 防御性检查
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

    nwqs <- as.matrix(valid_trans[, expected_cols, drop = FALSE]) %*% combined_coefs
    data_valid$nwqs <- as.vector(nwqs)

    # ──────────────────────────────────────────────────────────────────────
    # [NEW] 兼容 Clogit 和 GLM 的拟合与指标萃取
    # ──────────────────────────────────────────────────────────────────────
    if (family == "clogit") {
      fit <- survival::clogit(formula_final, data = data_valid)
      aic_val <- stats::extractAIC(fit)[2]
      null_dev <- -2 * fit$loglik[1]
      res_dev <- -2 * fit$loglik[2]
      df_n <- length(fit$coefficients)
      df_r <- fit$n - length(fit$coefficients)
    } else {
      fit <- glm(formula_final, data = data_valid, family = family)
      aic_val <- if (family == "quasipoisson") NA_real_ else AIC(fit)
      null_dev <- fit$null.deviance
      res_dev <- fit$deviance
      df_n <- fit$df.null
      df_r <- fit$df.residual
    }

    coefs_fit <- coef(fit)

    if (!("nwqs" %in% names(coefs_fit)) || !is.finite(unname(coefs_fit["nwqs"]))) {
      return(NULL)
    }

    if (!is.finite(stats::sd(data_valid$nwqs, na.rm = TRUE)) ||
      stats::sd(data_valid$nwqs, na.rm = TRUE) < 1e-8) {
      return(NULL)
    }

    list(
      fit_obj = fit, weights = final_weights_iter, shapes = normalized_shapes_iter,
      coefs = coefs_fit, aic = aic_val, null_dev = null_dev,
      res_dev = res_dev, df_null = df_n, df_res = df_r
    )
  }

  if (use_parallel) {
    rh_results <- future.apply::future_lapply(
      seq_len(rh), one_rh,
      future.seed = if (!is.null(seed)) seed else TRUE
    )
  } else {
    rh_results <- lapply(seq_len(rh), one_rh)
  }

  rh_results <- Filter(function(x) {
    !is.null(x) &&
      !is.null(x$coefs) &&
      "nwqs" %in% names(x$coefs) &&
      is.finite(unname(x$coefs["nwqs"]))
  }, rh_results)
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
  final_data$nwqs <- as.vector(
    as.matrix(full_trans[, expected_cols_full, drop = FALSE]) %*% combined_coefs_full
  )

  if (rh == 1) {
    single_res <- rh_results[[1]]
    final_obj <- single_res$fit_obj
    coef_summary <- as.data.frame(summary(final_obj)$coefficients)

    # ──────────────────────────────────────────────────────────────────────────
    # [NEW] Clogit 的 summary 结构矫正 (与 glm 统一命名)
    # ──────────────────────────────────────────────────────────────────────────
    if (family == "clogit") {
      if ("coef" %in% names(coef_summary)) names(coef_summary)[names(coef_summary) == "coef"] <- "Estimate"
      if ("se(coef)" %in% names(coef_summary)) names(coef_summary)[names(coef_summary) == "se(coef)"] <- "Std. Error"
      if ("z" %in% names(coef_summary)) names(coef_summary)[names(coef_summary) == "z"] <- "z value"
      if ("Pr(>|z|)" %in% names(coef_summary) == FALSE && "p" %in% names(coef_summary)) {
        names(coef_summary)[names(coef_summary) == "p"] <- "Pr(>|z|)"
      }
    }

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

  if (!quiet) {
    warning(
      "When rh > 1, the Standard Errors and 95% CIs in `fit$coefficients` ",
      "represent ONLY the data-splitting (algorithmic) variance across holdout ",
      "iterations, NOT true sampling variance. They are inherently too narrow ",
      "and should NOT be used for statistical inference. Please use `nwqs_boot()` ",
      "to obtain valid percentile bootstrap confidence intervals."
    )
  }

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


#' @title NWQS 模型的 Bootstrap 置信区间估计
#'
#' @description
#' `nwqs_boot` 为 NWQS 方法执行外部 Bootstrap 重抽样，旨在近似估计模型参数的真实抽样变异 (Sampling Variability)。
#' 与仅计算单一对比差异不同，此函数会同时提取所有模型项的效应，并提供具备学术出版级质量的置信区间表格。
#'
#' @details
#' \strong{重抽样策略与统计学注意事项:} \cr
#' \itemize{
#'   \item \strong{独立样本 vs. 匹配数据:} 算法根据所选的误差分布自适应抽样策略。对于标准 GLM 模型，执行简单的随机有放回抽样。
#'     而对于条件逻辑回归 (`family = "clogit"`)，算法执行\strong{整群重抽样 (Cluster Resampling)}。这意味着它会基于
#'     匹配变量 (`strata_col`) 作为一个整体进行抽样，并在内部动态重构匹配 ID。这在流行病学设计中极为关键，防止了同一病例-对照组被
#'     拆散或在多次抽中时混淆独立观察单位。
#'   \item \strong{内部 RH 参数 (`rh_inner`):} 此参数控制单个 Bootstrap 样本内部执行的重复保留次数。
#'     强烈建议在 Bootstrap 时将其设置为 1（默认值），以防止因“双重嵌套循环”引发的灾难性计算开销。外部 Bootstrap 足以提供
#'     推断所需的稳健分布。
#'   \item \strong{解释与选择偏倚:} Bootstrap 百分位法得出的 95\% CI 提供的是统计抽样方差的良好估计。然而，请注意它**无法纠正**
#'     原始设计中的系统性选择偏倚。如果原始队列/样本并不代表目标人群，或者失访率较高且与暴露结局均相关，则通过 Bootstrap 生成的
#'     置信区间会精确地“围绕”在这个有偏的估计量周围。
#'   \item \strong{流行病学效应量化:} 对于指数族模型（`binomial`, `poisson`, `quasipoisson`, `clogit`），算法会自动计算并
#'     报告指数化后的平均值与置信区间上限/下限，以便直接解释为比值比 (OR) 或相对危险度 (RR)。
#' }
#'
#' @param data data.frame。原始数据集。
#' @param mix_name Character vector。混合物组分列名。
#' @param covariates Character vector 或 \code{NULL}。需要调整的协变量/混杂因素。
#' @param outcome Character。结局变量列名，默认为 "y"。
#' @param strata_col Character。匹配数据中条件逻辑回归所需的层/组 ID 列名。
#' @param family Character。误差分布与连接函数，可选 "gaussian", "binomial", "poisson", "quasipoisson", "clogit"。
#' @param n_boot Integer。外部 Bootstrap 抽样重复的次数。考虑到百分位置信区间的稳定性，建议在最终发表前设置为 500 或 1000（默认为 100）。
#' @param rh_inner Integer。在每个 `nwqs()` 拟合内部使用的 RH 迭代次数。为避免计算量激增，默认为 1。
#' @param n_permutation Integer。传递给底层 `nwqs()` 函数的置换次数，用于计算变量重要性权重。默认为 100。
#' @param conf_level Numeric。所需置信区间的置信水平，默认为 0.95 (即 95\% 置信区间)。
#' @param seed Integer 或 \code{NULL}。随机数种子，用于保证抽样的可重复性。
#' @param keep_fits Logical。是否在内存中保留所有的 Bootstrap `nwqs` 模型对象。通常为了节省 RAM 设为 \code{FALSE}。
#' @param plan_strategy Character。外部 Bootstrap 所采用的并行计算策略。
#' @param n_workers Integer 或 \code{NULL}。并行核心数。
#' @param ... 传递给核心 \code{nwqs()} 函数的其他参数（例如分位点数量 \code{q}，样条自由度 \code{df_spline} 等）。
#'
#' @return 一个 \code{c("nwqs_boot", "list")} 类的列表，主要包含以下内容供下游统计报告使用：
#' \itemize{
#'   \item \code{ci_table}: 包含均值及 Bootstrap 百分位置信区间的长格式数据框。若适用，自动报告 OR/RR 及其置信区间。
#'   \item \code{formatted_table}: 可直接用于医学期刊附表的宽格式汇总表（带有权重降序排列）。
#'   \item \code{boot_table}: 包含所有成功 Bootstrap 迭代中提取的原始估计值的长格式数据。
#'   \item \code{final_weights}: 从原始点估计拟合中获取的混合物组分平均相对权重。
#'   \item \code{n_success}: 成功收敛的 Bootstrap 迭代次数。
#'   \item ... 及其他包含数据、公式调用和内部架构（样条结点）的辅助属性。
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
                      strata_col = NULL,
                      family = c("gaussian", "binomial", "poisson", "quasipoisson", "clogit"),
                      n_boot = 100,
                      rh_inner = 1,
                      n_permutation = 100,
                      conf_level = 0.95,
                      seed = 1234,
                      keep_fits = FALSE,
                      plan_strategy = c("sequential", "multisession", "multicore"),
                      n_workers = NULL,
                      q = 4,
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

  cols_to_keep <- c("Target", "Term", "Estimate")

  # ── 1. Parallel Configuration ────────────────────────────────────────────────
  old_plan <- tryCatch(
    configure_parallel_plan(
      loop_number = n_boot,
      strategy    = plan_strategy,
      n_workers   = n_workers,
      verbose     = FALSE
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
  if (!use_parallel && !is.null(seed)) set.seed(seed)

  n_obs <- nrow(data)
  alpha <- 1 - conf_level

  # ── 2. Single Bootstrap Iteration ───────────────────────────────────────────
  one_boot <- function(b) {
    # ────────────────────────────────────────────────────────────────────────
    # Clogit 的整群重抽样与独立 ID 重构
    # ────────────────────────────────────────────────────────────────────────
    if (family == "clogit") {
      unique_strata <- unique(data[[strata_col]])
      sampled_strata <- sample(unique_strata, size = length(unique_strata), replace = TRUE)

      boot_data_list <- lapply(seq_along(sampled_strata), function(idx) {
        sub_data <- data[data[[strata_col]] == sampled_strata[idx], , drop = FALSE]
        # 极其关键：防止同一匹配组被抽中多次而导致合并
        sub_data[[strata_col]] <- paste0(sampled_strata[idx], "_boot_", idx)
        return(sub_data)
      })
      data_boot <- do.call(rbind, boot_data_list)
    } else {
      idx_boot <- sample(seq_len(n_obs), size = n_obs, replace = TRUE)
      data_boot <- data[idx_boot, , drop = FALSE]
    }

    fit_b <- tryCatch(
      {
        nwqs(
          data = data_boot, mix_name = mix_name, covariates = covariates, outcome = outcome, q = q,
          strata_col = strata_col, family = family, plan_strategy = plan_strategy, rh = rh_inner,
          n_permutation = n_permutation, seed = NULL, quiet = TRUE, ...
        )
      },
      error = function(e) {
        message("Bootstrap ", b, " failed: ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(fit_b)) {
      return(list(Success = FALSE, Effects = NULL, Fit = NULL))
    }

    eff_b <- tryCatch(extract_nwqs_effects(fit_b), error = function(e) NULL)

    if (is.null(eff_b)) {
      return(list(Success = FALSE, Effects = NULL, Fit = if (keep_fits) fit_b else NULL))
    }

    eff_b_clean <- eff_b[, names(eff_b) %in% cols_to_keep, drop = FALSE]

    # 直接提取单次迭代的权重、形状和骨架参数
    list(
      Success = TRUE,
      Effects = eff_b_clean,
      Weights = fit_b$final_weights,
      Shapes = fit_b$mean_shapes,
      Struct = list(
        df_spline = fit_b$df_spline,
        spline_knots = fit_b$spline_knots,
        spline_boundary = fit_b$spline_boundary
      ),
      Fit = if (keep_fits) fit_b else NULL
    )
  }

  # ── 3. Execute Bootstrap Loop ────────────────────────────────────────────────
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

  # ── 4. Aggregate Results ─────────────────────────────────────────────────────
  valid_results <- boot_results[boot_success]

  # 获取基础结构信息（从第一个成功的迭代中提取）
  first_struct <- valid_results[[1]]$Struct

  # 提取所有成功的效应值
  valid_effs <- lapply(valid_results, function(x) x$Effects)
  success_indices <- which(boot_success)
  for (i in seq_along(valid_effs)) {
    valid_effs[[i]]$Boot_ID <- success_indices[i]
  }
  all_effs_df <- do.call(rbind, valid_effs)

  # 计算平均权重
  valid_weights <- lapply(valid_results, function(x) x$Weights)
  weights_mat <- do.call(rbind, valid_weights)
  avg_weights <- colMeans(weights_mat, na.rm = TRUE)
  avg_weights <- avg_weights / sum(avg_weights, na.rm = TRUE) # 归一化

  # 计算平均形状
  valid_shapes <- lapply(valid_results, function(x) x$Shapes)

  # 增加 is.numeric() 判断，因为 mean_shapes 通常是数值向量
  if (is.numeric(valid_shapes[[1]]) || is.data.frame(valid_shapes[[1]]) || is.matrix(valid_shapes[[1]])) {
    avg_shapes <- Reduce("+", valid_shapes) / length(valid_shapes)
  } else {
    avg_shapes <- NULL
    warning("Complex shape structure detected; could not average shapes.")
  }

  # ---- Bootstrap CI ----
  ci_lower <- aggregate(
    Estimate ~ Target + Term,
    data = all_effs_df,
    FUN = function(x) quantile(x, probs = alpha / 2, na.rm = TRUE)
  )
  names(ci_lower)[names(ci_lower) == "Estimate"] <- "Boot_CI_Lower"

  ci_upper <- aggregate(
    Estimate ~ Target + Term,
    data = all_effs_df,
    FUN = function(x) quantile(x, probs = 1 - alpha / 2, na.rm = TRUE)
  )
  names(ci_upper)[names(ci_upper) == "Estimate"] <- "Boot_CI_Upper"

  # ---- Bootstrap mean ----
  boot_mean <- aggregate(
    Estimate ~ Target + Term,
    data = all_effs_df,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  names(boot_mean)[names(boot_mean) == "Estimate"] <- "Boot_Mean"

  # ---- Merge CI + boot mean (彻底摒弃点估计) ----
  ci_table <- merge(boot_mean, ci_lower, by = c("Target", "Term"), all.x = TRUE)
  ci_table <- merge(ci_table, ci_upper, by = c("Target", "Term"), all.x = TRUE)
  ci_table$N_Success <- n_success

  # ──────────────────────────────────────────────────────────────────────────
  # Clogit 等指数族还原 (展示为 OR/RR)
  # ──────────────────────────────────────────────────────────────────────────
  is_exp_family <- family %in% c("binomial", "poisson", "quasipoisson", "clogit")

  disp_lcl <- ci_table$Boot_CI_Lower
  disp_ucl <- ci_table$Boot_CI_Upper
  disp_mean <- ci_table$Boot_Mean

  if (is_exp_family) {
    disp_lcl <- exp(disp_lcl)
    disp_ucl <- exp(disp_ucl)
    disp_mean <- exp(disp_mean)
  }

  ci_table$Formatted <- sprintf(
    "%.3f [%.3f, %.3f]",
    disp_mean, disp_lcl, disp_ucl
  )

  formatted_table <- reshape(
    ci_table[, c("Term", "Target", "Formatted")],
    idvar = "Term", timevar = "Target", direction = "wide"
  )
  names(formatted_table) <- gsub("^Formatted\\.", "", names(formatted_table))

  # 合并平均权重到格式化表格
  weight_df <- data.frame(
    Term = names(avg_weights),
    Weight = round(avg_weights, 3),
    stringsAsFactors = FALSE
  )
  weight_df <- rbind(
    data.frame(Term = "Overall", Weight = NA_real_, stringsAsFactors = FALSE),
    weight_df
  )

  formatted_table <- merge(weight_df, formatted_table, by = "Term", all.y = TRUE)
  formatted_table <- formatted_table[order(
    formatted_table$Term != "Overall",
    -formatted_table$Weight,
    na.last = TRUE
  ), ]
  rownames(formatted_table) <- NULL

  boot_raw <- all_effs_df
  if (is_exp_family) {
    boot_raw$Estimate <- exp(boot_raw$Estimate)
  }
  boot_raw$ColName <- paste(boot_raw$Term, boot_raw$Target, sep = "_")

  boot_contrast_matrix <- reshape(
    boot_raw[, c("Boot_ID", "ColName", "Estimate")],
    idvar = "Boot_ID", timevar = "ColName", direction = "wide"
  )
  names(boot_contrast_matrix) <- gsub("^Estimate\\.", "", names(boot_contrast_matrix))
  rownames(boot_contrast_matrix) <- NULL

  boot_fits <- if (keep_fits) lapply(valid_results, function(x) x$Fit) else NULL

  # ── 5. Output Construction ───────────────────────────────────────────────────
  out <- list(
    ci_table = ci_table,
    formatted_table = formatted_table,
    boot_table = all_effs_df,
    boot_contrast_mat = boot_contrast_matrix,
    boot_fits = boot_fits,
    conf_level = conf_level,
    n_boot = n_boot,
    n_success = n_success,
    rh_inner = rh_inner,
    n_permutation = n_permutation,
    final_weights = avg_weights,
    mean_shapes = avg_shapes,
    family = family,
    q = q,
    df_spline = first_struct$df_spline,
    spline_knots = first_struct$spline_knots,
    spline_boundary = first_struct$spline_boundary,
    data = data,
    call = match.call()
  )

  class(out) <- c("nwqs_boot", "list")

  message(sprintf(
    "NWQS Bootstrap completed: %d/%d successful fits. Time taken: %.2f mins",
    n_success, n_boot,
    as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  ))

  return(out)
}
