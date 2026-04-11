#' @title 运行袋外置换与并行重抽样调度器 (OOB Permutation Orchestrator)
#'
#' @description
#' `run_oob_permutation` 是 NWQS 架构中的核心调度函数。它负责准备非线性样条的设计矩阵（Design Matrix），
#' 管理多线程并行计算，并向底层的 `weight_engine`（通常为 `permutation_scorer`）分发带有严格随机种子的重抽样任务。
#'
#' @details
#' \strong{数据准备与防偏倚机制:} \cr
#' \itemize{
#'   \item \strong{安全族解析 (Safe Family Parsing):} 针对流行病学中常用的条件逻辑回归 (`clogit`) 进行了特殊拦截，
#'     避免常规 `get()` 函数在调用 `survival` 包对象时引发环境报错。
#'   \item \strong{匹配结构隔离:} 将 `strata_col` 从主数据集中剥离并作为独立向量 `strata_id` 传递。这一操作
#'     确保了在后续生成设计矩阵 (`model.matrix`) 时，匹配 ID 不会被错误地视为常规数值或因子型协变量，
#'     从而避免自由度无意义的消耗和潜在的完全分离 (Complete Separation) 偏倚。
#'   \item \strong{极速路径 (Fast-Path):} 当无协变量（粗模型，Crude Model）时，函数会绕过高开销的 `as.formula`
#'     和 `model.matrix` 解析，直接拼接截距列。当存在协变量（调整模型，Adjusted Model）时，由于需要控制混杂因素 (Confounding factors)，
#'     函数会严谨地生成包含虚拟变量编码 (Dummy Coding) 的完整设计矩阵。
#' }
#'
#' @param data \code{data.frame}。包含所需变量的原始数据集。
#' @param mix_name Character vector。混合物组分的列名。
#' @param outcome Character。结局变量的列名，默认为 \code{"y"}。
#' @param covariates Character vector 或 \code{NULL}。需要调整的协变量/混杂因素。
#' @param weight_engine Function。底层用于计算权重和形状的引擎函数，默认为 \code{\link{permutation_scorer}}。
#' @param n_permutation Integer。内部 Bootstrap 或置换迭代的次数，默认为 100。
#' @param seed Integer 或 \code{NULL}。并行计算的随机数种子，确保基于 future 的并行结果可重现。
#' @param boot_strategy Character。并行策略，可选 \code{"sequential"}, \code{"multicore"}, 或 \code{"multisession"}。
#' @param boot_n_workers Integer 或 \code{NULL}。并行工作节点数。
#' @param ... 传递给样条扩展函数 (\code{wqs_nonlinear_expand}) 或 \code{weight_engine} 的其他核心参数
#'   （如 \code{model_knots}, \code{model_boundary}, \code{q}, \code{df_spline} 等）。
#'
#' @return 一个列表，长度等于 \code{n_permutation}，包含每次迭代中由 \code{weight_engine} 返回的权重和形状向量。
#'
#' @importFrom stats model.matrix as.formula
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#' @export
run_oob_permutation <- function(data, mix_name, outcome = "y",
                                covariates = NULL,
                                weight_engine = permutation_scorer,
                                n_permutation = 100, seed = NULL,
                                boot_strategy = c("sequential", "multicore", "multisession"),
                                boot_n_workers = NULL, ...) {
  args <- list(...)

  # ──────────────────────────────────────────────────────────────────
  # [安全解析 family]：拦截 clogit，避免触发 get() 报错
  # ──────────────────────────────────────────────────────────────────
  fam_arg <- if (!is.null(args$family)) args$family else "gaussian"
  if (is.character(fam_arg) && fam_arg == "clogit") {
    fam_obj <- list(family = "clogit")
  } else if (is.character(fam_arg)) {
    fam_obj <- get(fam_arg, mode = "function")()
  } else {
    fam_obj <- fam_arg
  }
  args$family <- NULL

  # ──────────────────────────────────────────────────────────────────
  # [提取 strata_col]：单独保存为 strata_id
  # ──────────────────────────────────────────────────────────────────
  strata_col <- args$strata_col
  if (!is.null(strata_col)) {
    if (!(strata_col %in% names(data))) {
      stop("strata_col '", strata_col, "' not found in data.")
    }
    strata_id <- data[[strata_col]]
    args$strata_col <- NULL
  } else {
    strata_id <- NULL
  }

  boot_strategy <- match.arg(boot_strategy)

  old_plan <- configure_parallel_plan(
    loop_number = n_permutation,
    strategy = boot_strategy,
    n_workers = boot_n_workers
  )
  on.exit(future::plan(old_plan), add = TRUE)

  use.seed <- if (is.null(seed)) TRUE else seed
  expand_func <- if (!is.null(args$expand_func)) args$expand_func else wqs_nonlinear_expand

  q_val <- if (!is.null(args$q)) args$q else 4
  df_val <- if (!is.null(args$df_spline)) args$df_spline else 3

  if (is.null(args$model_knots) || is.null(args$model_boundary)) {
    stop("'model_knots' and 'model_boundary' must be provided.")
  }

  full_spline_data <- expand_func(
    data, mix_name,
    df_spline = df_val,
    knots = args$model_knots,
    boundary = args$model_boundary
  )
  spline_vars <- colnames(full_spline_data)

  if (is.null(covariates)) covariates <- character(0)

  missing_cov <- setdiff(covariates, names(data))
  if (length(missing_cov) > 0) stop("Missing covariates: ", paste(missing_cov, collapse = ", "))

  if (length(covariates) > 0) {
    # model.matrix is unavoidable here to handle potential dummy coding for factor covariates
    formula_str <- paste("~", paste(c(covariates, spline_vars), collapse = " + "))
    temp_data <- cbind(data[, covariates, drop = FALSE], full_spline_data)
    internal_formula <- as.formula(formula_str)
    X_matrix <- model.matrix(internal_formula, data = temp_data)
  } else {
    # CRUDE MODEL FAST-PATH:
    # Bypass data.frame coercion and formula parsing entirely for maximum speed.
    # We just manually append the intercept to match model.matrix output structure.
    X_matrix <- cbind("(Intercept)" = 1, full_spline_data)
  }

  y_raw <- data[[outcome]]


  if (fam_obj$family == "binomial") {
    if (is.factor(y_raw)) {
      if (nlevels(y_raw) != 2) stop("For binomial family, factor outcome must have exactly 2 levels.")
      y_vector <- as.numeric(y_raw) - 1
    } else {
      y_vector <- as.numeric(y_raw)
    }
  } else {
    y_vector <- as.numeric(y_raw)
  }

  results <- future.apply::future_lapply(
    seq_len(n_permutation),
    function(i) {
      tryCatch(
        {
          do.call(weight_engine, c(list(
            x = X_matrix,
            y = y_vector,
            strata_id = strata_id, # 将匹配组ID无损传给打分引擎
            mix_name = mix_name,
            spline_vars = spline_vars,
            family = fam_obj
          ), args))
        },
        error = function(e) {
          message("Permutation iteration ", i, " failed: ", conditionMessage(e))
          return(NULL)
        }
      )
    },
    future.seed = use.seed
  )

  names(results) <- paste0("B_", seq_len(n_permutation))
  results
}


# #' @title 极速广义线性与置换打分引擎 (Fast GLM-Permutation Scorer)
# #'
# #' @description
# #' `permutation_scorer` 是计算混合物相对重要性（权重）与非线性剂量反应形状的核心引擎。
# #' 该函数通过 Bootstrap 划分训练集与袋外集 (OOB)，在训练集上估计样条参数，随后在 OOB 集上
# #' 通过置换特定暴露的样条基函数来量化模型损失的增加量，从而推导出暴露的相对权重。
# #'
# #' @details
# #' \strong{核心统计算法与性能优化:} \cr
# #' \enumerate{
# #'   \item \strong{自适应抽样设计:}
# #'     \itemize{
# #'       \item \emph{IID 抽样:} 对于常规广义线性模型 (`gaussian`, `binomial`, `poisson`)，执行标准的有放回随机抽样。
# #'       \item \emph{整群重抽样 (Cluster Resampling):} 对于配对/分层数据 (`clogit`)，算法严格执行基于 `strata_id` 的整群有放回抽样，
# #'       并在内部动态重构匹配组 ID (`new_strata`)。这防止了同一匹配组被破坏，避免了由于不平衡抽样引入的结构性选择偏倚。
# #'     }
# #'   \item \strong{底层损失函数极速求值 (Fast-Path Evaluation):}
# #'     \itemize{
# #'       \item 常规模型使用向量化的偏差 (Deviance) 或均方误差 (MSE) 作为损失函数。
# #'       \item \emph{Clogit 极简评估:} 对于条件逻辑回归，在 OOB 集上验证特征重要性时，没有使用高开销的 `clogit()` 重新拟合。
# #'       而是直接调用了底层的 \code{survival::coxph.fit}，强行传入训练集的系数 (\code{init = coefs_init})，并将最大迭代次数设为 0
# #'       (\code{control = coxph.control(iter.max = 0)})，配合 Efron 近似。这使得算法可以在毫秒级提取出偏似然 (Partial Log-likelihood)，
# #'       大幅提升了多维暴露置换时的性能。
# #'     }
# #'   \item \strong{权重归一化:} 特征重要性被定义为：将某暴露的所有非线性样条基向量同时随机打乱后，OOB 损失相对于基线 OOB 损失的平均增加量。
# #'     最终权重与重要性得分的平方根成正比，并归一化为和为 1 的相对比例。
# #' }
# #'
# #' @param x Matrix 或 data.frame。包含截距（可选）、混杂因素和所有暴露样条基函数的设计矩阵。
# #' @param y Numeric vector。结局变量向量。对于二分类或计数模型，必须转换为数值型。
# #' @param mix_name Character vector。混合物组分的名称列表。
# #' @param spline_vars Character vector。设计矩阵中所有属于样条基函数的列名。
# #' @param family \code{family} 对象或包含 \code{$family} 属性的列表。用于指定误差分布。
# #' @param n_permutation Integer。针对单个组分在 OOB 集上执行置换混洗的次数，默认为 100。
# #' @param strata_id Vector 或 \code{NULL}。若 \code{family$family == "clogit"}，此向量提供每行数据所属的匹配组 ID。
# #' @param ... 接收其他潜在的控制参数。
# #'
# #' @return 返回一个包含两个元素的列表：
# #' \itemize{
# #'   \item \code{weights}: 命名数值向量，表示各混合物组分的相对重要性权重（基于置换损失的平方根进行归一化）。
# #'   \item \code{shapes}: 命名数值向量，包含从训练集中提取的样条基函数回归系数（不含截距）。若模型因共线性导致变量被剔除，对应系数设为 0。
# #' }
# #'
# #' @importFrom stats coef predict glm.fit as.formula
# #' @importFrom survival Surv clogit coxph.fit coxph.control
# #' @export
# permutation_scorer <- function(x, y, mix_name, spline_vars, family, n_permutation = 100, strata_id = NULL, ...) {
#   n_obs <- nrow(x)
#   fam_name <- family$family

#   if (fam_name != "clogit") {
#     linkinv <- family$linkinv
#   }

#   # ──────────────────────────────────────────────────────────────────────
#   # 1. 内部 Bootstrap: 区分 IID 抽样与整群 (Cluster) 抽样
#   # ──────────────────────────────────────────────────────────────────────
#   if (fam_name == "clogit") {
#     if (is.null(strata_id)) stop("strata_id must be provided for family = 'clogit'")

#     unique_strata <- unique(strata_id)
#     n_strata <- length(unique_strata)

#     sampled_strata <- sample(unique_strata, size = n_strata, replace = TRUE)

#     idx_list <- lapply(seq_along(sampled_strata), function(i) {
#       orig_idx <- which(strata_id == sampled_strata[i])
#       data.frame(
#         orig_row = orig_idx,
#         new_strata = paste0(sampled_strata[i], "_boot_", i),
#         stringsAsFactors = FALSE
#       )
#     })
#     train_map <- do.call(rbind, idx_list)

#     idx <- train_map$orig_row
#     strata_train <- train_map$new_strata

#     oob_strata <- setdiff(unique_strata, sampled_strata)
#     if (length(oob_strata) == 0) {
#       return(NULL)
#     }

#     oob_idx <- which(strata_id %in% oob_strata)
#     strata_oob <- strata_id[oob_idx]
#   } else {
#     idx <- sample(seq_len(n_obs), size = n_obs, replace = TRUE)
#     oob_idx <- setdiff(seq_len(n_obs), idx)
#     if (length(oob_idx) == 0) {
#       return(NULL)
#     }
#   }

#   x_train <- x[idx, , drop = FALSE]
#   y_train <- y[idx]
#   x_oob <- x[oob_idx, , drop = FALSE]
#   y_oob <- y[oob_idx]

#   int_col <- match("(Intercept)", colnames(x_train))
#   if (!is.na(int_col)) {
#     x_train_net <- x_train[, -int_col, drop = FALSE]
#     x_oob_net <- x_oob[, -int_col, drop = FALSE]
#   } else {
#     x_train_net <- x_train
#     x_oob_net <- x_oob
#   }

#   # ──────────────────────────────────────────────────────────────────────
#   # 2. 损失函数定义 (极速底层实现版 - 预编译处理)
#   # ──────────────────────────────────────────────────────────────────────
#   if (fam_name == "clogit") {
#     y_surv_oob <- survival::Surv(rep(1, length(y_oob)), y_oob)
#     strata_int_oob <- as.integer(as.factor(strata_oob))

#     fast_clogit_loss <- function(x_new, coefs_init) {
#       fit_eval <- tryCatch(
#         {
#           survival::coxph.fit(
#             x = as.matrix(x_new),
#             y = y_surv_oob,
#             strata = strata_int_oob,
#             init = coefs_init,
#             control = survival::coxph.control(iter.max = 0),
#             method = "efron", # efron 计算极快，非常适合数据量大的项目
#             rownames = NULL
#           )
#         },
#         error = function(e) NULL
#       )

#       if (is.null(fit_eval)) {
#         return(NA_real_)
#       }
#       return(-2 * fit_eval$loglik[1])
#     }
#   } else {
#     calc_loss <- function(y_true, mu_pred) {
#       if (fam_name == "gaussian") {
#         return(mean((y_true - mu_pred)^2))
#       }
#       if (fam_name == "binomial") {
#         mu_pred <- pmax(pmin(mu_pred, 1 - 1e-7), 1e-7)
#         return(-2 * mean(y_true * log(mu_pred) + (1 - y_true) * log(1 - mu_pred)))
#       }
#       if (fam_name %in% c("poisson", "quasipoisson")) {
#         mu_pred <- pmax(mu_pred, 1e-7)
#         term1 <- ifelse(y_true == 0, 0, y_true * log(y_true / mu_pred))
#         return(2 * mean(term1 - (y_true - mu_pred)))
#       }
#       return(mean((y_true - mu_pred)^2))
#     }
#   }

#   # ──────────────────────────────────────────────────────────────────────
#   # 3 & 4. 模型拟合与基线 OOB 损失计算
#   # ──────────────────────────────────────────────────────────────────────
#   if (fam_name == "clogit") {
#     df_train <- data.frame(y_event = y_train, strata_id = strata_train)
#     df_train <- cbind(df_train, as.data.frame(x_train_net))

#     x_cols <- colnames(x_train_net)
#     form_str <- paste0("y_event ~ ", paste(sprintf("`%s`", x_cols), collapse = " + "), " + strata(strata_id)")

#     fit <- tryCatch(
#       {
#         survival::clogit(as.formula(form_str), data = df_train)
#       },
#       error = function(e) NULL
#     )

#     if (is.null(fit)) {
#       return(NULL)
#     }

#     coef_all <- coef(fit)
#     coef_all[is.na(coef_all)] <- 0
#     intercept_val <- 0
#     coefs_no_int <- coef_all

#     base_loss <- fast_clogit_loss(x_oob_net, coefs_no_int)
#     if (is.na(base_loss)) {
#       return(NULL)
#     }
#   } else {
#     x_train_glm <- cbind(Intercept = 1, as.matrix(x_train_net))
#     x_oob_glm <- cbind(Intercept = 1, as.matrix(x_oob_net))

#     fit <- stats::glm.fit(x = x_train_glm, y = y_train, family = family)
#     coef_all <- fit$coefficients
#     coef_all[is.na(coef_all)] <- 0

#     intercept_val <- unname(coef_all[1])
#     coefs_no_int <- coef_all[-1]

#     eta_oob <- as.numeric(x_oob_glm %*% coef_all)
#     mu_oob <- linkinv(eta_oob)
#     base_loss <- calc_loss(y_oob, mu_oob)
#   }

#   # ──────────────────────────────────────────────────────────────────────
#   # 5. 置换重洗评估特征重要性
#   # ──────────────────────────────────────────────────────────────────────
#   importance_scores <- numeric(length(mix_name))
#   names(importance_scores) <- mix_name

#   x_oob_shuffled <- x_oob_net
#   n_oob <- length(oob_idx)

#   for (var in mix_name) {
#     target_cols <- grep(paste0("^", var, "_B"), colnames(x_oob_net))

#     if (length(target_cols) == 0) {
#       warning(paste("No spline basis columns found for mixture component:", var))
#       return(NULL)
#     }

#     shuffled_loss_list <- numeric(n_permutation)

#     for (k in seq_len(n_permutation)) {
#       shuffle_idx <- sample(n_oob)
#       x_oob_shuffled[, target_cols] <- x_oob_net[shuffle_idx, target_cols, drop = FALSE]

#       if (fam_name == "clogit") {
#         loss_val <- fast_clogit_loss(x_oob_shuffled, coefs_no_int)
#         shuffled_loss_list[k] <- if (is.na(loss_val)) base_loss else loss_val
#       } else {
#         eta_shuffled <- intercept_val + as.numeric(as.matrix(x_oob_shuffled) %*% coefs_no_int)
#         mu_shuffled <- linkinv(eta_shuffled)
#         shuffled_loss_list[k] <- calc_loss(y_oob, mu_shuffled)
#       }
#     }

#     x_oob_shuffled[, target_cols] <- x_oob_net[, target_cols, drop = FALSE]
#     importance_scores[var] <- max(0, mean(shuffled_loss_list) - base_loss)
#   }

#   # ──────────────────────────────────────────────────────────────────────
#   # 6. 权重归一化
#   # ──────────────────────────────────────────────────────────────────────
#   if (sum(importance_scores) <= 0) {
#     weights <- rep(NA_real_, length(mix_name))
#     names(weights) <- mix_name
#     shape_coefs <- rep(NA_real_, length(spline_vars))
#     names(shape_coefs) <- spline_vars
#   } else {
#     weights <- sqrt(importance_scores) / sum(sqrt(importance_scores))
#     shape_coefs <- coefs_no_int[spline_vars]
#     shape_coefs[is.na(shape_coefs)] <- 0
#   }

#   return(list(weights = weights, shapes = shape_coefs))
# }
