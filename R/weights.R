#' @title 极速广义线性与条件逻辑回归的置换打分引擎
#'
#' @description
#' `permutation_scorer` 是 NWQS 架构中负责单次内部迭代的核心计算引擎。该函数在袋内 (In-Bag)
#' Bootstrap 样本上拟合无惩罚的广义线性模型 (GLM) 或条件逻辑回归 (clogit) 模型，并在袋外 (OOB)
#' 样本上评估预测损失，最终通过随机置换推导出各个混合物组分的相对重要性。
#'
#' @details
#' \strong{工作流与算法提速:}
#' \enumerate{
#'   \item \strong{自适应抽样与拟合:} 在袋内数据上拟合标准 GLM 或 clogit 以初步估计样条基系数。
#'   \item \strong{极速基线损失评估:} 在 OOB 样本上计算基线损失（如均方误差或偏对数似然）。对于 `clogit`，
#'         为突破计算瓶颈，底层直接剥离了 formula 解析，调用 C++ 层的 \code{survival::coxph.fit} 
#'         并设定 \code{iter.max = 0}，实现了毫秒级的精准偏对数似然 (Partial Log-likelihood) 评估。
#'   \item \strong{分组置换重要性 (Grouped Permutation Importance):} 针对每个混合物组分，联合打乱其所有对应的
#'         非线性样条基函数列，并重新计算 OOB 损失的变化量。
#' }
#'
#' \strong{医学统计学严谨性与偏倚防范:}
#' \itemize{
#'   \item \strong{匹配设计的整群保护 (Cluster Preservation):} 当 \code{family = "clogit"} 时，算法严格执行
#'         基于匹配层 (\code{strata_id}) 的整群有放回抽样 (Cluster Resampling)，并在内部动态重命名匹配组 ID。
#'         这一机制完美保护了 1:N 配对病例对照研究 (Matched Case-Control) 的设计完整性。如果在此时错误地
#'         采用 IID 独立抽样而打乱了配对结构，将不可避免地引入严重的选择偏倚 (Selection Bias) 并歪曲推断。
#'   \item \strong{降低残余混杂 (Residual Confounding):} 通过置换整个样条基矩阵而不是单一线性系数，引擎能够
#'         全面捕捉高度非线性的剂量反应关系（如阈值效应或 U 型暴露风险），这有效降低了因参数化形式误设 
#'         (Model Misspecification) 而残留的混杂风险。
#' }
#'
#' \strong{权重推导公式:}
#' 组分的相对权重由其置换后造成的损失增加量 (Delta Loss) 决定，为避免负重要性，截断了小于 0 的噪声波动：
#' $$Weights_i = \frac{\sqrt{\max(0, \Delta Loss_i)}}{\sum \sqrt{\max(0, \Delta Loss)}}$$
#'
#' @param x Numeric matrix。包含样条基函数列和调整协变量的设计矩阵。
#' @param y Numeric vector。结局变量。
#' @param mix_name Character vector。原始混合物组分的名称列表。
#' @param spline_vars Character vector。设计矩阵 \code{x} 中所有属于样条基函数的列名。
#' @param family List。GLM 的误差分布族对象（或包含 \code{$family = "clogit"} 的伪族列表）。
#' @param n_permutation Integer。用于稳定重要性得分的 OOB 样本置换次数，默认为 100。
#' @param strata_id Vector 或 \code{NULL}。条件逻辑回归模型所需的匹配层/组别 ID。
#' @param ... 其他兼容性扩展参数。
#'
#' @return 返回一个包含以下两个元素的列表：
#' \itemize{
#'   \item \code{weights}: 基于置换损失增加量归一化后的各混合物组分相对权重。
#'   \item \code{shapes}: 在袋内 (In-bag) 样本上估计得到的所有混合物样条系数。若模型因完全分离或共线性剔除变量，对应系数设为 0。
#' }
#'
#' @importFrom stats coef predict glm.fit as.formula
#' @importFrom survival Surv clogit coxph.fit coxph.control
#' @export
permutation_scorer <- function(x, y, mix_name, spline_vars, family, n_permutation = 100, strata_id = NULL, ...) {
  n_obs <- nrow(x)
  fam_name <- family$family
  
  if (fam_name != "clogit") {
    linkinv <- family$linkinv
  }

  # ──────────────────────────────────────────────────────────────────────
  # 1. 内部 Bootstrap: 必须区分 IID 抽样与整群 (Cluster) 抽样
  # ──────────────────────────────────────────────────────────────────────
  if (fam_name == "clogit") {
    if (is.null(strata_id)) stop("strata_id must be provided for family = 'clogit'")
    
    unique_strata <- unique(strata_id)
    n_strata <- length(unique_strata)
    
    # 随机抽取匹配组
    sampled_strata <- sample(unique_strata, size = n_strata, replace = TRUE)
    
    # 构建训练集映射字典，防止组别重名
    idx_list <- lapply(seq_along(sampled_strata), function(i) {
      orig_idx <- which(strata_id == sampled_strata[i])
      data.frame(
        orig_row = orig_idx,
        new_strata = paste0(sampled_strata[i], "_boot_", i),
        stringsAsFactors = FALSE
      )
    })
    train_map <- do.call(rbind, idx_list)
    
    idx <- train_map$orig_row
    strata_train <- train_map$new_strata
    
    # OOB 样本: 未被抽中的组
    oob_strata <- setdiff(unique_strata, sampled_strata)
    if (length(oob_strata) == 0) return(NULL)
    
    oob_idx <- which(strata_id %in% oob_strata)
    strata_oob <- strata_id[oob_idx]
    
  } else {
    idx <- sample(seq_len(n_obs), size = n_obs, replace = TRUE)
    oob_idx <- setdiff(seq_len(n_obs), idx)
    if (length(oob_idx) == 0) return(NULL)
  }

  x_train <- x[idx, , drop = FALSE]
  y_train <- y[idx]
  x_oob <- x[oob_idx, , drop = FALSE]
  y_oob <- y[oob_idx]

  # 移除拦截项以配合后续逻辑
  int_col <- match("(Intercept)", colnames(x_train))
  if (!is.na(int_col)) {
    x_train_net <- x_train[, -int_col, drop = FALSE]
    x_oob_net <- x_oob[, -int_col, drop = FALSE]
  } else {
    x_train_net <- x_train
    x_oob_net <- x_oob
  }

  # ──────────────────────────────────────────────────────────────────────
  # 2. 损失函数定义 (针对大数据量进行的极速底层实现)
  # ──────────────────────────────────────────────────────────────────────
  if (fam_name == "clogit") {
    # 绕开缓慢的 formula 解析，直接调用 C++ 引擎求解偏对数似然
    fast_clogit_loss <- function(x_new, y_new, strata_new, coefs_init) {
      y_surv <- survival::Surv(rep(1, length(y_new)), y_new)
      strata_int <- as.integer(as.factor(strata_new))
      
      fit_eval <- tryCatch({
        survival::coxph.fit(
          x = as.matrix(x_new),
          y = y_surv,
          strata = strata_int,
          init = coefs_init,
          control = survival::coxph.control(iter.max = 0), # 迭代0次，仅算当前特征下的Deviance
          method = "exact",
          rownames = NULL
        )
      }, error = function(e) NULL)
      
      if (is.null(fit_eval)) return(NA_real_)
      return(-2 * fit_eval$loglik[1])
    }
  } else {
    calc_loss <- function(y_true, mu_pred) {
      if (fam_name == "gaussian") return(mean((y_true - mu_pred)^2))
      if (fam_name == "binomial") {
        mu_pred <- pmax(pmin(mu_pred, 1 - 1e-7), 1e-7)
        return(-2 * mean(y_true * log(mu_pred) + (1 - y_true) * log(1 - mu_pred)))
      }
      if (fam_name %in% c("poisson", "quasipoisson")) {
        mu_pred <- pmax(mu_pred, 1e-7)
        term1 <- ifelse(y_true == 0, 0, y_true * log(y_true / mu_pred))
        return(2 * mean(term1 - (y_true - mu_pred)))
      }
      return(mean((y_true - mu_pred)^2))
    }
  }

  # ──────────────────────────────────────────────────────────────────────
  # 3 & 4. 模型拟合与基线 OOB 损失计算
  # ──────────────────────────────────────────────────────────────────────
  if (fam_name == "clogit") {
    df_train <- data.frame(y_event = y_train, strata_id = strata_train)
    df_train <- cbind(df_train, as.data.frame(x_train_net))
    
    x_cols <- colnames(x_train_net)
    form_str <- paste0("y_event ~ ", paste(sprintf("`%s`", x_cols), collapse = " + "), " + strata(strata_id)")
    
    fit <- tryCatch({
      survival::clogit(as.formula(form_str), data = df_train)
    }, error = function(e) NULL)
    
    if (is.null(fit)) return(NULL)
    
    coef_all <- coef(fit)
    coef_all[is.na(coef_all)] <- 0
    intercept_val <- 0
    coefs_no_int <- coef_all
    
    # 算基线 Deviance 
    base_loss <- fast_clogit_loss(x_oob_net, y_oob, strata_oob, coefs_no_int)
    if (is.na(base_loss)) return(NULL)
    
  } else {
    x_train_glm <- cbind(Intercept = 1, as.matrix(x_train_net))
    x_oob_glm <- cbind(Intercept = 1, as.matrix(x_oob_net))

    fit <- stats::glm.fit(x = x_train_glm, y = y_train, family = family)
    coef_all <- fit$coefficients
    coef_all[is.na(coef_all)] <- 0

    intercept_val <- unname(coef_all[1])
    coefs_no_int <- coef_all[-1]

    eta_oob <- as.numeric(x_oob_glm %*% coef_all)
    mu_oob <- linkinv(eta_oob)
    base_loss <- calc_loss(y_oob, mu_oob)
  }

  # ──────────────────────────────────────────────────────────────────────
  # 5. 置换重洗评估特征重要性 (Permutation Importance)
  # ──────────────────────────────────────────────────────────────────────
  importance_scores <- numeric(length(mix_name))
  names(importance_scores) <- mix_name

  x_oob_shuffled <- x_oob_net
  n_oob <- length(oob_idx)

  for (var in mix_name) {
    target_cols <- grep(paste0("^", var, "_B"), colnames(x_oob_net))

    if (length(target_cols) == 0) {
      warning(paste("No spline basis columns found for mixture component:", var))
      return(NULL)
    }

    shuffled_loss_list <- numeric(n_permutation)

    for (k in seq_len(n_permutation)) {
      shuffle_idx <- sample(n_oob)
      x_oob_shuffled[, target_cols] <- x_oob_net[shuffle_idx, target_cols, drop = FALSE]

      if (fam_name == "clogit") {
        # 计算打乱后的 OOB Deviance
        loss_val <- fast_clogit_loss(x_oob_shuffled, y_oob, strata_oob, coefs_no_int)
        shuffled_loss_list[k] <- if (is.na(loss_val)) base_loss else loss_val
      } else {
        eta_shuffled <- intercept_val + as.numeric(as.matrix(x_oob_shuffled) %*% coefs_no_int)
        mu_shuffled <- linkinv(eta_shuffled)
        shuffled_loss_list[k] <- calc_loss(y_oob, mu_shuffled)
      }
    }

    x_oob_shuffled[, target_cols] <- x_oob_net[, target_cols, drop = FALSE]
    importance_scores[var] <- max(0, mean(shuffled_loss_list) - base_loss)
  }

  # ──────────────────────────────────────────────────────────────────────
  # 6. 权重归一化与形变特征提取
  # ──────────────────────────────────────────────────────────────────────
  if (sum(importance_scores) <= 0) {
    weights <- rep(NA_real_, length(mix_name))
    names(weights) <- mix_name
    shape_coefs <- rep(NA_real_, length(spline_vars))
    names(shape_coefs) <- spline_vars
  } else {
    weights <- sqrt(importance_scores) / sum(sqrt(importance_scores))
    shape_coefs <- coefs_no_int[spline_vars]
    shape_coefs[is.na(shape_coefs)] <- 0
  }

  return(list(weights = weights, shapes = shape_coefs))
}
