#' Fast Ridge-Permutation Scorer for Spline WQS
#' 
#' 极速岭回归置换评分器：基于 OOB 网格寻优的权重提取引擎
#'
#' @description
#' This function serves as the core computational engine for a single iteration of NWQS. 
#' It avoids the heavy computational cost of internal k-fold cross-validation by using 
#' a predefined Lambda grid and selecting the optimal penalty based on Out-of-Bag (OOB) 
#' prediction error. Once the optimal shape is locked, it calculates variable importance 
#' through grouped permutation of spline bases.
#' \cr
#' 该函数是 NWQS 单次迭代的核心计算引擎。它通过预设 Lambda 网格并在袋外（OOB）数据上
#' 寻找最小损失来锁定最优惩罚参数，避开了昂贵的内部 K 折交叉验证。锁定最优形状后，
#' 通过对样条基组进行“成组置换”来提取组分的相对重要性。
#'
#' @details
#' \strong{Tuning Workflow (寻优流程):}
#' \enumerate{
#'   \item \strong{Path Fitting:} Fits a Ridge GLM path across 15 log-spaced Lambda values using In-Bag data.
#'   \item \strong{OOB Selection:} Evaluates all 15 models on OOB data to identify the \eqn{\lambda} 
#'         that minimizes the specified loss (MSE for Gaussian, Deviance for Binomial).
#'   \item \strong{Grouped Permutation:} For each mixture component, its associated spline basis 
#'         columns are shuffled simultaneously to break the outcome relationship while preserving 
#'         within-component spline structure.
#' }
#' 
#' \strong{Weight Derivation:}
#' \deqn{Weights_i = \frac{\sqrt{\Delta Loss_i}}{\sum \sqrt{\Delta Loss}}}
#' 
#' @param X_matrix numeric matrix. Design matrix containing spline bases and covariates.
#'   包含样条基函数和协变量的设计矩阵。
#' @param y_vector numeric vector. Outcome variable.
#'   因变量向量。
#' @param mix_name character vector. Names of original mixture components.
#'   原始混合物组分名称。
#' @param spline_vars character vector. Names of all spline basis columns in \code{X_matrix}.
#'   设计矩阵中所有样条基函数列的名称。
#' @param fam_obj list. A GLM family object (e.g., \code{gaussian()}, \code{binomial()}).
#'   包含 linkinv 和 family 名称的列表。
#' @param shuffle integer. Number of permutations to stabilize importance scores. Default is 100.
#'   置换洗牌次数。默认为 100。
#' @param ... Additional arguments (currently unused, for compatibility).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{weights}: Normalized importance-based weights for each mixture component.
#'   \item \code{shapes}: Estimated spline coefficients (\eqn{\theta}) at the optimal \eqn{\lambda}.
#' }
#' 
#' @importFrom glmnet glmnet
#' @importFrom stats coef predict
#' @export
ridge_permutation_scorer <- function(X_matrix, y_vector, mix_name, spline_vars, fam_obj, shuffle = 100, ...) {
    
    n_obs <- nrow(X_matrix)
    # 1. Bootstrap 采样划分 In-Bag (训练) 与 OOB (验证/测试)
    idx <- sample(seq_len(n_obs), size = n_obs, replace = TRUE)
    oob_idx <- setdiff(seq_len(n_obs), idx)
    
    if (length(oob_idx) == 0) return(NULL)

    X_train <- X_matrix[idx, , drop = FALSE]
    y_train <- y_vector[idx]
    X_oob <- X_matrix[oob_idx, , drop = FALSE]
    y_oob <- y_vector[oob_idx]

    # 2. 极速 Lambda 网格搜索 (15个对数间隔点)
    lambda_grid <- exp(seq(log(1e-5), log(1), length.out = 15))

    # 一次性拟合整条 Ridge 路径 (Alpha = 0)
    fit <- suppressWarnings(glmnet::glmnet(
        x = X_train, y = y_train, family = fam_obj$family,
        alpha = 0, lambda = lambda_grid, standardize = FALSE, intercept = TRUE
    ))

    # 3. 计算 OOB 损失函数，寻找最优 Lambda
    pred_oob_all <- predict(fit, newx = X_oob, type = "response")

    calc_oob_loss <- function(y, mu_mat) {
        if (fam_obj$family == "gaussian") {
            return(colMeans((y - mu_mat)^2))
        }
        if (fam_obj$family == "binomial") {
            # 防止对数 log(0) 溢出
            mu_mat <- pmax(pmin(mu_mat, 1 - 1e-7), 1e-7)
            return(-colMeans(y * log(mu_mat) + (1 - y) * log(1 - mu_mat)))
        }
        # 兜底返回 MSE
        return(colMeans((y - mu_mat)^2))
    }

    oob_losses <- calc_oob_loss(y_oob, pred_oob_all)
    best_idx <- which.min(oob_losses)
    best_coefs <- as.matrix(coef(fit, s = lambda_grid[best_idx]))

    intercept_val <- best_coefs[1, 1]
    coefs_no_int <- best_coefs[-1, 1]
    base_loss <- oob_losses[best_idx]

    # 4. 成组置换重要性提取 (Grouped Permutation)
    importance_scores <- numeric(length(mix_name))
    names(importance_scores) <- mix_name
    X_oob_shuffled <- X_oob

    # 预计算 OOB 的线性预测部分，提速洗牌过程
    for (var in mix_name) {
        # 寻找属于当前组分的所有样条基函数列
        target_cols <- grep(paste0("^", var, "_B"), colnames(X_oob))
        shuffled_loss_list <- numeric(shuffle)
        
        for (k in seq_len(shuffle)) {
            # 同步洗牌该组分的所有基函数
            X_oob_shuffled[, target_cols] <- X_oob[sample(length(oob_idx)), target_cols]
            
            # 计算洗牌后的 Loss
            eta_shuffled <- intercept_val + as.numeric(X_oob_shuffled %*% coefs_no_int)
            mu_shuffled <- fam_obj$linkinv(eta_shuffled)
            shuffled_loss_list[k] <- mean(calc_oob_loss(y_oob, matrix(mu_shuffled, ncol = 1)))
        }
        # 还原 OOB 矩阵用于下一个变量
        X_oob_shuffled[, target_cols] <- X_oob[, target_cols] 
        
        # 计算 Loss 的增加量 (Importance)
        importance_scores[var] <- max(0, mean(shuffled_loss_list) - base_loss)
    }

    # 5. 返回归一化权重与形状系数
    return(list(
        weights = sqrt(importance_scores) / sum(sqrt(importance_scores)),
        shapes = coefs_no_int[spline_vars]
    ))
}