#' Calculate Spline-WQS Weights (Single Iteration)
#' 计算 Spline-WQS 权重 (单次迭代)
#'
#' @description
#' Performs a single iteration of the WQS regression process using spline expansion.
#' The workflow includes:
#' 1. Bootstrapping the dataset (separating In-Bag and Out-of-Bag samples).
#' 2. Non-linearly expanding mixture variables using splines.
#' 3. Fitting a Generalized Linear Model (GLM) on the training set.
#' 4. Calculating variable importance via **Grouped Permutation** on the OOB set.
#' \cr
#' 执行单次 Spline-WQS 迭代。工作流包括：
#' 1. 数据集重抽样（分为袋内和袋外数据）。
#' 2. 使用样条函数对混合变量进行非线性展开。
#' 3. 在训练集上拟合广义线性模型 (GLM)。
#' 4. 在 OOB 数据集上通过 **成组置换 (Grouped Permutation)** 计算变量重要性。
#'
#' @details
#' **Grouped Permutation Importance (成组置换重要性):**
#' Unlike standard permutation importance where columns are shuffled independently, this function
#' respects the structure of spline basis functions. When calculating the importance of a mixture component
#' (e.g., "Chemical_A"), all its corresponding spline basis columns (e.g., "Chemical_A_basis1", "Chemical_A_basis2")
#' are shuffled **synchronously** using the same random index. This preserves the internal correlation structure
#' of the spline expansion while breaking the relationship with the outcome.
#' \cr
#' 与标准置换重要性不同，本函数尊重样条基函数的结构。在计算某个混合组分（如 "Chemical_A"）的重要性时，
#' 它对应的所有样条基函数列（如 "Chemical_A_basis1", "Chemical_A_basis2"）会使用相同的随机索引**同步洗牌**。
#' 这样既打破了与因变量的联系，又保留了样条展开内部的相关结构。
#'
#' **Weight Calculation:**
#' \deqn{w_i = \frac{\sqrt{Imp_i}}{\sum \sqrt{Imp}}}
#' Where \eqn{Imp_i} is the increase in MSE when variable \eqn{i} is shuffled.
#'
#' @param data data.frame. The full dataset containing mixture and covariates.
#'   包含混合物和协变量的完整数据框。
#' @param mix_name character vector. Names of the mixture components (original variable names before expansion).
#'   混合物组分的变量名向量（展开前的原始变量名）。
#' @param dependent_var character. Name of the dependent variable (outcome). Defaults to "y".
#'   因变量名称。默认为 "y"。
#' @param expand_func function. Function used to non-linearly expand the mixture variables.
#'   Defaults to `wqs_nonlinear_expand`.
#'   用于非线性展开混合变量的函数。默认为 `wqs_nonlinear_expand`。
#' @param shuffle integer. Number of permutations for calculating variable importance.
#'   Higher values provide more stable importance scores but increase computation time. Default is 100.
#'   计算变量重要性时的置换次数。数值越高评分越稳定，但计算时间增加。默认为 100。
#' @param ... Additional arguments passed to `expand_func`.
#'   Critical arguments include:
#'   \itemize{
#'     \item `df_spline`: Degrees of freedom for the spline expansion (default 3).
#'     \item `transform_fun`: Quantile transformation function.
#'   }
#'   传递给 `expand_func` 的额外参数。
#'
#' @return numeric vector. A named vector of normalized weights for each mixture component.
#'   Returns `NA` if all importance scores are zero or negative.
#'   Returns `NULL` if OOB sample size is zero (edge case).
#'   返回每个混合物组分的归一化权重向量。如果重要性均为负或0返回 NA。
#'
#' @importFrom stats glm predict gaussian as.formula
#' @export
calc_spline_wqs_weights = function(data, mix_name, dependent_var = "y",
                                   expand_func = wqs_nonlinear_expand,
                                   shuffle = 100, ...) {

    args = list(...)

    n_obs = nrow(data)
    # 产生 Bootstrap 索引 (有放回)
    idx = sample(seq_len(n_obs), size = n_obs, replace = TRUE)
    oob_idx = setdiff(seq_len(n_obs), idx)

    if (length(oob_idx) == 0) return(NULL)

    # 准备训练和 OOB 数据
    train_raw = data[idx, , drop = FALSE]
    oob_raw = data[oob_idx, , drop = FALSE]

    # 非线性转换 (在循环外只做一次，保证效率)
    # 这里的 ... 将 df_spline 等参数传递给 wqs_nonlinear_expand
    q = args$q
    df = args$df_spline
    train_data_spline = expand_func(train_raw, mix_name, df_spline = df, q = q)
    oob_data_spline = expand_func(oob_raw, mix_name, df_spline = df, q = q)

    # # 仅测试
    # train_data_spline = expand_func(train_raw, mix_name)
    # oob_data_spline = expand_func(oob_raw, mix_name)

    # 合并最终数据集
    base_cols = !(names(data) %in% mix_name)
    train_final = cbind(train_raw[, base_cols, drop = FALSE], train_data_spline)
    oob_final = cbind(oob_raw[, base_cols, drop = FALSE], oob_data_spline)

    # 动态构造公式
    spline_vars = colnames(train_data_spline)
    covariates = setdiff(names(train_raw)[base_cols], dependent_var)
    formula_str = paste(dependent_var, "~", paste(c(spline_vars, covariates), collapse = " + "))
    internal_formula = as.formula(formula_str)

    # 拟合模型 (目前固定为 gaussian，后续可考虑通过 ... 或新参数传入 family)
    fit = glm(formula = internal_formula, data = train_final, family = gaussian())

    # 计算基础 MSE
    y_true = oob_raw[[dependent_var]]
    base_pred = predict(fit, newdata = oob_final)
    base_mse = mean((y_true - base_pred)^2)

    # 核心优化：成组洗牌计算重要性
    importance_scores = numeric(length(mix_name))
    names(importance_scores) = mix_name

    for (var in mix_name) {
        # 识别该变量对应的所有样条基函数列
        # 使用正则匹配保证准确性 (如 Component1_B1, Component1_B2...)
        target_cols = spline_vars[grep(paste0("^", var, "_B"), spline_vars)]

        shuffled_mse_list = numeric(shuffle)

        for (k in seq_len(shuffle)) {
            temp_oob_final = oob_final

            # 同步洗牌：同一组分的所有基函数共用同一个随机索引，保持组分内部结构
            shuffle_idx = sample(nrow(temp_oob_final))
            temp_oob_final[, target_cols] = temp_oob_final[shuffle_idx, target_cols]

            # 直接预测，不再重新 expand，速度极快
            shuffled_pred = predict(fit, newdata = temp_oob_final)
            shuffled_mse_list[k] = mean((y_true - shuffled_pred)^2)
        }

        # 计算该变量的平均 MSE 增量 (Permutation Importance)
        importance_scores[var] = max(0, mean(shuffled_mse_list) - base_mse)
    }

    # 计算归一化权重
    out = importance_scores
    if (sum(out) <= 0) {
        warning("所有变量的重要性得分为 0 或负数，返回 NA")
        weights = NA
    } else {
        weights = sqrt(out) / sum(sqrt(out))
    }

    return(weights)
}
