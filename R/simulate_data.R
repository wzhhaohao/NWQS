#' Add Gaussian noise given target SNR
#'
#' @description
#' Calculates signal power and adds white noise based on target SNR (dB).
#' 计算信号功率并根据目标信噪比（dB）添加白噪声。
#'
#' @details
#' Signal-to-Noise Ratio (SNR) formula: SNR = P_signal / P_noise
#' SNR_dB = 10 * log10(SNR)
#' Higher SNR indicates better signal quality relative to noise.
#' 信噪比越高，表示信号相对于噪声越强，信号质量越好。
#'
#' @param signal_vec numeric vector. The clean signal vector. / 原始纯净信号。
#' @param snr_db numeric. Target SNR in decibels (e.g., 0, 5, 10, 20). / 目标信噪比 (dB)。
#' @return numeric vector. The noisy signal. / 添加噪音后的信号。
#' @export
add_noise_by_snr = function(signal_vec, snr_db) {
    stopifnot(is.numeric(signal_vec), length(signal_vec) > 1)

    # Calculate signal power: P_signal = E[(x - mu)^2]
    power_signal = mean((signal_vec - mean(signal_vec))^2)

    if (power_signal == 0) {
        return(signal_vec)
    }

    # Calculate noise standard deviation based on SNR formula
    # SNR_dB = 10 * log10(P_signal / P_noise)
    sigma_noise = sqrt(power_signal / 10^(snr_db / 10))

    noise_vec = rnorm(length(signal_vec), mean = 0, sd = sigma_noise)

    return(signal_vec + noise_vec)
}

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

#' Generate Covariance/Correlation Matrix / 生成不同类型的协方差/相关矩阵
#'
#' @description
#' Generates a positive-definite symmetric covariance matrix based on the specified mode.
#' 根据指定模式生成正定对称的协方差矩阵。
#'
#' @param n_vars integer. Number of variables. / 变量数量。
#' @param mode character. Correlation mode: "low", "medium", "high", "mixed". / 相关性模式。
#' @param rho numeric. Base correlation strength (optional). / 基础相关系数强度。
#' @param seed integer. Random seed for reproducibility. / 随机种子。
#' @return matrix. Correlation matrix. / 相关系数矩阵。
#' @export
generate_sigma = function(n_vars, mode = c("medium", "low", "high", "mixed"), rho = 0.7, seed = NULL) {
    mode = match.arg(mode)
    if (!is.null(seed)) set.seed(seed)

    if (mode == "low") {
        # Identity matrix + small noise
        A = diag(n_vars) + matrix(runif(n_vars^2, -0.1, 0.1), nrow = n_vars)
        sigma = cov2cor(t(A) %*% A)
    } else if (mode == "medium") {
        # Random Gram matrix
        A = matrix(runif(n_vars^2, -1, 1), ncol = n_vars)
        sigma = cov2cor(t(A) %*% A)
    } else if (mode == "high") {
        # Factor model structure
        sigma = matrix(rho, nrow = n_vars, ncol = n_vars)
        diag(sigma) = 1

        # Add jitter
        noise = matrix(runif(n_vars^2, -0.05, 0.05), nrow = n_vars)
        sigma = sigma + (noise + t(noise)) / 2

        # Repair eigenvalues to ensure positive definiteness
        eig = eigen(sigma)
        val = pmax(eig$values, 0.01)
        sigma = cov2cor(eig$vectors %*% diag(val) %*% t(eig$vectors))
    } else if (mode == "mixed") {
        # Block diagonal structure
        split_idx = floor(n_vars / 2)
        s1 = split_idx
        s2 = n_vars - split_idx

        # Block 1: High correlation
        B1 = matrix(0.8, nrow = s1, ncol = s1)
        diag(B1) = 1

        # Block 2: Medium correlation (Random)
        A2 = matrix(runif(s2^2, -1, 1), ncol = s2)
        B2 = cov2cor(t(A2) %*% A2)

        # Combine blocks
        sigma = matrix(0, nrow = n_vars, ncol = n_vars)
        sigma[1:s1, 1:s1] = B1
        sigma[(s1 + 1):n_vars, (s1 + 1):n_vars] = B2

        # Add noise and repair eigenvalues
        noise = matrix(runif(n_vars^2, -0.1, 0.1), nrow = n_vars)
        sigma = sigma + (noise + t(noise)) / 2

        eig = eigen(sigma)
        val = pmax(eig$values, 0.01)
        sigma = cov2cor(eig$vectors %*% diag(val) %*% t(eig$vectors))
    }

    return(sigma)
}

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

#' Generate Covariates
#'
#' @description
#' Generates data with continuous, binary, and categorical variables and their linear effects.
#' 生成包含连续、二值和分类变量的数据及其线性效应。
#'
#' @param n_obs integer. Sample size. / 样本量。
#' @param beta_cont numeric. Coefficient for continuous variable. / 连续变量系数。
#' @param beta_bin numeric. Coefficient for binary variable. / 二值变量系数。
#' @param beta_cat numeric vector. Coefficients for categorical variable. / 分类变量系数。
#' @param prob_bin numeric. Probability for binary variable. / 二值变量概率。
#' @param prob_cat numeric vector. Probabilities for categorical levels. / 分类变量各水平概率。
#' @return list. Model matrix, raw data, and linear effects. / 包含模型矩阵、原始数据和线性效应。
#' @export
generate_covariates = function(n_obs = 1000,
                                beta_cont = 0.5,
                                beta_bin = -0.8,
                                beta_cat = c(0, -0.5, 0.7),
                                prob_bin = 0.5,
                                prob_cat = c(1 / 3, 1 / 3, 1 / 3),
                                Intercept = 0) {
    x_cont = rnorm(n_obs, 0, 1)
    x_bin_raw = rbinom(n_obs, 1, prob_bin)
    x_cat_raw = sample(1:3, n_obs, replace = TRUE, prob = prob_cat)

    x_bin = factor(x_bin_raw, levels = c(0, 1))
    x_cat = factor(x_cat_raw, levels = 1:3)

    # Calculate linear predictor: eta = X * beta
    eta_cov = beta_cont * x_cont + beta_bin * x_bin_raw + beta_cat[x_cat_raw] + Intercept

    df_raw = data.frame(x_cont, x_bin, x_cat)
    # mat_model = model.matrix(~ x_cont + x_bin + x_cat, data = df_raw)[, -1, drop = FALSE]

    df_result = as.data.frame(cbind(eta_cov = eta_cov, df_raw))

    list(mm = df_result, original = df_raw, eta_cov = eta_cov)
}

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

#' Generate Linear Model Data
#'
#' @description
#' Generates main predictors (multivariate normal) and covariates, adding noise based on SNR.
#' 生成主预测变量（多元正态）和协变量，并基于SNR添加噪音。
#'
#' @param n_obs integer. Sample size. / 样本量。
#' @param mu_preds numeric vector. Mean vector for predictors. / 预测变量均值。
#' @param sigma_preds matrix. Covariance matrix for predictors. / 预测变量协方差矩阵。
#' @param beta_wqs numeric. Weight for predictor coefficients. / 预测变量系数的权重，默认为1，不做缩放。
#' @param beta_preds numeric vector. Coefficients for predictors. / 预测变量回归系数。
#' @param snr_db numeric. Signal-to-Noise Ratio (dB). / 信噪比。
#' @param transform_fun function. Optional transformation function. / 可选的数据变换函数。
#' @param ... arguments passed to generate_covariates. 传递给 generate_covariates 的参数。
#' @return data.frame. The final synthetic dataset. / 最终合成数据。
#' @export
generate_linear_data = function(n_obs = 1000,
                                 mu_preds,
                                 sigma_preds,
                                 beta_wqs = 1,
                                 beta_preds,
                                 snr_db = 10,
                                 transform_fun = NULL,
                                 seed = NULL,
                                 ...) {
    if (!is.null(seed)) set.seed(seed)
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' required")

    # Generate Components
    preds_raw = MASS::mvrnorm(n_obs, mu = mu_preds, Sigma = sigma_preds)
    preds_scaled = as.data.frame(scale(preds_raw))

    # Quantile/Percentile Transformation
    if (!is.null(transform_fun) && is.function(transform_fun)) {
        preds_final = transform_fun(preds_scaled)
    } else {
        preds_final = preds_scaled
    }

    names(preds_scaled) = paste0("Component", 1:ncol(preds_scaled))

    # Generate Background Covariates (passing '...')
    cov_list = generate_covariates(n_obs = n_obs, ...)

    beta_preds = beta_wqs * beta_preds

    # Calculate Clean Signal (Y = X * beta + cov_eta)
    y_clean = as.matrix(preds_final) %*% beta_preds + cov_list$eta_cov

    # Add SNR noise
    y_observed = add_noise_by_snr(as.vector(y_clean), snr_db = snr_db)

    # Combine Results
    cols_cov = setdiff(names(cov_list$mm), "eta_cov")
    final_df = cbind(
        y = y_observed,
        preds_scaled,
        cov_list$mm[, cols_cov, drop = FALSE]
    )

    as.data.frame(final_df)
}

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

#' Generate Non-linear Model Data with Splines
#' 生成基于自然样条(Natural Splines)的非线性模型数据
#'
#' @description
#' Generates predictors (multivariate normal), transforms them (optional),
#' expands them using natural splines, adds covariates, and generates Y with SNR-based noise.
#' 生成预测变量，进行变换（可选），使用自然样条扩展，加入协变量，并根据信噪比生成 Y。
#'
#' @param n_obs integer. Sample size. / 样本量。
#' @param mu_preds numeric vector. Mean vector for predictors. / 预测变量均值。
#' @param sigma_preds matrix. Covariance matrix for predictors. / 预测变量协方差矩阵。
#' @param beta_wqs numeric. Scaling factor for coefficients (default 1). / 系数缩放因子。
#' @param beta_preds numeric vector. Coefficients for the predictors.
#'   Length must equal n_vars. / 预测变量的系数向量，长度必须等于变量数。
#'   Note: This coefficient is broadcasted to all spline basis functions of the variable.
#' @param snr_db numeric. Signal-to-Noise Ratio (dB). / 信噪比。
#' @param transform_fun function. Function to transform predictors (e.g., quantile binning). / 变换函数。 默认分位数转化
#' @param df_spline integer. Degrees of freedom for natural splines (default 3). / 自然样条自由度。
#' @param seed integer. Random seed. / 随机种子。
#' @param ... arguments passed to generate_covariates. 传递给 generate_covariates 的参数。
#'
#' @return data.frame. Contains Y, transformed predictors (before spline expansion), and covariates.
#' @export
gen_nonlinear_data = function(n_obs = 1000,
                               mu_preds,
                               sigma_preds,
                               beta_wqs = 1,
                               beta_preds,
                               snr_db = 10,
                               transform_fun = NULL,
                               q = 4, # 新增: 告知函数分位数的层数，用于求真值
                               df_spline = 3,
                               seed = NULL,
                               shape = "linear_like",
                               ...) {
    if (!requireNamespace("splines", quietly = TRUE)) stop("Package 'splines' required")
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' required")
    if (!is.null(seed)) set.seed(seed)

    preds_raw = MASS::mvrnorm(n_obs, mu = mu_preds, Sigma = sigma_preds)
    preds_scaled = as.data.frame(scale(preds_raw))
    n_vars = ncol(preds_scaled)
    names(preds_scaled) = paste0("Component", 1:n_vars)

    if (!is.null(transform_fun) && is.function(transform_fun)) {
        preds_trans = transform_fun(preds_scaled)
    } else {
        preds_trans = preds_scaled
    }

    mat_spline_list = lapply(preds_trans, function(x) splines::ns(x, df = df_spline))

    if (length(beta_preds) != n_vars) stop("Length of 'beta_preds' must match n_vars.")

    cov_list = generate_covariates(n_obs = n_obs, ...)

    # -----------------------------------------------------------
    # [架构升级]: 完美混合 纯线性(pure_linear) 与 非线性样条
    # -----------------------------------------------------------

    # 动态推导并锁定节点 (保证数据生成和模型拟合的尺子绝对一致)
    eval_pts = 0:(q - 1)

    # 利用 splines 内部算法自动计算最合理的 knots (支持任意的 q)
    temp_spline = splines::ns(eval_pts, df = df_spline)
    global_knots = attr(temp_spline, "knots")
    global_boundary = attr(temp_spline, "Boundary.knots")

    # 打印出来看看，如果是 q=4，它会自动算出 1 和 2 (极度聪明)
    # print(global_knots)

    # 数据生成侧：强制使用这把尺子
    mat_spline_list = lapply(preds_trans, function(x) {
        splines::ns(x, df = df_spline, knots = global_knots, Boundary.knots = global_boundary)
    })

    # 真值计算侧：也强制使用这把尺子 (无截距版)
    basis_std_true = splines::ns(eval_pts, df = df_spline, knots = global_knots, Boundary.knots = global_boundary, intercept = FALSE)

    if (length(shape) == 1) shape = rep(shape, n_vars)

    eta_components_raw = matrix(0, nrow = n_obs, ncol = n_vars)
    baseline_components = numeric(n_vars)

    true_eff_mat = matrix(0, nrow = n_vars + 1, ncol = q - 1)
    rownames(true_eff_mat) = c("Overall", names(preds_scaled))
    colnames(true_eff_mat) = paste0("Q", 2:q, "_vs_Q1")

    for (i in 1:n_vars) {
        current_shape = shape[i]
        b = beta_preds[i] * beta_wqs
        min_x = min(preds_trans[, i])

        if (current_shape == "pure_linear") {
            eta_components_raw[, i] = preds_trans[, i] * b
            baseline_components[i] = min_x * b
            for (k in 2:q) {
                true_eff_mat[i + 1, k - 1] = (eval_pts[k] - eval_pts[1]) * b
            }
        } else {
            if (current_shape == "linear_like") {
                pattern = c(1, 2, 3)
            } else if (current_shape == "neg_linear") {
                pattern = c(-1, -2, -3)
            } else if (current_shape == "u_shape") {
                pattern = c(1.5, -3.0, 1.5)
            } else if (current_shape == "inv_u_shape") {
                pattern = c(-1.5, 3.0, -1.5)
            } else if (current_shape == "s_shape") {
                pattern = c(1, -1, 1)
            } else if (current_shape == "threshold") {
                pattern = c(0, 0.5, 4.0)
            } else if (current_shape == "inv_threshold") {
                pattern = c(0, -0.5, -4.0)
            } else {
                pattern = rep(1, df_spline)
            }

            comp_beta = b * pattern
            eta_components_raw[, i] = as.vector(mat_spline_list[[i]] %*% comp_beta)
            baseline_components[i] = as.vector(predict(mat_spline_list[[i]], newx = min_x) %*% comp_beta)

            # [精准匹配计算真理值]：使用同一套 basis_std_true 计算
            for (k in 2:q) {
                b_diff = basis_std_true[k, ] - basis_std_true[1, ]
                true_eff_mat[i + 1, k - 1] = sum(b_diff * comp_beta)
            }
        }
    }

    # 合并 Overall 真理值
    for (k in 2:q) {
        true_eff_mat[1, k - 1] = sum(true_eff_mat[2:(n_vars + 1), k - 1])
    }

    # -----------------------------------------------------------
    # 汇总计算最终偏效应 (Anchoring to 0) 并生成 Y
    # -----------------------------------------------------------
    eta_spline_raw = rowSums(eta_components_raw)
    baseline_effect = sum(baseline_components)

    eta_spline_adjusted = eta_spline_raw - baseline_effect
    y_clean = eta_spline_adjusted + cov_list$eta_cov

    y_observed = add_noise_by_snr(as.vector(y_clean), snr_db = snr_db)

    cols_cov = setdiff(names(cov_list$mm), "eta_cov")
    final_df = cbind(y = y_observed, preds_scaled, cov_list$mm[, cols_cov, drop = FALSE])

    attr(final_df, "true_effect_mat") = true_eff_mat
    attr(final_df, "spline_knots") = global_knots # <--- 存起来
    attr(final_df, "spline_boundary") = global_boundary # <--- 存起来

    return(as.data.frame(final_df))
}

#' Generate Non-linear Binary Model Data (Auto-Balanced)
#' 生成基于自然样条的非线性二分类数据 (支持自动平衡 0/1 比例)
#'
#' @param n_obs integer. 样本量。
#' @param mu_preds numeric vector. 预测变量均值。
#' @param sigma_preds matrix. 预测变量协方差矩阵。
#' @param beta_wqs numeric. 混合物整体效应杠杆系数。
#' @param beta_preds numeric vector. 各物质权重分配的基础系数。
#' @param intercept numeric. 基础截距。如果 target_prop 不为 NULL，此参数将被覆盖。
#' @param target_prop numeric. 目标事件发生率 (如 0.5)。自动计算让 Y=1 比例达到目标的截距。
#' @param link character. 连接函数 ("logit", "probit", "cloglog")。
#' @param snr_db numeric. 信噪比 (添加到潜变量 eta 上)。Inf 表示无额外噪声。
#' @param transform_fun function. 预测变量的转换函数 (如分位数转换)。
#' @param q integer. 分位数层数，用于求真值。
#' @param df_spline integer. 样条自由度。
#' @param seed integer. 随机种子。
#' @param shape character. 混合物的非线性形状模式 ("linear_like", "u_shape", "s_shape")。
#' @param ... 传递给 generate_covariates 的参数。
#' @export
gen_nonlinear_bio_data = function(n_obs = 1000, mu_preds, sigma_preds, beta_wqs = 1, beta_preds,
                                   intercept = 0, target_prop = NULL, link = c("logit", "probit", "cloglog"),
                                   snr_db = Inf, transform_fun = NULL, q = 4, df_spline = 3, seed = NULL,
                                   shape = "linear_like", ...) {
    if (!requireNamespace("splines", quietly = TRUE)) stop("Package 'splines' required")
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' required")

    link = match.arg(link)
    if (!is.null(seed)) set.seed(seed)

    preds_raw = MASS::mvrnorm(n_obs, mu = mu_preds, Sigma = sigma_preds)
    preds_scaled = as.data.frame(scale(preds_raw))
    n_vars = ncol(preds_scaled)
    names(preds_scaled) = paste0("Component", 1:ncol(preds_scaled))

    if (!is.null(transform_fun) && is.function(transform_fun)) {
        preds_trans = transform_fun(preds_scaled)
    } else {
        preds_trans = preds_scaled
    }

    if (length(beta_preds) != n_vars) stop("Length of 'beta_preds' must match n_vars.")

    cov_list = generate_covariates(n_obs = n_obs, ...)

    # -----------------------------------------------------------
    # [架构升级]: 动态推导并锁定节点 (保证数据生成和模型拟合的尺子绝对一致)
    # -----------------------------------------------------------
    eval_pts = 0:(q - 1)
    temp_spline = splines::ns(eval_pts, df = df_spline)
    global_knots = attr(temp_spline, "knots")
    global_boundary = attr(temp_spline, "Boundary.knots")

    # 数据生成侧：强制使用这把尺子
    mat_spline_list = lapply(preds_trans, function(x) {
        splines::ns(x, df = df_spline, knots = global_knots, Boundary.knots = global_boundary)
    })

    # 真值计算侧：也强制使用这把尺子 (无截距版)
    basis_std_true = splines::ns(eval_pts, df = df_spline, knots = global_knots, Boundary.knots = global_boundary, intercept = FALSE)

    if (length(shape) == 1) shape = rep(shape, n_vars)

    eta_components_raw = matrix(0, nrow = n_obs, ncol = n_vars)
    baseline_components = numeric(n_vars)

    true_eff_mat = matrix(0, nrow = n_vars + 1, ncol = q - 1)
    rownames(true_eff_mat) = c("Overall", names(preds_scaled))
    colnames(true_eff_mat) = paste0("Q", 2:q, "_vs_Q1")

    for (i in 1:n_vars) {
        current_shape = shape[i]
        b = beta_preds[i] * beta_wqs
        min_x = min(preds_trans[, i])

        if (current_shape == "pure_linear") {
            # 纯粹线性效应 (Log-OR 尺度)
            eta_components_raw[, i] = preds_trans[, i] * b
            baseline_components[i] = min_x * b

            for (k in 2:q) {
                true_eff_mat[i + 1, k - 1] = (eval_pts[k] - eval_pts[1]) * b
            }
        } else {
            # 样条非线性效应
            if (current_shape == "linear_like") {
                pattern = c(1, 2, 3)
            } else if (current_shape == "neg_linear") {
                pattern = c(-1, -2, -3)
            } else if (current_shape == "u_shape") {
                pattern = c(1.5, -3.0, 1.5)
            } else if (current_shape == "inv_u_shape") {
                pattern = c(-1.5, 3.0, -1.5)
            } else if (current_shape == "s_shape") {
                pattern = c(1, -1, 1)
            } else if (current_shape == "threshold") {
                pattern = c(0, 0.5, 4.0)
            } else if (current_shape == "inv_threshold") {
                pattern = c(0, -0.5, -4.0)
            } else {
                pattern = rep(1, df_spline)
            }

            comp_beta = b * pattern

            eta_components_raw[, i] = as.vector(mat_spline_list[[i]] %*% comp_beta)
            baseline_components[i] = as.vector(predict(mat_spline_list[[i]], newx = min_x) %*% comp_beta)

            # [精准匹配计算真理值]：使用同一套 basis_std_true 计算
            for (k in 2:q) {
                b_diff = basis_std_true[k, ] - basis_std_true[1, ]
                true_eff_mat[i + 1, k - 1] = sum(b_diff * comp_beta)
            }
        }
    }

    for (k in 2:q) {
        true_eff_mat[1, k - 1] = sum(true_eff_mat[2:(n_vars + 1), k - 1])
    }

    # -----------------------------------------------------------
    # 汇总计算最终偏效应 (Anchoring to 0) 并生成 Binary Y
    # -----------------------------------------------------------
    eta_spline_raw = rowSums(eta_components_raw)
    baseline_effect = sum(baseline_components)

    eta_spline_adjusted = eta_spline_raw - baseline_effect
    eta_partial = eta_spline_adjusted + cov_list$eta_cov

    if (!is.null(snr_db) && is.finite(snr_db)) {
        eta_noisy_partial = add_noise_by_snr(eta_partial, snr_db)
    } else {
        eta_noisy_partial = eta_partial
    }

    final_intercept = intercept
    if (!is.null(target_prop)) {
        calc_mean_prob_diff = function(b0) {
            eta_temp = b0 + eta_noisy_partial
            if (link == "logit") {
                p = 1 / (1 + exp(-eta_temp))
            } else if (link == "probit") {
                p = pnorm(eta_temp)
            } else if (link == "cloglog") p = 1 - exp(-exp(eta_temp))
            return(mean(p) - target_prop)
        }
        tryCatch(
            {
                final_intercept = uniroot(calc_mean_prob_diff, interval = c(-50, 50))$root
            },
            error = function(e) {}
        )
    }

    eta_final = final_intercept + eta_noisy_partial

    if (link == "logit") {
        probs = 1 / (1 + exp(-eta_final))
    } else if (link == "probit") {
        probs = pnorm(eta_final)
    } else if (link == "cloglog") probs = 1 - exp(-exp(eta_final))

    y_binary = rbinom(n_obs, size = 1, prob = probs)

    cols_cov = setdiff(names(cov_list$mm), "eta_cov")
    final_df = cbind(y = y_binary, preds_scaled, cov_list$mm[, cols_cov, drop = FALSE])

    attr(final_df, "true_effect_mat") = true_eff_mat # Log-OR 尺度真值
    attr(final_df, "true_prob") = probs
    attr(final_df, "spline_knots") = global_knots # <--- 存起来
    attr(final_df, "spline_boundary") = global_boundary # <--- 存起来
    return(as.data.frame(final_df))
}


#' Generate Non-linear Count Model Data (Poisson)
#' 生成基于自然样条的非线性计数数据 (Poisson Regression)
#'
#' @description
#' Generates predictors, expands using splines, adds covariates, calculates expected counts (lambda)
#' via log-link, and samples count Y from Poisson distribution. Support explicit shape control.
#' 生成预测变量，使用样条扩展，加入协变量，通过对数连接函数计算期望次数，并从泊松分布中生成计数 Y。支持明确的非线性形状控制。
#'
#' @param n_obs integer. Sample size. / 样本量。
#' @param mu_preds numeric vector. Mean vector for predictors.
#' @param sigma_preds matrix. Covariance matrix.
#' @param beta_wqs numeric. Scaling factor.
#' @param beta_preds numeric vector. Coefficients.
#' @param intercept numeric. Intercept (log-scale).
#'   Controls the baseline event rate.
#'   e.g., intercept=0 -> mean count ~1; intercept=2 -> mean count ~7.4.
#'   / 模型截距（对数尺度）。控制基线事件发生率。
#' @param snr_db numeric. SNR for the linear predictor (latent noise).
#' @param transform_fun function.
#' @param q integer. 分位数层数，用于求真值。
#' @param df_spline integer.
#' @param seed integer.
#' @param shape character. 混合物的非线性形状模式 ("linear_like", "u_shape", "s_shape", "threshold")。
#' @param ... arguments passed to generate_covariates.
#'
#' @return data.frame. Contains Y (count), predictors, covariates.
#' @export
gen_nonlinear_count_data = function(n_obs = 1000, mu_preds, sigma_preds, beta_wqs = 1, beta_preds,
                                     intercept = 0, snr_db = Inf, transform_fun = NULL, q = 4, df_spline = 3,
                                     seed = NULL, shape = "linear_like", ...) {
    if (!requireNamespace("splines", quietly = TRUE)) stop("Package 'splines' required")
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' required")
    if (!is.null(seed)) set.seed(seed)

    preds_raw = MASS::mvrnorm(n_obs, mu = mu_preds, Sigma = sigma_preds)
    preds_scaled = as.data.frame(scale(preds_raw))
    n_vars = ncol(preds_scaled)
    names(preds_scaled) = paste0("Component", 1:n_vars)

    if (!is.null(transform_fun) && is.function(transform_fun)) {
        preds_trans = transform_fun(preds_scaled)
    } else {
        preds_trans = preds_scaled
    }

    if (length(beta_preds) != n_vars) stop("Length of 'beta_preds' must match n_vars.")

    cov_list = generate_covariates(n_obs = n_obs, ...)

    # -----------------------------------------------------------
    # [架构升级]: 动态推导并锁定节点 (保证数据生成和模型拟合的尺子绝对一致)
    # -----------------------------------------------------------
    eval_pts = 0:(q - 1)
    temp_spline = splines::ns(eval_pts, df = df_spline)
    global_knots = attr(temp_spline, "knots")
    global_boundary = attr(temp_spline, "Boundary.knots")

    # 数据生成侧：强制使用这把尺子
    mat_spline_list = lapply(preds_trans, function(x) {
        splines::ns(x, df = df_spline, knots = global_knots, Boundary.knots = global_boundary)
    })

    # 真值计算侧：也强制使用这把尺子 (无截距版)
    basis_std_true = splines::ns(eval_pts, df = df_spline, knots = global_knots, Boundary.knots = global_boundary, intercept = FALSE)

    if (length(shape) == 1) shape = rep(shape, n_vars)

    eta_components_raw = matrix(0, nrow = n_obs, ncol = n_vars)
    baseline_components = numeric(n_vars)

    true_eff_mat = matrix(0, nrow = n_vars + 1, ncol = q - 1)
    rownames(true_eff_mat) = c("Overall", names(preds_scaled))
    colnames(true_eff_mat) = paste0("Q", 2:q, "_vs_Q1")

    for (i in 1:n_vars) {
        current_shape = shape[i]
        b = beta_preds[i] * beta_wqs
        min_x = min(preds_trans[, i])

        if (current_shape == "pure_linear") {
            # 纯粹线性效应 (Log-RR 尺度)
            eta_components_raw[, i] = preds_trans[, i] * b
            baseline_components[i] = min_x * b

            for (k in 2:q) {
                true_eff_mat[i + 1, k - 1] = (eval_pts[k] - eval_pts[1]) * b
            }
        } else {
            # 样条非线性效应
            if (current_shape == "linear_like") {
                pattern = c(1, 2, 3)
            } else if (current_shape == "neg_linear") {
                pattern = c(-1, -2, -3)
            } else if (current_shape == "u_shape") {
                pattern = c(1.5, -3.0, 1.5)
            } else if (current_shape == "inv_u_shape") {
                pattern = c(-1.5, 3.0, -1.5)
            } else if (current_shape == "s_shape") {
                pattern = c(1, -1, 1)
            } else if (current_shape == "threshold") {
                pattern = c(0, 0.5, 4.0)
            } else if (current_shape == "inv_threshold") {
                pattern = c(0, -0.5, -4.0)
            } else {
                pattern = rep(1, df_spline)
            }

            comp_beta = b * pattern

            eta_components_raw[, i] = as.vector(mat_spline_list[[i]] %*% comp_beta)
            baseline_components[i] = as.vector(predict(mat_spline_list[[i]], newx = min_x) %*% comp_beta)

            # [精准匹配计算真理值]：使用同一套 basis_std_true 计算
            for (k in 2:q) {
                b_diff = basis_std_true[k, ] - basis_std_true[1, ]
                true_eff_mat[i + 1, k - 1] = sum(b_diff * comp_beta)
            }
        }
    }

    for (k in 2:q) {
        true_eff_mat[1, k - 1] = sum(true_eff_mat[2:(n_vars + 1), k - 1])
    }

    # -----------------------------------------------------------
    # 汇总计算最终偏效应 (Anchoring to 0) 并生成 Count Y
    # -----------------------------------------------------------
    eta_spline_raw = rowSums(eta_components_raw)
    baseline_effect = sum(baseline_components)

    eta_spline_adjusted = eta_spline_raw - baseline_effect
    eta_partial = eta_spline_adjusted + cov_list$eta_cov

    if (!is.null(snr_db) && is.finite(snr_db)) {
        eta_noisy_partial = add_noise_by_snr(eta_partial, snr_db)
    } else {
        eta_noisy_partial = eta_partial
    }

    eta_final = intercept + eta_noisy_partial
    lambda = exp(eta_final)

    if (any(lambda > 10000)) warning("Extremely high lambda values.")

    y_count = rpois(n_obs, lambda = lambda)

    cols_cov = setdiff(names(cov_list$mm), "eta_cov")
    final_df = cbind(y = y_count, preds_scaled, cov_list$mm[, cols_cov, drop = FALSE])

    attr(final_df, "true_effect_mat") = true_eff_mat # Log-RR 尺度真值
    attr(final_df, "spline_knots") = global_knots # <--- 存起来
    attr(final_df, "spline_boundary") = global_boundary # <--- 存起来
    return(as.data.frame(final_df))
}
