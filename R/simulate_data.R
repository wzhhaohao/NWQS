#' Add Gaussian noise given target SNR / 根据信噪比添加高斯噪音
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

    if (power_signal == 0) return(signal_vec)

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

#' Generate Covariates / 生成背景协变量
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
                               beta_cat = c(0.1, -0.4, 0.7),
                               prob_bin = 0.5,
                               prob_cat = c(1/3, 1/3, 1/3)) {
  
    x_cont = rnorm(n_obs, 0, 1)
    x_bin_raw = rbinom(n_obs, 1, prob_bin)
    x_cat_raw = sample(1:3, n_obs, replace = TRUE, prob = prob_cat)

    x_bin = factor(x_bin_raw, levels = c(0, 1))
    x_cat = factor(x_cat_raw, levels = 1:3)

    # Calculate linear predictor: eta = X * beta
    eta_cov = beta_cont * x_cont + beta_bin * x_bin_raw + beta_cat[x_cat_raw]

    df_raw = data.frame(x_cont, x_bin, x_cat)
    # mat_model = model.matrix(~ x_cont + x_bin + x_cat, data = df_raw)[, -1, drop = FALSE]
  
    df_result = as.data.frame(cbind(eta_cov = eta_cov, df_raw))

    list(mm = df_result, original = df_raw, eta_cov = eta_cov)
}

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

#' Generate Linear Model Data / 生成线性模型数据
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
    final_df = cbind(y = y_observed, 
                     preds_scaled, 
                     cov_list$mm[, cols_cov, drop = FALSE])

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
                              df_spline = 3,
                              seed = NULL,
                              ...) {
  
    if (!requireNamespace("splines", quietly = TRUE)) stop("Package 'splines' required")
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' required")
  
    if (!is.null(seed)) set.seed(seed)
  
    # Generate Raw Predictors (Multivariate Normal)
    preds_raw = MASS::mvrnorm(n_obs, mu = mu_preds, Sigma = sigma_preds)
    preds_scaled = as.data.frame(scale(preds_raw)) # Z-score standardization
    n_vars = ncol(preds_scaled)

    names(preds_scaled) = paste0("Component", 1:ncol(preds_scaled))
  
    # Transform Predictors (Optional)
    # FIXME: 关于 transform_fun 中，对分位数化后-1了，所以要保证所有逻辑都一致，其次就是到底要不要把transform_fun放在外面
    if (!is.null(transform_fun) && is.function(transform_fun)) {
        preds_trans = transform_fun(preds_scaled)
    } else {
        preds_trans = preds_scaled
    }
  
    # Create Spline Basis Matrix (Non-linear expansion)
    # Matrix dimension: n_obs x (n_vars * df_spline)
    mat_spline_list = lapply(preds_trans, function(x) splines::ns(x, df = df_spline))
    mat_spline_full = do.call(cbind, mat_spline_list)

    # Dimension Check
    # Check if beta_preds matches the number of variables (since we repeat it)
    if (length(beta_preds) != n_vars) {
        stop(sprintf(
            "Length of 'beta_preds' (%d) must match the number of variables n_vars (%d).", 
            length(beta_preds), n_vars
        ))
    }
  
    # Generate Covariates & Linear Effects
    cov_list = generate_covariates(n_obs = n_obs, ...)
  
    # Calculate Clean Signal (Y_clean)
    # Logic: Repeat the single coefficient for a variable across all its spline basis functions
    # eta = Spline_Matrix * (beta_expanded)
    beta_expanded = rep(beta_preds * beta_wqs, each = df_spline)
    eta_spline = as.matrix(mat_spline_full) %*% beta_expanded
    
    y_clean = eta_spline + cov_list$eta_cov
  
    # Add Noise based on SNR
    y_observed = add_noise_by_snr(as.vector(y_clean), snr_db = snr_db)
  
    # Combine Results
    cols_cov = setdiff(names(cov_list$mm), "eta_cov")
  
    final_df = cbind(y = y_observed, 
                     preds_scaled, 
                     cov_list$mm[, cols_cov, drop = FALSE])
  
    return(as.data.frame(final_df))
}


# TODO：gen_nonlinear_bio_data

# TODO：gen_nonlinear_multinonial_data

# TODO:设计一个重复测量的数据生成函数