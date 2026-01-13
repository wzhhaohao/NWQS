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
                               prob_cat = c(1/3, 1/3, 1/3),
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
#' Generate Non-linear Binary Model Data with Splines, Link Function, and SNR
#' 生成基于自然样条的非线性二分类数据 (支持 Logit/Probit/Cloglog 及 SNR 加噪)
#'
#' @description
#' Generates predictors, expands using splines, adds covariates, adds Gaussian noise 
#' to the linear predictor based on SNR, and samples binary Y using the specified link function.
#' 生成预测变量，使用样条扩展，加入协变量，根据信噪比向线性预测值添加高斯噪声，
#' 最后通过指定的连接函数生成二分类 Y。
#'
#' @param n_obs integer. Sample size. / 样本量。
#' @param mu_preds numeric vector. Mean vector for predictors. / 预测变量均值。
#' @param sigma_preds matrix. Covariance matrix for predictors. / 预测变量协方差矩阵。
#' @param beta_wqs numeric. Scaling factor for coefficients. / 系数缩放因子。
#' @param beta_preds numeric vector. Coefficients for the predictors. Length must equal n_vars.
#' @param intercept numeric. Intercept of the model. / 模型截距。
#' @param link character. Link function: "logit", "probit", "cloglog".
#' @param snr_db numeric. Signal-to-Noise Ratio (dB) for the linear predictor. 
#'   Default is Inf (no noise). / 线性预测值的信噪比。默认为 Inf（无额外噪声）。
#' @param transform_fun function. Function to transform predictors. / 变换函数。
#' @param df_spline integer. Degrees of freedom for natural splines. / 自然样条自由度。
#' @param seed integer. Random seed. / 随机种子。
#' @param ... arguments passed to generate_covariates.
#' 
#' @return data.frame. Contains Y (0/1), predictors, covariates.
#' @export
# gen_nonlinear_bio_data = function(n_obs = 1000, 
#                                   mu_preds, 
#                                   sigma_preds, 
#                                   beta_wqs = 1, 
#                                   beta_preds,
#                                   intercept = 0,
#                                   link = c("logit", "probit", "cloglog"), 
#                                   snr_db = Inf,  # 新增参数，避免传入 ... 导致报错
#                                   transform_fun = NULL, 
#                                   df_spline = 3,
#                                   seed = NULL,
#                                   ...) {
  
#     if (!requireNamespace("splines", quietly = TRUE)) stop("Package 'splines' required")
#     if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' required")
  
#     # Validate link argument
#     link = match.arg(link)

#     if (!is.null(seed)) set.seed(seed)
  
#     # 1. Generate Raw Predictors
#     preds_raw = MASS::mvrnorm(n_obs, mu = mu_preds, Sigma = sigma_preds)
#     preds_scaled = as.data.frame(scale(preds_raw)) 
#     n_vars = ncol(preds_scaled)
#     names(preds_scaled) = paste0("Component", 1:ncol(preds_scaled))
  
#     # 2. Transform Predictors
#     if (!is.null(transform_fun) && is.function(transform_fun)) {
#         preds_trans = transform_fun(preds_scaled)
#     } else {
#         preds_trans = preds_scaled
#     }
  
#     # 3. Create Spline Basis Matrix
#     mat_spline_list = lapply(preds_trans, function(x) splines::ns(x, df = df_spline))
#     mat_spline_full = do.call(cbind, mat_spline_list)

#     if (length(beta_preds) != n_vars) {
#         stop(sprintf("Length of 'beta_preds' (%d) must match n_vars (%d).", length(beta_preds), n_vars))
#     }
  
#     # 4. Generate Covariates
#     # 注意：此时 snr_db 已经被上面的参数捕获，不会进入 ...，所以不会报错
#     cov_list = generate_covariates(n_obs = n_obs, ...)
  
#     # 5. Calculate Linear Predictor (Eta) - Pure Signal
#     beta_expanded = rep(beta_preds * beta_wqs, each = df_spline)
    
#     eta_clean = intercept + 
#                 (as.matrix(mat_spline_full) %*% beta_expanded) + 
#                 cov_list$eta_cov
    
#     eta_clean = as.vector(eta_clean)

#     # 6. Add Noise to Linear Predictor based on SNR
#     # 使用您提供的 add_noise_by_snr 函数
#     # 如果 snr_db 是 Inf，则不加噪
#     if (!is.null(snr_db) && is.finite(snr_db)) {
#         eta_final = add_noise_by_snr(eta_clean, snr_db)
#     } else {
#         eta_final = eta_clean
#     }

#     # 7. Calculate Probability (Inverse Link Function) using Noisy Eta
#     if (link == "logit") {
#         probs = 1 / (1 + exp(-eta_final))
#     } else if (link == "probit") {
#         probs = pnorm(eta_final)
#     } else if (link == "cloglog") {
#         probs = 1 - exp(-exp(eta_final))
#     }

#     # 8. Generate Binary Outcome
#     y_binary = rbinom(n_obs, size = 1, prob = probs)
  
#     # 9. Combine Results
#     cols_cov = setdiff(names(cov_list$mm), "eta_cov")
  
#     final_df = cbind(y = y_binary, 
#                      preds_scaled, 
#                      cov_list$mm[, cols_cov, drop = FALSE])
  
#     # Attributes for debugging
#     attr(final_df, "true_prob") = probs
#     attr(final_df, "link_used") = link
#     attr(final_df, "snr_db") = snr_db
    
#     return(as.data.frame(final_df))
# }




#' Generate Non-linear Binary Model Data (Auto-Balanced)
#' 生成基于自然样条的非线性二分类数据 (支持自动平衡 0/1 比例)
#'
#' @param target_prop numeric. Target proportion of Y=1 (e.g., 0.5). 
#'   If NULL, uses the provided 'intercept'.
#'   目标事件发生率。如果提供（如 0.5），程序会自动忽略 'intercept' 参数，
#'   计算出能让 Y=1 比例达到该目标的截距。
#' @param ... (其他参数同前)
#' @export
gen_nonlinear_bio_data = function(n_obs = 1000, 
                                  mu_preds, 
                                  sigma_preds, 
                                  beta_wqs = 1, 
                                  beta_preds,
                                  intercept = 0, # 如果 target_prop 不为 NULL，此参数将被覆盖
                                  target_prop = NULL, # <--- 新增核心参数
                                  link = c("logit", "probit", "cloglog"), 
                                  snr_db = Inf,
                                  transform_fun = NULL, 
                                  df_spline = 3,
                                  seed = NULL,
                                  ...) {
  
    if (!requireNamespace("splines", quietly = TRUE)) stop("Package 'splines' required")
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' required")
  
    link = match.arg(link)
    if (!is.null(seed)) set.seed(seed)
  
    # 1. Generate Raw Predictors
    preds_raw = MASS::mvrnorm(n_obs, mu = mu_preds, Sigma = sigma_preds)
    preds_scaled = as.data.frame(scale(preds_raw)) 
    n_vars = ncol(preds_scaled)
    names(preds_scaled) = paste0("Component", 1:ncol(preds_scaled))
  
    # 2. Transform Predictors
    if (!is.null(transform_fun) && is.function(transform_fun)) {
        preds_trans = transform_fun(preds_scaled)
    } else {
        preds_trans = preds_scaled
    }
  
    # 3. Create Spline Basis Matrix
    mat_spline_list = lapply(preds_trans, function(x) splines::ns(x, df = df_spline))
    mat_spline_full = do.call(cbind, mat_spline_list)

    if (length(beta_preds) != n_vars) {
        stop(sprintf("Length of 'beta_preds' (%d) must match n_vars (%d).", length(beta_preds), n_vars))
    }
  
    # 4. Generate Covariates
    cov_list = generate_covariates(n_obs = n_obs, ...)
  
    # 5. Calculate Partial Linear Predictor (Excluding Intercept)
    # -----------------------------------------------------------
    beta_expanded = rep(beta_preds * beta_wqs, each = df_spline)
    
    # eta_partial 包含了样条效应 + 协变量效应 (cov_list$eta_cov 里通常包含了一个默认为0的intercept，这没关系，我们会在后面调整)
    eta_partial = (as.matrix(mat_spline_full) %*% beta_expanded) + cov_list$eta_cov
    eta_partial = as.vector(eta_partial)

    # 6. Add Noise to Linear Predictor (Latent Variable Noise)
    # -----------------------------------------------------------
    # 我们先加噪，再找截距。这样能保证最终观测到的 Y 是平衡的。
    if (!is.null(snr_db) && is.finite(snr_db)) {
        eta_noisy_partial = add_noise_by_snr(eta_partial, snr_db)
    } else {
        eta_noisy_partial = eta_partial
    }

    # 7. Auto-Calculate Intercept if target_prop is set
    # -----------------------------------------------------------
    final_intercept = intercept
    
    if (!is.null(target_prop)) {
        if (target_prop <= 0 || target_prop >= 1) stop("target_prop must be between 0 and 1")
        
        # 定义目标函数：给定截距 b0，计算平均概率与目标的差值
        calc_mean_prob_diff = function(b0) {
            eta_temp = b0 + eta_noisy_partial
            
            if (link == "logit") {
                p = 1 / (1 + exp(-eta_temp))
            } else if (link == "probit") {
                p = pnorm(eta_temp)
            } else if (link == "cloglog") {
                p = 1 - exp(-exp(eta_temp))
            }
            return(mean(p) - target_prop)
        }
        
        # 使用 uniroot 求解让差值为 0 的截距
        # 搜索范围 -50 到 50 通常足够覆盖绝大多数情况
        tryCatch({
            root_res = uniroot(calc_mean_prob_diff, interval = c(-50, 50), extendInt = "yes")
            final_intercept = root_res$root
        }, error = function(e) {
            warning("Could not find optimal intercept. Using default.")
        })
    }

    # 8. Calculate Final Probabilities with optimized intercept
    # -----------------------------------------------------------
    eta_final = final_intercept + eta_noisy_partial
    
    if (link == "logit") {
        probs = 1 / (1 + exp(-eta_final))
    } else if (link == "probit") {
        probs = pnorm(eta_final)
    } else if (link == "cloglog") {
        probs = 1 - exp(-exp(eta_final))
    }

    # 9. Generate Binary Outcome
    y_binary = rbinom(n_obs, size = 1, prob = probs)
  
    # 10. Combine Results
    cols_cov = setdiff(names(cov_list$mm), "eta_cov")
  
    final_df = cbind(y = y_binary, 
                     preds_scaled, 
                     cov_list$mm[, cols_cov, drop = FALSE])
  
    # Attributes for debugging
    attr(final_df, "true_prob") = probs
    attr(final_df, "link_used") = link
    attr(final_df, "snr_db") = snr_db
    attr(final_df, "intercept_used") = final_intercept # 记录自动计算的截距
    attr(final_df, "target_prop") = target_prop
    
    return(as.data.frame(final_df))
}
# TODO：gen_nonlinear_multinonial_data
# 看看文章，暂时还不认为需要些

# TODO: 设计一个poisson分布的数据生成函数
#' Generate Non-linear Count Model Data (Poisson)
#' 生成基于自然样条的非线性计数数据 (Poisson Regression)
#'
#' @description
#' Generates predictors, expands using splines, adds covariates, calculates expected counts (lambda)
#' via log-link, and samples count Y from Poisson distribution.
#' 生成预测变量，使用样条扩展，加入协变量，通过对数连接函数计算期望次数，并从泊松分布中生成计数 Y。
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
#' @param df_spline integer.
#' @param seed integer.
#' @param ... arguments passed to generate_covariates.
#'
#' @return data.frame. Contains Y (count), predictors, covariates.
#' @export
gen_nonlinear_count_data = function(n_obs = 1000, 
                                mu_preds, 
                                sigma_preds, 
                                beta_wqs = 1, 
                                beta_preds,
                                intercept = 0, # Log-scale intercept
                                snr_db = Inf, 
                                transform_fun = NULL, 
                                df_spline = 3,
                                seed = NULL,
                                ...) {

if (!is.null(seed)) set.seed(seed)

# 1. Generate & Transform Predictors
preds_raw = MASS::mvrnorm(n_obs, mu = mu_preds, Sigma = sigma_preds)
preds_scaled = as.data.frame(scale(preds_raw)) 
names(preds_scaled) = paste0("Component", 1:ncol(preds_scaled))

if (!is.null(transform_fun) && is.function(transform_fun)) {
    preds_trans = transform_fun(preds_scaled)
} else {
    preds_trans = preds_scaled
}

# 2. Spline Expansion
mat_spline_list = lapply(preds_trans, function(x) splines::ns(x, df = df_spline))
mat_spline_full = do.call(cbind, mat_spline_list)

if (length(beta_preds) != ncol(preds_scaled)) {
    stop("Length of beta_preds must match number of predictors.")
}

# 3. Covariates
cov_list = generate_covariates(n_obs = n_obs, ...)

# 4. Calculate Linear Predictor (Eta)
# Eta is on the LOG scale
beta_expanded = rep(beta_preds * beta_wqs, each = df_spline)

eta_clean = intercept + 
            (as.matrix(mat_spline_full) %*% beta_expanded) + 
            cov_list$eta_cov

eta_clean = as.vector(eta_clean)

# 5. Add Noise to Eta (Optional Overdispersion source)
if (!is.null(snr_db) && is.finite(snr_db)) {
    eta_final = add_noise_by_snr(eta_clean, snr_db)
} else {
    eta_final = eta_clean
}

# 6. Inverse Link: Log -> Count Mean (Lambda)
# lambda = exp(eta)
lambda = exp(eta_final)

# [安全检查] 防止 lambda 过大导致 rpois 溢出或产生不切实际的数据
# 在 Poisson 回归中，系数稍微大一点（比如 5），exp(5) 就变成 148 了，exp(10) 就是 22000
if (any(lambda > 10000)) {
    warning("Some lambda values are extremely high (>10000). Consider reducing betas or intercept.")
}

# 7. Sample Outcome (Poisson)
y_count = rpois(n_obs, lambda = lambda)

# 8. Combine Results
cols_cov = setdiff(names(cov_list$mm), "eta_cov")
final_df = cbind(y = y_count, 
                    preds_scaled, 
                    cov_list$mm[, cols_cov, drop = FALSE])

attr(final_df, "true_lambda") = lambda
attr(final_df, "snr_db") = snr_db

return(as.data.frame(final_df))
}


# TODO:设计一个重复测量的数据生成函数