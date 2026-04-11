#' @title 基于目标信噪比注入高斯白噪声 (Add Gaussian Noise given Target SNR)
#'
#' @description
#' 计算纯净信号的功率，并根据指定的目标信噪比 (Signal-to-Noise Ratio, SNR) 注入正态分布的白噪声。
#'
#' @details
#' \strong{流行病学模拟意义:} \cr
#' 在方法学研究中，该函数用于模拟**未测量的混杂因素 (Unmeasured Confounding)** 或 **暴露评估的测量误差 (Measurement Error)**。
#' 信噪比的计算公式为：
#' $$SNR_{dB} = 10 \log_{10}\left(\frac{P_{signal}}{P_{noise}}\right)$$
#' \code{snr_db} 越低，表明临床或环境数据中的噪声占比越大，这被用来严苛地测试 NWQS 在低信噪比下提取真实暴露权重的稳健性。
#'
#' @param signal_vec Numeric vector。无噪声的纯净信号（如真实的潜在线性预测子 \eqn{\eta}）。
#' @param snr_db Numeric。目标信噪比 (dB)。值越大代表信号质量越好，\code{Inf} 代表不添加噪声。
#' @return Numeric vector。叠加了测量误差/环境噪声后的观测信号。
#' @export
add_noise_by_snr <- function(signal_vec, snr_db) {
    stopifnot(is.numeric(signal_vec), length(signal_vec) > 1)

    # Calculate signal power: P_signal = E[(x - mu)^2]
    power_signal <- mean((signal_vec - mean(signal_vec))^2)

    if (power_signal == 0) {
        return(signal_vec)
    }

    # Calculate noise standard deviation based on SNR formula
    # SNR_dB = 10 * log10(P_signal / P_noise)
    sigma_noise <- sqrt(power_signal / 10^(snr_db / 10))

    noise_vec <- rnorm(length(signal_vec), mean = 0, sd = sigma_noise)

    return(signal_vec + noise_vec)
}

#' @title 生成具有特定相关性结构的协方差矩阵 (Generate Covariance Matrix)
#'
#' @description
#' 根据指定的模式生成正定对称的相关系数/协方差矩阵，专为模拟高维共线性暴露数据而设计。
#'
#' @details
#' \strong{环境暴露多重共线性模拟:} \cr
#' 公共卫生领域中的混合物数据（如全氟化合物 PFAS、多环芳烃 PAHs）通常具有高度相关性。
#' 此函数提供的模式完美契合不同的暴露场景：
#' \itemize{
#'   \item \code{"high"}: 模拟同源暴露（如同一污染源排放的多种同系物），所有变量间具有高强度基础相关性。
#'   \item \code{"mixed"}: 模拟真实世界中具有“区块对角 (Block Diagonal)”特征的暴露族群（部分高度相关，部分独立）。
#' }
#' 算法内部强制进行特征值修复 (Eigenvalue Repair)，确保生成的矩阵在数学上绝对正定 (Positive-Definite)。
#'
#' @param n_vars Integer。混合物组分（变量）的数量。
#' @param mode Character。相关性模式，可选 \code{"low"}, \code{"medium"}, \code{"high"}, \code{"mixed"}。
#' @param rho Numeric。基础相关系数强度（仅针对特定模式生效），默认为 0.7。
#' @param seed Integer 或 \code{NULL}。随机种子，用于控制 Monte Carlo 模拟的可重复性。
#' @return Matrix。正定对称的相关系数矩阵。
#' @export
generate_sigma <- function(n_vars, mode = c("medium", "low", "high", "mixed"), rho = 0.7, seed = NULL) {
    mode <- match.arg(mode)
    if (!is.null(seed)) set.seed(seed)

    if (mode == "low") {
        # Identity matrix + small noise
        A <- diag(n_vars) + matrix(runif(n_vars^2, -0.1, 0.1), nrow = n_vars)
        sigma <- cov2cor(t(A) %*% A)
    } else if (mode == "medium") {
        # Random Gram matrix
        A <- matrix(runif(n_vars^2, -1, 1), ncol = n_vars)
        sigma <- cov2cor(t(A) %*% A)
    } else if (mode == "high") {
        # Factor model structure
        sigma <- matrix(rho, nrow = n_vars, ncol = n_vars)
        diag(sigma) <- 1

        # Add jitter
        noise <- matrix(runif(n_vars^2, -0.05, 0.05), nrow = n_vars)
        sigma <- sigma + (noise + t(noise)) / 2

        # Repair eigenvalues to ensure positive definiteness
        eig <- eigen(sigma)
        val <- pmax(eig$values, 0.01)
        sigma <- cov2cor(eig$vectors %*% diag(val) %*% t(eig$vectors))
    } else if (mode == "mixed") {
        # Block diagonal structure
        split_idx <- floor(n_vars / 2)
        s1 <- split_idx
        s2 <- n_vars - split_idx

        # Block 1: High correlation
        B1 <- matrix(0.8, nrow = s1, ncol = s1)
        diag(B1) <- 1

        # Block 2: Medium correlation (Random)
        A2 <- matrix(runif(s2^2, -1, 1), ncol = s2)
        B2 <- cov2cor(t(A2) %*% A2)

        # Combine blocks
        sigma <- matrix(0, nrow = n_vars, ncol = n_vars)
        sigma[1:s1, 1:s1] <- B1
        sigma[(s1 + 1):n_vars, (s1 + 1):n_vars] <- B2

        # Add noise and repair eigenvalues
        noise <- matrix(runif(n_vars^2, -0.1, 0.1), nrow = n_vars)
        sigma <- sigma + (noise + t(noise)) / 2

        eig <- eigen(sigma)
        val <- pmax(eig$values, 0.01)
        sigma <- cov2cor(eig$vectors %*% diag(val) %*% t(eig$vectors))
    }

    return(sigma)
}

#' @title 生成流行病学测量的混杂因素集 (Generate Covariates)
#'
#' @description
#' 随机生成包含连续型、二分类和多分类变量的混杂因素数据集，并计算它们的真实线性效应。
#'
#' @details
#' 此函数用于模拟流行病学研究中常规收集的、需要被模型控制的协变量（如年龄、性别、BMI 分组或吸烟状态）。
#' 将这些协变量注入生成数据，旨在验证模型在存在可测量混杂 (Measured Confounding) 时，
#' 对混合物主效应进行无偏估计 (Unbiased Estimation) 的能力。
#'
#' @param n_obs Integer。模拟队列的样本量。
#' @param beta_cont Numeric。连续型变量（如年龄的标准化值）的真实回归系数。
#' @param beta_bin Numeric。二分类变量（如性别）的真实回归系数。
#' @param beta_cat Numeric vector。多分类变量（如吸烟状态类别）的真实回归系数向量。
#' @param prob_bin Numeric。二分类变量发生率。
#' @param prob_cat Numeric vector。多分类变量各水平的分布概率。
#' @param Intercept Numeric。基线截距。
#' @return List。包含生成的协变量数据框 \code{original}、设计矩阵形式的 \code{mm}，以及真实的线性预测子贡献 \code{eta_cov}。
#' @export
generate_covariates <- function(n_obs = 1000,
                                beta_cont = 0.5,
                                beta_bin = -0.8,
                                beta_cat = c(0, -0.5, 0.7),
                                prob_bin = 0.5,
                                prob_cat = c(1 / 3, 1 / 3, 1 / 3),
                                Intercept = 0) {
    x_cont <- rnorm(n_obs, 0, 1)
    x_bin_raw <- rbinom(n_obs, 1, prob_bin)
    x_cat_raw <- sample(1:3, n_obs, replace = TRUE, prob = prob_cat)

    x_bin <- factor(x_bin_raw, levels = c(0, 1))
    x_cat <- factor(x_cat_raw, levels = 1:3)

    # Calculate linear predictor: eta = X * beta
    eta_cov <- beta_cont * x_cont + beta_bin * x_bin_raw + beta_cat[x_cat_raw] + Intercept

    df_raw <- data.frame(x_cont, x_bin, x_cat)
    # mat_model = model.matrix(~ x_cont + x_bin + x_cat, data = df_raw)[, -1, drop = FALSE]

    df_result <- as.data.frame(cbind(eta_cov = eta_cov, df_raw))

    list(mm = df_result, original = df_raw, eta_cov = eta_cov)
}

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

#' @title 生成线性混合物效应的连续型数据 (Generate Linear Model Data)
#'
#' @description
#' 为 Monte Carlo 模拟生成标准线性响应数据。暴露组分通过多元正态分布生成，
#' 可选经过分位数转换，最后与协变量的线性预测子组合，并根据信噪比引入测量误差。
#'
#' @export
generate_linear_data <- function(n_obs = 1000,
                                 mu_preds,
                                 sigma_preds,
                                 beta_wqs = 1,
                                 beta_preds,
                                 snr_db = 10,
                                 transform_fun = NULL,
                                 seed = NULL,
                                 ...) {
    if (!is.null(seed)) set.seed(seed)
    if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("Package 'MASS' required")
    }

    preds_raw <- MASS::mvrnorm(
        n_obs,
        mu = mu_preds,
        Sigma = sigma_preds
    )

    preds_scaled <- as.data.frame(scale(preds_raw))
    names(preds_scaled) <- paste0("Component", seq_len(ncol(preds_scaled)))

    if (!is.null(transform_fun) && is.function(transform_fun)) {
        preds_final <- transform_fun(preds_scaled)
    } else {
        preds_final <- preds_scaled
    }

    preds_final <- as.data.frame(preds_final)
    names(preds_final) <- names(preds_scaled)

    cov_list <- generate_covariates(n_obs = n_obs, ...)

    beta_preds <- beta_wqs * beta_preds

    y_clean <- as.matrix(preds_final) %*% beta_preds + cov_list$eta_cov
    y_observed <- add_noise_by_snr(as.vector(y_clean), snr_db = snr_db)

    cols_cov <- setdiff(names(cov_list$mm), "eta_cov")

    final_df <- cbind(
        y = y_observed,
        preds_scaled,
        cov_list$mm[, cols_cov, drop = FALSE]
    )

    as.data.frame(final_df)
}

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

#' @title 生成非线性自然样条的连续型剂量反应数据 (Generate Non-linear Data)
#'
#' @description
#' 高级数据生成函数。此机制允许为每种暴露成分设定特定的非线性剂量反应轨迹（如 U 型、S 型或阈值效应），
#' 随后利用自然三次样条 (Natural Cubic Splines) 精确映射这些形状，并合成最终的连续型结局变量。
#'
#' @details
#' \strong{因果结构与基准事实 (Ground Truth) 的严格对齐:} \cr
#' 为了在 Monte Carlo 模拟中严谨地评估 NWQS 方法捕捉非线性的能力，本函数在生成数据时直接在底层推导出
#' 样条的基础节点 (\code{knots})，并将其作为\strong{不可变的度量尺}贯穿始终。函数同时会在属性 \code{true_effect_mat}
#' 中输出真正的偏效应对比值（如真正的 Q4 vs Q1 效应量）。这种将“因果形状生成”与“理论效应量计算”完美绑定的设计，
#' 彻底消除了评估模型性能时的基准误差。
#'
#' @param shape Character 或 Character vector。控制各组分真实因果关系曲线的形态。
#'   可选 \code{"linear_like"}（准线性）, \code{"u_shape"}（非单调 U 型）, \code{"inv_u_shape"}（倒 U 型）,
#'   \code{"s_shape"}（S 型阈值）或 \code{"threshold"}（硬阈值安全水平）。
#' @export
gen_nonlinear_data <- function(n_obs = 1000,
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

    preds_raw <- MASS::mvrnorm(n_obs, mu = mu_preds, Sigma = sigma_preds)
    preds_scaled <- as.data.frame(scale(preds_raw))
    n_vars <- ncol(preds_scaled)
    names(preds_scaled) <- paste0("Component", 1:n_vars)

    if (!is.null(transform_fun) && is.function(transform_fun)) {
        preds_trans <- transform_fun(preds_scaled)
    } else {
        preds_trans <- preds_scaled
    }

    mat_spline_list <- lapply(preds_trans, function(x) splines::ns(x, df = df_spline))

    if (length(beta_preds) != n_vars) stop("Length of 'beta_preds' must match n_vars.")

    cov_list <- generate_covariates(n_obs = n_obs, ...)

    # -----------------------------------------------------------
    # [架构升级]: 完美混合 纯线性(pure_linear) 与 非线性样条
    # -----------------------------------------------------------

    # 动态推导并锁定节点 (保证数据生成和模型拟合的尺子绝对一致)
    eval_pts <- 0:(q - 1)

    # 利用 splines 内部算法自动计算最合理的 knots (支持任意的 q)
    temp_spline <- splines::ns(eval_pts, df = df_spline)
    global_knots <- attr(temp_spline, "knots")
    global_boundary <- attr(temp_spline, "Boundary.knots")

    # 打印出来看看，如果是 q=4，它会自动算出 1 和 2 (极度聪明)
    # print(global_knots)

    # 数据生成侧：强制使用这把尺子
    mat_spline_list <- lapply(preds_trans, function(x) {
        splines::ns(x, df = df_spline, knots = global_knots, Boundary.knots = global_boundary)
    })

    # 真值计算侧：也强制使用这把尺子 (无截距版)
    basis_std_true <- splines::ns(eval_pts, df = df_spline, knots = global_knots, Boundary.knots = global_boundary, intercept = FALSE)

    if (length(shape) == 1) shape <- rep(shape, n_vars)

    eta_components_raw <- matrix(0, nrow = n_obs, ncol = n_vars)
    baseline_components <- numeric(n_vars)

    true_eff_mat <- matrix(0, nrow = n_vars + 1, ncol = q - 1)
    rownames(true_eff_mat) <- c("Overall", names(preds_scaled))
    colnames(true_eff_mat) <- paste0("Q", 2:q, "_vs_Q1")

    for (i in 1:n_vars) {
        current_shape <- shape[i]
        b <- beta_preds[i] * beta_wqs
        min_x <- min(preds_trans[, i])

        if (current_shape == "linear_like") {
            pattern <- c(1, 2, 3)
        } else if (current_shape == "neg_linear") {
            pattern <- c(-1, -2, -3)
        } else if (current_shape == "u_shape") {
            pattern <- c(1.5, -3.0, 1.5)
        } else if (current_shape == "inv_u_shape") {
            pattern <- c(-1.5, 3.0, -1.5)
        } else if (current_shape == "s_shape") {
            pattern <- c(1, -1, 1)
        } else if (current_shape == "threshold") {
            pattern <- c(0, 0.5, 4.0)
        } else if (current_shape == "inv_threshold") {
            pattern <- c(0, -0.5, -4.0)
        } else {
            pattern <- rep(1, df_spline)
        }

        comp_beta <- b * pattern
        eta_components_raw[, i] <- as.vector(mat_spline_list[[i]] %*% comp_beta)
        baseline_components[i] <- as.vector(predict(mat_spline_list[[i]], newx = min_x) %*% comp_beta)

        # [精准匹配计算真理值]：使用同一套 basis_std_true 计算
        for (k in 2:q) {
            b_diff <- basis_std_true[k, ] - basis_std_true[1, ]
            true_eff_mat[i + 1, k - 1] <- sum(b_diff * comp_beta)
        }
    }

    # 合并 Overall 真理值
    for (k in 2:q) {
        true_eff_mat[1, k - 1] <- sum(true_eff_mat[2:(n_vars + 1), k - 1])
    }

    # -----------------------------------------------------------
    # 汇总计算最终偏效应 (Anchoring to 0) 并生成 Y
    # -----------------------------------------------------------
    eta_spline_raw <- rowSums(eta_components_raw)
    baseline_effect <- sum(baseline_components)

    eta_spline_adjusted <- eta_spline_raw - baseline_effect
    y_clean <- eta_spline_adjusted + cov_list$eta_cov

    y_observed <- add_noise_by_snr(as.vector(y_clean), snr_db = snr_db)

    cols_cov <- setdiff(names(cov_list$mm), "eta_cov")
    final_df <- cbind(y = y_observed, preds_scaled, cov_list$mm[, cols_cov, drop = FALSE])

    attr(final_df, "true_effect_mat") <- true_eff_mat
    attr(final_df, "spline_knots") <- global_knots # <--- 存起来
    attr(final_df, "spline_boundary") <- global_boundary # <--- 存起来

    return(as.data.frame(final_df))
}

#' @title 生成非线性剂量反应的二分类结局数据 (Generate Non-linear Binary Data)
#'
#' @description
#' 基于特定连接函数（Logit, Probit 或 Cloglog）生成非线性二分类结局数据，支持自动截距校准以控制罕见事件发生率。
#'
#' @details
#' \strong{罕见病与基线风险校准 (Auto-Balanced Incidence Rate):} \cr
#' 在病例对照研究模拟中，任意设置回归系数会导致合成队列中事件发生率（Incidence Rate）极端失衡（如全为 0 或全为 1）。
#' 如果提供了 \code{target_prop}（例如 0.05 代表 5\% 罕见病），算法将通过 \code{uniroot} 动态逆向求解截距，
#' 确保最终生成的二分类响应精确匹配目标疾病流行率。属性 \code{true_effect_mat} 中保存的即为真实的**对数比值比 (Log-OR)**。
#'
#' @param target_prop Numeric (0, 1) 或 \code{NULL}。目标疾病发生率。若提供，算法将自动搜索截距以平衡队列中病例的比例。
#' @param link Character。广义线性模型的连接函数，可选 \code{"logit"}（对应 OR）, \code{"probit"}, 或 \code{"cloglog"}。
#' @export
gen_nonlinear_bio_data <- function(n_obs = 1000, mu_preds, sigma_preds, beta_wqs = 1, beta_preds,
                                   intercept = 0, target_prop = NULL, link = c("logit", "probit", "cloglog"),
                                   snr_db = Inf, transform_fun = NULL, q = 4, df_spline = 3, seed = NULL,
                                   shape = "linear_like", ...) {
    if (!requireNamespace("splines", quietly = TRUE)) stop("Package 'splines' required")
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' required")

    link <- match.arg(link)
    if (!is.null(seed)) set.seed(seed)

    preds_raw <- MASS::mvrnorm(n_obs, mu = mu_preds, Sigma = sigma_preds)
    preds_scaled <- as.data.frame(scale(preds_raw))
    n_vars <- ncol(preds_scaled)
    names(preds_scaled) <- paste0("Component", 1:ncol(preds_scaled))

    if (!is.null(transform_fun) && is.function(transform_fun)) {
        preds_trans <- transform_fun(preds_scaled)
    } else {
        preds_trans <- preds_scaled
    }

    if (length(beta_preds) != n_vars) stop("Length of 'beta_preds' must match n_vars.")

    cov_list <- generate_covariates(n_obs = n_obs, ...)

    # -----------------------------------------------------------
    # [架构升级]: 动态推导并锁定节点 (保证数据生成和模型拟合的尺子绝对一致)
    # -----------------------------------------------------------
    eval_pts <- 0:(q - 1)
    temp_spline <- splines::ns(eval_pts, df = df_spline)
    global_knots <- attr(temp_spline, "knots")
    global_boundary <- attr(temp_spline, "Boundary.knots")

    # 数据生成侧：强制使用这把尺子
    mat_spline_list <- lapply(preds_trans, function(x) {
        splines::ns(x, df = df_spline, knots = global_knots, Boundary.knots = global_boundary)
    })

    # 真值计算侧：也强制使用这把尺子 (无截距版)
    basis_std_true <- splines::ns(eval_pts, df = df_spline, knots = global_knots, Boundary.knots = global_boundary, intercept = FALSE)

    if (length(shape) == 1) shape <- rep(shape, n_vars)

    eta_components_raw <- matrix(0, nrow = n_obs, ncol = n_vars)
    baseline_components <- numeric(n_vars)

    true_eff_mat <- matrix(0, nrow = n_vars + 1, ncol = q - 1)
    rownames(true_eff_mat) <- c("Overall", names(preds_scaled))
    colnames(true_eff_mat) <- paste0("Q", 2:q, "_vs_Q1")

    for (i in 1:n_vars) {
        current_shape <- shape[i]
        b <- beta_preds[i] * beta_wqs
        min_x <- min(preds_trans[, i])

        if (current_shape == "pure_linear") {
            # 纯粹线性效应 (Log-OR 尺度)
            eta_components_raw[, i] <- preds_trans[, i] * b
            baseline_components[i] <- min_x * b

            for (k in 2:q) {
                true_eff_mat[i + 1, k - 1] <- (eval_pts[k] - eval_pts[1]) * b
            }
        } else {
            # 样条非线性效应
            if (current_shape == "linear_like") {
                pattern <- c(1, 2, 3)
            } else if (current_shape == "neg_linear") {
                pattern <- c(-1, -2, -3)
            } else if (current_shape == "u_shape") {
                pattern <- c(1.5, -3.0, 1.5)
            } else if (current_shape == "inv_u_shape") {
                pattern <- c(-1.5, 3.0, -1.5)
            } else if (current_shape == "s_shape") {
                pattern <- c(1, -1, 1)
            } else if (current_shape == "threshold") {
                pattern <- c(0, 0.5, 4.0)
            } else if (current_shape == "inv_threshold") {
                pattern <- c(0, -0.5, -4.0)
            } else {
                pattern <- rep(1, df_spline)
            }

            comp_beta <- b * pattern

            eta_components_raw[, i] <- as.vector(mat_spline_list[[i]] %*% comp_beta)
            baseline_components[i] <- as.vector(predict(mat_spline_list[[i]], newx = min_x) %*% comp_beta)

            # [精准匹配计算真理值]：使用同一套 basis_std_true 计算
            for (k in 2:q) {
                b_diff <- basis_std_true[k, ] - basis_std_true[1, ]
                true_eff_mat[i + 1, k - 1] <- sum(b_diff * comp_beta)
            }
        }
    }

    for (k in 2:q) {
        true_eff_mat[1, k - 1] <- sum(true_eff_mat[2:(n_vars + 1), k - 1])
    }

    # -----------------------------------------------------------
    # 汇总计算最终偏效应 (Anchoring to 0) 并生成 Binary Y
    # -----------------------------------------------------------
    eta_spline_raw <- rowSums(eta_components_raw)
    baseline_effect <- sum(baseline_components)

    eta_spline_adjusted <- eta_spline_raw - baseline_effect
    eta_partial <- eta_spline_adjusted + cov_list$eta_cov

    if (!is.null(snr_db) && is.finite(snr_db)) {
        eta_noisy_partial <- add_noise_by_snr(eta_partial, snr_db)
    } else {
        eta_noisy_partial <- eta_partial
    }

    final_intercept <- intercept
    if (!is.null(target_prop)) {
        calc_mean_prob_diff <- function(b0) {
            eta_temp <- b0 + eta_noisy_partial
            if (link == "logit") {
                p <- 1 / (1 + exp(-eta_temp))
            } else if (link == "probit") {
                p <- pnorm(eta_temp)
            } else if (link == "cloglog") p <- 1 - exp(-exp(eta_temp))
            return(mean(p) - target_prop)
        }
        tryCatch(
            {
                final_intercept <- uniroot(calc_mean_prob_diff, interval = c(-50, 50))$root
            },
            error = function(e) {}
        )
    }

    eta_final <- final_intercept + eta_noisy_partial

    if (link == "logit") {
        probs <- 1 / (1 + exp(-eta_final))
    } else if (link == "probit") {
        probs <- pnorm(eta_final)
    } else if (link == "cloglog") probs <- 1 - exp(-exp(eta_final))

    y_binary <- rbinom(n_obs, size = 1, prob = probs)

    cols_cov <- setdiff(names(cov_list$mm), "eta_cov")
    final_df <- cbind(y = y_binary, preds_scaled, cov_list$mm[, cols_cov, drop = FALSE])

    attr(final_df, "true_effect_mat") <- true_eff_mat # Log-OR 尺度真值
    attr(final_df, "true_prob") <- probs
    attr(final_df, "spline_knots") <- global_knots # <--- 存起来
    attr(final_df, "spline_boundary") <- global_boundary # <--- 存起来
    return(as.data.frame(final_df))
}


#' @title 生成非线性暴露的计数结局数据 (Generate Non-linear Poisson Count Data)
#'
#' @description
#' 通过对数连接 (Log-link) 将混合物的非线性样条特征映射为事件发生的期望发生率 (\eqn{\lambda})，
#' 并从泊松过程 (Poisson Process) 中生成离散的计数结局。
#'
#' @details
#' \strong{截距与泊松基线风险:} \cr
#' 截距参数在对数尺度上直接控制基线事件发生率。例如 \code{intercept=0} 意味着协变量和暴露处于基准时，期望计数约为 1；
#' 而 \code{intercept=2} 则期望计数提升至 \eqn{e^2 \approx 7.4}。属性 \code{true_effect_mat} 提供的是纯正的**对数相对危险度 (Log-RR)**。
#'
#' @param intercept Numeric。模型截距（对数尺度），直接控制泊松过程的基础 \eqn{\lambda}。
#' @export
gen_nonlinear_count_data <- function(n_obs = 1000, mu_preds, sigma_preds, beta_wqs = 1, beta_preds,
                                     intercept = 0, snr_db = Inf, transform_fun = NULL, q = 4, df_spline = 3,
                                     seed = NULL, shape = "linear_like", ...) {
    if (!requireNamespace("splines", quietly = TRUE)) stop("Package 'splines' required")
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' required")
    if (!is.null(seed)) set.seed(seed)

    preds_raw <- MASS::mvrnorm(n_obs, mu = mu_preds, Sigma = sigma_preds)
    preds_scaled <- as.data.frame(scale(preds_raw))
    n_vars <- ncol(preds_scaled)
    names(preds_scaled) <- paste0("Component", 1:n_vars)

    if (!is.null(transform_fun) && is.function(transform_fun)) {
        preds_trans <- transform_fun(preds_scaled)
    } else {
        preds_trans <- preds_scaled
    }

    if (length(beta_preds) != n_vars) stop("Length of 'beta_preds' must match n_vars.")

    cov_list <- generate_covariates(n_obs = n_obs, ...)

    # -----------------------------------------------------------
    # [架构升级]: 动态推导并锁定节点 (保证数据生成和模型拟合的尺子绝对一致)
    # -----------------------------------------------------------
    eval_pts <- 0:(q - 1)
    temp_spline <- splines::ns(eval_pts, df = df_spline)
    global_knots <- attr(temp_spline, "knots")
    global_boundary <- attr(temp_spline, "Boundary.knots")

    # 数据生成侧：强制使用这把尺子
    mat_spline_list <- lapply(preds_trans, function(x) {
        splines::ns(x, df = df_spline, knots = global_knots, Boundary.knots = global_boundary)
    })

    # 真值计算侧：也强制使用这把尺子 (无截距版)
    basis_std_true <- splines::ns(eval_pts, df = df_spline, knots = global_knots, Boundary.knots = global_boundary, intercept = FALSE)

    if (length(shape) == 1) shape <- rep(shape, n_vars)

    eta_components_raw <- matrix(0, nrow = n_obs, ncol = n_vars)
    baseline_components <- numeric(n_vars)

    true_eff_mat <- matrix(0, nrow = n_vars + 1, ncol = q - 1)
    rownames(true_eff_mat) <- c("Overall", names(preds_scaled))
    colnames(true_eff_mat) <- paste0("Q", 2:q, "_vs_Q1")

    for (i in 1:n_vars) {
        current_shape <- shape[i]
        b <- beta_preds[i] * beta_wqs
        min_x <- min(preds_trans[, i])

        if (current_shape == "pure_linear") {
            # 纯粹线性效应 (Log-RR 尺度)
            eta_components_raw[, i] <- preds_trans[, i] * b
            baseline_components[i] <- min_x * b

            for (k in 2:q) {
                true_eff_mat[i + 1, k - 1] <- (eval_pts[k] - eval_pts[1]) * b
            }
        } else {
            # 样条非线性效应
            if (current_shape == "linear_like") {
                pattern <- c(1, 2, 3)
            } else if (current_shape == "neg_linear") {
                pattern <- c(-1, -2, -3)
            } else if (current_shape == "u_shape") {
                pattern <- c(1.5, -3.0, 1.5)
            } else if (current_shape == "inv_u_shape") {
                pattern <- c(-1.5, 3.0, -1.5)
            } else if (current_shape == "s_shape") {
                pattern <- c(1, -1, 1)
            } else if (current_shape == "threshold") {
                pattern <- c(0, 0.5, 4.0)
            } else if (current_shape == "inv_threshold") {
                pattern <- c(0, -0.5, -4.0)
            } else {
                pattern <- rep(1, df_spline)
            }

            comp_beta <- b * pattern

            eta_components_raw[, i] <- as.vector(mat_spline_list[[i]] %*% comp_beta)
            baseline_components[i] <- as.vector(predict(mat_spline_list[[i]], newx = min_x) %*% comp_beta)

            # [精准匹配计算真理值]：使用同一套 basis_std_true 计算
            for (k in 2:q) {
                b_diff <- basis_std_true[k, ] - basis_std_true[1, ]
                true_eff_mat[i + 1, k - 1] <- sum(b_diff * comp_beta)
            }
        }
    }

    for (k in 2:q) {
        true_eff_mat[1, k - 1] <- sum(true_eff_mat[2:(n_vars + 1), k - 1])
    }

    # -----------------------------------------------------------
    # 汇总计算最终偏效应 (Anchoring to 0) 并生成 Count Y
    # -----------------------------------------------------------
    eta_spline_raw <- rowSums(eta_components_raw)
    baseline_effect <- sum(baseline_components)

    eta_spline_adjusted <- eta_spline_raw - baseline_effect
    eta_partial <- eta_spline_adjusted + cov_list$eta_cov

    if (!is.null(snr_db) && is.finite(snr_db)) {
        eta_noisy_partial <- add_noise_by_snr(eta_partial, snr_db)
    } else {
        eta_noisy_partial <- eta_partial
    }

    eta_final <- intercept + eta_noisy_partial
    lambda <- exp(eta_final)

    if (any(lambda > 10000)) warning("Extremely high lambda values.")

    y_count <- rpois(n_obs, lambda = lambda)

    cols_cov <- setdiff(names(cov_list$mm), "eta_cov")
    final_df <- cbind(y = y_count, preds_scaled, cov_list$mm[, cols_cov, drop = FALSE])

    attr(final_df, "true_effect_mat") <- true_eff_mat # Log-RR 尺度真值
    attr(final_df, "spline_knots") <- global_knots # <--- 存起来
    attr(final_df, "spline_boundary") <- global_boundary # <--- 存起来
    return(as.data.frame(final_df))
}
