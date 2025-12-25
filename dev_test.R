# dev_test.R

# 1. 清理环境 (可选，保持干净)
rm(list = ls())

# 2. 核心步骤：加载当前正在开发的包
# 这相当于 library(NWQS)，但它加载的是你本地 R/ 文件夹里的最新代码
# 每次修改完代码后，只需要重新运行这一行即可，无需 install
devtools::load_all()

# -------------------------------------------------------------------------
# 下面是你的测试流程
# -------------------------------------------------------------------------

# 1. 设置参数
n_vars = 4
mix_name = paste0("Component", 1:n_vars)
mu_preds = rep(0, n_vars)
sigma_preds = generate_sigma(n_vars = n_vars, mode = "mixed", seed = 525)
beta_preds = c(0.1, 0.2, 0.3, 0.4)
beta_wqs = 3
transform_fun = function(x) trans_quantile(x, q = 4)

# 2. 生成模拟数据
data = gen_nonlinear_data(n_obs = 1000, 
                            mu_preds = mu_preds, 
                            sigma_preds = sigma_preds, 
                            beta_preds = beta_preds,
                            beta_wqs = beta_wqs, 
                            snr_db = 10, 
                            transform_fun = transform_fun,
                            df_spline = 3, 
                            seed = 525)


# 3. 运行 NWQS 模型
message("正在运行 NWQS 模型...")
nwqs_result = nwqs(data = data, 
                    mix_name = mix_name,
                    covariates = c("x_cont", "x_bin", "x_cat"),
                    dependent_var = "y",
                    model_func = calc_spline_wqs_weights,
                    q = 4,
                    split_prop = 0.6,
                    seed = 525,
                    rh = 100,
                    family = "gaussian",
                    transform_fun = transform_fun,
                    plan_strategy = "multicore")
nwqs_result

# 5. 简单验证
# 检查权重是否计算出来了
print(nwqs_result$final_weights)


# TODO:测试误差
# 检查是否和真实值有相关性
if (!is.null(nwqs_result$final_weights)) {
    diff_sum <- mean(abs(nwqs_result$final_weights - beta_preds))
    cat("与真实权重的平均绝对差异:", diff_sum, "\n")

    beta_wqs_diff = mean(nwqs_result$rh_coefs[, "wqs_score"] - beta_wqs)
    beta_wqs_diff_mse = mean((nwqs_result$rh_coefs[, "wqs_score"] - beta_wqs)^2)
    cat("与真实 beta_wqs 的平均差异:", beta_wqs_diff, "\n")
    cat("与真实 beta_wqs 的均方误差 (MSE):", beta_wqs_diff_mse, "\n")
}



sum(nwqs_result$rh_coefs[, "wqs_score"] > beta_wqs)
sum(nwqs_result$rh_coefs[, "(Intercept)"] > 0)
mean(nwqs_result$rh_coefs[, "(Intercept)"])





library(gWQS)

a = gwqs(formula = y ~ wqs + x_cont + x_bin + x_cat, data = data,
      mix_name = mix_name,
      y = "y",
      q = 4,
      validation = 0.6,
      b = 100,
      rh = 100,
      plan_strategy = "multicore",
      family = "gaussian",
      seed = 525)


# 进阶测试法(如何去测试某个函数或者脚本)
usethis::use_testthat()
usethis::use_test("nwqs.R")
