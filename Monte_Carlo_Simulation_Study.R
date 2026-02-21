# =========================================================================
# 蒙特卡洛模拟研究 (Monte Carlo Simulation Study)
# 脚本名称: simulation.R
# 目标: 模拟 100 个独立数据集，全面覆盖 连续型、二分类、计数型 结局，
#       对比 NWQS 与 gWQS 在非线性(U型/阈值)场景下的无偏性与稳健性
# =========================================================================

# -------------------------------------------------------------------------
# 0. 环境准备与全局设置
# -------------------------------------------------------------------------
rm(list = ls())
devtools::load_all() 
library(gWQS)
library(dplyr)
library(tidyr)
library(ggplot2)

N_SIMULATIONS <- 10
RH_SIM <- 5  # 为了节约模拟时间，模拟研究中的 RH 通常设为 5-10 即可

# 设定基础参数 
n_vars <- 4
mix_name <- paste0("Component", 1:n_vars)
mu_preds <- rep(0, n_vars)
set.seed(525) 
sigma_preds <- generate_sigma(n_vars = n_vars, mode = "mixed", seed = 525)
beta_preds <- c(0.1, 0.2, 0.3, 0.4)
w_true <- c(Component1=0.1, Component2=0.2, Component3=0.3, Component4=0.4)
beta_wqs <- 3
transform_fun <- function(x) trans_quantile(x, q = 4)

# -------------------------------------------------------------------------
# 1. 结果存储容器初始化
# -------------------------------------------------------------------------
results_gauss <- list()
results_bin   <- list()
results_count <- list()
temp_save_file <- "simulation_temp_results.rds"

message(sprintf("=================================================="))
message(sprintf("  开始 Monte Carlo 全场景模拟研究 (总计 %d 次)", N_SIMULATIONS))
message(sprintf("  包含: Gaussian, Binomial, Quasi-Poisson"))
message(sprintf("=================================================="))

start_time <- Sys.time()

# -------------------------------------------------------------------------
# 2. 开启 100 次的模拟主循环
# -------------------------------------------------------------------------
for (i in 1:N_SIMULATIONS) {
  message(sprintf("\n>>> [主循环] 正在运行模拟迭代: %d / %d ...", i, N_SIMULATIONS))
  sim_seed <- 10000 + i  # 保证每次循环数据完全独立
  
  # =========================================================
  # 模块 A: 连续型 (Gaussian) - 使用 U 型效应
  # =========================================================
  message("    -> 正在测试 Gaussian...")
  sim_data_gauss <- gen_nonlinear_data(
    n_obs = 1000, mu_preds = mu_preds, sigma_preds = sigma_preds, 
    beta_preds = beta_preds, beta_wqs = beta_wqs, snr_db = 10, 
    transform_fun = transform_fun, df_spline = 3, shape = "u_shape", seed = sim_seed
  )
  
  nwqs_gauss <- nwqs(
    data = sim_data_gauss, mix_name = mix_name, covariates = c("x_cont", "x_bin", "x_cat"),
    dependent_var = "y", model_func = calc_spline_wqs_weights, q = 4, split_prop = 0.6, 
    seed = sim_seed, rh = RH_SIM, family = "gaussian", transform_fun = transform_fun, plan_strategy = "multicore", n_workers = 8
  )
  gwqs_gauss <- gwqs(
    formula = y ~ wqs + x_cont + x_bin + x_cat, data = sim_data_gauss, mix_name = mix_name,
    y = "y", q = 4, validation = 0.6, b = 100, rh = RH_SIM, plan_strategy = "multicore", family = "gaussian", seed = sim_seed
  )
  
  w_gwqs_g <- gwqs_gauss$final_weights$Estimate; names(w_gwqs_g) <- gwqs_gauss$final_weights$mix_name
  res_gauss_df <- data.frame(
    Iteration = i, 
    NWQS_AIC = nwqs_gauss$mean_aic, gWQS_AIC = mean(gwqs_gauss$fit$aic),
    NWQS_Dev = nwqs_gauss$mean_res_dev, gWQS_Dev = mean(gwqs_gauss$fit$deviance),
    NWQS_SAE = calc_weight_error(nwqs_gauss$final_weights, w_true)$SAE,
    gWQS_SAE = calc_weight_error(w_gwqs_g, w_true)$SAE
  )
  w_nwqs_df_g <- as.data.frame(t(nwqs_gauss$final_weights)); colnames(w_nwqs_df_g) <- paste0("NWQS_", colnames(w_nwqs_df_g))
  w_gwqs_df_g <- as.data.frame(t(w_gwqs_g[names(w_true)])); colnames(w_gwqs_df_g) <- paste0("gWQS_", colnames(w_gwqs_df_g))
  results_gauss[[i]] <- cbind(res_gauss_df, w_nwqs_df_g, w_gwqs_df_g)

  # =========================================================
  # 模块 B: 二分类 (Binomial) - 使用 阈值 (Threshold) 效应
  # =========================================================
  message("    -> 正在测试 Binomial...")
  sim_data_bin <- gen_nonlinear_bio_data(
    n_obs = 1000, mu_preds = mu_preds, sigma_preds = sigma_preds, 
    beta_preds = beta_preds, beta_wqs = 1, target_prop = 0.3, link = "logit",
    snr_db = 10, transform_fun = transform_fun, df_spline = 3, shape = "threshold", seed = sim_seed
  )
  
  nwqs_bin <- nwqs(
    data = sim_data_bin, mix_name = mix_name, covariates = c("x_cont", "x_bin", "x_cat"),
    dependent_var = "y", model_func = calc_spline_wqs_weights, q = 4, split_prop = 0.6, 
    seed = sim_seed, rh = RH_SIM, family = "binomial", transform_fun = transform_fun, plan_strategy = "multicore", n_workers = 8
  )
  gwqs_bin <- gwqs(
    formula = y ~ wqs + x_cont + x_bin + x_cat, data = sim_data_bin, mix_name = mix_name,
    y = "y", q = 4, validation = 0.6, b = 100, rh = RH_SIM, plan_strategy = "multicore", family = "binomial", seed = sim_seed
  )
  
  w_gwqs_b <- gwqs_bin$final_weights$Estimate; names(w_gwqs_b) <- gwqs_bin$final_weights$mix_name
  res_bin_df <- data.frame(
    Iteration = i, 
    NWQS_AIC = nwqs_bin$mean_aic, gWQS_AIC = mean(gwqs_bin$fit$aic),
    NWQS_Dev = nwqs_bin$mean_res_dev, gWQS_Dev = mean(gwqs_bin$fit$deviance),
    NWQS_SAE = calc_weight_error(nwqs_bin$final_weights, w_true)$SAE,
    gWQS_SAE = calc_weight_error(w_gwqs_b, w_true)$SAE
  )
  w_nwqs_df_b <- as.data.frame(t(nwqs_bin$final_weights)); colnames(w_nwqs_df_b) <- paste0("NWQS_", colnames(w_nwqs_df_b))
  w_gwqs_df_b <- as.data.frame(t(w_gwqs_b[names(w_true)])); colnames(w_gwqs_df_b) <- paste0("gWQS_", colnames(w_gwqs_df_b))
  results_bin[[i]] <- cbind(res_bin_df, w_nwqs_df_b, w_gwqs_df_b)

  # =========================================================
  # 模块 C: 计数型 (Quasi-Poisson) - 使用 阈值 (Threshold) 效应
  # =========================================================
  message("    -> 正在测试 Quasi-Poisson...")
  sim_data_count <- gen_nonlinear_count_data(
    n_obs = 1000, mu_preds = mu_preds, sigma_preds = sigma_preds, 
    beta_preds = beta_preds, beta_wqs = 1, intercept = 0, snr_db = 10, 
    transform_fun = transform_fun, df_spline = 3, shape = "threshold", seed = sim_seed
  )
  
  nwqs_count <- nwqs(
    data = sim_data_count, mix_name = mix_name, covariates = c("x_cont", "x_bin", "x_cat"),
    dependent_var = "y", model_func = calc_spline_wqs_weights, q = 4, split_prop = 0.6, 
    seed = sim_seed, rh = RH_SIM, family = "quasipoisson", transform_fun = transform_fun, plan_strategy = "multicore", n_workers = 8
  )
  gwqs_count <- gwqs(
    formula = y ~ wqs + x_cont + x_bin + x_cat, data = sim_data_count, mix_name = mix_name,
    y = "y", q = 4, validation = 0.6, b = 100, rh = RH_SIM, plan_strategy = "multicore", family = "quasipoisson", seed = sim_seed
  )
  
  w_gwqs_c <- gwqs_count$final_weights$Estimate; names(w_gwqs_c) <- gwqs_count$final_weights$mix_name
  res_count_df <- data.frame(
    Iteration = i, 
    NWQS_AIC = NA, gWQS_AIC = NA, # Quasi-Poisson 没有真实的 AIC
    NWQS_Dev = nwqs_count$mean_res_dev, gWQS_Dev = mean(gwqs_count$fit$deviance),
    NWQS_SAE = calc_weight_error(nwqs_count$final_weights, w_true)$SAE,
    gWQS_SAE = calc_weight_error(w_gwqs_c, w_true)$SAE
  )
  w_nwqs_df_c <- as.data.frame(t(nwqs_count$final_weights)); colnames(w_nwqs_df_c) <- paste0("NWQS_", colnames(w_nwqs_df_c))
  w_gwqs_df_c <- as.data.frame(t(w_gwqs_c[names(w_true)])); colnames(w_gwqs_df_c) <- paste0("gWQS_", colnames(w_gwqs_df_c))
  results_count[[i]] <- cbind(res_count_df, w_nwqs_df_c, w_gwqs_df_c)

  # =========================================================
  # 安全机制：每 10 次存盘
  # =========================================================
  if (i %% 10 == 0) {
    saveRDS(list(Gaussian = results_gauss, Binomial = results_bin, Poisson = results_count), file = temp_save_file)
    message(sprintf("   [存盘] 已完成 %d 次，全场景中间结果已安全备份。", i))
  }
}

future::plan(future::sequential)
end_time <- Sys.time()
message(sprintf("\n=================================================="))
message(sprintf("  Monte Carlo 模拟完成！总耗时: %.2f mins", as.numeric(difftime(end_time, start_time, units="mins"))))
message(sprintf("=================================================="))


# -------------------------------------------------------------------------
# 3. 汇总性能指标并绘制 100 次权重分配的并排箱线图
# -------------------------------------------------------------------------

# 合并三个 DataFrame
df_gauss <- do.call(rbind, results_gauss)
df_bin   <- do.call(rbind, results_bin)
df_count <- do.call(rbind, results_count)

# 保存最终数据
saveRDS(list(Gaussian = df_gauss, Binomial = df_bin, Poisson = df_count), file = "Final_Simulation_100_Runs.rds")

# 提取并打印整体性能差异表格
print_summary <- function(df, name) {
  stats <- df %>% summarise(
    NWQS_Dev = mean(NWQS_Dev), gWQS_Dev = mean(gWQS_Dev),
    NWQS_SAE = mean(NWQS_SAE), gWQS_SAE = mean(gWQS_SAE)
  )
  cat(sprintf("\n--- %s 结局 100次平均表现 ---\n", name))
  cat(sprintf("Deviance: NWQS = %.2f, gWQS = %.2f\n", stats$NWQS_Dev, stats$gWQS_Dev))
  cat(sprintf("Weight SAE: NWQS = %.4f, gWQS = %.4f\n", stats$NWQS_SAE, stats$gWQS_SAE))
}

print_summary(df_gauss, "Gaussian")
print_summary(df_bin, "Binomial")
print_summary(df_count, "Quasi-Poisson")


# -------------------------------------------------------------------------
# 4. 自动化制图核心函数：绘制 100 次模拟的权重箱线图 (降维打击审稿人)
# -------------------------------------------------------------------------
plot_simulation_weights <- function(sim_df, title, true_w = w_true) {
  
  # 提取 NWQS 权重
  nwqs_w <- sim_df %>% select(Iteration, starts_with("NWQS_Component")) %>%
    pivot_longer(-Iteration, names_to = "Component", values_to = "Weight") %>%
    mutate(Model = "NWQS", Component = gsub("NWQS_", "", Component))
  
  # 提取 gWQS 权重
  gwqs_w <- sim_df %>% select(Iteration, starts_with("gWQS_Component")) %>%
    pivot_longer(-Iteration, names_to = "Component", values_to = "Weight") %>%
    mutate(Model = "gWQS", Component = gsub("gWQS_", "", Component))
  
  # 合并绘图数据
  plot_df <- bind_rows(nwqs_w, gwqs_w)
  
  # 生成真实权重的参考数据
  true_df <- data.frame(
    Component = names(true_w),
    True_Weight = as.numeric(true_w)
  )
  
  # 画图
  p <- ggplot(plot_df, aes(x = Component, y = Weight, fill = Model)) +
    geom_boxplot(alpha = 0.8, outlier.size = 0.5, position = position_dodge(0.8)) +
    # 添加红色虚线段代表各个成分真实的客观权重
    geom_segment(data = true_df, aes(x = as.numeric(as.factor(Component)) - 0.4, 
                                     xend = as.numeric(as.factor(Component)) + 0.4, 
                                     y = True_Weight, yend = True_Weight),
                 color = "red", linetype = "dashed", linewidth = 1, inherit.aes = FALSE) +
    scale_fill_manual(values = c("NWQS" = "#E74C3C", "gWQS" = "#95A5A6")) +
    theme_bw(base_size = 14) +
    theme(legend.position = "top", legend.title = element_blank()) +
    labs(title = title,
         subtitle = "Red dashed lines indicate True Weight Generation Parameters",
         y = "Estimated Weight across 100 Simulations", x = "Mixture Components")
  
  return(p)
}

# 批量生成三张图并保存
p_w_gauss <- plot_simulation_weights(df_gauss, "Weight Distribution (100 Runs): Gaussian Outcome")
p_w_bin   <- plot_simulation_weights(df_bin,   "Weight Distribution (100 Runs): Binomial Outcome")
p_w_count <- plot_simulation_weights(df_count, "Weight Distribution (100 Runs): Quasi-Poisson Outcome")

ggsave("Sim_Weights_Gaussian.png", plot = p_w_gauss, width = 10, height = 6, dpi = 300)
ggsave("Sim_Weights_Binomial.png", plot = p_w_bin, width = 10, height = 6, dpi = 300)
ggsave("Sim_Weights_Poisson.png", plot = p_w_count, width = 10, height = 6, dpi = 300)

message("\n>>> 所有模拟结果已汇总，三张权重稳定性箱线图已自动保存至本地！")