# =========================================================================
# NWQS 超参数与稳定性敏感性分析 (终极全自动流水线版)
# 特性: 智能OS识别、自动依赖加载、完美色彩映射、NWQS智能负载均衡并发
# =========================================================================

# -------------------------------------------------------------------------
# 0. 环境准备与全自动包管理
# -------------------------------------------------------------------------
rm(list = ls())

message("\n📦 正在检查并加载所需的 R 包...")
required_packages <- c("dplyr", "tidyr", "ggplot2", "patchwork", "future")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("   [-] 正在安装缺失的包: %s", pkg))
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}
message("✅ 所有依赖包加载完毕！")

# 加载本地开发的 NWQS 包
# 如果你是通过 RStudio Build 运行，可注释掉这行；若是独立脚本则保留
devtools::load_all() 

# -------------------------------------------------------------------------
# 1. 智能 OS 嗅探与并行策略动态配置
# -------------------------------------------------------------------------
os_type <- .Platform$OS.type

if (os_type == "unix") {
  auto_plan_strategy <- "multicore"
  message("🖥️  系统识别: Unix/macOS/Linux")
} else {
  auto_plan_strategy <- "multisession"
  message("🖥️  系统识别: Windows")
}

# 将核心数设为 NULL，彻底激活你包内 configure_parallel_plan 的“智能负载均衡”
auto_workers <- NULL 
message(sprintf("⚡ 并行策略已锁定 [%s]。核心数分配交由 NWQS 智能负载均衡引擎处理！", auto_plan_strategy))

# -------------------------------------------------------------------------
# 2. 全局开关与参数网格生成
# -------------------------------------------------------------------------
TARGET_FAMILY <- "gaussian"
N_SIMULATIONS <- 2 # ⚠️ 首次试跑设为 2，正式跑图请设为 100

base_out_dir <- file.path("results", "Sensitivity_Analysis", toupper(TARGET_FAMILY))
if (!dir.exists(base_out_dir)) dir.create(base_out_dir, recursive = TRUE)

# 定义 P 对应的真实权重字典
w_true_list <- list(
  "4"  = c(0.10, 0.20, 0.30, 0.40),
  "8"  = c(0.04, 0.06, 0.08, 0.10, 0.14, 0.16, 0.18, 0.24),
  "12" = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.11, 0.12, 0.13, 0.13, 0.13)
)

# 生成 5 x 4 x 3 = 60 种组合的网格
param_grid <- expand.grid(
  rh = c(10, 30, 50, 100, 200),
  N  = c(200, 500, 1000, 2000),
  P  = c(4, 8, 12)
)

# -------------------------------------------------------------------------
# 3. 核心单次测试引擎
# -------------------------------------------------------------------------
run_nwqs_sensitivity <- function(N_val, P_val, rh_val, sim_rounds, plan_strat, n_work) {
  options(future.globals.maxSize = 4000 * 1024^2)
  options(future.rng.onMisuse = "ignore") # 忽略并行的种子警告
  
  mix_name <- paste0("Component", 1:P_val)
  mu_preds <- rep(0, P_val)
  sigma_preds <- generate_sigma(n_vars = P_val, mode = "mixed", seed = 525)
  w_true <- w_true_list[[as.character(P_val)]]
  names(w_true) <- mix_name

  results_list <- list()

  for (i in 1:sim_rounds) {
    current_iter_seed <- as.integer(10000 + i + (N_val + P_val + rh_val)) 
    set.seed(current_iter_seed)
    
    sim_data <- gen_nonlinear_data(
      n_obs = N_val, mu_preds = mu_preds, sigma_preds = sigma_preds,
      beta_preds = w_true, beta_wqs = 3, snr_db = 10,
      transform_fun = function(x) trans_quantile(x, q = 4), 
      q = 4, df_spline = 3, shape = "threshold", seed = current_iter_seed
    )

    # 🚀 n_workers 传入 NULL，触发底层 configure_parallel_plan()
    nwqs_fit <- nwqs(
      data = sim_data, mix_name = mix_name, covariates = c("x_cont", "x_bin", "x_cat"),
      dependent_var = "y", model_func = ridge_permutation_scorer,
      q = 4, split_prop = 0.6, seed = current_iter_seed, 
      rh = rh_val, # 这里 rh_val 也就是你的 loop_number
      family = TARGET_FAMILY, 
      plan_strategy = plan_strat, 
      n_workers = n_work, # NULL 
      shuffle = 50
    )

    sae_val <- calc_weight_error(nwqs_fit$final_weights, w_true)$SAE
    dev_val <- nwqs_fit$fit$deviance

    results_list[[i]] <- data.frame(
      N = N_val, P = P_val, rh = rh_val, Iteration = i, SAE = sae_val, Deviance = dev_val
    )
  }
  
  return(do.call(rbind, results_list))
}

# -------------------------------------------------------------------------
# 4. 全量自动化执行流水线
# -------------------------------------------------------------------------
total_tasks <- nrow(param_grid)
all_results <- list()
start_time_global <- Sys.time()

message(sprintf("\n🚀 引擎已就绪！即将执行 %d 个全量敏感性测试任务...", total_tasks))

for (task_idx in 1:total_tasks) {
  curr_N  <- param_grid$N[task_idx]
  curr_P  <- param_grid$P[task_idx]
  curr_rh <- param_grid$rh[task_idx]
  
  message(sprintf("▶ [%d/%d] 正在测算: N=%-4d | P=%-2d | rh=%-3d", 
                  task_idx, total_tasks, curr_N, curr_P, curr_rh))
  
  tryCatch({
    res_df <- run_nwqs_sensitivity(
      N_val = curr_N, P_val = curr_P, rh_val = curr_rh, 
      sim_rounds = N_SIMULATIONS, 
      plan_strat = auto_plan_strategy, 
      n_work = auto_workers # 传入 NULL 激活智能算法
    )
    all_results[[task_idx]] <- res_df
  }, error = function(e) {
    message(sprintf("❌ 出错跳过! 错误: %s", e$message))
  })
  
  gc(verbose = FALSE)
}

final_df <- do.call(rbind, all_results)

# 保存数据
write.csv(final_df, file.path(base_out_dir, "Sensitivity_Raw_Details.csv"), row.names = FALSE)
summary_df <- final_df %>%
  group_by(N, P, rh) %>%
  summarize(
    Mean_SAE = mean(SAE, na.rm = TRUE), SD_SAE = sd(SAE, na.rm = TRUE),
    Mean_Dev = mean(Deviance, na.rm = TRUE), .groups = 'drop'
  )
write.csv(summary_df, file.path(base_out_dir, "Sensitivity_Summary.csv"), row.names = FALSE)

message(sprintf("\n✅ 模拟计算完毕！总耗时: %.2f mins\n", 
                as.numeric(difftime(Sys.time(), start_time_global, units = "mins"))))

# -------------------------------------------------------------------------
# 5. 三维联动敏感性绘图引擎 (精准色彩映射版)
# -------------------------------------------------------------------------
plot_sensitivity_split <- function(raw_df, out_dir) {
  message("🎨 正在生成学术级精美报表...")
  
  pal_3 <- c("#4A90C8", "#D92828", "#6EC44A") 
  pal_4 <- c("#8B6FB8", "#00B4D8", "#006B3C", "#A8D8EA") 
  pal_5 <- c("#F4B6B6", "#5BA3D0", "#E03030", "#7AD450", "#9B7FC0") 

  plot_df <- raw_df %>%
    mutate(
      N_factor  = factor(N, levels = c(200, 500, 1000, 2000), labels = paste0("N = ", c(200, 500, 1000, 2000))),
      P_factor  = factor(P, levels = c(4, 8, 12), labels = paste0("P = ", c(4, 8, 12))),
      rh_factor = factor(rh, levels = c(10, 30, 50, 100, 200), labels = paste0("rh = ", c(10, 30, 50, 100, 200))),
      MSE_Dev   = Deviance / N 
    )

  base_theme <- theme_bw(base_size = 13) +
    theme(
      legend.position = "bottom", strip.background = element_rect(fill = "#2C3E50"),
      strip.text = element_text(color = "white", face = "bold", size = 12), panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 14), axis.title = element_text(face = "bold")
    )
  dodge_width <- 0.75

  # 图 1：rh
  p1_sae <- ggplot(plot_df, aes(x = rh_factor, y = SAE, fill = N_factor, color = N_factor)) +
    geom_boxplot(position = position_dodge(width = dodge_width), alpha = 0.5, outlier.size = 0.5) +
    stat_summary(fun = mean, geom = "line", aes(group = N_factor), position = position_dodge(width = dodge_width), linewidth = 0.8) +
    stat_summary(fun = mean, geom = "point", aes(group = N_factor), position = position_dodge(width = dodge_width), size = 2, color = "black", shape = 21, fill = "white") +
    facet_wrap(~P_factor) + scale_fill_manual(values = pal_4) + scale_color_manual(values = pal_4) +
    labs(title = "A. Impact of Bootstrap Iterations on Weight Extraction Error", x = "rh Iterations", y = "SAE", fill = "Sample Size (N)", color = "Sample Size (N)") + base_theme
  
  p1_dev <- ggplot(plot_df, aes(x = rh_factor, y = MSE_Dev, fill = N_factor, color = N_factor)) +
    geom_boxplot(position = position_dodge(width = dodge_width), alpha = 0.5, outlier.size = 0.5) +
    stat_summary(fun = mean, geom = "line", aes(group = N_factor), position = position_dodge(width = dodge_width), linewidth = 0.8) +
    stat_summary(fun = mean, geom = "point", aes(group = N_factor), position = position_dodge(width = dodge_width), size = 2, color = "black", shape = 21, fill = "white") +
    facet_wrap(~P_factor, scales = "free_y") + scale_fill_manual(values = pal_4) + scale_color_manual(values = pal_4) +
    labs(title = "B. Impact of Bootstrap Iterations on Model Fit", x = "rh Iterations", y = "Mean Squared Error (Dev/N)", fill = "Sample Size (N)", color = "Sample Size (N)") + base_theme
  plot_rh <- (p1_sae / p1_dev) + patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom")

  # 图 2：N
  p2_sae <- ggplot(plot_df, aes(x = N_factor, y = SAE, fill = rh_factor, color = rh_factor)) +
    geom_boxplot(position = position_dodge(width = dodge_width), alpha = 0.5, outlier.size = 0.5) +
    stat_summary(fun = mean, geom = "line", aes(group = rh_factor), position = position_dodge(width = dodge_width), linewidth = 0.8) +
    stat_summary(fun = mean, geom = "point", aes(group = rh_factor), position = position_dodge(width = dodge_width), size = 2, color = "black", shape = 21, fill = "white") +
    facet_wrap(~P_factor) + scale_fill_manual(values = pal_5) + scale_color_manual(values = pal_5) +
    labs(title = "A. Impact of Sample Size on Weight Extraction Error", x = "Sample Size (N)", y = "SAE", fill = "rh Iterations", color = "rh Iterations") + base_theme
  
  p2_dev <- ggplot(plot_df, aes(x = N_factor, y = MSE_Dev, fill = rh_factor, color = rh_factor)) +
    geom_boxplot(position = position_dodge(width = dodge_width), alpha = 0.5, outlier.size = 0.5) +
    stat_summary(fun = mean, geom = "line", aes(group = rh_factor), position = position_dodge(width = dodge_width), linewidth = 0.8) +
    stat_summary(fun = mean, geom = "point", aes(group = rh_factor), position = position_dodge(width = dodge_width), size = 2, color = "black", shape = 21, fill = "white") +
    facet_wrap(~P_factor, scales = "free_y") + scale_fill_manual(values = pal_5) + scale_color_manual(values = pal_5) +
    labs(title = "B. Impact of Sample Size on Model Fit", x = "Sample Size (N)", y = "Mean Squared Error (Dev/N)", fill = "rh Iterations", color = "rh Iterations") + base_theme
  plot_N <- (p2_sae / p2_dev) + patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom")

  # 图 3：P
  p3_sae <- ggplot(plot_df, aes(x = N_factor, y = SAE, fill = P_factor, color = P_factor)) +
    geom_boxplot(position = position_dodge(width = dodge_width), alpha = 0.5, outlier.size = 0.5) +
    stat_summary(fun = mean, geom = "line", aes(group = P_factor), position = position_dodge(width = dodge_width), linewidth = 0.8) +
    stat_summary(fun = mean, geom = "point", aes(group = P_factor), position = position_dodge(width = dodge_width), size = 2, color = "black", shape = 21, fill = "white") +
    facet_wrap(~rh_factor, nrow = 1) + scale_fill_manual(values = pal_3) + scale_color_manual(values = pal_3) +
    labs(title = "A. Impact of Dimensionality on Weight Extraction Error", x = "Sample Size (N)", y = "SAE", fill = "Dimensionality (P)", color = "Dimensionality (P)") + base_theme
  
  p3_dev <- ggplot(plot_df, aes(x = N_factor, y = MSE_Dev, fill = P_factor, color = P_factor)) +
    geom_boxplot(position = position_dodge(width = dodge_width), alpha = 0.5, outlier.size = 0.5) +
    stat_summary(fun = mean, geom = "line", aes(group = P_factor), position = position_dodge(width = dodge_width), linewidth = 0.8) +
    stat_summary(fun = mean, geom = "point", aes(group = P_factor), position = position_dodge(width = dodge_width), size = 2, color = "black", shape = 21, fill = "white") +
    facet_wrap(~rh_factor, scales = "free_y", nrow = 1) + scale_fill_manual(values = pal_3) + scale_color_manual(values = pal_3) +
    labs(title = "B. Impact of Dimensionality on Model Fit", x = "Sample Size (N)", y = "Mean Squared Error (Dev/N)", fill = "Dimensionality (P)", color = "Dimensionality (P)") + base_theme
  plot_P <- (p3_sae / p3_dev) + patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom")

  ggsave(file.path(out_dir, "Sensitivity_Fig1_Impact_of_RH.pdf"), plot = plot_rh, width = 12, height = 9, device = "pdf")
  ggsave(file.path(out_dir, "Sensitivity_Fig2_Impact_of_N.pdf"), plot = plot_N, width = 12, height = 9, device = "pdf")
  ggsave(file.path(out_dir, "Sensitivity_Fig3_Impact_of_P.pdf"), plot = plot_P, width = 16, height = 9, device = "pdf") 
  
  message("✅ 所有出图流程完毕！请前往 results 文件夹查看结果。")
  return(list(plot_rh = plot_rh, plot_N = plot_N, plot_P = plot_P))
}

# -------------------------------------------------------------------------
# 6. 一键点火出图
# -------------------------------------------------------------------------
final_plots <- plot_sensitivity_split(final_df, base_out_dir)