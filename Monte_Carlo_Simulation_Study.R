# =========================================================================
# NWQS 蒙特卡洛模拟研究: 全量基准对比引擎 (Mass Benchmark Engine)
# 包含: NWQS, gWQS, QGcomp, Ridge, Lasso, ElasticNet, RandomForest
# =========================================================================

# -------------------------------------------------------------------------
# 0. 全局配置与环境准备
# -------------------------------------------------------------------------
rm(list = ls())
devtools::load_all() 
library(dplyr)
library(tidyr)
library(ggplot2)
library(gWQS)
library(qgcomp)
library(glmnet)
library(randomForest)

# 【控制开关】 
TARGET_FAMILY <- "gaussian" 
N_SIMULATIONS <- 100  # 蒙特卡洛外层大循环次数
RH_SIM        <- 100  # NWQS/gWQS 的内层自举次数

base_out_dir <- file.path("results", "Monte_Carlo_Results", toupper(TARGET_FAMILY))
if (!dir.exists(base_out_dir)) dir.create(base_out_dir, recursive = TRUE)

message(sprintf("\n=================================================="))
message(sprintf(" 🚀 初始化全场景模拟流水线 | 分布族: [%s]", toupper(TARGET_FAMILY)))
message(sprintf("==================================================\n"))

# -------------------------------------------------------------------------
run_simulation_benchmark_engine <- function(family, scen_name, params) {
  options(future.globals.maxSize = 4000 * 1024^2)
  message(sprintf("\n▶ 正在执行: %s | 参数: N=%d, P=%d", scen_name, params$N, params$P))

  n_vars <- params$P
  mix_name <- paste0("Component", 1:n_vars)
  mu_preds <- rep(0, n_vars)
  sigma_preds <- generate_sigma(n_vars = n_vars, mode = params$corr, seed = 525)
  w_true <- params$w
  names(w_true) <- mix_name

  list_dev <- list()
  list_sae <- list()
  list_weights <- list()

  # 偏差计算器
  calc_dev <- function(y, pred, fam) {
    if (fam == "gaussian") return(sum((y - pred)^2))
    if (fam == "binomial") {
      p <- pmax(pmin(pred, 1 - 1e-7), 1e-7)
      return(-2 * sum(y * log(p) + (1 - y) * log(1 - p)))
    }
    if (fam %in% c("poisson", "quasipoisson")) {
      pred <- pmax(pred, 1e-7)
      return(2 * sum(ifelse(y == 0, 0, y * log(y / pred)) - (y - pred)))
    }
  }

  fam_obj_glmnet <- ifelse(family == "quasipoisson", "poisson", family)
  fam_obj_qgcomp <- if (family == "gaussian") gaussian() else if (family == "binomial") binomial() else poisson()

  start_time <- Sys.time()

  for (i in 1:N_SIMULATIONS) {
    if (i %% 10 == 0) message(sprintf("  -> 进度: %d / %d", i, N_SIMULATIONS))
    sim_seed <- 10000 + i
    transform_fun_sim <- function(x) trans_quantile(x, q = 4)

    # 1. 数据生成
    if (family == "gaussian") {
      sim_data <- gen_nonlinear_data(n_obs = params$N, mu_preds = mu_preds, sigma_preds = sigma_preds, beta_preds = w_true, beta_wqs = 3, snr_db = params$snr_db, transform_fun = transform_fun_sim, q = 4, df_spline = 3, shape = params$shape, seed = sim_seed)
    } else if (family == "binomial") {
      sim_data <- gen_nonlinear_bio_data(n_obs = params$N, mu_preds = mu_preds, sigma_preds = sigma_preds, beta_preds = w_true, beta_wqs = 1, target_prop = 0.3, link = "logit", snr_db = params$snr_db, transform_fun = transform_fun_sim, q = 4, df_spline = 3, shape = params$shape, seed = sim_seed)
    } else if (family == "quasipoisson") {
      sim_data <- gen_nonlinear_count_data(n_obs = params$N, mu_preds = mu_preds, sigma_preds = sigma_preds, beta_preds = w_true, beta_wqs = 1, intercept = 0, snr_db = params$snr_db, transform_fun = transform_fun_sim, q = 4, df_spline = 3, shape = params$shape, seed = sim_seed)
    }
    Y_target <- sim_data$y

    # 2. 拟合竞技场
    nwqs_fit <- nwqs(data = sim_data, mix_name = mix_name, covariates = c("x_cont", "x_bin", "x_cat"), dependent_var = "y", model_func = ridge_permutation_scorer, q = 4, split_prop = 0.6, seed = sim_seed, rh = RH_SIM, family = family, plan_strategy = "multicore", n_workers = 8)
    gwqs_fit <- suppressWarnings(gWQS::gwqs(formula = y ~ wqs + x_cont + x_bin + x_cat, data = sim_data, mix_name = mix_name, q = 4, validation = 0.6, b = 50, rh = RH_SIM, plan_strategy = "multicore", family = family, seed = sim_seed))
    
    data_qgcomp <- sim_data
    data_qgcomp[mix_name] <- transform_fun_sim(data_qgcomp[mix_name])
    qgcomp_fit <- qgcomp::qgcomp.noboot(f = y ~ ., expnms = mix_name, data = data_qgcomp[, c("y", mix_name, "x_cont", "x_bin", "x_cat")], family = fam_obj_qgcomp, q = NULL)

    df_spline <- 3
    temp_sp <- splines::ns(0:3, df = df_spline)
    sp_mat <- wqs_nonlinear_expand(data_qgcomp, mix_name, knots = attr(temp_sp, "knots"), boundary = attr(temp_sp, "Boundary.knots"))
    X_sp <- cbind(sp_mat, model.matrix(~ x_cont + x_bin + x_cat - 1, data = sim_data))

    cv_ridge <- glmnet::cv.glmnet(X_sp, Y_target, alpha = 0, family = fam_obj_glmnet)
    cv_lasso <- glmnet::cv.glmnet(X_sp, Y_target, alpha = 1, family = fam_obj_glmnet)
    cv_enet  <- glmnet::cv.glmnet(X_sp, Y_target, alpha = 0.5, family = fam_obj_glmnet)
    rf_fit <- randomForest::randomForest(x = data_qgcomp[, c(mix_name, "x_cont", "x_bin", "x_cat")], y = if (family == "binomial") as.factor(Y_target) else Y_target, importance = TRUE, ntree = 500)

    # 3. 结果权重提取与标准化 (✅ 修复了 w_ridge 缺失的问题)
    w_nwqs <- nwqs_fit$final_weights
    w_gwqs <- gwqs_fit$final_weights$Estimate
    names(w_gwqs) <- gwqs_fit$final_weights$mix_name
    
    w_qg <- qgcomp_fit$pos.weights
    missing_names <- setdiff(mix_name, names(w_qg))
    w_qg <- c(w_qg, setNames(rep(0, length(missing_names)), missing_names))[mix_name]

    # --- 新增 Ridge 权重计算 ---
    ridge_cf <- as.matrix(coef(cv_ridge, s = "lambda.min"))[-1, 1]
    w_ridge <- sapply(mix_name, function(c) sum(abs(ridge_cf[grep(paste0("^", c, "_B"), names(ridge_cf))])))
    w_ridge <- w_ridge / sum(w_ridge)
    # ---------------------------

    lasso_cf <- as.matrix(coef(cv_lasso, s = "lambda.min"))[-1, 1]
    w_lasso <- sapply(mix_name, function(c) sum(abs(lasso_cf[grep(paste0("^", c, "_B"), names(lasso_cf))])))
    w_lasso <- w_lasso / sum(w_lasso)

    enet_cf <- as.matrix(coef(cv_enet, s = "lambda.min"))[-1, 1]
    w_enet <- sapply(mix_name, function(c) sum(abs(enet_cf[grep(paste0("^", c, "_B"), names(enet_cf))])))
    w_enet <- w_enet / sum(w_enet)

    rf_imp <- randomForest::importance(rf_fit)[mix_name, if (family == "binomial") "MeanDecreaseGini" else "%IncMSE"]
    w_rf <- pmax(rf_imp, 0) / sum(pmax(rf_imp, 0))

    # 4. 填充数据
    list_dev[[i]] <- data.frame(
      Iteration  = i, NWQS = nwqs_fit$fit$deviance, gWQS = mean(gwqs_fit$fit$deviance), QG_Comp = calc_dev(Y_target, predict(qgcomp_fit), family),
      Ridge = calc_dev(Y_target, predict(cv_ridge, X_sp, s = "lambda.min", type = "response"), family),
      Lasso = calc_dev(Y_target, predict(cv_lasso, X_sp, s = "lambda.min", type = "response"), family),
      ElasticNet = calc_dev(Y_target, predict(cv_enet, X_sp, s = "lambda.min", type = "response"), family),
      RandomForest = calc_dev(Y_target, if (family == "binomial") predict(rf_fit, type = "prob")[, 2] else predict(rf_fit), family)
    ) %>% pivot_longer(-Iteration, names_to = "Model", values_to = "Deviance")

    list_sae[[i]] <- data.frame(
      Iteration  = i, NWQS = calc_weight_error(w_nwqs, w_true)$SAE, gWQS = calc_weight_error(w_gwqs, w_true)$SAE, QG_Comp = calc_weight_error(w_qg, w_true)$SAE,
      Ridge = calc_weight_error(w_ridge, w_true)$SAE, Lasso = calc_weight_error(w_lasso, w_true)$SAE, ElasticNet = calc_weight_error(w_enet, w_true)$SAE,
      RandomForest = calc_weight_error(w_rf, w_true)$SAE
    ) %>% pivot_longer(-Iteration, names_to = "Model", values_to = "SAE")

    list_weights[[i]] <- data.frame(
      Iteration  = i, Component = mix_name, True_Value = w_true, NWQS = w_nwqs[mix_name], gWQS = w_gwqs[mix_name], QG_Comp = w_qg[mix_name],
      Ridge = w_ridge[mix_name], Lasso = w_lasso[mix_name], ElasticNet = w_enet[mix_name], RandomForest = w_rf[mix_name]
    ) %>% pivot_longer(NWQS:RandomForest, names_to = "Model", values_to = "Estimated_Weight")

    rm(sim_data, nwqs_fit, gwqs_fit, qgcomp_fit, cv_ridge, cv_lasso, cv_enet, rf_fit); gc(verbose = FALSE)
  }

  # 5. 数据保存
  df_all_dev <- do.call(rbind, list_dev)
  df_all_sae <- do.call(rbind, list_sae)
  df_all_weights <- do.call(rbind, list_weights)

  base_test_content <- gsub("_N[0-9]+$", "", scen_name)
  final_out_dir <- file.path(base_out_dir, base_test_content)
  if (!dir.exists(final_out_dir)) dir.create(final_out_dir, recursive = TRUE)

  summary_all <- df_all_dev %>% group_by(Model) %>% summarize(Mean_Dev = mean(Deviance), SD_Dev = sd(Deviance)) %>%
    left_join(df_all_sae %>% group_by(Model) %>% summarize(Mean_SAE = mean(SAE)), by = "Model")
  write.csv(summary_all, file.path(final_out_dir, paste0("Summary_", scen_name, ".csv")), row.names = FALSE)

  img_base <- file.path(final_out_dir, sprintf("MC_Benchmark_%s_%s", family, scen_name))
  final_mc_plot <- plot_monte_carlo_benchmark(dev_data = df_all_dev, sae_data = df_all_sae, weight_data = df_all_weights, save_path = paste0(img_base, ".png"))
  
  dynamic_nrow <- ceiling(params$P / 7)
  ggplot2::ggsave(filename = paste0(img_base, ".pdf"), plot = final_mc_plot, width = 16, height = 11 + (dynamic_nrow - 1) * 3.5, device = "pdf")

  cat(sprintf("\n✅ 场景 %s 完成。耗时: %.2f mins\n", scen_name, as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
  return(invisible(list(plot = final_mc_plot, summary = summary_all)))
}



# -------------------------------------------------------------------------
# 3. 全量场景字典生成 (10大场景 x 3种样本量 = 30个组合)
# -------------------------------------------------------------------------
w_norm_4    <- c(0.10, 0.20, 0.30, 0.40)
w_sparse_4  <- c(0.60, 0.40, 0.00, 0.00)
w_norm_8    <- c(0.04, 0.06, 0.08, 0.10, 0.14, 0.16, 0.18, 0.24)
w_norm_12   <- c(0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.11, 0.12, 0.13, 0.13, 0.13)
w_sparse_12 <- c(0.40, 0.30, 0.20, 0.10, rep(0, 8))

base_settings <- list(
  "S1_Base"        = list(P = 4,  corr = "mixed", shape = "threshold", w = w_norm_4,   snr_db = 10),
  "S3_Corr_Low"    = list(P = 4,  corr = "low",   shape = "threshold", w = w_norm_4,   snr_db = 10),
  "S3_Corr_High"   = list(P = 4,  corr = "high",  shape = "threshold", w = w_norm_4,   snr_db = 10),
  "S4_Dim_P8"      = list(P = 8,  corr = "mixed", shape = "threshold", w = w_norm_8,   snr_db = 10),
  "S4_Dim_P12"     = list(P = 12, corr = "mixed", shape = "threshold", w = w_norm_12,  snr_db = 10),
  "S5_Sparse_P4"   = list(P = 4,  corr = "mixed", shape = "threshold", w = w_sparse_4, snr_db = 10),
  "S5_Sparse_P12"  = list(P = 12, corr = "high",  shape = "threshold", w = w_sparse_12,snr_db = 10),
  "S6_Hetero_P4"   = list(P = 4,  corr = "mixed", w = w_norm_4, snr_db = 10,
                          shape = c("threshold", "inv_threshold", "neg_linear", "pure_linear")),
  "S7_Linear_P4"   = list(P = 4,  corr = "mixed", shape = "pure_linear", w = w_norm_4, snr_db = 10),
  "S8_ShapeMix_P4" = list(P = 4,  corr = "mixed", w = w_norm_4, snr_db = 10,
                          shape = c("u_shape", "threshold", "pure_linear", "s_shape"))
)

N_values <- c(200, 500, 1000) 
scenarios <- list()

for (b_name in names(base_settings)) {
  for (n in N_values) {
    full_name <- sprintf("%s_N%d", b_name, n)
    curr <- base_settings[[b_name]]
    curr$N <- n
    scenarios[[full_name]] <- curr
  }
}

# -------------------------------------------------------------------------
# 4. 全自动化点火循环 (无人值守挂机版)
# -------------------------------------------------------------------------
total_tasks <- length(names(scenarios))
current_task <- 0

message(sprintf("\n🚀 引擎已就绪！即将执行 %d 个全量测试任务...", total_tasks))

for (scen_name in names(scenarios)) {
  current_task <- current_task + 1
  message(sprintf("\n=================================================="))
  message(sprintf("🚀 [AUTO-PILOT] 任务 %d/%d 启动场景: %s", current_task, total_tasks, scen_name))
  message(sprintf("=================================================="))
  
  tryCatch({
    run_simulation_benchmark_engine(
      family    = TARGET_FAMILY,
      scen_name = scen_name,
      params    = scenarios[[scen_name]]
    )
  }, error = function(e) {
    message(sprintf("❌ 出错跳过: %s | 错误信息: %s", scen_name, e$message))
  })
  
  gc(verbose = FALSE)
}

message("\n🎉 满配全流水线任务已收割完毕！请在 results/ 目录下查看所有 CSV 和 高清 PDF。")