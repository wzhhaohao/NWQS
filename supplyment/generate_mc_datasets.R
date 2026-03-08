# # ============================================================================
# # generate_mc_datasets.R
# # ---------------------------------------------------------------------------
# # 目的:
# #   生成 Monte Carlo 模拟数据库，并将每套模拟数据保存为 .rds 文件
# #
# # 输出目录结构:
# # results/
# #   Monte_Carlo_DB/
# #     GAUSSIAN/
# #       S1_Base_N200/
# #         sim_001.rds
# #         sim_002.rds
# #         ...
# #         sim_100.rds
# #       S1_Base_N500/
# #       ...
# #
# # 说明:
# #   1. 本脚本只负责“数据生成与落盘”
# #   2. 不做任何模型拟合
# #   3. 每个 .rds 文件保存为一个 list，包含:
# #      - sim_id
# #      - scen_name
# #      - family
# #      - seed
# #      - data
# #      - true_weight
# #      - true_effect_mat
# #      - params
# #      - meta
# # ============================================================================

# rm(list = ls())

# # ---------------------------------------------------------------------------
# # 0. 加载依赖
# # ---------------------------------------------------------------------------
# library(NWQS)
# library(future)
# library(future.apply)

# # ---------------------------------------------------------------------------
# # 1. 全局配置
# # ---------------------------------------------------------------------------
# TARGET_FAMILY <- "gaussian"
# N_SIMULATIONS <- 100

# message(sprintf("📦 Target family: %s", toupper(TARGET_FAMILY)))
# message(sprintf("🔁 Number of simulations per scenario: %d", N_SIMULATIONS))

# # 根目录
# base_result_dir <- file.path("results")
# db_root_dir <- file.path(base_result_dir, "Monte_Carlo_DB", toupper(TARGET_FAMILY))

# if (!dir.exists(db_root_dir)) {
#     dir.create(db_root_dir, recursive = TRUE)
# }

# # ---------------------------------------------------------------------------
# # 2. 定义场景参数
# # ---------------------------------------------------------------------------
# w_norm_4 <- c(0.10, 0.20, 0.30, 0.40)
# w_sparse_4 <- c(0.60, 0.40, 0.00, 0.00)
# w_norm_8 <- c(0.04, 0.06, 0.08, 0.10, 0.14, 0.16, 0.18, 0.24)
# w_norm_12 <- c(
#     0.02, 0.03, 0.04, 0.05, 0.06, 0.08,
#     0.10, 0.11, 0.12, 0.13, 0.13, 0.13
# )
# w_sparse_12 <- c(0.40, 0.30, 0.20, 0.10, rep(0, 8))

# base_settings <- list(
#     "S1_Base" = list(
#         P = 4, corr = "mixed", shape = "threshold",
#         w = w_norm_4, snr_db = 10
#     ),
#     "S2_Linear" = list(
#         P = 4, corr = "mixed", shape = "pure_linear",
#         w = w_norm_4, snr_db = 10
#     ),
#     "S3_HighCorr" = list(
#         P = 4, corr = "high", shape = "threshold",
#         w = w_norm_4, snr_db = 10
#     ),
#     "S4_ComplexShape" = list(
#         P = 4, corr = "mixed",
#         shape = c("u_shape", "inv_threshold", "pure_linear", "s_shape"),
#         w = w_norm_4, snr_db = 10
#     ),
#     "S7_LowSNR" = list(
#         P = 4, corr = "mixed", shape = "threshold",
#         w = w_norm_4, snr_db = 2
#     ),
#     "S5_HighDim" = list(
#         P = 12, corr = "mixed", shape = "threshold",
#         w = w_norm_12, snr_db = 10
#     ),
#     "S6_Sparse_HighDim" = list(
#         P = 12, corr = "high", shape = "threshold",
#         w = w_sparse_12, snr_db = 10
#     )
# )

# N_values <- c(200, 500, 1000)

# scenarios <- list()
# for (b_name in names(base_settings)) {
#     for (n in N_values) {
#         full_name <- sprintf("%s_N%d", b_name, n)
#         curr <- base_settings[[b_name]]
#         curr$N <- n
#         scenarios[[full_name]] <- curr
#     }
# }

# message(sprintf("📚 Total scenarios to generate: %d", length(scenarios)))

# # ---------------------------------------------------------------------------
# # 3. 读取固定 Sigma 矩阵字典
# # ---------------------------------------------------------------------------
# message("\n📂 Loading fixed Sigma matrices from local files ...")

# sigma_dir <- file.path("data", "Sigma_Matrices")
# sigma_dict <- list()

# unique_combos <- unique(data.frame(
#     P = sapply(base_settings, `[[`, "P"),
#     corr = sapply(base_settings, `[[`, "corr"),
#     stringsAsFactors = FALSE
# ))

# for (i in seq_len(nrow(unique_combos))) {
#     p_val <- unique_combos$P[i]
#     corr_val <- unique_combos$corr[i]
#     dict_key <- paste0("P", p_val, "_", corr_val)

#     file_name <- sprintf("True_Sigma_Matrix_P%d_%s.csv", p_val, corr_val)
#     file_path <- file.path(sigma_dir, file_name)

#     if (!file.exists(file_path)) {
#         stop(sprintf(
#             "❌ Missing Sigma matrix file: %s\nPlease generate it first.",
#             file_path
#         ))
#     }

#     mat_df <- read.csv(file_path, row.names = 1)
#     sigma_dict[[dict_key]] <- as.matrix(mat_df)

#     message(sprintf("   ✅ Loaded: %s", dict_key))
# }

# message("✅ Sigma dictionary loaded.\n")

# # ---------------------------------------------------------------------------
# # 4. 数据生成辅助函数
# # ---------------------------------------------------------------------------
# generate_one_dataset <- function(sim_id,
#                                  scen_name,
#                                  params,
#                                  family,
#                                  sigma_dict_ref) {
#     n_vars <- params$P
#     mix_name <- paste0("Component", seq_len(n_vars))
#     mu_preds <- rep(0, n_vars)

#     dict_key <- paste0("P", params$P, "_", params$corr)
#     sigma_preds <- sigma_dict_ref[[dict_key]]

#     if (is.null(sigma_preds)) {
#         stop(sprintf("Sigma matrix not found for key: %s", dict_key))
#     }

#     w_true <- params$w
#     names(w_true) <- mix_name

#     scen_offset <- sum(utf8ToInt(scen_name))
#     sim_seed <- 10000 + scen_offset + sim_id

#     transform_fun_sim <- function(x) {
#         trans_quantile(x, q = 4)
#     }

#     # -------------------------------------------------------------------------
#     # 根据 family 调用不同生成器
#     # -------------------------------------------------------------------------
#     is_pure_linear <- is.character(params$shape) &&
#         length(params$shape) == 1 &&
#         identical(params$shape, "pure_linear")

#     if (family == "gaussian") {
#         if (is_pure_linear) {
#             # ===============================================================
#             # [修正]
#             # 纯线性场景单独使用线性数据生成器，避免混入非线性 shape 机制
#             # ===============================================================
#             sim_data <- generate_linear_data(
#                 n_obs = params$N,
#                 mu_preds = mu_preds,
#                 sigma_preds = sigma_preds,
#                 beta_preds = w_true,
#                 beta_wqs = 1,
#                 snr_db = params$snr_db,
#                 transform_fun = transform_fun_sim,
#                 seed = sim_seed
#             )

#             # ---------------------------------------------------------------
#             # 为了和其他生成器保持一致，这里补充 true_effect_mat
#             # 线性场景下，Qk vs Q1 的 true overall effect 可直接由
#             # beta_wqs * sum(w_true * (q_k - q_1)) 得到。
#             # 在 q = 4 且分位编码为 0,1,2,3 时：
#             # Q2 vs Q1 = 1 * beta_wqs * sum(w_true)
#             # Q3 vs Q1 = 2 * beta_wqs * sum(w_true)
#             # Q4 vs Q1 = 3 * beta_wqs * sum(w_true)
#             # ---------------------------------------------------------------
#             beta_eff <- w_true
#             comp_names <- names(w_true)

#             true_effect_mat <- matrix(
#                 0,
#                 nrow = length(c("Overall", comp_names)),
#                 ncol = 3,
#                 dimnames = list(
#                     c("Overall", comp_names),
#                     c("Q2_vs_Q1", "Q3_vs_Q1", "Q4_vs_Q1")
#                 )
#             )

#             for (k in 1:3) {
#                 comp_eff <- beta_eff * k
#                 true_effect_mat[comp_names, k] <- comp_eff
#                 true_effect_mat["Overall", k] <- sum(comp_eff)
#             }

#             attr(sim_data, "true_effect_mat") <- true_effect_mat
#         } else {
#             sim_data <- gen_nonlinear_data(
#                 n_obs = params$N,
#                 mu_preds = mu_preds,
#                 sigma_preds = sigma_preds,
#                 beta_preds = w_true,
#                 beta_wqs = 1,
#                 snr_db = params$snr_db,
#                 transform_fun = transform_fun_sim,
#                 q = 4,
#                 df_spline = 3,
#                 shape = params$shape,
#                 seed = sim_seed
#             )
#         }
#     } else if (family == "binomial") {
#         sim_data <- gen_nonlinear_bio_data(
#             n_obs = params$N,
#             mu_preds = mu_preds,
#             sigma_preds = sigma_preds,
#             beta_preds = w_true,
#             beta_wqs = 1,
#             target_prop = 0.3,
#             link = "logit",
#             snr_db = params$snr_db,
#             transform_fun = transform_fun_sim,
#             q = 4,
#             df_spline = 3,
#             shape = params$shape,
#             seed = sim_seed
#         )
#     } else if (family %in% c("poisson", "quasipoisson")) {
#         sim_data <- gen_nonlinear_count_data(
#             n_obs = params$N,
#             mu_preds = mu_preds,
#             sigma_preds = sigma_preds,
#             beta_preds = w_true,
#             beta_wqs = 1,
#             intercept = 0,
#             snr_db = params$snr_db,
#             transform_fun = transform_fun_sim,
#             q = 4,
#             df_spline = 3,
#             shape = params$shape,
#             seed = sim_seed
#         )
#     } else {
#         stop(sprintf("Unsupported family: %s", family))
#     }

#     # -------------------------------------------------------------------------
#     # 组织输出对象
#     # -------------------------------------------------------------------------
#     out_obj <- list(
#         sim_id = sim_id,
#         scen_name = scen_name,
#         family = family,
#         seed = sim_seed,
#         data = sim_data,
#         true_weight = w_true,
#         true_effect_mat = attr(sim_data, "true_effect_mat"),
#         params = params,
#         meta = list(
#             P = params$P,
#             N = params$N,
#             corr = params$corr,
#             shape = params$shape,
#             snr_db = params$snr_db,
#             mix_name = mix_name,
#             sigma_key = dict_key,
#             created_at = Sys.time()
#         )
#     )

#     return(out_obj)
# }

# # ---------------------------------------------------------------------------
# # 5. 单场景数据库生成函数
# # ---------------------------------------------------------------------------
# generate_scenario_db <- function(scen_name,
#                                  params,
#                                  family,
#                                  sigma_dict_ref,
#                                  n_sim = 100,
#                                  db_root_dir) {
#     scen_dir <- file.path(db_root_dir, scen_name)

#     if (!dir.exists(scen_dir)) {
#         dir.create(scen_dir, recursive = TRUE)
#     }

#     message(sprintf("\n=================================================="))
#     message(sprintf("🧪 Generating datasets for scenario: %s", scen_name))
#     message(sprintf("   N = %d | P = %d | corr = %s", params$N, params$P, params$corr))
#     message(sprintf("=================================================="))

#     # ---------------------------------------------------------------
#     # 断点续跑: 只生成缺失的 sim_XXX.rds
#     # ---------------------------------------------------------------
#     target_files <- sprintf("sim_%03d.rds", seq_len(n_sim))
#     target_paths <- file.path(scen_dir, target_files)
#     need_idx <- which(!file.exists(target_paths))

#     if (length(need_idx) == 0) {
#         message("⏭️ All datasets already exist. Skip.")
#         return(invisible(NULL))
#     }

#     message(sprintf("📌 Missing datasets to generate: %d / %d", length(need_idx), n_sim))

#     # ---------------------------------------------------------------
#     # 并行生成缺失数据
#     # ---------------------------------------------------------------
#     gen_results <- future.apply::future_lapply(
#         need_idx,
#         function(i) {
#             tryCatch(
#                 {
#                     obj <- generate_one_dataset(
#                         sim_id = i,
#                         scen_name = scen_name,
#                         params = params,
#                         family = family,
#                         sigma_dict_ref = sigma_dict_ref
#                     )

#                     out_file <- file.path(scen_dir, sprintf("sim_%03d.rds", i))
#                     saveRDS(obj, out_file, compress = "xz")

#                     return(list(success = TRUE, sim_id = i, file = out_file))
#                 },
#                 error = function(e) {
#                     return(list(
#                         success = FALSE,
#                         sim_id = i,
#                         error_msg = e$message
#                     ))
#                 }
#             )
#         },
#         future.seed = TRUE
#     )

#     # ---------------------------------------------------------------
#     # 汇总日志
#     # ---------------------------------------------------------------
#     ok_res <- Filter(function(x) isTRUE(x$success), gen_results)
#     bad_res <- Filter(function(x) !isTRUE(x$success), gen_results)

#     message(sprintf("✅ Successfully generated: %d", length(ok_res)))

#     if (length(bad_res) > 0) {
#         message(sprintf("⚠️ Failed datasets: %d", length(bad_res)))
#         for (k in seq_len(min(5, length(bad_res)))) {
#             message(sprintf(
#                 "   -> sim_%03d failed: %s",
#                 bad_res[[k]]$sim_id,
#                 bad_res[[k]]$error_msg
#             ))
#         }
#     }

#     # ---------------------------------------------------------------
#     # 保存场景级元信息
#     # ---------------------------------------------------------------
#     scen_meta <- list(
#         scen_name = scen_name,
#         family = family,
#         n_simulations = n_sim,
#         params = params,
#         generated_at = Sys.time()
#     )

#     saveRDS(
#         scen_meta,
#         file.path(scen_dir, "scenario_meta.rds"),
#         compress = "xz"
#     )

#     invisible(gen_results)
# }

# # ---------------------------------------------------------------------------
# # 6. 主循环：生成全部场景数据库
# # ---------------------------------------------------------------------------
# start_time <- Sys.time()
# total_tasks <- length(names(scenarios))
# current_task <- 0

# message("\n🚀 Starting Monte Carlo dataset generation pipeline ...")

# for (scen_name in names(scenarios)) {
#     current_task <- current_task + 1
#     params <- scenarios[[scen_name]]

#     message(sprintf(
#         "\n[%d/%d] Processing scenario: %s",
#         current_task, total_tasks, scen_name
#     ))

#     tryCatch(
#         {
#             generate_scenario_db(
#                 scen_name = scen_name,
#                 params = params,
#                 family = TARGET_FAMILY,
#                 sigma_dict_ref = sigma_dict,
#                 n_sim = N_SIMULATIONS,
#                 db_root_dir = db_root_dir
#             )
#         },
#         error = function(e) {
#             message(sprintf(
#                 "❌ Scenario failed: %s | Error: %s",
#                 scen_name, e$message
#             ))
#         }
#     )

#     gc(verbose = FALSE)
# }

# elapsed_min <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

# message("\n==================================================")
# message("🎉 Monte Carlo dataset generation completed!")
# message(sprintf("📁 Output root: %s", db_root_dir))
# message(sprintf("⏱️ Total elapsed time: %.2f minutes", elapsed_min))
# message("==================================================\n")

# future::plan(future::sequential)
