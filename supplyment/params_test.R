# =========================================================================
# NWQS 超参数敏感性分析 (nwqs_boot 版)
#
# 维度: N (样本量) × P (混合物数量) × SNR (信噪比)
# 数据: 使用 gen_nonlinear_data() 即时生成，引用预存 Sigma 矩阵 (mixed)
# 核心指标: SAE, R², Pearson, Spearman, β_wqs, Coverage
# 并行策略: multisession (跨平台兼容)
# =========================================================================

rm(list = ls())

library(NWQS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(future)
library(future.apply)

# ── Configuration ────────────────────────────────────────────────────────
PROJECT_ROOT <- '/Users/wangzhehao/temporary/packages/NWQS result'
SIGMA_DIR    <- file.path(PROJECT_ROOT, "data", "Sigma_Matrices")

N_SIMS_PER_CELL <- 20    # 每个网格 20 次模拟
N_BOOT          <- 100    # nwqs_boot 外层 bootstrap
N_PERM          <- 50     # 内部置换次数
WORKER_CORES    <- 8

OUT_DIR <- file.path(PROJECT_ROOT, "results", "Sensitivity_Analysis")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ── Weight dictionaries ──────────────────────────────────────────────────
w_dict <- list(
  "4"  = c(0.10, 0.20, 0.30, 0.40),
  "8"  = c(0.04, 0.06, 0.08, 0.10, 0.14, 0.16, 0.18, 0.24),
  "12" = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.11, 0.12, 0.13, 0.13, 0.13)
)

# ── Pre-load Sigma matrices (mixed correlation only) ─────────────────────
sigma_dict <- list()
for (p in c(4, 8, 12)) {
  sigma_file <- file.path(SIGMA_DIR, sprintf("True_Sigma_Matrix_P%d_mixed.csv", p))
  if (!file.exists(sigma_file)) stop(sprintf("Sigma matrix not found: %s", sigma_file))
  sigma_dict[[as.character(p)]] <- as.matrix(read.csv(sigma_file, row.names = 1))
  message(sprintf("  Loaded Sigma P=%d (%dx%d)", p, p, p))
}

# ── Parameter grid: N × P × SNR ─────────────────────────────────────────
param_grid <- expand.grid(
  N   = c(200, 500, 1000),
  P   = c(4, 8, 12),
  SNR = c(2, 5, 10, 20),
  stringsAsFactors = FALSE
)

# ── Parallel setup ───────────────────────────────────────────────────────
plan(multisession, workers = WORKER_CORES)

total_fits <- nrow(param_grid) * N_SIMS_PER_CELL
cat(sprintf("=== Sensitivity Analysis ===\n"))
cat(sprintf("Grid: %d cells (N: %s × P: %s × SNR: %s)\n",
            nrow(param_grid),
            paste(unique(param_grid$N), collapse = ", "),
            paste(unique(param_grid$P), collapse = ", "),
            paste(unique(param_grid$SNR), collapse = ", ")))
cat(sprintf("Sims per cell: %d | Total fits: %d\n", N_SIMS_PER_CELL, total_fits))
cat(sprintf("Parallel: multisession, %d workers\n\n", WORKER_CORES))

# ── Helpers ──────────────────────────────────────────────────────────────
eval_full <- function(wqs_score, data, covariates) {
  data$wqs_score <- wqs_score
  fm <- as.formula(paste("y ~ wqs_score +", paste(covariates, collapse = " + ")))
  fit <- glm(fm, data = data, family = "gaussian")
  y <- data$y; yh <- fitted(fit); r <- y - yh
  list(R2 = 1 - sum(r^2) / sum((y - mean(y))^2), RMSE = sqrt(mean(r^2)))
}

safe_normalize <- function(x, nm) {
  x <- as.numeric(x); x[!is.finite(x) | x < 0] <- 0
  out <- if (sum(x) <= 0) rep(1 / length(x), length(x)) else x / sum(x)
  names(out) <- nm; out
}

# =========================================================================
# Core: Single Simulation (self-contained, safe for future workers)
# =========================================================================
run_one_sim <- function(N_val, P_val, SNR_val, sim_id,
                        sigma_mat, beta_preds, mix_name) {

  covariates    <- c("x_cont", "x_bin", "x_cat")
  mu_preds      <- rep(0, P_val)
  transform_fun <- function(x) trans_quantile(x, q = 4)

  # Deterministic seed per (N, P, SNR, sim_id)
  iter_seed <- as.integer(10000 + sim_id + N_val * 7 + P_val * 113 + SNR_val * 31)

  # Generate data on the fly
  sim_data <- gen_nonlinear_data(
    n_obs = N_val, mu_preds = mu_preds, sigma_preds = sigma_mat,
    beta_preds = beta_preds, beta_wqs = 1, snr_db = SNR_val,
    transform_fun = transform_fun, q = 4, df_spline = 3,
    shape = "threshold", seed = iter_seed
  )
  true_eff_mat <- attr(sim_data, "true_effect_mat")

  # True importance from actual spline effects (not input weights)
  w_true <- calc_true_importance(true_eff_mat, mix_name)

  # Fit nwqs_boot (rh_inner=1, inner sequential)
  fit <- nwqs_boot(
    data = sim_data, mix_name = mix_name, covariates = covariates,
    outcome = "y", family = "gaussian",
    n_boot = N_BOOT, rh_inner = 1, n_permutation = N_PERM,
    q = 4, split_prop = 0.6, seed = NULL,
    transform_fun = transform_fun,
    plan_strategy = "sequential",
    keep_fits = FALSE
  )

  w_est <- safe_normalize(fit$final_weights[mix_name], mix_name)
  err   <- calc_weight_error(w_est, w_true)
  perf  <- eval_full(fit$data$wqs_score, sim_data, covariates)

  # Bootstrap coverage: Overall highest-Q vs Q1
  q_max <- fit$q
  target_col <- paste0("Q", q_max, "_vs_Q1")
  ci_row <- fit$ci_table[fit$ci_table$Target == target_col &
                          fit$ci_table$Term == "Overall", ]
  true_overall <- true_eff_mat["Overall", target_col]
  covered <- if (nrow(ci_row) > 0 && !is.na(ci_row$Boot_CI_Lower[1])) {
    ci_row$Boot_CI_Lower[1] <= true_overall & ci_row$Boot_CI_Upper[1] >= true_overall
  } else NA

  data.frame(
    N = N_val, P = P_val, SNR = SNR_val, Sim_ID = sim_id,
    SAE = err$SAE, Pearson = err$Pearson, Spearman = err$Spearman,
    R2 = perf$R2, RMSE = perf$RMSE,
    Beta_WQS = fit$mean_coefs["wqs_score"],
    Covered = as.logical(covered)
  )
}

# =========================================================================
# Execute Grid
# =========================================================================
all_results <- list()
start_global <- Sys.time()

for (task_idx in seq_len(nrow(param_grid))) {
  curr_N   <- param_grid$N[task_idx]
  curr_P   <- param_grid$P[task_idx]
  curr_SNR <- param_grid$SNR[task_idx]

  mix_name   <- paste0("Component", 1:curr_P)
  sigma_mat  <- sigma_dict[[as.character(curr_P)]]
  beta_preds <- w_dict[[as.character(curr_P)]]

  cat(sprintf("[%02d/%d] N=%4d, P=%2d, SNR=%2d ...",
              task_idx, nrow(param_grid), curr_N, curr_P, curr_SNR))
  t0 <- Sys.time()

  # Parallelize the 20 sims within each cell
  cell_results <- future.apply::future_lapply(
    seq_len(N_SIMS_PER_CELL),
    function(sid) {
      tryCatch(
        run_one_sim(curr_N, curr_P, curr_SNR, sid,
                    sigma_mat, beta_preds, mix_name),
        error = function(e) NULL
      )
    },
    future.seed = TRUE,
    future.packages = c("NWQS", "dplyr")
  )

  cell_df <- do.call(rbind, Filter(Negate(is.null), cell_results))
  all_results[[task_idx]] <- cell_df

  n_ok <- if (!is.null(cell_df)) nrow(cell_df) else 0
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  cat(sprintf(" %d/%d ok, %.1fs\n", n_ok, N_SIMS_PER_CELL, elapsed))

  gc(verbose = FALSE)
}

final_df <- do.call(rbind, all_results)
rownames(final_df) <- NULL
write.csv(final_df, file.path(OUT_DIR, "Sensitivity_Raw.csv"), row.names = FALSE)

# =========================================================================
# Summary
# =========================================================================
options(max.print = 9999)

summary_df <- final_df %>%
  group_by(N, P, SNR) %>%
  summarise(
    N_Valid       = n(),
    Mean_SAE      = mean(SAE, na.rm = TRUE),
    SD_SAE        = sd(SAE, na.rm = TRUE),
    Mean_R2       = mean(R2, na.rm = TRUE),
    SD_R2         = sd(R2, na.rm = TRUE),
    Mean_Pearson  = mean(Pearson, na.rm = TRUE),
    Pearson_Pass  = mean(Pearson >= 0.8, na.rm = TRUE),
    Mean_Spearman = mean(Spearman, na.rm = TRUE),
    Spearman_Pass = mean(Spearman >= 0.8, na.rm = TRUE),
    Mean_Beta     = mean(Beta_WQS, na.rm = TRUE),
    SD_Beta       = sd(Beta_WQS, na.rm = TRUE),
    Coverage      = mean(Covered, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(P, N, SNR)

write.csv(summary_df, file.path(OUT_DIR, "Sensitivity_Summary.csv"), row.names = FALSE)

cat("\n========== Summary ==========\n")
print(as.data.frame(summary_df), digits = 4)

# =========================================================================
# Visualization
# =========================================================================
plot_df <- final_df %>%
  mutate(
    N_label   = factor(sprintf("N = %d", N), levels = sprintf("N = %d", c(200, 500, 1000))),
    P_label   = factor(sprintf("P = %d", P), levels = sprintf("P = %d", c(4, 8, 12))),
    SNR_label = factor(sprintf("SNR = %d dB", SNR),
                       levels = sprintf("SNR = %d dB", c(2, 5, 10, 20)))
  )

base_theme <- theme_bw(base_size = 13) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "#ECF0F1"),
        strip.text = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold"))

pal_N   <- c("N = 200" = "#D92828", "N = 500" = "#8B6FB8", "N = 1000" = "#4A90C8")
pal_SNR <- c("SNR = 2 dB" = "#D92828", "SNR = 5 dB" = "#8B6FB8",
             "SNR = 10 dB" = "#4A90C8", "SNR = 20 dB" = "#6EC44A")

# ── Fig 1: Impact of N (faceted by P, colored by SNR) ────────────────────
p1_sae <- ggplot(plot_df, aes(x = N_label, y = SAE, fill = SNR_label)) +
  geom_boxplot(position = position_dodge(0.75), alpha = 0.6, outlier.size = 0.5) +
  facet_wrap(~P_label) + scale_fill_manual(values = pal_SNR) + base_theme +
  labs(title = "A. Weight Recovery Error (SAE) by Sample Size",
       x = "Sample Size", y = "SAE", fill = "Signal-to-Noise Ratio")

p1_r2 <- ggplot(plot_df, aes(x = N_label, y = R2, fill = SNR_label)) +
  geom_boxplot(position = position_dodge(0.75), alpha = 0.6, outlier.size = 0.5) +
  facet_wrap(~P_label, scales = "free_y") + scale_fill_manual(values = pal_SNR) + base_theme +
  labs(title = "B. Prediction Accuracy (R\u00B2) by Sample Size",
       x = "Sample Size", y = "R\u00B2", fill = "Signal-to-Noise Ratio")

fig1 <- (p1_sae / p1_r2) +
  patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "Sensitivity_Fig1_Impact_N.png"), fig1, width = 14, height = 10, dpi = 300)
ggsave(file.path(OUT_DIR, "Sensitivity_Fig1_Impact_N.pdf"), fig1, width = 14, height = 10, device = "pdf")

# ── Fig 2: Impact of P (faceted by SNR, colored by N) ────────────────────
p2_sae <- ggplot(plot_df, aes(x = P_label, y = SAE, fill = N_label)) +
  geom_boxplot(position = position_dodge(0.75), alpha = 0.6, outlier.size = 0.5) +
  facet_wrap(~SNR_label) + scale_fill_manual(values = pal_N) + base_theme +
  labs(title = "A. Weight Recovery Error (SAE) by Dimensionality",
       x = "Dimensionality (P)", y = "SAE", fill = "Sample Size")

p2_pearson <- ggplot(plot_df, aes(x = P_label, y = Pearson, fill = N_label)) +
  geom_boxplot(position = position_dodge(0.75), alpha = 0.6, outlier.size = 0.5) +
  geom_hline(yintercept = 0.8, linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~SNR_label) + scale_fill_manual(values = pal_N) + base_theme +
  labs(title = "B. Weight-Truth Pearson Correlation by Dimensionality",
       x = "Dimensionality (P)", y = "Pearson r", fill = "Sample Size")

fig2 <- (p2_sae / p2_pearson) +
  patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "Sensitivity_Fig2_Impact_P.png"), fig2, width = 14, height = 10, dpi = 300)
ggsave(file.path(OUT_DIR, "Sensitivity_Fig2_Impact_P.pdf"), fig2, width = 14, height = 10, device = "pdf")

# ── Fig 3: Impact of SNR (faceted by P, colored by N) ────────────────────
p3_sae <- ggplot(plot_df, aes(x = SNR_label, y = SAE, fill = N_label)) +
  geom_boxplot(position = position_dodge(0.75), alpha = 0.6, outlier.size = 0.5) +
  facet_wrap(~P_label) + scale_fill_manual(values = pal_N) + base_theme +
  labs(title = "A. Weight Recovery Error (SAE) by Signal-to-Noise Ratio",
       x = "SNR (dB)", y = "SAE", fill = "Sample Size")

p3_beta <- ggplot(plot_df, aes(x = SNR_label, y = Beta_WQS, fill = N_label)) +
  geom_boxplot(position = position_dodge(0.75), alpha = 0.6, outlier.size = 0.5) +
  facet_wrap(~P_label, scales = "free_y") + scale_fill_manual(values = pal_N) + base_theme +
  labs(title = "B. Effect Estimate Stability (\u03B2_wqs) by SNR",
       x = "SNR (dB)", y = "\u03B2_wqs", fill = "Sample Size")

fig3 <- (p3_sae / p3_beta) +
  patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "Sensitivity_Fig3_Impact_SNR.png"), fig3, width = 14, height = 10, dpi = 300)
ggsave(file.path(OUT_DIR, "Sensitivity_Fig3_Impact_SNR.pdf"), fig3, width = 14, height = 10, device = "pdf")

# ── Cleanup ──────────────────────────────────────────────────────────────
future::plan(future::sequential)

elapsed_total <- as.numeric(difftime(Sys.time(), start_global, units = "mins"))
cat(sprintf("\n=== Sensitivity Analysis Complete! %.1f mins ===\n", elapsed_total))
cat(sprintf("Results: %s\n", OUT_DIR))
