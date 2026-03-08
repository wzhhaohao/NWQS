# =========================================================================
# Appendix S1: Split Ratio Sensitivity Analysis (Publication Code)
#
# Scenarios: S1_Base + S3_HighCorr (robustness check)
# Sample sizes: N = 200, 500, 1000
# Split proportions: 0.3, 0.4, 0.5, 0.6, 0.7, 0.8
# Replicates per cell: 50 simulation datasets, rh = 30 each
# =========================================================================

library(NWQS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Load MC utilities (calc_true_importance, calc_weight_error)
# These are in NWQS/R/monte_carlo.R — loaded automatically via library(NWQS)
# If running without installing the package, uncomment:
# source('/Users/wangzhehao/temporary/packages/NWQS/R/monte_carlo.R')

# ── Configuration ────────────────────────────────────────────────────────
PROJECT_ROOT <- '/Users/wangzhehao/temporary/packages/NWQS result'
DB_ROOT      <- file.path(PROJECT_ROOT, "results/Monte_Carlo_DB/GAUSSIAN")

SPLIT_PROPS   <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
SAMPLE_SIZES  <- c(200, 500, 1000)
SCENARIOS     <- c("S1_Base", "S3_HighCorr")  # Main + Robustness
SIM_IDS       <- 1:50                          # 50 replicates per cell
RH_EACH       <- 30
N_PERM        <- 100
COVARIATES    <- c("x_cont", "x_bin", "x_cat")

OUT_DIR <- file.path(PROJECT_ROOT, "results", "Split_Ratio_Sensitivity")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

options(max.print = 9999)  # Prevent R from truncating summary table output
set.seed(2025)

# ── Helper ───────────────────────────────────────────────────────────────
eval_full <- function(wqs_score, data, outcome = "y", covariates, family = "gaussian") {
  data$wqs_score <- wqs_score
  fm <- as.formula(paste(outcome, "~ wqs_score +", paste(covariates, collapse = " + ")))
  fit <- glm(fm, data = data, family = family)
  y <- data[[outcome]]; yh <- fitted(fit); r <- y - yh
  list(R2 = 1 - sum(r^2) / sum((y - mean(y))^2), RMSE = sqrt(mean(r^2)))
}

# =========================================================================
# Main Experiment
# =========================================================================
results_all <- list()
total_cells <- length(SCENARIOS) * length(SAMPLE_SIZES) * length(SPLIT_PROPS) * length(SIM_IDS)
counter <- 0

for (scen_base in SCENARIOS) {
  for (N in SAMPLE_SIZES) {
    scen_name <- sprintf("%s_N%d", scen_base, N)
    scen_dir  <- file.path(DB_ROOT, scen_name)

    if (!dir.exists(scen_dir)) {
      message(sprintf("  [SKIP] %s: directory not found", scen_name))
      next
    }

    # Get mix_name from first sim
    first_sim <- readRDS(file.path(scen_dir, "sim_001.rds"))
    mix_name  <- first_sim$meta$mix_name

    for (sim_id in SIM_IDS) {
      sim_path <- file.path(scen_dir, sprintf("sim_%03d.rds", sim_id))
      if (!file.exists(sim_path)) next

      sim_obj  <- readRDS(sim_path)
      sim_data <- sim_obj$data
      w_true   <- calc_true_importance(sim_obj$true_effect_mat, mix_name)

      for (sp in SPLIT_PROPS) {
        counter <- counter + 1

        fit <- tryCatch(
          nwqs(
            data = sim_data, mix_name = mix_name, covariates = COVARIATES,
            outcome = "y", q = 4, split_prop = sp, rh = RH_EACH,
            n_permutation = N_PERM, family = "gaussian", seed = NULL,
            transform_fun = function(x) trans_quantile(x, q = 4),
            plan_strategy = "multisession", n_workers = 8, quiet = TRUE
          ),
          error = function(e) NULL
        )
        if (is.null(fit)) next

        w_est <- fit$final_weights[mix_name]
        err   <- calc_weight_error(w_est, w_true)
        perf  <- eval_full(fit$data$wqs_score, sim_data, covariates = COVARIATES)

        results_all[[length(results_all) + 1]] <- data.frame(
          Scenario   = scen_base,
          N          = N,
          split_prop = sp,
          Sim_ID     = sim_id,
          SAE        = err$SAE,
          Pearson    = err$Pearson,
          Spearman   = err$Spearman,
          R2         = perf$R2,
          RMSE       = perf$RMSE,
          Beta_WQS   = fit$mean_coefs["wqs_score"]
        )
      }

      if (counter %% 150 == 0)
        cat(sprintf("  Progress: %d/%d (%s, N=%d, sim_%03d)\n",
                    counter, total_cells, scen_base, N, sim_id))
    }
    n_done <- sum(sapply(results_all, function(x) x$Scenario[1] == scen_base & x$N[1] == N))
    cat(sprintf("  Completed: %s (%d fits)\n", scen_name, n_done))
  }
}

df_all <- do.call(rbind, results_all)
rownames(df_all) <- NULL
write.csv(df_all, file.path(OUT_DIR, "Raw_Split_Sensitivity.csv"), row.names = FALSE)

# =========================================================================
# Summary Table
# =========================================================================
summary_table <- df_all %>%
  group_by(Scenario, N, split_prop) %>%
  summarise(
    N_Sims        = n(),
    Mean_SAE      = mean(SAE, na.rm = TRUE),
    SD_SAE        = sd(SAE, na.rm = TRUE),
    Mean_Pearson  = mean(Pearson, na.rm = TRUE),
    Mean_Spearman = mean(Spearman, na.rm = TRUE),
    Pearson_GE08  = mean(Pearson >= 0.95, na.rm = TRUE),
    Spearman_GE08 = mean(Spearman >= 0.95, na.rm = TRUE),
    Mean_R2       = mean(R2, na.rm = TRUE),
    SD_R2         = sd(R2, na.rm = TRUE),
    Mean_RMSE     = mean(RMSE, na.rm = TRUE),
    Mean_Beta     = mean(Beta_WQS, na.rm = TRUE),
    SD_Beta       = sd(Beta_WQS, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Scenario, N, split_prop)

write.csv(summary_table, file.path(OUT_DIR, "Summary_Split_Sensitivity.csv"), row.names = FALSE)

# Print per-scenario summary for readability
for (scen in unique(summary_table$Scenario)) {
  cat(sprintf("\n========== %s ==========\n", scen))
  sub <- summary_table %>% filter(Scenario == scen) %>% select(-Scenario)
  print(as.data.frame(sub), digits = 4, row.names = FALSE)
}

# =========================================================================
# Optimal split per (Scenario, N) — by minimum Mean_SAE
# Also show runner-up (split_prop = 0.6) for comparison
# =========================================================================
optimal <- summary_table %>%
  group_by(Scenario, N) %>%
  slice_min(Mean_SAE, n = 1) %>%
  ungroup() %>%
  select(Scenario, N, split_prop, Mean_SAE, Mean_R2, SD_Beta, Pearson_GE08, Spearman_GE08)

cat("\n========== Optimal Split Proportion (by min Mean_SAE) ==========\n")
print(as.data.frame(optimal), row.names = FALSE)

# Show split_prop = 0.6 row for comparison
ref_06 <- summary_table %>%
  filter(split_prop == 0.6) %>%
  select(Scenario, N, split_prop, Mean_SAE, Mean_R2, SD_Beta, Pearson_GE08, Spearman_GE08)

cat("\n========== Reference: split_prop = 0.6 (proposed default) ==========\n")
print(as.data.frame(ref_06), row.names = FALSE)

# Compute the SAE gap between optimal and 0.6
comparison <- merge(
  optimal %>% select(Scenario, N, Opt_split = split_prop, Opt_SAE = Mean_SAE, Opt_SD_Beta = SD_Beta),
  ref_06 %>% select(Scenario, N, Ref_SAE = Mean_SAE, Ref_SD_Beta = SD_Beta),
  by = c("Scenario", "N")
)
comparison$SAE_Gap <- comparison$Ref_SAE - comparison$Opt_SAE
comparison$Beta_SD_Gap <- comparison$Ref_SD_Beta - comparison$Opt_SD_Beta

cat("\n========== Gap Analysis: Optimal vs split_prop=0.6 ==========\n")
print(as.data.frame(comparison), digits = 4, row.names = FALSE)

# =========================================================================
# Visualization — Faceted by Scenario (rows) × N (columns)
# =========================================================================
df_plot <- df_all
df_plot$split_label <- sprintf("%.0f/%.0f",
  df_plot$split_prop * 100, (1 - df_plot$split_prop) * 100)
df_plot$split_label <- factor(df_plot$split_label,
  levels = sprintf("%.0f/%.0f", SPLIT_PROPS * 100, (1 - SPLIT_PROPS) * 100))
df_plot$N_label <- factor(sprintf("N = %d", df_plot$N),
  levels = sprintf("N = %d", SAMPLE_SIZES))

# ── Figure S1a: SAE ─────────────────────────────────────────────────────
p_sae <- ggplot(df_plot, aes(x = split_label, y = SAE)) +
  geom_boxplot(fill = "#4A90C8", alpha = 0.7, outlier.size = 0.3, width = 0.6) +
  facet_grid(Scenario ~ N_label, scales = "free_y") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#ECF0F1"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "A. Weight Recovery Error (SAE)",
       subtitle = "Lower SAE indicates more accurate weight estimation",
       x = "Train / Validate (%)", y = "SAE")

# ── Figure S1b: R² ──────────────────────────────────────────────────────
p_r2 <- ggplot(df_plot, aes(x = split_label, y = R2)) +
  geom_boxplot(fill = "#6EC44A", alpha = 0.7, outlier.size = 0.3, width = 0.6) +
  facet_grid(Scenario ~ N_label, scales = "free_y") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#ECF0F1"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "B. Full-Data Prediction Accuracy (R\u00B2)",
       subtitle = "Higher R\u00B2 indicates better overall model fit",
       x = "Train / Validate (%)", y = "R\u00B2")

# ── Figure S1c: Pearson ─────────────────────────────────────────────────
p_pear <- ggplot(df_plot, aes(x = split_label, y = Pearson)) +
  geom_boxplot(fill = "#D92828", alpha = 0.7, outlier.size = 0.3, width = 0.6) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black", linewidth = 0.4) +
  facet_grid(Scenario ~ N_label, scales = "free_y") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#ECF0F1"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "C. Weight-Truth Pearson Correlation",
       subtitle = "Dashed line: r = 0.95 threshold",
       x = "Train / Validate (%)", y = "Pearson r")

# ── Figure S1d: β stability ─────────────────────────────────────────────
p_beta <- ggplot(df_plot, aes(x = split_label, y = Beta_WQS)) +
  geom_boxplot(fill = "#8B6FB8", alpha = 0.7, outlier.size = 0.3, width = 0.6) +
  facet_grid(Scenario ~ N_label, scales = "free_y") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#ECF0F1"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "D. Mixture Effect Coefficient Stability (\u03B2_wqs)",
       subtitle = "Tighter distribution indicates more stable estimation",
       x = "Train / Validate (%)", y = "\u03B2_wqs")

# ── Combine ──────────────────────────────────────────────────────────────
final_fig <- (p_sae / p_r2 / p_pear / p_beta) +
  patchwork::plot_annotation(
    title    = "Figure S1. Split Ratio Sensitivity Analysis",
    subtitle = sprintf("Gaussian family | %d replicates per cell | rh = %d per fit", length(SIM_IDS), RH_EACH),
    theme    = theme(
      plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#7F8C8D")
    )
  )

fig_path_png <- file.path(OUT_DIR, "Figure_S1_Split_Sensitivity.png")
fig_path_pdf <- file.path(OUT_DIR, "Figure_S1_Split_Sensitivity.pdf")
ggsave(fig_path_png, final_fig, width = 16, height = 22, dpi = 300)
ggsave(fig_path_pdf, final_fig, width = 16, height = 22, device = "pdf")

cat(sprintf("\nFigures saved:\n  %s\n  %s\n", fig_path_png, fig_path_pdf))
cat(sprintf("Tables saved in: %s\n", OUT_DIR))
message("\n=== Split Ratio Sensitivity Analysis Complete ===")
