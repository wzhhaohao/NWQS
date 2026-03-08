# =========================================================================
# NWQS Monte Carlo Simulation Study — Main Analysis
# Benchmark: NWQS Boot vs gWQS Boot vs QGcomp Boot
#
# Evaluation Metrics:
#   1. R² and RMSE (full-data refit, sample-size invariant)
#   2. SAE (weight recovery error)
#   3. Weight correlation success rate (Pearson ≥ 0.8, Spearman ≥ 0.8)
#   4. Weight recovery distribution (per-component)
#   5. NWQS Bootstrap CI coverage (self-calibration only)
#
# Data: Pre-generated in Monte_Carlo_DB/{FAMILY}/{Scenario}/sim_{001..100}.rds
# =========================================================================

# -------------------------------------------------------------------------
# 0. Configuration
# -------------------------------------------------------------------------
rm(list = ls())

library(NWQS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gWQS)
library(qgcomp)
library(future)
library(future.apply)

# ── User Settings ────────────────────────────────────────────────────────
TARGET_FAMILY <- "gaussian" # "gaussian", "binomial", "quasipoisson"
N_SIMULATIONS <- 100
N_BOOT_NWQS <- 100 # nwqs_boot outer bootstrap
N_BOOT_GWQS <- 100 # gWQS bootstrap samples
N_BOOT_QGCOMP <- 100 # qgcomp.boot B
N_PERM_NWQS <- 100 # nwqs internal permutations
WORKER_CORES <- 8

# ── Paths ────────────────────────────────────────────────────────────────
PROJECT_ROOT <- "/Users/wangzhehao/temporary/packages/NWQS result"
db_root_dir <- file.path(PROJECT_ROOT, "results", "Monte_Carlo_DB", toupper(TARGET_FAMILY))
out_root_dir <- file.path(PROJECT_ROOT, "results", "Monte_Carlo_Results", toupper(TARGET_FAMILY))

if (!dir.exists(db_root_dir)) stop(sprintf("DB not found: %s", db_root_dir))
if (!dir.exists(out_root_dir)) dir.create(out_root_dir, recursive = TRUE)

plan(multisession, workers = WORKER_CORES)
message(sprintf("Parallel: %d workers | Family: %s", WORKER_CORES, toupper(TARGET_FAMILY)))

# -------------------------------------------------------------------------
# 1. Utility Functions
# -------------------------------------------------------------------------

safe_normalize <- function(x, nm = names(x)) {
  x <- as.numeric(x)
  x[!is.finite(x) | x < 0] <- 0
  out <- if (sum(x) <= 0) rep(1 / length(x), length(x)) else x / sum(x)
  if (!is.null(nm)) names(out) <- nm
  out
}

#' Refit a GLM on full data using a pre-computed wqs_score vector
#' Returns R² and RMSE (or pseudo-R² for non-gaussian)
refit_and_evaluate <- function(wqs_score, data, outcome, covariates, family) {
  data$wqs_score <- wqs_score
  formula <- as.formula(paste(outcome, "~ wqs_score +", paste(covariates, collapse = " + ")))
  fit <- glm(formula, data = data, family = family)

  y <- data[[outcome]]
  y_hat <- fitted(fit)
  residuals <- y - y_hat

  rmse <- sqrt(mean(residuals^2))
  ss_res <- sum(residuals^2)
  ss_tot <- sum((y - mean(y))^2)
  r2 <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_

  list(R2 = r2, RMSE = rmse)
}

# -------------------------------------------------------------------------
# 2. Core Engine: Single Simulation
# -------------------------------------------------------------------------
run_one_sim <- function(sim_file, family, mix_name) {
  sim_obj <- readRDS(sim_file)
  sim_data <- sim_obj$data
  sim_id <- sim_obj$sim_id
  true_eff_mat <- sim_obj$true_effect_mat
  Y <- sim_data$y
  covariates <- c("x_cont", "x_bin", "x_cat")
  transform_fun <- function(x) trans_quantile(x, q = 4)

  w_true <- calc_true_importance(true_eff_mat, mix_name)

  fam_qgcomp <- switch(family,
    gaussian = gaussian(),
    binomial = binomial(),
    poisson()
  )

  # ── 1. NWQS Boot (rh_inner=1 → valid sampling CI) ─────────────────
  nwqs_res <- nwqs_boot(
    data = sim_data, mix_name = mix_name, covariates = covariates,
    outcome = "y", family = family,
    n_boot = N_BOOT_NWQS, rh_inner = 1, n_permutation = N_PERM_NWQS,
    q = 4, split_prop = 0.6, seed = NULL,
    transform_fun = transform_fun,
    plan_strategy = "sequential", keep_fits = FALSE
  )
  w_nwqs <- safe_normalize(nwqs_res$final_weights[mix_name], mix_name)

  # NWQS wqs_score is already on the full data
  eval_nwqs <- refit_and_evaluate(
    nwqs_res$data$wqs_score, sim_data, "y", covariates, family
  )

  # ── 2. gWQS Boot (rh=1, b=100) ────────────────────────────────────
  gwqs_res <- suppressWarnings(gWQS::gwqs(
    formula = y ~ wqs + x_cont + x_bin + x_cat,
    data = sim_data, mix_name = mix_name,
    q = 4, validation = 0.6, b = N_BOOT_GWQS, rh = 1,
    plan_strategy = "sequential", family = family, seed = NULL
  ))
  w_gwqs_raw <- gwqs_res$final_weights$mean_weight
  names(w_gwqs_raw) <- gwqs_res$final_weights$mix_name
  w_gwqs <- safe_normalize(w_gwqs_raw[mix_name], mix_name)

  # gWQS: reconstruct wqs_score on full data, then refit
  data_q <- as.data.frame(transform_fun(sim_data[mix_name]))
  wqs_gwqs_full <- as.vector(as.matrix(data_q[mix_name]) %*% w_gwqs[mix_name])
  eval_gwqs <- refit_and_evaluate(wqs_gwqs_full, sim_data, "y", covariates, family)

  # ── 3. QGcomp Boot ────────────────────────────────────────────────
  data_q_full <- sim_data
  data_q_full[mix_name] <- transform_fun(data_q_full[mix_name])
  formula_qg <- as.formula(paste("y ~", paste(c(mix_name, covariates), collapse = " + ")))

  qgcomp_res <- suppressWarnings(qgcomp::qgcomp.boot(
    f = formula_qg, expnms = mix_name, data = data_q_full,
    family = fam_qgcomp, q = NULL, B = N_BOOT_QGCOMP,
    degree = 2, seed = NULL
  ))

  # QGcomp weights: combine pos + neg absolute contributions
  w_qg_all <- rep(0, length(mix_name))
  names(w_qg_all) <- mix_name
  if (!is.null(qgcomp_res$pos.weights)) {
    for (nm in names(qgcomp_res$pos.weights)) w_qg_all[nm] <- w_qg_all[nm] + abs(qgcomp_res$pos.weights[nm])
  }
  if (!is.null(qgcomp_res$neg.weights)) {
    for (nm in names(qgcomp_res$neg.weights)) w_qg_all[nm] <- w_qg_all[nm] + abs(qgcomp_res$neg.weights[nm])
  }
  w_qgcomp <- safe_normalize(w_qg_all, mix_name)

  # QGcomp: predict on full data
  qg_pred <- predict(qgcomp_res)
  qg_resid <- Y - qg_pred
  eval_qgcomp <- list(
    R2   = if (sum((Y - mean(Y))^2) > 0) 1 - sum(qg_resid^2) / sum((Y - mean(Y))^2) else NA_real_,
    RMSE = sqrt(mean(qg_resid^2))
  )

  # ── Weight Error Metrics ──────────────────────────────────────────
  err_nwqs <- calc_weight_error(w_nwqs, w_true)
  err_gwqs <- calc_weight_error(w_gwqs, w_true)
  err_qgcomp <- calc_weight_error(w_qgcomp, w_true)

  # ── NWQS Bootstrap CI Coverage (self-calibration only) ────────────
  cov_df <- tryCatch(
    {
      ci_sub <- nwqs_res$ci_table
      true_long <- as.data.frame(as.table(true_eff_mat))
      colnames(true_long) <- c("Term", "Target", "True_Value")
      true_long$Term <- as.character(true_long$Term)
      true_long$Target <- as.character(true_long$Target)

      merged <- merge(ci_sub, true_long, by = c("Target", "Term"), all.x = TRUE)
      merged$Bias <- merged$Estimate - merged$True_Value
      merged$Covered <- !is.na(merged$Boot_CI_Lower) & !is.na(merged$Boot_CI_Upper) &
        merged$Boot_CI_Lower <= merged$True_Value & merged$Boot_CI_Upper >= merged$True_Value
      merged$Sim_ID <- sim_id
      merged
    },
    error = function(e) NULL
  )

  # ── Assemble All Results ──────────────────────────────────────────
  models <- c("NWQS", "gWQS", "QGcomp")

  # Per-model performance
  df_perf <- data.frame(
    Sim_ID = sim_id,
    Model = models,
    R2 = c(eval_nwqs$R2, eval_gwqs$R2, eval_qgcomp$R2),
    RMSE = c(eval_nwqs$RMSE, eval_gwqs$RMSE, eval_qgcomp$RMSE),
    SAE = c(err_nwqs$SAE, err_gwqs$SAE, err_qgcomp$SAE),
    Pearson = c(err_nwqs$Pearson, err_gwqs$Pearson, err_qgcomp$Pearson),
    Spearman = c(err_nwqs$Spearman, err_gwqs$Spearman, err_qgcomp$Spearman)
  )

  # Per-component weights
  df_weights <- do.call(rbind, lapply(seq_along(models), function(k) {
    w_est <- list(w_nwqs, w_gwqs, w_qgcomp)[[k]]
    data.frame(
      Sim_ID = sim_id, Model = models[k], Component = mix_name,
      True_Value = as.numeric(w_true[mix_name]),
      Estimated_Weight = as.numeric(w_est[mix_name])
    )
  }))

  list(perf = df_perf, weights = df_weights, coverage = cov_df)
}


# -------------------------------------------------------------------------
# 3. Scenario Definitions
# -------------------------------------------------------------------------
w_norm_4 <- c(0.10, 0.20, 0.30, 0.40)
w_sparse_4 <- c(0.60, 0.40, 0.00, 0.00)
w_norm_8 <- c(0.04, 0.06, 0.08, 0.10, 0.14, 0.16, 0.18, 0.24)
w_norm_12 <- c(0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.11, 0.12, 0.13, 0.13, 0.13)
w_sparse_12 <- c(0.40, 0.30, 0.20, 0.10, rep(0, 8))

base_settings <- list(
  "S1_Base" = list(P = 4, corr = "mixed", shape = "threshold", w = w_norm_4, snr_db = 10),
  "S2_Linear" = list(P = 4, corr = "mixed", shape = "pure_linear", w = w_norm_4, snr_db = 10),
  "S3_HighCorr" = list(P = 4, corr = "high", shape = "threshold", w = w_norm_4, snr_db = 10),
  "S4_ComplexShape" = list(
    P = 4, corr = "mixed",
    shape = c("u_shape", "inv_threshold", "pure_linear", "s_shape"),
    w = w_norm_4, snr_db = 10
  ),
  "S5_HighDim" = list(P = 12, corr = "mixed", shape = "threshold", w = w_norm_12, snr_db = 10),
  "S6_Sparse_HighDim" = list(P = 12, corr = "high", shape = "threshold", w = w_sparse_12, snr_db = 10),
  "S7_LowSNR" = list(P = 4, corr = "mixed", shape = "threshold", w = w_norm_4, snr_db = 2)
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
# 4. Main Loop
# -------------------------------------------------------------------------
total_tasks <- length(scenarios)
current_task <- 0

message(sprintf(
  "\n=== MC Study: %s | %d scenarios x %d sims ===\n",
  toupper(TARGET_FAMILY), total_tasks, N_SIMULATIONS
))

for (scen_name in names(scenarios)) {
  current_task <- current_task + 1
  params <- scenarios[[scen_name]]
  mix_name <- paste0("Component", seq_len(params$P))

  base_scen <- gsub("_N[0-9]+$", "", scen_name)
  final_out_dir <- file.path(out_root_dir, base_scen, scen_name)
  if (!dir.exists(final_out_dir)) dir.create(final_out_dir, recursive = TRUE)

  # Skip if done (both PNG + PDF exist)
  target_png <- file.path(final_out_dir, sprintf("MC_Benchmark_%s_%s.png", TARGET_FAMILY, scen_name))
  target_pdf <- file.path(final_out_dir, sprintf("MC_Benchmark_%s_%s.pdf", TARGET_FAMILY, scen_name))
  if (file.exists(target_png) && file.exists(target_pdf)) {
    message(sprintf("[SKIP] %d/%d %s", current_task, total_tasks, scen_name))
    next
  }

  message(sprintf(
    "\n[%d/%d] %s (N=%d, P=%d) ...",
    current_task, total_tasks, scen_name, params$N, params$P
  ))

  scen_db_dir <- file.path(db_root_dir, scen_name)
  if (!dir.exists(scen_db_dir)) {
    message("  DB not found, skip.")
    next
  }

  sim_files <- file.path(scen_db_dir, sprintf("sim_%03d.rds", 1:N_SIMULATIONS))
  sim_files <- sim_files[file.exists(sim_files)]
  if (length(sim_files) == 0) {
    message("  No sim files, skip.")
    next
  }

  start_time <- Sys.time()

  # ── Parallel execution ─────────────────────────────────────────────
  sim_results <- future.apply::future_lapply(
    sim_files,
    function(f) {
      tryCatch(
        run_one_sim(f, family = TARGET_FAMILY, mix_name = mix_name),
        error = function(e) list(IS_ERROR = TRUE, msg = conditionMessage(e))
      )
    },
    future.seed = TRUE,
    future.packages = c("NWQS", "gWQS", "qgcomp", "dplyr", "tidyr")
  )

  valid <- Filter(function(x) is.null(x$IS_ERROR), sim_results)
  n_fail <- length(sim_files) - length(valid)
  if (n_fail > 0) message(sprintf("  %d/%d failed.", n_fail, length(sim_files)))
  if (length(valid) == 0) {
    message("  All failed, skip.")
    next
  }

  # ── Aggregate ──────────────────────────────────────────────────────
  df_perf <- do.call(rbind, lapply(valid, `[[`, "perf"))
  df_weights <- do.call(rbind, lapply(valid, `[[`, "weights"))
  df_coverage <- do.call(rbind, Filter(Negate(is.null), lapply(valid, `[[`, "coverage")))

  # ── Save Raw CSVs ──────────────────────────────────────────────────
  write.csv(df_perf, file.path(final_out_dir, paste0("Raw_Performance_", scen_name, ".csv")), row.names = FALSE)
  write.csv(df_weights, file.path(final_out_dir, paste0("Raw_Weights_", scen_name, ".csv")), row.names = FALSE)
  if (!is.null(df_coverage) && nrow(df_coverage) > 0) {
    write.csv(df_coverage, file.path(final_out_dir, paste0("Raw_Coverage_NWQS_", scen_name, ".csv")), row.names = FALSE)
  }

  # ── Summary Table 1: Model Performance ─────────────────────────────
  summary_perf <- df_perf %>%
    group_by(Model) %>%
    summarise(
      Mean_R2 = mean(R2, na.rm = TRUE),
      SD_R2 = sd(R2, na.rm = TRUE),
      Mean_RMSE = mean(RMSE, na.rm = TRUE),
      SD_RMSE = sd(RMSE, na.rm = TRUE),
      Mean_SAE = mean(SAE, na.rm = TRUE),
      SD_SAE = sd(SAE, na.rm = TRUE),
      Pearson_Pass = mean(Pearson >= 0.8, na.rm = TRUE), # proportion with r ≥ 0.8
      Spearman_Pass = mean(Spearman >= 0.8, na.rm = TRUE), # proportion with ρ ≥ 0.8
      Mean_Pearson = mean(Pearson, na.rm = TRUE),
      Mean_Spearman = mean(Spearman, na.rm = TRUE),
      N_Valid = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(Mean_R2))

  write.csv(summary_perf, file.path(final_out_dir, paste0("Summary_Performance_", scen_name, ".csv")), row.names = FALSE)

  cat(sprintf("\n  --- %s Performance Summary ---\n", scen_name))
  print(as.data.frame(summary_perf), row.names = FALSE)

  # ── Summary Table 2: Weight Recovery per Component ─────────────────
  summary_weights <- df_weights %>%
    group_by(Model, Component) %>%
    summarise(
      True_Value = mean(True_Value, na.rm = TRUE),
      Mean_Weight = mean(Estimated_Weight, na.rm = TRUE),
      SD_Weight = sd(Estimated_Weight, na.rm = TRUE),
      Bias = mean(Estimated_Weight - True_Value, na.rm = TRUE),
      CI_Lower = quantile(Estimated_Weight, 0.025, na.rm = TRUE),
      CI_Upper = quantile(Estimated_Weight, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  write.csv(summary_weights, file.path(final_out_dir, paste0("Summary_Weights_", scen_name, ".csv")), row.names = FALSE)

  # ── Summary Table 3: NWQS Bootstrap CI Coverage (self-calibration) ──
  if (!is.null(df_coverage) && nrow(df_coverage) > 0) {
    summary_cp <- df_coverage %>%
      group_by(Target, Term) %>%
      summarise(
        N_Valid = n(),
        Coverage_Prob = mean(Covered, na.rm = TRUE),
        Mean_Bias = mean(Bias, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(Target, factor(Term, levels = c("Overall", mix_name)))
    write.csv(summary_cp, file.path(final_out_dir, paste0("Summary_Coverage_NWQS_", scen_name, ".csv")), row.names = FALSE)

    overall_cp <- summary_cp %>% filter(Term == "Overall", grepl("Q[0-9]+_vs_Q1", Target))
    if (nrow(overall_cp) > 0) {
      message(sprintf(
        "  NWQS Coverage (Overall, max Q vs Q1): %.1f%%",
        100 * tail(overall_cp$Coverage_Prob, 1)
      ))
    }
  }

  # ── Benchmark Plot (PNG + PDF) ─────────────────────────────────────
  tryCatch(
    {
      dynamic_nrow <- ceiling(params$P / 7)
      custom_palette <- c("NWQS" = "#4A90C8", "gWQS" = "#D92828", "QGcomp" = "#6EC44A")
      df_perf$Model <- factor(df_perf$Model, levels = c("NWQS", "gWQS", "QGcomp"))
      df_weights$Model <- factor(df_weights$Model, levels = c("NWQS", "gWQS", "QGcomp"))

      # Shared raincloud builder
      make_rain <- function(data, yvar, title, subtitle, ylab) {
        ggplot(data, aes(x = Model, y = .data[[yvar]], fill = Model, color = Model)) +
          ggdist::stat_halfeye(
            adjust = 0.5, width = 0.6, .width = 0,
            justification = -0.3, point_colour = NA, alpha = 0.7
          ) +
          geom_boxplot(
            width = 0.2, outlier.shape = NA, alpha = 0.5, color = "black",
            position = position_nudge(x = -0.1)
          ) +
          geom_point(size = 0.6, alpha = 0.25, position = position_jitter(width = 0.04)) +
          scale_fill_manual(values = custom_palette) +
          scale_color_manual(values = custom_palette) +
          theme_bw(base_size = 14) +
          theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 20, hjust = 1, face = "bold"),
            panel.grid.minor = element_blank()
          ) +
          labs(title = title, subtitle = subtitle, x = "", y = ylab)
      }

      # Row 1: Three panels
      p_r2 <- make_rain(df_perf, "R2", "A. Prediction Accuracy (R\u00B2)", "Higher is better", "R\u00B2")
      p_sae <- make_rain(df_perf, "SAE", "B. Weight Recovery Error (SAE)", "Lower is better", "SAE")

      # Panel C: Pearson + Spearman faceted
      df_corr <- df_perf %>%
        select(Sim_ID, Model, Pearson, Spearman) %>%
        tidyr::pivot_longer(cols = c(Pearson, Spearman), names_to = "Metric", values_to = "Correlation") %>%
        mutate(Metric = factor(Metric, levels = c("Pearson", "Spearman")))

      p_corr <- ggplot(df_corr, aes(x = Model, y = Correlation, fill = Model)) +
        geom_boxplot(alpha = 0.8, outlier.size = 0.5, color = "black", width = 0.6) +
        geom_jitter(size = 0.4, alpha = 0.2, width = 0.12, color = "gray30") +
        geom_hline(yintercept = 0.95, linetype = "dashed", color = "#2C3E50", linewidth = 0.5) +
        facet_wrap(~Metric) +
        scale_fill_manual(values = custom_palette) +
        theme_bw(base_size = 14) +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 20, hjust = 1, face = "bold"),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "#ECF0F1"),
          strip.text = element_text(face = "bold", size = 12)
        ) +
        labs(
          title = "C. Weight-Truth Correlation",
          subtitle = "Pearson: absolute | Spearman: ranking | Dashed = 0.95",
          x = "", y = "Correlation"
        )

      # Row 2: Weight recovery faceted boxplot
      p_weights <- ggplot(df_weights, aes(x = Model, y = Estimated_Weight, fill = Model)) +
        geom_boxplot(alpha = 0.8, outlier.size = 0.5, color = "black", width = 0.6) +
        geom_hline(aes(yintercept = True_Value), linetype = "dashed", color = "black", linewidth = 1) +
        facet_wrap(~Component, scales = "free_y", nrow = dynamic_nrow) +
        scale_fill_manual(values = custom_palette) +
        theme_bw(base_size = 14) +
        theme(
          legend.position = "bottom", legend.title = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.background = element_rect(fill = "#ECF0F1"),
          strip.text = element_text(face = "bold")
        ) +
        labs(
          title = "D. Component-Specific Weight Recovery",
          subtitle = "Dashed line = true weight", x = "", y = "Estimated Weight"
        )

      dynamic_height <- 0.8 * dynamic_nrow
      row1 <- p_r2 + p_sae + p_corr + patchwork::plot_layout(widths = c(1, 1, 1.4))
      mc_plot <- row1 / p_weights +
        patchwork::plot_layout(heights = c(1, dynamic_height)) +
        patchwork::plot_annotation(
          title = sprintf("Monte Carlo Benchmark: %s | %s", scen_name, toupper(TARGET_FAMILY)),
          theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
        )

      ggsave(target_png, mc_plot, width = 18, height = 11 + (dynamic_nrow - 1) * 3.5, dpi = 300)
      ggsave(target_pdf, mc_plot, width = 18, height = 11 + (dynamic_nrow - 1) * 3.5, device = "pdf")
    },
    error = function(e) message(sprintf("  Plot failed: %s", e$message))
  )

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  message(sprintf("  Done: %.1f mins (%d/%d valid)", elapsed, length(valid), length(sim_files)))
  gc(verbose = FALSE)
}

# -------------------------------------------------------------------------
# 5. Cleanup
# -------------------------------------------------------------------------
future::plan(future::sequential)
message(sprintf("\n=== MC Study Complete! Results: %s ===", out_root_dir))
