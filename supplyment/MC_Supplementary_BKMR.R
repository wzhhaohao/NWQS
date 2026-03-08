# =========================================================================
# NWQS Monte Carlo — Supplementary: BKMR Comparison (Appendix)
#
# BKMR is run separately due to extreme computational cost.
# Only runs on scenarios with N ≤ 500 and P ≤ 8.
# Results are saved independently and can be merged with main analysis.
# =========================================================================

rm(list = ls())

library(NWQS)
library(bkmr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(future)
library(future.apply)

# ── Settings ─────────────────────────────────────────────────────────────
TARGET_FAMILY  <- "gaussian"
N_SIMULATIONS  <- 100
BKMR_ITER      <- 10000             # MCMC iterations (5000-10000 recommended)
BKMR_MAX_N     <- 500               # Skip if N > this
BKMR_MAX_P     <- 8                 # Skip if P > this
WORKER_CORES   <- 4                 # Fewer workers (BKMR is memory-hungry)

PROJECT_ROOT <- '/Users/wangzhehao/temporary/packages/NWQS result'
db_root_dir  <- file.path(PROJECT_ROOT, "results", "Monte_Carlo_DB", toupper(TARGET_FAMILY))
out_root_dir <- file.path(PROJECT_ROOT, "results", "Monte_Carlo_Results_BKMR", toupper(TARGET_FAMILY))

if (!dir.exists(db_root_dir)) stop("DB not found")
if (!dir.exists(out_root_dir)) dir.create(out_root_dir, recursive = TRUE)

plan(multisession, workers = WORKER_CORES)

# ── Utility ──────────────────────────────────────────────────────────────
safe_normalize <- function(x, nm = names(x)) {
  x <- as.numeric(x); x[!is.finite(x) | x < 0] <- 0
  out <- if (sum(x) <= 0) rep(1/length(x), length(x)) else x/sum(x)
  if (!is.null(nm)) names(out) <- nm; out
}

# ── Single BKMR Simulation ──────────────────────────────────────────────
run_bkmr_sim <- function(sim_file, family, mix_name) {
  sim_obj  <- readRDS(sim_file)
  sim_data <- sim_obj$data
  sim_id   <- sim_obj$sim_id
  true_eff <- sim_obj$true_effect_mat
  Y        <- sim_data$y
  w_true   <- calc_true_importance(true_eff, mix_name)

  transform_fun <- function(x) trans_quantile(x, q = 4)
  data_q <- sim_data
  data_q[mix_name] <- transform_fun(data_q[mix_name])

  Z <- as.matrix(data_q[, mix_name])
  X <- model.matrix(~ x_cont + x_bin + x_cat, data = sim_data)[, -1]

  # Fit BKMR
  bkmr_fit <- bkmr::kmbayes(
    y = Y, Z = Z, X = X,
    iter = BKMR_ITER, verbose = FALSE, varsel = TRUE
  )

  # Extract PIPs as importance weights
  pips <- bkmr::ExtractPIPs(bkmr_fit)
  w_bkmr_raw <- pips$PIP
  names(w_bkmr_raw) <- pips$variable
  w_bkmr <- safe_normalize(w_bkmr_raw[mix_name], mix_name)

  # Prediction performance
  bkmr_pred <- rowMeans(bkmr_fit$ypred)  # posterior mean predictions
  resid <- Y - bkmr_pred
  rmse <- sqrt(mean(resid^2))
  r2 <- if (sum((Y - mean(Y))^2) > 0) 1 - sum(resid^2) / sum((Y - mean(Y))^2) else NA

  err <- calc_weight_error(w_bkmr, w_true)

  df_perf <- data.frame(
    Sim_ID = sim_id, Model = "BKMR",
    R2 = r2, RMSE = rmse, SAE = err$SAE,
    Pearson = err$Pearson, Spearman = err$Spearman
  )

  df_weights <- data.frame(
    Sim_ID = sim_id, Model = "BKMR", Component = mix_name,
    True_Value = as.numeric(w_true[mix_name]),
    Estimated_Weight = as.numeric(w_bkmr[mix_name])
  )

  list(perf = df_perf, weights = df_weights)
}

# ── Scenario Definitions (subset eligible for BKMR) ─────────────────────
w_norm_4  <- c(0.10, 0.20, 0.30, 0.40)
w_norm_8  <- c(0.04, 0.06, 0.08, 0.10, 0.14, 0.16, 0.18, 0.24)

base_settings <- list(
  "S1_Base"         = list(P=4, corr="mixed", shape="threshold", w=w_norm_4, snr_db=10),
  "S2_Linear"       = list(P=4, corr="mixed", shape="pure_linear", w=w_norm_4, snr_db=10),
  "S3_HighCorr"     = list(P=4, corr="high", shape="threshold", w=w_norm_4, snr_db=10),
  "S4_ComplexShape"  = list(P=4, corr="mixed",
                            shape=c("u_shape","inv_threshold","pure_linear","s_shape"),
                            w=w_norm_4, snr_db=10),
  "S7_LowSNR"       = list(P=4, corr="mixed", shape="threshold", w=w_norm_4, snr_db=2)
)

N_values <- c(200, 500)  # BKMR only on small/medium N
scenarios <- list()
for (b_name in names(base_settings)) {
  for (n in N_values) {
    full_name <- sprintf("%s_N%d", b_name, n)
    curr <- base_settings[[b_name]]
    curr$N <- n
    if (curr$N <= BKMR_MAX_N && curr$P <= BKMR_MAX_P)
      scenarios[[full_name]] <- curr
  }
}

# ── Main Loop ────────────────────────────────────────────────────────────
total <- length(scenarios); idx <- 0

for (scen_name in names(scenarios)) {
  idx <- idx + 1
  params <- scenarios[[scen_name]]
  mix_name <- paste0("Component", seq_len(params$P))

  final_out_dir <- file.path(out_root_dir, gsub("_N[0-9]+$", "", scen_name), scen_name)
  if (!dir.exists(final_out_dir)) dir.create(final_out_dir, recursive = TRUE)

  done_file <- file.path(final_out_dir, paste0("BKMR_Performance_", scen_name, ".csv"))
  if (file.exists(done_file)) { message(sprintf("[SKIP] %d/%d %s", idx, total, scen_name)); next }

  message(sprintf("\n[%d/%d] BKMR: %s (N=%d, P=%d, iter=%d) ...",
                  idx, total, scen_name, params$N, params$P, BKMR_ITER))

  scen_db_dir <- file.path(db_root_dir, scen_name)
  sim_files <- file.path(scen_db_dir, sprintf("sim_%03d.rds", 1:N_SIMULATIONS))
  sim_files <- sim_files[file.exists(sim_files)]
  if (length(sim_files) == 0) { message("  No data, skip."); next }

  start_time <- Sys.time()

  results <- future.apply::future_lapply(
    sim_files,
    function(f) tryCatch(run_bkmr_sim(f, TARGET_FAMILY, mix_name),
                         error = function(e) list(IS_ERROR = TRUE)),
    future.seed = TRUE, future.packages = c("NWQS", "bkmr", "dplyr")
  )

  valid <- Filter(function(x) is.null(x$IS_ERROR), results)
  message(sprintf("  %d/%d succeeded.", length(valid), length(sim_files)))
  if (length(valid) == 0) next

  df_perf    <- do.call(rbind, lapply(valid, `[[`, "perf"))
  df_weights <- do.call(rbind, lapply(valid, `[[`, "weights"))

  write.csv(df_perf,    done_file, row.names = FALSE)
  write.csv(df_weights, file.path(final_out_dir, paste0("BKMR_Weights_", scen_name, ".csv")), row.names = FALSE)

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  message(sprintf("  Done: %.1f mins", elapsed))
  gc(verbose = FALSE)
}

future::plan(future::sequential)
message("\n=== BKMR Supplementary Analysis Complete ===")
