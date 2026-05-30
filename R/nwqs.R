#' @title Fit a Non-Linear Weighted Quantile Sum (NWQS) Regression Model
#'
#' @description
#' Core function for fitting the Non-linear Weighted Quantile Sum (NWQS)
#' regression model. Assesses the overall joint effect and relative importance
#' of highly correlated mixture exposures while flexibly accommodating
#' non-linear dose-response relationships (e.g., threshold, U-shaped, S-shaped
#' curves). The algorithm strictly separates penalized representation learning
#' (weight and shape discovery) from unpenalized final effect estimation.
#'
#' @details
#' The function implements a robust Repeated Holdout (RH) architecture with
#' shape decoupling:
#' \enumerate{
#'   \item \strong{Outer RH splitting:} Data is randomly split into training
#'     (weight/shape discovery) and validation (effect estimation) sets.
#'   \item \strong{Weight and shape discovery:} A multi-parameter engine
#'     estimates spline basis coefficients (shapes) and derives relative
#'     importance (weights) via OOB permutation.
#'   \item \strong{Shape normalization:} Shapes are standardized to unit
#'     variance on the training set to decouple shape from effect size.
#'   \item \strong{1-DoF effect estimation:} Normalized shapes and weights are
#'     projected onto the validation set to construct a single \code{nwqs}
#'     index, then a GLM estimates the overall mixture effect.
#' }
#'
#' \strong{Important warnings:}
#' \itemize{
#'   \item When \code{rh > 1}, the standard errors and CIs in
#'     \code{fit$coefficients} reflect only data-splitting (algorithmic)
#'     variance, NOT true sampling variance. Use \code{\link{nwqs_boot}} for
#'     valid inference.
#'   \item Ensure adequate confounders are specified via \code{covariates} to
#'     minimize residual confounding.
#' }
#'
#' @param data \code{data.frame}. Contains mixture components, covariates,
#'   matching variables, and the outcome.
#' @param mix_name Character vector. Column names of mixture components.
#' @param covariates Character vector. Covariate/confounder column names. If
#'   \code{NULL}, fits an unadjusted model.
#' @param outcome Character. Outcome variable column name. Default is
#'   \code{"y"}.
#' @param weight_engine Function. Engine for weight and shape discovery.
#'   Default is \code{permutation_scorer}.
#' @param min_shape_sd Numeric. When a component's training-set partial
#'   linear-predictor standard deviation falls below this threshold, the
#'   per-component shape normalization is bypassed (the component is
#'   treated as carrying no information that iteration). Defaults to
#'   \code{1e-8}. When the threshold fires and \code{quiet = FALSE}, a
#'   one-line \code{message()} reports which component and which RH
#'   iteration was affected.
#' @param transform_type Character. Either \code{"percentile_rank"} (default;
#'   continuous empirical CDF on each mixture column) or \code{"q_bin"}
#'   (legacy 0.1.x behavior: discrete quantile bins).
#' @param q Integer. With \code{transform_type = "q_bin"} this is the number
#'   of discrete bins; with \code{transform_type = "percentile_rank"} it
#'   controls the number of contrast points used downstream by
#'   \code{extract_nwqs_effects()} and \code{nwqs_contrast()}. Default is 4.
#' @param ties Character. Tie-handling rule for \code{"percentile_rank"};
#'   passed to \code{rank(ties.method = )}. One of \code{"average"} (default),
#'   \code{"min"}, \code{"max"}, \code{"random"}.
#' @param df_spline Integer. Degrees of freedom for natural cubic splines.
#'   Default is 3.
#' @param transform_fun Function. Custom transformation for mixture components.
#'   If \code{NULL}, uses \code{trans_quantile}.
#' @param train_prop Numeric in (0, 1). Training set proportion per RH
#'   iteration. Default is 0.6.
#' @param rh Integer. Number of repeated holdout iterations. Default is 10.
#' @param seed Integer. Random seed for reproducible splits. Default is 1234.
#' @param n_permutation Integer. Number of internal permutations for variable
#'   importance. Default is 30 (raised from 10 in 0.2.0 for more stable
#'   weight estimates on small samples).
#' @param family Character. Error distribution: \code{"gaussian"},
#'   \code{"binomial"}, \code{"poisson"}, \code{"quasipoisson"}, or
#'   \code{"negbin"}. \code{"negbin"} fits the final NWQS regression
#'   with \code{MASS::glm.nb()}; the inner OOB permutation loss uses a
#'   Poisson surrogate for speed (the importance ranking is unchanged).
#' @param plan_strategy Character. Parallel strategy: \code{"sequential"},
#'   \code{"multicore"}, or \code{"multisession"}.
#' @param n_workers Integer. Number of parallel workers. If \code{NULL},
#'   auto-detected.
#' @param control Object of class \code{nwqs_control} returned by
#'   \code{\link{nwqs_control}()}. Bundles advanced knobs that do not live
#'   on this signature: \code{custom_knots}, \code{custom_boundary},
#'   \code{zero_weight_action}.
#' @param quiet Logical. If \code{TRUE}, suppresses RH variance warnings.
#'   Typically used when called internally by \code{nwqs_boot()}.
#' @param ... Additional arguments passed to \code{run_oob_permutation} or
#'   \code{weight_engine}.
#'
#' @return An object of class \code{c("nwqs", "list")} containing:
#' \itemize{
#'   \item \code{fit}: GLM summary (coefficients, deviance, AIC).
#'   \item \code{final_weights}: Ensemble-averaged relative weights.
#'   \item \code{mean_shapes}: Ensemble-averaged normalized spline coefficients.
#'   \item \code{rh_coefs}: Matrix of GLM coefficients across RH iterations.
#'   \item \code{rh_weights}: Matrix of weights across RH iterations.
#'   \item \code{data}: Original data with the computed \code{nwqs} index.
#' }
#'
#' @seealso \code{\link{nwqs_boot}}, \code{\link{plot.nwqs}},
#'   \code{\link{summary.nwqs}}
#'
#' @importFrom stats glm coef AIC pnorm sd median as.formula
#' @importFrom splines ns
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#' @export
nwqs <- function(data, mix_name, covariates = NULL, outcome = "y",
                 weight_engine = permutation_scorer,
                 transform_type = c("percentile_rank", "q_bin"),
                 q = NWQS_DEFAULTS$q,
                 ties = c("average", "min", "max", "random"),
                 df_spline = NWQS_DEFAULTS$df_spline,
                 min_shape_sd = NWQS_DEFAULTS$min_shape_sd,
                 transform_fun = NULL,
                 train_prop = NWQS_DEFAULTS$train_prop,
                 rh = NWQS_DEFAULTS$rh,
                 seed = NWQS_DEFAULTS$seed,
                 n_permutation = NWQS_DEFAULTS$n_permutation,
                 family = c("gaussian", "binomial", "poisson", "quasipoisson", "negbin"),
                 plan_strategy = c("sequential", "multisession", "multicore"),
                 n_workers = NULL,
                 control = nwqs_control(),
                 quiet = FALSE, ...) {
  family <- match.arg(family)
  plan_strategy <- match.arg(plan_strategy)
  transform_type <- match.arg(transform_type)
  ties <- match.arg(ties)
  if (length(covariates) == 0) covariates <- NULL

  t_start <- Sys.time()
  args <- list(...)

  if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
    stop("Please install 'future' and 'future.apply' packages.")
  }
  if (rh < 1) stop("'rh' must be at least 1.")
  if (train_prop <= 0 || train_prop >= 1) stop("'train_prop' must be in (0, 1).")

  current_reserve <- if (!is.null(args$cpu_reserve)) args$cpu_reserve else 0.2

  old_plan <- tryCatch(
    configure_parallel_plan(
      loop_number = rh, strategy = plan_strategy,
      n_workers = n_workers, reserve_cpu = current_reserve,
      verbose = !isTRUE(quiet)
    ),
    error = function(e) {
      warning(
        "Failed to configure parallel plan: ", conditionMessage(e),
        ". Falling back to sequential."
      )
      NULL
    }
  )
  on.exit(
    {
      if (!is.null(old_plan)) future::plan(old_plan)
    },
    add = TRUE
  )

  use_parallel <- !inherits(future::plan(), "sequential")

  if (is.null(transform_fun)) {
    transform_fun <- function(x) {
      trans_quantile(x, type = transform_type, q = q, ties = ties)
    }
  }

  train_components_sorted <- lapply(mix_name, function(comp) {
    sort(data[[comp]])
  })
  names(train_components_sorted) <- mix_name

  data_Q <- data
  data_Q[mix_name] <- transform_fun(data[mix_name])

  basis_info <- build_spline_basis_knots(
    transform_type  = transform_type,
    q               = q,
    df_spline       = df_spline,
    custom_knots    = control$custom_knots,
    custom_boundary = control$custom_boundary
  )
  model_knots <- basis_info$knots
  model_boundary <- basis_info$boundary

  if (!use_parallel && !is.null(seed)) set.seed(seed)

  n_obs <- nrow(data)

  if (is.null(covariates)) {
    formula_str <- paste(outcome, "~ nwqs")
  } else {
    missing_cov <- setdiff(covariates, names(data))
    if (length(missing_cov) > 0) {
      stop(paste("Missing covariates:", paste(missing_cov, collapse = ", ")))
    }
    formula_str <- paste(outcome, "~ nwqs +", paste(covariates, collapse = " + "))
  }
  formula_final <- as.formula(formula_str)

  one_rh <- function(i) {
    idx_all <- sample(seq_len(n_obs))
    n_train <- floor(n_obs * train_prop)
    train_idx <- idx_all[seq_len(n_train)]
    valid_idx <- idx_all[(n_train + 1):n_obs]

    data_train <- data_Q[train_idx, , drop = FALSE]
    data_valid <- data_Q[valid_idx, , drop = FALSE]

    vars_needed <- unique(c(mix_name, outcome, covariates))

    boot_res <- run_oob_permutation(
      data = data_train[, vars_needed, drop = FALSE],
      mix_name = mix_name,
      outcome = outcome,
      covariates = covariates,
      weight_engine = weight_engine,
      n_permutation = n_permutation,
      q = q,
      df_spline = df_spline,
      model_knots = model_knots,
      model_boundary = model_boundary,
      family = family,
      boot_strategy = "sequential",
      ...
    )
    valid_res <- Filter(Negate(is.null), boot_res)
    if (length(valid_res) == 0) {
      return(NULL)
    }

    w_matrix_iter <- do.call(rbind, lapply(valid_res, function(x) x$weights))
    s_matrix_iter <- do.call(rbind, lapply(valid_res, function(x) x$shapes))

    mean_weights_iter <- colMeans(w_matrix_iter, na.rm = TRUE)
    mean_shapes_iter <- colMeans(s_matrix_iter, na.rm = TRUE)

    if (!all(is.finite(mean_weights_iter)) || sum(mean_weights_iter, na.rm = TRUE) <= 0) {
      return(NULL)
    }

    final_weights_iter <- mean_weights_iter / sum(mean_weights_iter)

    train_trans <- wqs_nonlinear_expand(
      data_train, mix_name,
      knots = model_knots, boundary = model_boundary
    )
    valid_trans <- wqs_nonlinear_expand(
      data_valid, mix_name,
      knots = model_knots, boundary = model_boundary
    )

    expected_cols <- names(mean_shapes_iter)
    actual_train_cols <- colnames(train_trans)
    actual_valid_cols <- colnames(valid_trans)
    if (!all(expected_cols %in% actual_train_cols) ||
      !all(expected_cols %in% actual_valid_cols)) {
      warning(
        "Column mismatch in wqs_nonlinear_expand() output at RH iteration ", i,
        ". Skipping this iteration."
      )
      return(NULL)
    }

    normalized_shapes_iter <- numeric(length(mean_shapes_iter))
    names(normalized_shapes_iter) <- names(mean_shapes_iter)

    combined_coefs <- numeric(length(mean_shapes_iter))
    names(combined_coefs) <- names(mean_shapes_iter)

    for (comp in mix_name) {
      comp_cols <- paste0(comp, "_B", 1:df_spline)
      theta_raw <- mean_shapes_iter[comp_cols]

      partial_eta <- as.vector(as.matrix(train_trans[, comp_cols, drop = FALSE]) %*% theta_raw)
      sd_eta <- sd(partial_eta, na.rm = TRUE)
      if (is.na(sd_eta) || sd_eta < min_shape_sd) {
        if (!isTRUE(quiet)) {
          message(sprintf(
            "nwqs: component '%s' degenerate at RH iteration %d (sd_eta < %g); shape normalization skipped.",
            comp, i, min_shape_sd
          ))
        }
        sd_eta <- 1
      }

      theta_norm <- theta_raw / sd_eta
      normalized_shapes_iter[comp_cols] <- theta_norm
      combined_coefs[comp_cols] <- theta_norm * final_weights_iter[comp]
    }

    nwqs <- as.matrix(valid_trans[, expected_cols, drop = FALSE]) %*% combined_coefs
    data_valid$nwqs <- as.vector(nwqs)

    if (family == "negbin") {
      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("Please install the 'MASS' package to use family = 'negbin'.")
      }
      fit <- tryCatch(
        suppressWarnings(MASS::glm.nb(formula_final, data = data_valid)),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NULL)
    } else {
      fit <- glm(formula_final, data = data_valid, family = family)
    }
    aic_val <- if (family == "quasipoisson") NA_real_ else AIC(fit)
    null_dev <- fit$null.deviance
    res_dev <- fit$deviance
    df_n <- fit$df.null
    df_r <- fit$df.residual

    coefs_fit <- coef(fit)

    if (!("nwqs" %in% names(coefs_fit)) || !is.finite(unname(coefs_fit["nwqs"]))) {
      return(NULL)
    }

    if (!is.finite(stats::sd(data_valid$nwqs, na.rm = TRUE)) ||
      stats::sd(data_valid$nwqs, na.rm = TRUE) < 1e-8) {
      return(NULL)
    }

    list(
      fit_obj = fit, weights = final_weights_iter, shapes = normalized_shapes_iter,
      coefs = coefs_fit, aic = aic_val, null_dev = null_dev,
      res_dev = res_dev, df_null = df_n, df_res = df_r
    )
  }

  if (use_parallel) {
    rh_results <- future.apply::future_lapply(
      seq_len(rh), one_rh,
      future.seed = if (!is.null(seed)) seed else TRUE
    )
  } else {
    rh_results <- lapply(seq_len(rh), one_rh)
  }

  rh_results <- Filter(function(x) {
    !is.null(x) &&
      !is.null(x$coefs) &&
      "nwqs" %in% names(x$coefs) &&
      is.finite(unname(x$coefs["nwqs"]))
  }, rh_results)
  if (length(rh_results) == 0) stop("All iterations failed.")

  if (rh == 1) {
    final_w_global <- rh_results[[1]]$weights
    final_s_global <- rh_results[[1]]$shapes
  } else {
    weight_mat_temp <- do.call(rbind, lapply(rh_results, function(x) x$weights))
    mean_weights_temp <- colMeans(weight_mat_temp, na.rm = TRUE)
    final_w_global <- mean_weights_temp / sum(mean_weights_temp)

    shape_mat_temp <- do.call(rbind, lapply(rh_results, function(x) x$shapes))
    final_s_global <- colMeans(shape_mat_temp, na.rm = TRUE)
  }

  full_trans <- wqs_nonlinear_expand(data_Q, mix_name, knots = model_knots, boundary = model_boundary)

  expected_cols_full <- names(final_s_global)
  if (!all(expected_cols_full %in% colnames(full_trans))) {
    stop("Column mismatch between final shapes and wqs_nonlinear_expand() output on full data.")
  }

  combined_coefs_full <- numeric(length(final_s_global))
  names(combined_coefs_full) <- names(final_s_global)

  for (comp in mix_name) {
    comp_cols <- paste0(comp, "_B", 1:df_spline)
    combined_coefs_full[comp_cols] <- final_s_global[comp_cols] * final_w_global[comp]
  }

  final_data <- data
  final_data$nwqs <- as.vector(
    as.matrix(full_trans[, expected_cols_full, drop = FALSE]) %*% combined_coefs_full
  )

  if (rh == 1) {
    single_res <- rh_results[[1]]
    final_obj <- single_res$fit_obj
    coef_summary <- as.data.frame(summary(final_obj)$coefficients)

    fit_obj <- list(
      coefficients = coef_summary, aic = single_res$aic, deviance = single_res$res_dev,
      null.deviance = single_res$null_dev, df.residual = single_res$df_res, df.null = single_res$df_null
    )

    results <- list(
      fit = fit_obj, final_weights = single_res$weights, mean_coefs = single_res$coefs,
      mean_shapes = single_res$shapes, rh_coefs = t(as.matrix(single_res$coefs)),
      rh_weights = t(as.matrix(single_res$weights)), rh_shapes = t(as.matrix(single_res$shapes)),
      rh = 1, b = n_permutation, q = q, df_spline = df_spline, family = family,
      transform_type = transform_type, ties = ties,
      train_components_sorted = train_components_sorted,
      spline_knots = model_knots, spline_boundary = model_boundary,
      formula = formula_final, model_obj = single_res$fit_obj,
      call = match.call(), data = final_data
    )

    class(results) <- c("nwqs", "list")
    return(results)
  }

  coef_mat <- do.call(rbind, lapply(rh_results, function(x) x$coefs))
  weight_mat <- do.call(rbind, lapply(rh_results, function(x) x$weights))
  shape_mat <- do.call(rbind, lapply(rh_results, function(x) x$shapes))

  mean_coefs <- colMeans(coef_mat, na.rm = TRUE)
  mean_weights <- colMeans(weight_mat, na.rm = TRUE)
  mean_weights <- mean_weights / sum(mean_weights)
  mean_shapes <- colMeans(shape_mat, na.rm = TRUE)

  if (family == "quasipoisson") {
    mean_aic <- NA_real_
  } else {
    mean_aic <- mean(vapply(rh_results, function(x) x$aic, numeric(1)), na.rm = TRUE)
  }

  mean_null_dev <- mean(vapply(rh_results, function(x) x$null_dev, numeric(1)), na.rm = TRUE)
  mean_res_dev <- mean(vapply(rh_results, function(x) x$res_dev, numeric(1)), na.rm = TRUE)

  df_null <- rh_results[[1]]$df_null
  df_res <- rh_results[[1]]$df_res

  coef_mean <- colMeans(coef_mat, na.rm = TRUE)
  coef_sd <- apply(coef_mat, 2, sd, na.rm = TRUE)

  if (!quiet) {
    warning(
      "When rh > 1, the Standard Errors and 95% CIs in `fit$coefficients` ",
      "represent ONLY the data-splitting (algorithmic) variance across holdout ",
      "iterations, NOT true sampling variance. They are inherently too narrow ",
      "and should NOT be used for statistical inference. Please use `nwqs_boot()` ",
      "to obtain valid percentile bootstrap confidence intervals."
    )
  }

  z_value <- ifelse(coef_sd > 0, coef_mean / coef_sd, NA_real_)
  p_value <- ifelse(is.na(z_value), NA_real_, 2 * pnorm(-abs(z_value)))

  coef_summary <- data.frame(
    Estimate     = coef_mean,
    `Std. Error` = coef_sd,
    `z value`    = z_value,
    `Pr(>|z|)`   = p_value,
    `2.5 %`      = coef_mean - 1.96 * coef_sd,
    `97.5 %`     = coef_mean + 1.96 * coef_sd,
    check.names  = FALSE
  )

  fit_obj <- list(
    coefficients = coef_summary, aic = mean_aic, deviance = mean_res_dev,
    null.deviance = mean_null_dev, df.residual = df_res, df.null = df_null
  )

  results <- list(
    fit = fit_obj, final_weights = mean_weights, mean_coefs = mean_coefs, mean_shapes = mean_shapes,
    rh_coefs = coef_mat, rh_weights = weight_mat, rh_shapes = shape_mat,
    rh = rh, b = n_permutation, q = q, df_spline = df_spline, family = family, call = match.call(),
    transform_type = transform_type, ties = ties,
    train_components_sorted = train_components_sorted,
    transform_fun = transform_fun, data = final_data,
    spline_knots = model_knots, spline_boundary = model_boundary,
    formula = formula_final, model_obj = NULL
  )

  class(results) <- c("nwqs", "list")
  return(results)
}


#' @title Bootstrap Confidence Interval Estimation for NWQS
#'
#' @description
#' Performs external bootstrap resampling for the NWQS method to approximate
#' the true sampling variability of model parameters. Extracts effects for
#' all model terms and provides publication-quality confidence interval tables.
#'
#' @details
#' Bootstrap resampling uses simple random sampling with replacement at the
#' observation level.
#'
#' For exponential family models (\code{binomial}, \code{poisson},
#' \code{quasipoisson}), results are automatically reported as exponentiated
#' values (OR or RR) with corresponding CIs.
#'
#' @param data data.frame. Original dataset.
#' @param mix_name Character vector. Mixture component column names.
#' @param covariates Character vector or \code{NULL}. Covariates/confounders.
#' @param outcome Character. Outcome column name. Default is \code{"y"}.
#' @param family Character. Error distribution: \code{"gaussian"},
#'   \code{"binomial"}, \code{"poisson"}, \code{"quasipoisson"}, or
#'   \code{"negbin"}.
#' @param n_boot Integer. Number of bootstrap replicates. Default is 100.
#' @param rh_inner Integer. RH iterations per bootstrap replicate. Default
#'   is 1.
#' @param n_permutation Integer. Permutation count for variable importance.
#'   Default is 30 (raised from 10 in 0.2.0).
#' @param conf_level Numeric. Confidence level. Default is 0.95.
#' @param seed Integer or \code{NULL}. Random seed.
#' @param keep_fits Logical. Whether to retain all bootstrap model objects.
#'   Default is \code{FALSE}.
#' @param plan_strategy Character. Parallel strategy for the outer bootstrap
#'   loop.
#' @param n_workers Integer or \code{NULL}. Number of parallel workers.
#' @param control Object of class \code{nwqs_control} returned by
#'   \code{\link{nwqs_control}()}. Forwarded to the inner \code{nwqs()}
#'   call on every bootstrap iteration.
#' @param transform_type Character. Either \code{"percentile_rank"} (default)
#'   or \code{"q_bin"}. Forwarded to \code{nwqs()}; see its documentation for
#'   the precise mathematical contract.
#' @param q Integer. Number of quantile bins (when
#'   \code{transform_type = "q_bin"}) or contrast points (when
#'   \code{transform_type = "percentile_rank"}). Default is 4.
#' @param ties Character. Tie-handling rule for \code{"percentile_rank"};
#'   passed to \code{rank(ties.method = )}. One of \code{"average"} (default),
#'   \code{"min"}, \code{"max"}, \code{"random"}.
#' @param quiet Logical. If \code{TRUE}, suppresses verbose output. Default
#'   is \code{TRUE}.
#' @param ... Additional arguments passed to \code{nwqs()}.
#'
#' @return An object of class \code{c("nwqs_boot", "list")} containing:
#' \itemize{
#'   \item \code{ci_table}: Long-format CI table with bootstrap percentile
#'     intervals.
#'   \item \code{formatted_table}: Wide-format summary table for publication.
#'   \item \code{boot_table}: Raw estimates from all successful bootstrap
#'     iterations.
#'   \item \code{final_weights}: Ensemble-averaged component weights.
#'   \item \code{n_success}: Number of successful bootstrap iterations.
#' }
#'
#' @seealso \code{\link{nwqs}}, \code{\link{plot.nwqs_boot}},
#'   \code{\link{extract_nwqs_effects}}
#'
#' @importFrom stats aggregate quantile
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#' @export
nwqs_boot <- function(data,
                      mix_name,
                      covariates = NULL,
                      outcome = "y",
                      family = c("gaussian", "binomial", "poisson", "quasipoisson", "negbin"),
                      n_boot = NWQS_DEFAULTS$n_boot,
                      rh_inner = NWQS_DEFAULTS$rh_inner,
                      n_permutation = NWQS_DEFAULTS$n_permutation,
                      conf_level = NWQS_DEFAULTS$conf_level,
                      seed = NWQS_DEFAULTS$seed,
                      keep_fits = FALSE,
                      plan_strategy = c("sequential", "multisession", "multicore"),
                      n_workers = NULL,
                      transform_type = c("percentile_rank", "q_bin"),
                      q = NWQS_DEFAULTS$q,
                      ties = c("average", "min", "max", "random"),
                      control = nwqs_control(),
                      quiet = TRUE,
                      ...) {
  start_time <- Sys.time()

  family <- match.arg(family)
  plan_strategy <- match.arg(plan_strategy)
  transform_type <- match.arg(transform_type)
  ties <- match.arg(ties)

  if (n_boot < 20) {
    warning("'n_boot' is quite small; bootstrap percentile CI may be unstable.")
  }
  if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
    stop("Please install 'future' and 'future.apply' packages.")
  }

  cols_to_keep <- c("Target", "Term", "Estimate")

  old_plan <- tryCatch(
    configure_parallel_plan(
      loop_number = n_boot,
      strategy    = plan_strategy,
      n_workers   = n_workers,
      verbose     = FALSE
    ),
    error = function(e) {
      warning(
        "Failed to configure parallel plan: ", conditionMessage(e),
        ". Falling back to sequential."
      )
      NULL
    }
  )
  on.exit(
    {
      if (!is.null(old_plan)) future::plan(old_plan)
    },
    add = TRUE
  )

  use_parallel <- !inherits(future::plan(), "sequential")
  if (!use_parallel && !is.null(seed)) set.seed(seed)

  n_obs <- nrow(data)
  alpha <- 1 - conf_level

  one_boot <- function(b) {
    idx_boot <- sample(seq_len(n_obs), size = n_obs, replace = TRUE)
    data_boot <- data[idx_boot, , drop = FALSE]

    fit_b <- tryCatch(
      {
        nwqs(
          data = data_boot, mix_name = mix_name, covariates = covariates,
          outcome = outcome,
          transform_type = transform_type, q = q, ties = ties,
          family = family,
          plan_strategy = plan_strategy, rh = rh_inner,
          n_permutation = n_permutation, seed = NULL,
          control = control, quiet = TRUE, ...
        )
      },
      error = function(e) {
        message("Bootstrap ", b, " failed: ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(fit_b)) {
      return(list(Success = FALSE, Effects = NULL, Fit = NULL))
    }

    eff_b <- tryCatch(extract_nwqs_effects(fit_b), error = function(e) NULL)

    if (is.null(eff_b)) {
      return(list(Success = FALSE, Effects = NULL, Fit = if (keep_fits) fit_b else NULL))
    }

    eff_b_clean <- eff_b[, names(eff_b) %in% cols_to_keep, drop = FALSE]

    list(
      Success = TRUE,
      Effects = eff_b_clean,
      Weights = fit_b$final_weights,
      Shapes = fit_b$mean_shapes,
      Coefs = fit_b$mean_coefs,
      Struct = list(
        df_spline = fit_b$df_spline,
        spline_knots = fit_b$spline_knots,
        spline_boundary = fit_b$spline_boundary,
        formula = fit_b$formula
      ),
      Fit = if (keep_fits) fit_b else NULL
    )
  }

  if (use_parallel) {
    boot_results <- future.apply::future_lapply(
      seq_len(n_boot), one_boot,
      future.seed = if (!is.null(seed)) seed else TRUE
    )
  } else {
    boot_results <- lapply(seq_len(n_boot), one_boot)
  }

  boot_success <- vapply(boot_results, function(x) x$Success, logical(1))
  n_success <- sum(boot_success)

  if (n_success < max(20, ceiling(0.5 * n_boot))) {
    warning("A large proportion of bootstrap fits failed. Bootstrap CI may be unstable.")
  }
  if (n_success == 0) {
    stop("All bootstrap replicates failed.")
  }

  valid_results <- boot_results[boot_success]

  first_struct <- valid_results[[1]]$Struct

  valid_effs <- lapply(valid_results, function(x) x$Effects)
  success_indices <- which(boot_success)
  for (i in seq_along(valid_effs)) {
    valid_effs[[i]]$Boot_ID <- success_indices[i]
  }
  all_effs_df <- do.call(rbind, valid_effs)

  valid_weights <- lapply(valid_results, function(x) x$Weights)
  weights_mat <- do.call(rbind, valid_weights)
  avg_weights <- colMeans(weights_mat, na.rm = TRUE)
  avg_weights <- avg_weights / sum(avg_weights, na.rm = TRUE)

  valid_shapes <- lapply(valid_results, function(x) x$Shapes)

  if (is.numeric(valid_shapes[[1]]) || is.data.frame(valid_shapes[[1]]) || is.matrix(valid_shapes[[1]])) {
    avg_shapes <- Reduce("+", valid_shapes) / length(valid_shapes)
  } else {
    avg_shapes <- NULL
    warning("Complex shape structure detected; could not average shapes.")
  }

  valid_coefs <- lapply(valid_results, function(x) x$Coefs)
  coefs_mat <- do.call(rbind, valid_coefs)
  avg_coefs <- colMeans(coefs_mat, na.rm = TRUE)

  ci_lower <- aggregate(
    Estimate ~ Target + Term,
    data = all_effs_df,
    FUN = function(x) quantile(x, probs = alpha / 2, na.rm = TRUE)
  )
  names(ci_lower)[names(ci_lower) == "Estimate"] <- "Boot_CI_Lower"

  ci_upper <- aggregate(
    Estimate ~ Target + Term,
    data = all_effs_df,
    FUN = function(x) quantile(x, probs = 1 - alpha / 2, na.rm = TRUE)
  )
  names(ci_upper)[names(ci_upper) == "Estimate"] <- "Boot_CI_Upper"

  boot_mean <- aggregate(
    Estimate ~ Target + Term,
    data = all_effs_df,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  names(boot_mean)[names(boot_mean) == "Estimate"] <- "Boot_Mean"

  ci_table <- merge(boot_mean, ci_lower, by = c("Target", "Term"), all.x = TRUE)
  ci_table <- merge(ci_table, ci_upper, by = c("Target", "Term"), all.x = TRUE)
  ci_table$N_Success <- n_success

  is_exp_family <- family %in% c("binomial", "poisson", "quasipoisson", "negbin")

  disp_lcl <- ci_table$Boot_CI_Lower
  disp_ucl <- ci_table$Boot_CI_Upper
  disp_mean <- ci_table$Boot_Mean

  if (is_exp_family) {
    disp_lcl <- exp(disp_lcl)
    disp_ucl <- exp(disp_ucl)
    disp_mean <- exp(disp_mean)
  }

  ci_table$Formatted <- sprintf(
    "%.3f [%.3f, %.3f]",
    disp_mean, disp_lcl, disp_ucl
  )

  formatted_table <- reshape(
    ci_table[, c("Term", "Target", "Formatted")],
    idvar = "Term", timevar = "Target", direction = "wide"
  )
  names(formatted_table) <- gsub("^Formatted\\.", "", names(formatted_table))

  weight_df <- data.frame(
    Term = names(avg_weights),
    Weight = round(avg_weights, 3),
    stringsAsFactors = FALSE
  )
  weight_df <- rbind(
    data.frame(Term = "Overall", Weight = NA_real_, stringsAsFactors = FALSE),
    weight_df
  )

  formatted_table <- merge(weight_df, formatted_table, by = "Term", all.y = TRUE)
  formatted_table <- formatted_table[order(
    formatted_table$Term != "Overall",
    -formatted_table$Weight,
    na.last = TRUE
  ), ]
  rownames(formatted_table) <- NULL

  boot_raw <- all_effs_df
  if (is_exp_family) {
    boot_raw$Estimate <- exp(boot_raw$Estimate)
  }
  boot_raw$ColName <- paste(boot_raw$Term, boot_raw$Target, sep = "_")

  boot_contrast_matrix <- reshape(
    boot_raw[, c("Boot_ID", "ColName", "Estimate")],
    idvar = "Boot_ID", timevar = "ColName", direction = "wide"
  )
  names(boot_contrast_matrix) <- gsub("^Estimate\\.", "", names(boot_contrast_matrix))
  rownames(boot_contrast_matrix) <- NULL

  boot_fits <- if (keep_fits) lapply(valid_results, function(x) x$Fit) else NULL

  train_components_sorted <- lapply(mix_name, function(comp) sort(data[[comp]]))
  names(train_components_sorted) <- mix_name

  out <- list(
    ci_table = ci_table,
    formatted_table = formatted_table,
    boot_table = all_effs_df,
    boot_contrast_mat = boot_contrast_matrix,
    boot_fits = boot_fits,
    conf_level = conf_level,
    n_boot = n_boot,
    n_success = n_success,
    rh_inner = rh_inner,
    n_permutation = n_permutation,
    final_weights = avg_weights,
    mean_shapes = avg_shapes,
    family = family,
    transform_type = transform_type,
    q = q,
    ties = ties,
    df_spline = first_struct$df_spline,
    spline_knots = first_struct$spline_knots,
    spline_boundary = first_struct$spline_boundary,
    train_components_sorted = train_components_sorted,
    mix_name = mix_name,
    covariates = covariates,
    outcome = outcome,
    formula = first_struct$formula,
    mean_coefs = avg_coefs,
    rh_coefs_boot = coefs_mat,
    rh_weights_boot = weights_mat,
    data = data,
    call = match.call()
  )

  class(out) <- c("nwqs_boot", "list")

  message(sprintf(
    "NWQS Bootstrap completed: %d/%d successful fits. Time taken: %.2f mins",
    n_success, n_boot,
    as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  ))

  return(out)
}
