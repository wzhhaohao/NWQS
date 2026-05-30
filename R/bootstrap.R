#' @title OOB Permutation and Parallel Resampling Orchestrator
#'
#' @description
#' Core dispatching function in the NWQS framework. Prepares the non-linear
#' spline design matrix, manages multi-threaded parallel computation, and
#' distributes resampling tasks with strict random seeds to the underlying
#' weight engine (typically \code{\link{permutation_scorer}}).
#'
#' @param data \code{data.frame}. The original dataset containing all required
#'   variables.
#' @param mix_name Character vector. Column names of mixture components.
#' @param outcome Character. Column name of the outcome variable. Default is
#'   \code{"y"}.
#' @param covariates Character vector or \code{NULL}. Covariates/confounders to
#'   adjust for.
#' @param weight_engine Function. The engine function for computing weights and
#'   shapes. Default is \code{\link{permutation_scorer}}.
#' @param n_permutation Integer. Number of internal bootstrap or permutation
#'   iterations. Default is 30 (raised from 10 in 0.2.0).
#' @param seed Integer or \code{NULL}. Random seed for reproducible parallel
#'   computation.
#' @param boot_strategy Character. Parallel strategy: \code{"sequential"},
#'   \code{"multicore"}, or \code{"multisession"}.
#' @param boot_n_workers Integer or \code{NULL}. Number of parallel workers.
#' @param ... Additional arguments passed to \code{wqs_nonlinear_expand} or
#'   \code{weight_engine} (e.g., \code{model_knots}, \code{model_boundary},
#'   \code{q}, \code{df_spline}).
#'
#' @return A list of length \code{n_permutation}, each element containing the
#'   weights and shapes returned by \code{weight_engine}.
#'
#' @importFrom stats model.matrix as.formula
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#' @export
run_oob_permutation <- function(data, mix_name, outcome = "y",
                                covariates = NULL,
                                weight_engine = permutation_scorer,
                                n_permutation = 30, seed = NULL,
                                boot_strategy = c("sequential", "multicore", "multisession"),
                                boot_n_workers = NULL, ...) {
  args <- list(...)

  fam_arg <- if (!is.null(args$family)) args$family else "gaussian"
  if (is.character(fam_arg)) {
    fam_obj <- get(fam_arg, mode = "function")()
  } else {
    fam_obj <- fam_arg
  }
  args$family <- NULL

  boot_strategy <- match.arg(boot_strategy)

  old_plan <- configure_parallel_plan(
    loop_number = n_permutation,
    strategy = boot_strategy,
    n_workers = boot_n_workers
  )
  on.exit(future::plan(old_plan), add = TRUE)

  use.seed <- if (is.null(seed)) TRUE else seed
  expand_func <- if (!is.null(args$expand_func)) args$expand_func else wqs_nonlinear_expand

  q_val <- if (!is.null(args$q)) args$q else 4
  df_val <- if (!is.null(args$df_spline)) args$df_spline else 3

  if (is.null(args$model_knots) || is.null(args$model_boundary)) {
    stop("'model_knots' and 'model_boundary' must be provided.")
  }

  full_spline_data <- expand_func(
    data, mix_name,
    df_spline = df_val,
    knots = args$model_knots,
    boundary = args$model_boundary
  )
  spline_vars <- colnames(full_spline_data)

  if (is.null(covariates)) covariates <- character(0)

  missing_cov <- setdiff(covariates, names(data))
  if (length(missing_cov) > 0) stop("Missing covariates: ", paste(missing_cov, collapse = ", "))

  if (length(covariates) > 0) {
    formula_str <- paste("~", paste(c(covariates, spline_vars), collapse = " + "))
    temp_data <- cbind(data[, covariates, drop = FALSE], full_spline_data)
    internal_formula <- as.formula(formula_str)
    X_matrix <- model.matrix(internal_formula, data = temp_data)
  } else {
    X_matrix <- cbind("(Intercept)" = 1, full_spline_data)
  }

  y_raw <- data[[outcome]]

  if (fam_obj$family == "binomial") {
    if (is.factor(y_raw)) {
      if (nlevels(y_raw) != 2) stop("For binomial family, factor outcome must have exactly 2 levels.")
      y_vector <- as.numeric(y_raw) - 1
    } else {
      y_vector <- as.numeric(y_raw)
    }
  } else {
    y_vector <- as.numeric(y_raw)
  }

  results <- future.apply::future_lapply(
    seq_len(n_permutation),
    function(i) {
      tryCatch(
        {
          do.call(weight_engine, c(list(
            x = X_matrix,
            y = y_vector,
            mix_name = mix_name,
            spline_vars = spline_vars,
            family = fam_obj
          ), args))
        },
        error = function(e) {
          message("Permutation iteration ", i, " failed: ", conditionMessage(e))
          return(NULL)
        }
      )
    },
    future.seed = use.seed
  )

  names(results) <- paste0("B_", seq_len(n_permutation))
  results
}
