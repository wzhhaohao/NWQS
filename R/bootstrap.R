#' Bootstrap General Runner (Parallel)
#'
#' @description
#' Orchestrates repeated bootstrap / OOB-based internal model fitting for
#' NWQS regression. It automatically handles parallel environment
#' configuration (via `configure_parallel_plan`) and error handling for
#' each iteration.
#' \cr
#' 协调 NWQS 回归内部 Bootstrap / OOB 过程的并行执行。
#' 它自动处理并行环境配置（通过 `configure_parallel_plan`）
#' 以及每次迭代的错误处理。
#'
#' @details
#' This function acts as a high-level wrapper that:
#' \enumerate{
#'   \item Configures parallelism for the internal bootstrap loop.
#'   \item Pre-computes the spline-expanded design matrix once, using the
#'         globally fixed spline knots and boundaries passed from the
#'         outer `nwqs()` routine.
#'   \item Executes `B` iterations using `future_lapply`. Inside each
#'         iteration, the `weight_engine` is called on the full precomputed
#'         design matrix, while resampling logic is handled inside
#'         `weight_engine` (e.g., `permutation_scorer`).
#'   \item Uses `tryCatch` so that a failure in a single iteration does
#'         not crash the entire process.
#' }
#'
#' @param data data.frame. The original dataset containing mixture and
#'   covariates.
#' @param mix_name character vector. Names of the mixture components.
#' @param outcome character. Name of the outcome variable.
#'   Defaults to "y".
#' @param weight_engine function. The internal modeling function to be
#'   executed for each bootstrap sample. Defaults to
#'   `permutation_scorer`.
#' @param n_permutation integer. Number of permutation iterations. Default is 100.
#' @param seed integer/logical. Random seed for parallel reproducibility.
#'   If NULL (default), it uses `TRUE` to generate a stream of random
#'   seeds automatically.
#' @param boot_strategy character. Parallel strategy for this bootstrap
#'   run: "sequential", "multicore", or "multisession".
#'   If running inside an already parallelized outer RH loop, this should
#'   typically be "sequential".
#' @param boot_n_workers integer. Number of cores/workers to use.
#'   If NULL, the optimal number is calculated by
#'   `configure_parallel_plan`.
#' @param ... Additional arguments passed directly to `weight_engine`.
#'
#' @return list. A named list of length `B` containing results from all
#'   iterations. Failed iterations contain `NULL`.
#'
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#' @export
run_oob_permutation <- function(data, mix_name, outcome = "y",
                                weight_engine = permutation_scorer,
                                n_permutation = 100, seed = NULL,
                                boot_strategy = c("sequential","multicore","multisession"),
                                boot_n_workers = NULL, ...) {
    args <- list(...)
    boot_strategy <- match.arg(boot_strategy)

    old_plan <- configure_parallel_plan(
        loop_number = n_permutation,
        strategy = boot_strategy,
        n_workers = boot_n_workers
    )
    on.exit(future::plan(old_plan), add = TRUE)

    use.seed <- if (is.null(seed)) TRUE else seed

    expand_func <- if (!is.null(args$expand_func)) {
        args$expand_func
    } else {
        wqs_nonlinear_expand
    }

    q_val <- if (!is.null(args$q)) args$q else 4
    df_val <- if (!is.null(args$df_spline)) args$df_spline else 3

    if (is.null(args$model_knots) || is.null(args$model_boundary)) {
        stop(
            "Critical Error: 'model_knots' and 'model_boundary' must be ",
            "passed to run_oob_permutation to maintain a fixed spline basis."
        )
    }

    fam_arg <- if (!is.null(args$family)) args$family else "gaussian"
    fam_obj <- if (is.character(fam_arg)) {
        get(fam_arg, mode = "function")()
    } else {
        fam_arg
    }

    full_spline_data <- expand_func(
        data, mix_name,
        df_spline = df_val,
        knots = args$model_knots,
        boundary = args$model_boundary
    )
    spline_vars <- colnames(full_spline_data)

    covariates <- setdiff(names(data), c(mix_name, outcome))
    if (length(covariates) > 0) {
        formula_str <- paste("~", paste(c(covariates, spline_vars), collapse = " + "))
    } else {
        formula_str <- paste("~", paste(spline_vars, collapse = " + "))
    }

    internal_formula <- as.formula(formula_str)
    temp_data <- cbind(data[, covariates, drop = FALSE], full_spline_data)
    X_matrix <- model.matrix(internal_formula, data = temp_data)

    y_raw <- data[[outcome]]
    if (fam_obj$family == "binomial") {
        if (is.factor(y_raw)) {
            if (nlevels(y_raw) != 2) {
                stop("For binomial family, factor outcome must have exactly 2 levels.")
            }
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
                    warning(paste("Bootstrap iteration", i, "failed with error:", e$message))
                    return(NULL)
                }
            )
        },
        future.seed = use.seed
    )

    names(results) <- paste0("B_", seq_len(n_permutation))
    return(results)
}
