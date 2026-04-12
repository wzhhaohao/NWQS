#' @title Intelligent Parallel Plan Configuration with Load Balancing
#'
#' @description
#' Automatically configures a parallel processing plan using the \pkg{future}
#' framework. Employs smart load balancing to safely and efficiently allocate
#' computing resources based on the total number of tasks and available cores.
#'
#' @details
#' Key design features:
#' \itemize{
#'   \item \strong{Safety guard:} Reserves a proportion of system cores
#'     (\code{reserve_cpu}) to prevent system freeze during large-scale
#'     bootstrap or Monte Carlo simulations.
#'   \item \strong{Load balancing:} Avoids uneven tail-end allocation. For
#'     example, running 10 tasks on 8 cores would naively produce 1 round of
#'     8 + 1 round of 2; this algorithm smooths it to 2 rounds of 5.
#'   \item \strong{Non-invasive:} If an external parallel plan is already
#'     active, the function silently returns without overriding it.
#' }
#'
#' @param loop_number Integer. Total number of parallel tasks (e.g., bootstrap
#'   replicates or repeated holdout iterations).
#' @param strategy Character. Parallel strategy. Default is
#'   \code{"multisession"} (cross-platform). Linux/macOS users may use
#'   \code{"multicore"} for better performance.
#' @param n_workers Integer or \code{NULL}. Manual core count override. If
#'   \code{NULL} (recommended), the automatic optimization algorithm is used.
#' @param reserve_cpu Numeric in (0, 1). Proportion of CPU cores reserved for
#'   the operating system. Default is 0.2 (20 percent).
#' @param verbose Logical. Whether to print load balancing information.
#'   Default is \code{TRUE}.
#' @param ... Additional arguments passed to \code{future::plan()}.
#'
#' @return Invisibly returns the \emph{previous} future plan. Should be used
#'   with \code{on.exit()} to restore the user's environment on function exit.
#'
#' @importFrom future plan availableCores
#' @export
configure_parallel_plan <- function(loop_number, strategy = "multisession",
                                    n_workers = NULL, reserve_cpu = 0.2,
                                    verbose = TRUE, ...) {
  current_plan <- future::plan()
  is_already_parallel <- !inherits(current_plan, "sequential")

  if (is_already_parallel) {
    return(invisible(current_plan))
  }

  if (strategy == "sequential" || loop_number <= 1) {
    old_plan <- future::plan("sequential")
    return(invisible(old_plan))
  }

  if (is.null(n_workers)) {
    total_cores <- future::availableCores()
    safe_limit <- floor(total_cores * (1 - reserve_cpu))
    if (safe_limit < 1) safe_limit <- 1L

    min_batches <- ceiling(loop_number / safe_limit)
    workers_final <- ceiling(loop_number / min_batches)
    workers_final <- as.integer(workers_final)

    if (isTRUE(verbose)) {
      message(sprintf(
        "Auto-Parallel: %d Cores Available (Limit %d). %d Loops split into %d rounds x %d workers.",
        total_cores, safe_limit, loop_number, min_batches, workers_final
      ))
    }
  } else {
    workers_final <- max(1L, as.integer(n_workers))
  }

  old_plan <- future::plan(strategy, workers = workers_final)

  return(invisible(old_plan))
}
