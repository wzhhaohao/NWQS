#' @title 自动配置带负载均衡的智能并行计划 (Intelligent Parallel Plan Configuration)
#'
#' @description
#' 利用 \pkg{future} 包自动配置并行处理计划。该函数采用“智能负载均衡 (Smart Load Balancing)”策略，
#' 能够根据总任务数和系统可用物理核心，安全且高效地优化计算资源的分配。
#'
#' @details
#' \strong{系统稳定性与统计重现性保障:}
#' \enumerate{
#'   \item \strong{防死机机制 (Safety Guard):} 强制保留 \code{reserve_cpu} 比例的系统核心，避免在进行
#'         大规模 Bootstrap 或 Monte Carlo 模拟时耗尽计算资源导致系统崩溃。
#'   \item \strong{负载均衡 (Load Balancing):} 避免了原生的末端分配不均问题。例如，在 8 个可用核心上运行 10 个任务，
#'         常规分配可能导致第一轮跑 8 个，第二轮仅跑 2 个（引发内存突刺和 CPU 争用）。本算法会将其智能平滑为 2 轮各 5 个工作节点。
#'   \item \strong{无侵入设计:} 如果外部环境（如用户的 \code{.Rprofile} 或更高层脚本）已经激活了并行计划，
#'         本函数将静默退出，绝不覆盖用户的全局设定。
#' }
#'
#' @param loop_number Integer。需要并行执行的循环/任务总数（例如 Bootstrap 的重抽样次数 \code{B} 或 \code{rh}）。
#' @param strategy Character。并行策略。默认为 \code{"multisession"}（对 Windows 友好且兼容性高）。Linux/macOS 用户为追求极致性能可指定 \code{"multicore"}。
#' @param n_workers Integer 或 \code{NULL}。手动指定的核心数。若为 \code{NULL}（推荐），则触发自动优化算法。
#' @param reserve_cpu Numeric (0, 1)。保留给操作系统的 CPU 核心比例，默认为 0.2 (20\%)。
#' @param verbose Logical。是否打印自动并行负载均衡信息。默认为 \code{TRUE}。
#' @param ... 传递给 \code{future::plan()} 的额外参数。
#'
#' @return 隐式返回 \emph{先前} 的 future 并行计划。必须配合 \code{on.exit()} 使用，以确保在函数退出时还原用户环境。
#'
#' @importFrom future plan availableCores
#' @export
configure_parallel_plan <- function(loop_number, strategy = "multisession", n_workers = NULL, reserve_cpu = 0.2, verbose = TRUE, ...) {
    # 0. 获取当前计划
    current_plan <- future::plan()
    is_already_parallel <- !inherits(current_plan, "sequential")

    # A. 如果外部已经设置了并行，则不做干扰
    if (is_already_parallel) {
        return(invisible(current_plan))
    }

    # B. 如果策略是串行，或者任务数只有 1，直接设为串行
    if (strategy == "sequential" || loop_number <= 1) {
        old_plan <- future::plan("sequential")
        return(invisible(old_plan))
    }

    # C. 智能核心数计算 (Smart Load Balancing)
    if (is.null(n_workers)) {
        # 1. 获取物理核心总数
        total_cores <- future::availableCores()

        # 2. 计算安全上限
        safe_limit <- floor(total_cores * (1 - reserve_cpu))
        if (safe_limit < 1) safe_limit <- 1L

        # 3. 负载均衡算法
        # 计算跑完所有循环所需的最少轮次 (Batches)
        # 这里用 loop_number
        min_batches <- ceiling(loop_number / safe_limit)

        # 倒推该轮次下的均摊核心数
        workers_final <- ceiling(loop_number / min_batches)
        workers_final <- as.integer(workers_final)

        if (isTRUE(verbose)) {
            message(sprintf(
                "Auto-Parallel: %d Cores Available (Limit %d). %d Loops split into %d rounds x %d workers.",
                total_cores, safe_limit, loop_number, min_batches, workers_final
            ))
        }
    } else {
        # 用户强制指定
        workers_final <- max(1L, as.integer(n_workers))
    }

    # D. 应用新计划并返回旧计划
    old_plan <- future::plan(strategy, workers = workers_final)

    return(invisible(old_plan))
}
