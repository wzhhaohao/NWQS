#' @importFrom future plan availableCores

#' Configure Intelligent Parallel Plan with Load Balancing
#' 智能负载均衡并行策略配置
#'
#' @description
#' Automatically configures the parallel processing plan using the `future` package.
#' It implements a "Smart Load Balancing" strategy to optimize resource usage based on
#' the total number of tasks (`loop_number`) and available CPU cores.
#' \cr
#' 该函数利用 `future` 包自动配置并行处理计划。它采用“智能负载均衡”策略，根据任务总数和可用 CPU 核心优化资源使用。
#'
#' @details
#' This function performs the following steps to ensure efficient and safe parallelization:
#' \enumerate{
#'   \item **Check Existing Plan**: If a parallel plan is already active (set externally), it respects the current settings and returns immediately to avoid conflicts.
#'   \item **Sequential Mode**: If `strategy` is "sequential" or `loop_number` <= 1, it sets a sequential plan.
#'   \item **Smart Load Balancing**: If `n_workers` is NULL, it calculates the optimal number of workers:
#'   \itemize{
#'     \item **Safety**: Limits usage to `(1 - reserve_cpu)` of total cores to prevent system freeze.
#'     \item **Efficiency**: Calculates the minimum number of batches required to finish all loops.
#'     \item **Balance**: Distributes workers evenly across batches to minimize memory spikes and CPU contention.
#'     \item *Example*: Running 10 tasks on a machine with 8 safe cores. Instead of running 8 workers then 2 workers (unbalanced), it runs 2 rounds of 5 workers.
#'   }
#' }
#'
#' @param loop_number integer. The total number of iterations or tasks to be executed (e.g., bootstrap samples `B` or repeated holdouts `rh`).
#'   需要并行执行的任务/循环总数。
#' @param strategy character. The parallel strategy to use. Defaults to "multicore" (efficient for Linux/macOS).
#'   **Note**: Windows users should strictly use "multisession".
#'   并行策略，默认为 "multicore"。Windows 用户请务必使用 "multisession"。
#' @param n_workers integer. Optional. Manually specify the number of workers.
#'   If NULL (default), the function calculates the optimal number automatically.
#'   手动指定核心数。若为 NULL 则触发自动优化算法。
#' @param reserve_cpu numeric. The proportion of system CPU to reserve (0 to 1). Default is 0.2 (20%).
#'   保留给系统的 CPU 比例，防止死机。默认为 0.2。
#' @param ... Additional arguments passed to `future::plan()`.
#'
#' @return The *previous* future plan (invisibly). This is intended to be used with `on.exit()` to restore the environment state after execution.
#'   返回旧的并行计划（不可见对象），必须配合 `on.exit()` 使用以在函数结束时还原环境。
#'
#' @export
configure_parallel_plan = function(loop_number, strategy = "multicore", n_workers = NULL, reserve_cpu = 0.2, ...) {
    
    # 0. 获取当前计划
    current_plan = future::plan()
    is_already_parallel = !inherits(current_plan, "sequential")
    
    # A. 如果外部已经设置了并行，则不做干扰
    if (is_already_parallel) {
        return(invisible(current_plan))
    }
    
    # B. 如果策略是串行，或者任务数只有 1，直接设为串行
    if (strategy == "sequential" || loop_number <= 1) {
        old_plan = future::plan("sequential")
        return(invisible(old_plan))
    }
    
    # C. 智能核心数计算 (Smart Load Balancing)
    if (is.null(n_workers)) {
        # 1. 获取物理核心总数
        total_cores = future::availableCores()
        
        # 2. 计算安全上限
        safe_limit = floor(total_cores * (1 - reserve_cpu))
        if (safe_limit < 1) safe_limit = 1L
        
        # 3. 负载均衡算法
        # 计算跑完所有循环所需的最少轮次 (Batches)
        # 这里用 loop_number 
        min_batches = ceiling(loop_number / safe_limit)      
        
        # 倒推该轮次下的均摊核心数
        workers_final = ceiling(loop_number / min_batches)   
        workers_final = as.integer(workers_final)
        
        message(sprintf("Auto-Parallel: %d Cores Available (Limit %d). %d Loops split into %d rounds x %d workers.", 
                        total_cores, safe_limit, loop_number, min_batches, workers_final))
        
    } else {
        # 用户强制指定
        workers_final = max(1L, as.integer(n_workers))
    }
    
    # D. 应用新计划并返回旧计划
    old_plan = future::plan(strategy, workers = workers_final)
    
    return(invisible(old_plan))
}
