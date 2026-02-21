# #' Bootstrap General Runner (Parallel)
# #'
# #' @description
# #' Orchestrates the parallel execution of the bootstrap process for WQS regression.
# #' It automatically handles data resampling, parallel environment configuration (via `configure_parallel_plan`),
# #' and error handling for each iteration.
# #' \cr
# #' 协调 WQS 回归 Bootstrap 过程的并行执行。
# #' 它自动处理数据重抽样、并行环境配置（通过 `configure_parallel_plan`）以及每次迭代的错误处理。
# #'
# #' @details
# #' This function acts as a high-level wrapper that:
# #' \enumerate{
# #'   \item **Configures Parallelism**: Calls `configure_parallel_plan` to set up an optimal `future` plan based on the number of bootstrap samples (`B`) and available cores. It ensures that any temporary parallel plan is restored to the previous state upon exit.
# #'   \item **Executes Bootstrap**: Runs `B` iterations using `future_lapply`. Inside each iteration, the `model_func` is called on the full dataset (resampling logic is typically handled within `model_func`, e.g., `calc_spline_wqs_weights`).
# #'   \item **Handles Errors**: Uses `tryCatch` to ensure that a failure in a single iteration does not crash the entire process. Failed iterations return `NULL` with a warning.
# #' }
# #' \cr
# #' 该函数作为一个高级封装器，执行以下操作：
# #' \enumerate{
# #'   \item **配置并行**：调用 `configure_parallel_plan`，根据 Bootstrap 样本数 (`B`) 和可用核心设置最优的 `future` 计划。它确保任何临时并行计划在退出时恢复到之前的状态。
# #'   \item **执行 Bootstrap**：使用 `future_lapply` 运行 `B` 次迭代。在每次迭代中，调用 `model_func`（重抽样逻辑通常在 `model_func` 内部处理，例如 `calc_spline_wqs_weights`）。
# #'   \item **错误处理**：使用 `tryCatch` 确保单次迭代失败不会导致整个进程崩溃。失败的迭代将返回 `NULL` 并伴随警告。
# #' }
# #'
# #' @param data data.frame. The original dataset containing mixture and covariates.
# #'   包含混合物和协变量的原始数据框。
# #' @param mix_name character vector. Names of the mixture components.
# #'   混合物组分的变量名向量。
# #' @param dependent_var character. Name of the outcome variable. Defaults to "y".
# #'   因变量（结局变量）名称。默认为 "y"。
# #' @param model_func function. The internal modeling function to be executed for each bootstrap sample.
# #'   Defaults to `calc_spline_wqs_weights`.
# #'   针对每个 Bootstrap 样本执行的内部模型函数。默认为 `calc_spline_wqs_weights`。
# #' @param B integer. Number of bootstrap iterations. Default is 100.
# #'   Bootstrap 迭代次数。
# #' @param seed integer/logical. Random seed for parallel reproducibility.
# #'   If NULL (default), it uses `TRUE` to generate a stream of random seeds automatically.
# #'   并行计算的随机种子。如果为 NULL，则使用 `TRUE` 自动生成随机种子流。
# #' @param boot_strategy character. Parallel strategy for this bootstrap run: "sequential", "multicore", or "multicore".
# #'   Note: If running inside an already parallelized loop (e.g., outer Repeated Holdout), this should typically be "sequential".
# #'   本次 Bootstrap 的并行策略。注意：如果是在已并行的循环（如外层 RH）内部运行，此处通常应设为 "sequential"。
# #' @param boot_n_workers integer. Number of cores/workers to use.
# #'   If NULL, the optimal number is calculated by `configure_parallel_plan`.
# #'   使用的核心数。如果为 NULL，由 `configure_parallel_plan` 自动计算最优数量。
# #' @param ... Additional arguments passed directly to `model_func`.
# #'   Examples include `df_spline`, `shuffle`, `transform_fun`, etc.
# #'   直接传递给 `model_func` 的额外参数。例如 `df_spline`, `shuffle`, `transform_fun` 等。
# #'
# #' @return list. A named list of length `B` containing results from all iterations.
# #'   Failed iterations contain `NULL`.
# #'   返回一个长度为 `B` 的命名列表，包含所有迭代的结果。失败的迭代包含 `NULL`。
# #'
# #' @importFrom future plan
# #' @importFrom future.apply future_lapply
# #' @export
# run_bootstrap = function(data, 
#                          mix_name,
#                          dependent_var = "y",
#                          model_func = calc_spline_wqs_weights, 
#                          B = B, 
#                          seed = NULL,
#                          boot_strategy = c("sequential", "multicore", "multicore"),
#                          boot_n_workers = NULL, 
#                          ...) {
#     args = list(...)
#     boot_strategy = match.arg(boot_strategy)

#     # ---------- 并行策略配置 (调用独立函数) ----------
#     # 将 B (Bootstrap次数) 传给 loop_number (任务总数) 以进行负载均衡计算
#     old_plan = configure_parallel_plan(loop_number = B, 
#                                        strategy = boot_strategy, 
#                                        n_workers = boot_n_workers)
    
#     # 确保退出时还原旧计划
#     on.exit(future::plan(old_plan), add = TRUE)
#     # ------------------------------------------------

#     # 处理随机种子
#     use.seed = if (is.null(seed)) TRUE else seed

#     results = future.apply::future_lapply(seq_len(B), function(i) {

#         tryCatch({
#             do.call(
#               model_func,
#               c(list(data = data,
#                      mix_name = mix_name,
#                      dependent_var = dependent_var),
#                 args) 
#             )
#         }, error = function(e) {
#             warning(paste("Bootstrap iteration", i, "failed with error:", e$message))
#             return(NULL)
#         })

#     }, future.seed = use.seed)

#     names(results) = paste0("B_", seq_len(B))
#     return(results)
# }


#' Bootstrap General Runner (Parallel & Matrix Optimized)
#' 极限提速版：在外层预计算设计矩阵，避免在 Bootstrap 内部重复造轮子
#'
#' @export
run_bootstrap = function(data, 
                         mix_name,
                         dependent_var = "y",
                         model_func = calc_spline_wqs_weights, 
                         B = 100, 
                         seed = NULL,
                         boot_strategy = c("sequential", "multicore", "multisession"),
                         boot_n_workers = NULL, 
                         ...) {
    args = list(...)
    boot_strategy = match.arg(boot_strategy)

    # ---------- 1. 并行策略配置 ----------
    old_plan = configure_parallel_plan(loop_number = B, strategy = boot_strategy, n_workers = boot_n_workers)
    on.exit(future::plan(old_plan), add = TRUE)

    use.seed = if (is.null(seed)) TRUE else seed

    # ---------- 2. 极限优化：全局预计算设计矩阵 (Pre-computing) ----------
    # 提取转换参数与 family
    expand_func = if(!is.null(args$expand_func)) args$expand_func else wqs_nonlinear_expand
    q_val = if(!is.null(args$q)) args$q else 4
    df_val = if(!is.null(args$df_spline)) args$df_spline else 3
    
    fam_arg = if(!is.null(args$family)) args$family else "gaussian"
    fam_obj = if(is.character(fam_arg)) get(fam_arg, mode = "function")() else fam_arg

    # 仅执行一次非线性基函数展开
    full_spline_data = expand_func(data, mix_name, df_spline = df_val, q = q_val)
    spline_vars = colnames(full_spline_data)

    # 获取除了结局变量和毒物之外的所有协变量
    covariates = setdiff(names(data), c(mix_name, dependent_var))
    
    # 构建纯粹的设计矩阵模型 (完美处理分类变量的 dummy 编码)
    if (length(covariates) > 0) {
        formula_str = paste("~", paste(c(covariates, spline_vars), collapse = " + "))
    } else {
        formula_str = paste("~", paste(spline_vars, collapse = " + "))
    }
    internal_formula = as.formula(formula_str)
    
    # 组装用于 model.matrix 的临时 dataframe
    temp_data = cbind(data[, covariates, drop = FALSE], full_spline_data)
    
    # ★ 终极设计矩阵与结局向量 ★
    X_matrix = model.matrix(internal_formula, data = temp_data)
    y_vector = as.numeric(data[[dependent_var]])

    # ---------- 3. 执行 Bootstrap ----------
    results = future.apply::future_lapply(seq_len(B), function(i) {
        tryCatch({
            # 直接将矩阵和 fam_obj 传递给底层的计算函数，零数据处理损耗！
            do.call(
              model_func,
              c(list(X_matrix = X_matrix,
                     y_vector = y_vector,
                     mix_name = mix_name,
                     spline_vars = spline_vars,
                     fam_obj = fam_obj),
                args) 
            )
        }, error = function(e) {
            warning(paste("Bootstrap iteration", i, "failed with error:", e$message))
            return(NULL)
        })
    }, future.seed = use.seed)

    names(results) = paste0("B_", seq_len(B))
    return(results)
}
