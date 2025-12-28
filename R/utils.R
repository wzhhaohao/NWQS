#' @importFrom splines ns
#' @importFrom stats quantile

#' Quantile or Percentile Transformation / 分位数或百分位数变换
#'
#' @description
#' A hybrid function combining flexible ranking methods:
#' 1. "quantile": gWQS-style integer binning (handles ties and boundaries robustly).
#' 2. "percentile": Continuous percentile ranking (0 to 1).
#' 这是一个融合函数，结合了灵活的排名方法：
#' 1. "quantile": gWQS 风格的整数分箱（稳健处理重复值和边界）。
#' 2. "percentile": 连续百分位数排名（0 到 1）。
#'
#' @param data data.frame. Input data. / 输入数据。
#' @param method character. "quantile" (default) or "percentile". / 变换方法。
#' @param q integer. Number of quantiles (only used if method = "quantile"). / 分位数数量。
#' @return data.frame. Transformed data. / 变换后的数据。
#' @export
trans_quantile = function(data, method = c("quantile", "percentile"), q = 4) {
    data = as.data.frame(data)
    method = match.arg(method)
  
    if (method == "percentile") {
        transform_func = function(x) {
            rank(x) / (length(x) + 1)
        }
    
        # Apply to all columns
        res_list = lapply(data, transform_func)
    
    } else {
        if (!is.numeric(q) || q < 1) stop("'q' must be a positive number")

        transform_func = function(x) {
            breaks = unique(quantile(x, probs = seq(0, 1, by = 1/q), na.rm = TRUE))
      
            # Handle Boundaries (gWQS style: Force -Inf / Inf for robustness)
            if (length(breaks) == 1) {
                breaks = c(-Inf, breaks)
            } else {
                breaks[1] = -Inf
                breaks[length(breaks)] = Inf
            }
      
            # Cut (Binning): Returns integer values from 0 to q-1
            as.numeric(cut(x, breaks = breaks, labels = FALSE, include.lowest = TRUE)) - 1
        }
    
        # Apply to all columns
        res_list = lapply(data, transform_func)
    }
  
    # Convert list back to data.frame and preserve names
    res_df = as.data.frame(res_list)
    names(res_df) = names(data)
  
    return(res_df)
}

# -------------------------------------------------------------------------
#' Nonlinear Expansion for WQS (Natural Splines) / WQS 非线性展开 (自然样条)
#'
#' @description
#' Transforms mixture variables into natural cubic spline bases to capture nonlinear effects.
#' 将混合物变量转换为自然三次样条基函数，以捕捉非线性效应。
#'
#' @details
#' By default, this function performs a quantile transformation (quartiles) before spline expansion
#' if no `transform_fun` is provided. It uses `splines::ns` for the basis expansion.
#' 默认情况下，如果未提供 `transform_fun`，该函数会在样条展开前执行四分位数转换。
#' 它使用 `splines::ns` 生成样条基底。
#'
#' @param data data.frame. The dataset containing the mixture variables.
#'   包含混合物变量的数据集。
#' @param mix_name character vector. Names of the mixture components to be expanded.
#'   需要展开的混合物组分名称。
#' @param transform_fun function. Optional custom transformation function applied before spline expansion.
#'   If NULL, applies a default quantile transformation (q=4).
#'   可选的自定义转换函数。如果为 NULL，则应用默认的四分位数转换。
#' @param df_spline integer. Degrees of freedom for the natural spline. Default is 3.
#'   自然样条的自由度。默认为 3。
#'
#' @return matrix. A matrix containing the spline basis functions for all mixture components.
#'   Column names are formatted as `{Component}_B{BasisIndex}`.
#'   返回包含所有混合物组分样条基函数的矩阵。
#' @export
#' FIXME q 参数没有传递过去,之后默认4，先测试传入
wqs_nonlinear_expand = function(data, mix_name, df_spline = 3, q = 4) {

    # FIXME：可能是这里的q没有传递过去
    trans_data = data[, mix_name] 

    X = 0:(q - 1)
    temp_spline = splines::ns(X, df = df_spline)
    temp_knots = attr(temp_spline, "knots")
    temp_boundary = attr(temp_spline, "Boundary.knots")

    # 对每一列应用自然样条展开
    mat_spline_list = lapply(trans_data, function(x) splines::ns(x, df = df_spline, knots = temp_knots, Boundary.knots = temp_boundary))
    mat_spline_full = do.call(cbind, mat_spline_list)

    # 生成列名
    total_cols = ncol(mat_spline_full)
    cols_per_mix = total_cols / length(mix_name)

    colnames(mat_spline_full) = paste0(
        rep(mix_name, each = cols_per_mix), 
        "_B", 
        rep(1:cols_per_mix, times = length(mix_name)))
    
    return(mat_spline_full)
}