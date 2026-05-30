# NWQS: 非线性加权分位数和回归

<!-- badges: start -->
<!-- badges: end -->

[English](README.md) | [中文](README_CN.md)

**NWQS** 在经典加权分位数和（WQS）回归的基础上引入自然三次样条，能够捕获环境混合暴露分析中的非线性剂量-反应关系（如阈值效应、U型、S型曲线等）。

## 核心特性

- **非线性剂量-反应建模**：通过自然三次样条基展开，支持阈值效应、U型、倒U型、S型等复杂暴露-反应模式。
- **置换重要性变量筛选**：基于袋外（OOB）损失变化推导各成分的相对权重，避免正则化偏差。
- **重复留出（Repeated Holdout）架构**：通过多次随机数据划分实现稳健的权重与形状估计。
- **外部 Bootstrap 推断**（`nwqs_boot()`）：提供基于百分位法的有效置信区间，反映真实抽样变异性。
- **多种 GLM 族**：支持高斯、二项、泊松、准泊松、负二项（`MASS::glm.nb`）。（0.1.x 曾支持条件 Logistic 回归 `clogit`，在 0.2.0 中已暂时移除；有序回归 `MASS::polr` 计划于后续版本加入。详见 `NEWS.md`。）
- **出版级可视化**：包含剂量-反应曲线、成分权重柱状图和 Bootstrap 对比箱线图。
- **蒙特卡洛仿真工具集**：用于在受控场景下对模型性能进行基准测试和验证。
- **自动并行计算**：基于 `future` 框架实现智能负载均衡。

## 安装

从 GitHub 安装开发版本：

```r
# install.packages("devtools")
devtools::install_github("wzhhaohao/NWQS")
```

## 快速上手

### 1. 模拟数据

```r
library(NWQS)

n_vars <- 7
mu <- rep(0, n_vars)
sigma <- generate_sigma(n_vars, mode = "medium", seed = 42)

dat <- gen_nonlinear_data(
  n_obs       = 500,
  mu_preds    = mu,
  sigma_preds = sigma,
  beta_wqs    = 1,
  beta_preds  = c(0.4, 0.3, 0.2, 0.1, 0, 0, 0),
  snr_db      = 10,
  shape       = c("linear_like", "threshold", "u_shape",
                   "s_shape", "linear_like", "linear_like", "linear_like"),
  seed        = 123
)
```

### 2. 拟合 NWQS 模型

```r
mix_name   <- paste0("Component", 1:n_vars)
covariates <- c("x_cont", "x_bin", "x_cat")

fit <- nwqs(
  data           = dat,
  mix_name       = mix_name,
  covariates     = covariates,
  outcome        = "y",
  family         = "gaussian",
  transform_type = "percentile_rank",  # 0.2.0 默认；传 "q_bin" 复现 0.1.x 行为
  q              = 4,                  # percentile_rank 模式下表示对比点数量
  df_spline      = 3,
  rh             = 10,
  n_permutation  = 50,
  plan_strategy  = "sequential",
  seed           = 1234
)

print(fit)
plot(fit)
```

### 3. Bootstrap 推断

```r
boot_res <- nwqs_boot(
  data           = dat,
  mix_name       = mix_name,
  covariates     = covariates,
  outcome        = "y",
  family         = "gaussian",
  n_boot         = 100,
  rh_inner       = 1,
  n_permutation  = 50,
  plan_strategy  = "multisession",
  seed           = 42
)

print(boot_res)
plot(boot_res)
```

## 支持的模型族

| 模型族 | 结局类型 | 效应量度 |
|---|---|---|
| `gaussian` | 连续型 | 线性效应 (Delta Y) |
| `binomial` | 二分类 (0/1) | 对数比值比 / 比值比 (OR) |
| `poisson` | 计数型 | 对数率比 / 率比 (RR) |
| `quasipoisson` | 过离散计数 | 对数率比 / 率比 (RR) |
| `negbin` | 过离散计数（负二项）| 对数率比 / 率比 (RR) |

> `clogit`（匹配病例对照）在 0.1.x 中受支持，在 0.2.0 中为配合主路径
> 重构而暂时移除。后续小版本会以独立函数 `nwqs_clogit()` 的形式重新
> 引入条件 Logistic 回归。

## 核心函数

| 函数 | 说明 |
|---|---|
| `nwqs()` | 拟合带有重复留出的非线性 WQS 模型 |
| `nwqs_boot()` | Bootstrap 置信区间（有效推断） |
| `nwqs_contrast()` | 联合暴露分位数对比显著性检验 |
| `extract_nwqs_effects()` | 提取各成分特异性分位数对比效应 |
| `plot.nwqs()` | 剂量-反应曲线及权重诊断图 |
| `plot.nwqs_boot()` | Bootstrap 对比箱线图 |

## 数据模拟函数

| 函数 | 说明 |
|---|---|
| `gen_nonlinear_data()` | 连续结局 + 非线性样条效应 |
| `gen_nonlinear_bio_data()` | 二分类结局 + 自动患病率校准 |
| `gen_nonlinear_count_data()` | 泊松计数结局 |
| `generate_linear_data()` | 线性混合效应（基线对照） |
| `generate_sigma()` | 相关矩阵（低/中/高/混合相关） |

## 方法概述

NWQS 算法采用三阶段架构：

1. **权重与形状发现**（训练集）：OOB 置换重要性推导各成分权重；样条基系数捕获非线性形状。
2. **形状标准化**：将形状标准化至单位方差，从而解耦形状与效应量。
3. **单自由度效应估计**（验证集）：将标准化形状和权重投影为单一 `nwqs` 指数，再通过标准 GLM 估计混合暴露的整体效应。

当 `rh > 1` 时，`nwqs()` 输出的标准误仅反映数据划分变异性，**不可用于统计推断**。请使用 `nwqs_boot()` 获取有效的百分位 Bootstrap 置信区间。

## 许可证

GPL-3 许可证。详见 [LICENSE](LICENSE)。

## 作者

1. **Wang Zhehao**（[wangzhehao_bill@foxmail.com](mailto:wangzhehao_bill@foxmail.com)）— 暨南大学基础医学与公共卫生学院，广州 510632
2. **Chen Shirui**（[chenshr7@bjmu.edu.cn](mailto:chenshr7@bjmu.edu.cn)）— 北京大学公众健康与重大疫情防控战略研究中心，北京 100191
3. **Lin Ziqiang**（[linziqiang0314@gmail.com](mailto:linziqiang0314@gmail.com)）— 暨南大学基础医学与公共卫生学院，广州 510632（通讯作者）
