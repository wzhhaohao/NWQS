# =========================================================================
# 独立脚本：预生成并永久锁定 Sigma 协方差矩阵 (DGP Parameters)
# =========================================================================
rm(list = ls())
library(NWQS) # 确保加载你的包

# 创建专属的存放目录
out_dir <- file.path("data", "Sigma_Matrices")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 定义需要生成的所有维度和相关性结构的穷举网格
P_values <- c(4, 8, 12)
corr_modes <- c("mixed", "low", "high")

message("🔒 开始生成全场景 Sigma 矩阵并保存为 CSV...")

for (p in P_values) {
  for (mode in corr_modes) {
    dict_key <- paste0("P", p, "_", mode)
    
    # 强行按住种子，确保绝对可复现
    set.seed(525)
    mat <- generate_sigma(n_vars = p, mode = mode, seed = 525)
    
    # 添加严谨的行列名称
    comp_names <- paste0("Component", 1:p)
    rownames(mat) <- comp_names
    colnames(mat) <- comp_names
    
    # 保存为 CSV (保留行名)
    file_path <- file.path(out_dir, paste0("True_Sigma_Matrix_", dict_key, ".csv"))
    write.csv(mat, file_path, row.names = TRUE)
    
    message(sprintf("   -> [OK] 维度: P=%-2d | 相关性: %-5s | 保存至: %s", p, mode, file_path))
  }
}

message("\n✅ 所有 9 个矩阵已成功锁死并存盘！以后主脚本直接读取即可。")