#!/usr/bin/env Rscript
# install_r_github_pkgs.R — 安装 scrna-integration-r 环境的重型 R 包
#
# 使用方式:
#   conda activate scrna-integration-r
#   Rscript scripts/install_r_github_pkgs.R
#
# 也可以在 Python 环境中通过 subprocess 调用:
#   subprocess.run([RSCRIPT_BIN, "--vanilla", "scripts/install_r_github_pkgs.R"])
#
# 依赖: scrna-integration-r conda 环境已激活, 编译器工具链就绪
# 时间: 首次运行 20-40 分钟（源码编译）
#
# 2026-06-15: 从 docs/cross-platform-exceptions.md 提取安装命令固化

Sys.setenv(R_INSTALL_STAGED = "FALSE")
options(timeout = 900)
options(repos = c(CRAN = "https://cloud.r-project.org"))

cat("=== Step 1/4: WGCNA (CRAN) ===\n")
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  install.packages("WGCNA", dependencies = TRUE)
}
cat("  WGCNA:", as.character(packageVersion("WGCNA")), "\n")

cat("\n=== Step 2/4: monocle3 (Bioconductor) ===\n")
if (!requireNamespace("monocle3", quietly = TRUE)) {
  # 确保 BiocManager 可用
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  # 先尝试安装缺失的 Bioconductor 依赖
  if (!requireNamespace("batchelor", quietly = TRUE))
    BiocManager::install("batchelor", update = FALSE, ask = FALSE)
  if (!requireNamespace("BPCells", quietly = TRUE))
    BiocManager::install("BPCells", update = FALSE, ask = FALSE)
  if (!requireNamespace("DelayedMatrixStats", quietly = TRUE))
    BiocManager::install("DelayedMatrixStats", update = FALSE, ask = FALSE)

  # 安装关键数据包 (如果下载超时，参见 docs/cross-platform-exceptions.md)
  if (!requireNamespace("GenomeInfoDbData", quietly = TRUE)) {
    cat("  下载 GenomeInfoDbData (12.3MB)...\n")
    BiocManager::install("GenomeInfoDbData", update = FALSE, ask = FALSE)
  }

  BiocManager::install("monocle3", update = FALSE, ask = FALSE)
}
cat("  monocle3:", as.character(packageVersion("monocle3")), "\n")

cat("\n=== Step 3/4: hdWGCNA (GitHub) ===\n")
if (!requireNamespace("hdWGCNA", quietly = TRUE)) {
  # 确保 UCell 可用
  if (!requireNamespace("UCell", quietly = TRUE))
    BiocManager::install("UCell", update = FALSE, ask = FALSE)
  remotes::install_github("smorabit/hdWGCNA", upgrade = FALSE)
}
cat("  hdWGCNA:", as.character(packageVersion("hdWGCNA")), "\n")

cat("\n=== Step 4/4: CellChat (GitHub) ===\n")
if (!requireNamespace("CellChat", quietly = TRUE)) {
  remotes::install_github("jinworks/CellChat", upgrade = FALSE)
}
cat("  CellChat:", as.character(packageVersion("CellChat")), "\n")

cat("\n=== 安装完成，最终状态 ===\n")
for (pkg in c("WGCNA", "monocle3", "hdWGCNA", "CellChat")) {
  cat(sprintf("  %-20s: %s\n", pkg,
    tryCatch(as.character(packageVersion(pkg)), error = function(e) "FAILED")))
}
