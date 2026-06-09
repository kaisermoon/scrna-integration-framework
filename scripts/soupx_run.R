#!/usr/bin/env Rscript --vanilla
# soupx_run.R — SoupX 环境 RNA 校正（去除 ambient RNA 污染）
#
# 用法（由 stage2 notebook subprocess 调用）:
#   Rscript --vanilla scripts/soupx_run.R <work_dir> <filtered_mtx_dir> <raw_mtx_dir> <sample_id>
#
# 参数:
#   work_dir         — 工作目录（在此产出 corrected_counts.mtx + barcodes.tsv + features.tsv）
#   filtered_mtx_dir — Python 侧导出的过滤后计数矩阵目录（含 matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz）
#   raw_mtx_dir      — cellranger raw_feature_bc_matrix 目录（含 matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz）
#   sample_id        — 样本标识（用于日志输出）
#
# 输出:
#   {work_dir}/corrected_counts.mtx  — SoupX 校正后的稀疏 counts 矩阵 (MatrixMarket format)
#   {work_dir}/barcodes.tsv          — 细胞条形码（与输入 filtered 一致）
#   {work_dir}/features.tsv          — 基因信息（与输入 filtered 一致）
#
# SoupX 工作原理:
#   - 利用 raw（含空液滴/碎片）和 filtered（只含真实细胞）两套计数之间的差异
#   - 估计每个基因的环境污染比例（ambient contamination fraction）
#   - 从 filtered 计数中减去估计的污染部分，产出校正后计数
#   - 需要两套矩阵来自同一个 cellranger 运行（同一套 barcodes 取子集）

# ---- 入参解析 ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("usage: Rscript --vanilla soupx_run.R <work_dir> <filtered_mtx_dir> <raw_mtx_dir> <sample_id>")
}
work_dir         <- args[1]
filtered_mtx_dir <- args[2]
raw_mtx_dir      <- args[3]
sample_id        <- args[4]

cat(sprintf("[SoupX] sample=%s  work_dir=%s\n", sample_id, work_dir))
cat(sprintf("[SoupX] filtered_mtx_dir=%s\n", filtered_mtx_dir))
cat(sprintf("[SoupX] raw_mtx_dir=%s\n", raw_mtx_dir))

# ---- 静默加载包 ----
suppressPackageStartupMessages({
  library(Matrix)
  library(SoupX)
})

# ---- 辅助函数：读 10x 式 mtx 目录（兼容 .gz 和未压缩） ----
read_10x_like <- function(mtx_dir, label) {
  # 三个文件: matrix.mtx[.gz], barcodes.tsv[.gz], features.tsv[.gz]
  find_file <- function(base) {
    fp_gz <- file.path(mtx_dir, paste0(base, ".gz"))
    fp    <- file.path(mtx_dir, base)
    if (file.exists(fp_gz)) return(fp_gz)
    if (file.exists(fp))    return(fp)
    stop(sprintf("[SoupX] %s: '%s' or '%s.gz' not found in %s", label, base, base, mtx_dir))
  }

  mtx_fp  <- find_file("matrix.mtx")
  bc_fp   <- find_file("barcodes.tsv")
  feat_fp <- find_file("features.tsv")

  counts <- Matrix::readMM(mtx_fp)
  barcodes <- readLines(bc_fp)
  features <- read.table(feat_fp, sep = "\t", header = FALSE,
                         col.names = c("feature_id", "feature_name", "feature_type"),
                         stringsAsFactors = FALSE, comment.char = "")

  # 行 = 基因 (features)，列 = 细胞 (barcodes)
  # Matrix::readMM 输出是 dgTMatrix，转为 dgCMatrix 以便后续操作
  counts <- as(counts, "dgCMatrix")
  rownames(counts) <- features$feature_name
  colnames(counts) <- barcodes

  cat(sprintf("[SoupX] %s: loaded %d genes x %d barcodes\n", label, nrow(counts), ncol(counts)))
  return(list(counts = counts, features = features, barcodes = barcodes))
}

# ---- 读入 filtered 和 raw 矩阵 ----
filtered <- read_10x_like(filtered_mtx_dir, "filtered")
raw      <- read_10x_like(raw_mtx_dir, "raw")

# ---- 基因交集断言 ----
# SoupX 要求 filtered 和 raw 的基因集一致（两个矩阵来自同一个 cellranger 运行的
# raw_feature_bc_matrix 和 filtered_feature_bc_matrix，基因集应完全相同）。
common_genes <- intersect(rownames(filtered$counts), rownames(raw$counts))
if (length(common_genes) < nrow(filtered$counts)) {
  cat(sprintf("[SoupX] WARNING: gene mismatch — filtered=%d, raw=%d, common=%d. Using common set.\n",
              nrow(filtered$counts), nrow(raw$counts), length(common_genes)))
  filtered$counts <- filtered$counts[common_genes, , drop = FALSE]
  raw$counts      <- raw$counts[common_genes, , drop = FALSE]
} else {
  cat(sprintf("[SoupX] gene intersection OK: %d genes\n", length(common_genes)))
}

# ---- 条形码匹配：filtered barcodes ⊆ raw barcodes ----
# cellranger 的 filtered 是 raw 的子集（只保留真实细胞）。如果 Python 侧导出的是
# QC 后进一步缩小的细胞集，barcodes 仍然是 raw 的子集。
filtered_bc <- colnames(filtered$counts)
raw_bc      <- colnames(raw$counts)

matched_bc <- intersect(filtered_bc, raw_bc)
if (length(matched_bc) == 0) {
  stop(sprintf("[SoupX] FATAL: no overlapping barcodes between filtered (%d) and raw (%d)",
               length(filtered_bc), length(raw_bc)))
}
if (length(matched_bc) < length(filtered_bc)) {
  cat(sprintf("[SoupX] barcode match: %d/%d filtered cells found in raw (%d missing)\n",
              length(matched_bc), length(filtered_bc), length(filtered_bc) - length(matched_bc)))
} else {
  cat(sprintf("[SoupX] barcode match OK: all %d filtered cells found in raw (%d total raw barcodes)\n",
              length(filtered_bc), length(raw_bc)))
}

# 对齐矩阵——只用匹配到的 cells
filtered_counts <- filtered$counts[, matched_bc, drop = FALSE]
raw_counts      <- raw$counts[, matched_bc, drop = FALSE]

# ---- 快速聚类（为 SoupX 提供 cluster labels） ----
# SoupX 利用 cluster 信息更好地区分"细胞类型特异性表达"与"环境 RNA 污染"。
# 用简单的 PCA + 图聚类给 filtered counts 做快速分群（不需精细生物学聚类）。
set.seed(42)
tryCatch({
  # log-normalize
  norm_counts <- log1p(t(t(filtered_counts) / Matrix::colSums(filtered_counts)) * 1e4)
  # 选 top 2000 高变基因（只选方差 >0 的基因，避免全零基因导致 PCA 失败）
  gene_vars <- apply(norm_counts, 1, var)
  hvg_candidates <- which(gene_vars > 0)
  if (length(hvg_candidates) < 10) stop("too few genes with non-zero variance")
  n_hvg <- min(2000, length(hvg_candidates))
  hvg <- names(sort(gene_vars[hvg_candidates], decreasing = TRUE))[1:n_hvg]
  # PCA
  pca_data <- t(norm_counts[hvg, , drop = FALSE])
  if (ncol(pca_data) == 0 || nrow(pca_data) == 0) stop("PCA input dimensions are zero")
  rank_val <- min(20, nrow(pca_data), ncol(pca_data))
  if (rank_val < 1) stop(sprintf("invalid PCA rank: %d (cells=%d, genes=%d)", rank_val, nrow(pca_data), ncol(pca_data)))
  pca <- prcomp(as.matrix(pca_data), center = TRUE, scale. = TRUE, rank. = rank_val)
  # kNN + 图聚类
  k_knn <- min(10, ncol(filtered_counts) - 1)
  knn <- RANN::nn2(pca$x, k = k_knn)
  adj <- matrix(0, ncol(filtered_counts), ncol(filtered_counts))
  for (i in seq_len(ncol(filtered_counts))) {
    adj[i, knn$nn.idx[i, ]] <- 1
  }
  adj <- (adj + t(adj)) / 2
  adj[adj < 1] <- 0
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
  clusters <- igraph::cluster_louvain(g)$membership
  names(clusters) <- colnames(filtered_counts)
  cat(sprintf("[SoupX] quick clustering done: %d clusters\n", length(unique(clusters))))
}, error = function(e) {
  cat(sprintf("[SoupX] quick clustering failed (%s), using uniform clusters\n", e$message))
  clusters <<- NULL
})

# ---- SoupX 核心流程 ----
cat("[SoupX] creating SoupChannel...\n")
sc <- SoupChannel(raw_counts, filtered_counts)

if (!is.null(clusters) && length(clusters) > 0) {
  sc <- setClusters(sc, clusters)
  cat(sprintf("[SoupX] setClusters: %d clusters\n", length(unique(clusters))))
} else {
  cat("[SoupX] setClusters: skipped (clustering unavailable), autoEstCont will estimate without cluster info\n")
}

# autoEstCont 容错：先用默认参数，若失败则逐步放宽条件
# SoupX 在某些数据上（基因极稀疏、大量全零基因）会因 soupProf$est 产生
# NaN 而报错 "missing values and NaN's not allowed"。渐进式回退策略：
#   1. 默认参数
#   2. 降低 tfidfMin（SoupX vignette 推荐放宽值）+ soupQuantile（放宽标记基因筛选门槛）
#   3. forceAccept=TRUE, priorRho=0.05（10x 数据典型污染比例经验值，跳过自动判断直接用 priorRho）
cat("[SoupX] running autoEstCont (estimating contamination fraction)...\n")
sc <- tryCatch({
  autoEstCont(sc)
}, error = function(e1) {
  cat(sprintf("[SoupX] default autoEstCont failed: %s\n", e1$message))
  cat("[SoupX] retrying with tfidfMin=0.5, soupQuantile=0.5...\n")
  tryCatch({
    autoEstCont(sc, tfidfMin = 0.5, soupQuantile = 0.5)
  }, error = function(e2) {
    cat(sprintf("[SoupX] retry failed: %s\n", e2$message))
    cat("[SoupX] retrying with forceAccept=TRUE, priorRho=0.05...\n")
    tryCatch({
      autoEstCont(sc, forceAccept = TRUE, priorRho = 0.05, doPlot = FALSE)
    }, error = function(e3) {
      cat(sprintf("[SoupX] all autoEstCont attempts failed: %s\n", e3$message))
      cat("[SoupX] WARNING: cannot estimate contamination, returning uncorrected counts\n")
      return(NULL)
    })
  })
})

if (is.null(sc)) {
  cat("[SoupX] autoEstCont failed, skipping correction\n")
  # 输出原始 filtered counts 作为"校正后"（即不做校正，但保持流程完整）
  corrected <- filtered_counts
} else {
  # 输出污染比例摘要
  rho <- sc$metaData$rho
  if (!is.null(rho)) {
    cat(sprintf("[SoupX] estimated contamination (rho): median=%.4f, mean=%.4f, max=%.4f\n",
                median(rho, na.rm = TRUE), mean(rho, na.rm = TRUE), max(rho, na.rm = TRUE)))
  }

  cat("[SoupX] running adjustCounts (correcting counts)...\n")
  corrected <- tryCatch({
    adjustCounts(sc)
  }, error = function(e) {
    cat(sprintf("[SoupX] adjustCounts failed: %s, returning filtered counts\n", e$message))
    return(filtered_counts)
  })
}

cat(sprintf("[SoupX] corrected counts: %d genes x %d cells\n", nrow(corrected), ncol(corrected)))

# ---- 写出校正矩阵 ----
# 写入 work_dir 下的三个文件，以便 Python 侧读回
if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)

corrected_mtx_fp <- file.path(work_dir, "corrected_counts.mtx")
bc_fp_out        <- file.path(work_dir, "barcodes.tsv")
feat_fp_out      <- file.path(work_dir, "features.tsv")

Matrix::writeMM(corrected, corrected_mtx_fp)
writeLines(colnames(corrected), bc_fp_out)

# features.tsv: 用 filtered 的 feature 信息（只保留基因名）
# 行顺序跟随 corrected 矩阵（即 common_genes 过滤后的 filtered features 顺序）
feat_out <- filtered$features[match(rownames(corrected), filtered$features$feature_name), , drop = FALSE]
write.table(feat_out, feat_fp_out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat(sprintf("[SoupX] wrote: %s, %s, %s\n", corrected_mtx_fp, bc_fp_out, feat_fp_out))

# ---- 快速摘要 ----
n_cells <- ncol(corrected)
n_genes <- nrow(corrected)
total_umi_before <- sum(filtered_counts)
total_umi_after  <- sum(corrected)
umi_retained <- 100 * total_umi_after / total_umi_before
cat(sprintf("[SoupX] summary: %d cells, %d genes, UMI retained: %.1f%% (before=%.0f, after=%.0f)\n",
            n_cells, n_genes, umi_retained, total_umi_before, total_umi_after))
