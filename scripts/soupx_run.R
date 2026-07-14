#!/usr/bin/env Rscript --vanilla
# soupx_run.R — SoupX 环境 RNA 校正（去除 ambient RNA 污染）
#
# 用法（由 01_nancang notebook cell 49ea7e09 subprocess 调用）:
#   Rscript --vanilla scripts/soupx_run.R <work_dir> <filtered_mtx_dir> <raw_mtx_dir> <sample_id>
#
# 参数:
#   work_dir         — 产物目录（在此写出 corrected_counts.mtx + barcodes.tsv + features.tsv
#                      + rho.txt + soupx_status.json）
#   filtered_mtx_dir — Python 侧导出的过滤后计数矩阵目录
#                      （含 matrix.mtx, barcodes.tsv, features.tsv，非压缩）
#   raw_mtx_dir      — cellranger raw_feature_bc_matrix 目录
#                      （含空液滴，兼容 .gz 与非 .gz）
#   sample_id        — 样本标识（用于日志输出）
#
# 输出（成功路径，全部写入 work_dir）:
#   corrected_counts.mtx  — SoupX adjustCounts 后稀疏矩阵，行=基因、列=cell
#                           （行列顺序 == 输入 filtered）
#   barcodes.tsv          — 细胞条形码（顺序 == 输入 filtered barcode）
#   features.tsv          — 基因信息（三列无表头 tab，第 2 列 = 输入 var_names 顺序）
#   rho.txt               — 全局污染分数
#   soupx_status.json     — 状态文件（供 P1-b 判定成功/失败并写契约）
#
# 失败路径：不写 corrected_counts.mtx，写 soupx_status.json(status=failed)，
# stderr 打印 reason，quit(status=1)。
#
# 红线（reviewer 卡）:
#   1. tod（raw droplets）必须保全完整 barcode，绝不裁到 filtered barcode 交集
#   2. 估计/校正任何一步失败 → exit(1)，绝不把未校正矩阵当校正结果写出
#   3. 内存纪律：dgCMatrix 稀疏保持，不 as.matrix 全量 raw droplet
#   4. 注释中文讲 why

# ===========================================================================
# 段 A：入参解析 + fail() 辅助函数
# ===========================================================================

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

# 尽早建 work_dir，确保 fail() 能将 status.json 落盘
if (!dir.exists(work_dir)) {
  dir.create(work_dir, recursive = TRUE)
}

# fail() —— 统一失败出口
# 所有失败路径必须走这里：写 soupx_status.json(status=failed)、
# stderr 打印原因、quit(status=1)。绝不产出 corrected_counts.mtx。
# 为什么必须硬 fail（exit 1）而非静默跳过？
#   — 若 R 端失败却 exit 0，Python 侧靠"文件不存在"跳过只是次级保险；
#     主保险是 exit code，供 subprocess.check_call / returncode 判定。
#   — returned filtered_counts 冒充 corrected 在 I/O 形态上无差异
#     （都是 genes×cells matrix），Python 侧读回后无法自动识别，
#     导致下游分析建立在未校正数据上而毫无察觉。
fail <- function(reason) {
  status_file <- file.path(work_dir, "soupx_status.json")
  # 手工写 JSON（不依赖 jsonlite 包，降低跨环境故障面）
  # reason 中 \ 和 " 需转义保证 JSON 合法
  reason_escaped <- gsub("\\", "\\\\", reason, fixed = TRUE)
  reason_escaped <- gsub('"', '\\"', reason_escaped)
  json_lines <- c(
    "{",
    '  "status": "failed",',
    sprintf('  "reason": "%s",', reason_escaped),
    '  "rho_global": null,',
    '  "n_droplets": null,',
    '  "n_genes": null,',
    '  "n_common_genes": null,',
    '  "method_params": null',
    "}"
  )
  writeLines(json_lines, status_file)
  message(sprintf("[SoupX] FAIL: %s", reason))
  quit(save = "no", status = 1)
}

# ===========================================================================
# 段 B：包加载
# ===========================================================================

if (!requireNamespace("Matrix", quietly = TRUE)) {
  fail("Matrix package unavailable: install.packages('Matrix')")
}
if (!requireNamespace("SoupX", quietly = TRUE)) {
  fail("SoupX package unavailable: install.packages('SoupX')")
}
suppressPackageStartupMessages({
  library(Matrix)
  library(SoupX)
})

# ===========================================================================
# 段 C：read_10x_like 辅助函数
# ===========================================================================
# 读 10x 式 mtx 目录，兼容 .gz 与未压缩文件。
# 返回 list(counts=dgCMatrix, features=data.frame, barcodes=character)。
# 文件缺失 → fail()（统一走失败出口，不再 stop）。
read_10x_like <- function(mtx_dir, label) {
  find_file <- function(base) {
    fp_gz <- file.path(mtx_dir, paste0(base, ".gz"))
    fp    <- file.path(mtx_dir, base)
    if (file.exists(fp_gz)) return(fp_gz)
    if (file.exists(fp))    return(fp)
    fail(sprintf("%s: '%s' or '%s.gz' not found in %s", label, base, base, mtx_dir))
  }

  mtx_fp  <- find_file("matrix.mtx")
  bc_fp   <- find_file("barcodes.tsv")
  feat_fp <- find_file("features.tsv")

  counts <- Matrix::readMM(mtx_fp)
  barcodes <- readLines(bc_fp)
  features <- read.table(feat_fp, sep = "\t", header = FALSE,
                         col.names = c("feature_id", "feature_name", "feature_type"),
                         stringsAsFactors = FALSE, comment.char = "")

  # 行 = 基因 (features)，列 = 细胞/液滴 (barcodes)
  # readMM 输出 dgTMatrix，转为 dgCMatrix 以便后续稀疏操作
  # 绝不 as.matrix——raw 可达百万 barcode，稠密化必炸内存
  counts <- as(counts, "dgCMatrix")
  rownames(counts) <- features$feature_name
  colnames(counts) <- barcodes

  cat(sprintf("[SoupX] %s: loaded %d genes x %d barcodes\n",
              label, nrow(counts), ncol(counts)))
  return(list(counts = counts, features = features, barcodes = barcodes))
}

# ===========================================================================
# 段 D：读入 filtered 和 raw 矩阵
# ===========================================================================

filtered <- read_10x_like(filtered_mtx_dir, "filtered")
raw      <- read_10x_like(raw_mtx_dir, "raw")

# ===========================================================================
# 段 E：基因对齐（仅对齐基因维度，不对齐 barcode）
# ===========================================================================
# SoupX 要求 filtered 和 raw 的基因集一致。
# 只取共同基因交集对齐行维度（基因方向），columns（barcode）各保留各的。
# 为什么只对齐基因、不对齐 barcode？
#   — SoupX 靠 raw 中的空液滴（不含真实细胞的 droplet）估计环境 RNA 谱。
#     若把 raw 的 barcode 裁到与 filtered 交集，空液滴被丢弃，
#     污染估计的数据基础被破坏，校正结果不可靠。
#   — cellranger 语义：filtered barcodes ⊆ raw barcodes 天然成立；
#     只需检查 filtered cells 是否在 raw 全集中存在。

common_genes <- intersect(rownames(filtered$counts), rownames(raw$counts))
n_common_genes <- length(common_genes)

if (n_common_genes < 10) {
  fail(sprintf("too few common genes between raw and filtered: %d (need >=10)", n_common_genes))
}

if (n_common_genes < nrow(filtered$counts)) {
  cat(sprintf("[SoupX] gene mismatch — filtered=%d, raw=%d, common=%d. Using common set.\n",
              nrow(filtered$counts), nrow(raw$counts), n_common_genes))
} else {
  cat(sprintf("[SoupX] gene intersection OK: %d genes\n", n_common_genes))
}

# 两矩阵均按 common_genes 取子集（只对齐基因方向）
filtered$counts <- filtered$counts[common_genes, , drop = FALSE]
raw$counts      <- raw$counts[common_genes, , drop = FALSE]

# ===========================================================================
# 段 F：barcode 健全性检查 + 构建 toc/tod
# ===========================================================================
# 核心修复（P0 bug）：删除原 92-113 行把 raw 裁到 matched_bc 的逻辑。
# tod（SoupChannel 第一参数）必须保留完整 droplet barcode（含空液滴）。
# 只做健全性检查：filtered barcode 是否 ⊆ raw barcode（cellranger 语义）。

filtered_bc <- colnames(filtered$counts)
raw_bc      <- colnames(raw$counts)
n_droplets  <- length(raw_bc)       # 完整 raw droplet 数
n_cells_filtered <- length(filtered_bc)

cat(sprintf("[SoupX] barcode counts: raw=%d droplets, filtered=%d cells\n",
            n_droplets, n_cells_filtered))

# 检查 filtered 是否 ⊆ raw
overlap_bc <- intersect(filtered_bc, raw_bc)
n_overlap <- length(overlap_bc)

if (n_overlap == 0) {
  fail(sprintf("no overlapping barcodes: filtered cells not found in raw droplets (filtered=%d, raw=%d)",
               n_cells_filtered, n_droplets))
}

if (n_overlap < n_cells_filtered) {
  cat(sprintf("[SoupX] WARNING: %d/%d filtered cells not found in raw droplets\n",
              n_cells_filtered - n_overlap, n_cells_filtered))
  # 不 fail：Python 侧可能已 QC 过滤部分 cell，但剩余 cell 应在 raw 中以空液滴形式存在
}

# 构建 toc 和 tod
# toc = filtered cells 矩阵（genes=common_genes × cells=完整 filtered barcode）
# tod = 完整 raw droplet 矩阵（genes=common_genes × 完整 raw barcode，含空液滴）
# 关键：tod 不裁 barcode——空液滴是 SoupX 估计环境 RNA 谱的数据基础，
#       裁掉等于破坏估计前提。
toc <- filtered$counts
tod <- raw$counts

cat(sprintf("[SoupX] tod=%d genes x %d droplets (完整 raw，未裁 barcode)\n",
            nrow(tod), ncol(tod)))
cat(sprintf("[SoupX] toc=%d genes x %d cells (完整 filtered)\n",
            nrow(toc), ncol(toc)))

# ===========================================================================
# 段 G：quick clustering（仅在 toc / filtered cells 上做）
# ===========================================================================
# SoupX 利用 cluster 信息更好地区分"细胞类型特异性表达"与"环境 RNA 污染"。
# PCA + kNN + Louvain 快速分群（不需精细生物学聚类）。
# 为什么只在 toc 上做 PCA？
#   — tod 含百万级空液滴，稠密化 PCA 必炸内存；且空液滴本身无生物学分组意义。
#   — 聚类失败（如 RANN/igraph 未安装）可降级：SoupX 可无 cluster 估计。

set.seed(42)
tryCatch({
  norm_counts <- log1p(t(t(toc) / Matrix::colSums(toc)) * 1e4)
  gene_vars <- apply(norm_counts, 1, var)
  hvg_candidates <- which(gene_vars > 0)
  if (length(hvg_candidates) < 10) stop("too few genes with non-zero variance")
  n_hvg <- min(2000, length(hvg_candidates))
  hvg <- names(sort(gene_vars[hvg_candidates], decreasing = TRUE))[1:n_hvg]
  pca_data <- t(norm_counts[hvg, , drop = FALSE])
  if (ncol(pca_data) == 0 || nrow(pca_data) == 0) stop("PCA input dimensions are zero")
  rank_val <- min(20, nrow(pca_data), ncol(pca_data))
  if (rank_val < 1) stop(sprintf("invalid PCA rank: %d (cells=%d, genes=%d)",
                                   rank_val, nrow(pca_data), ncol(pca_data)))
  pca <- prcomp(as.matrix(pca_data), center = TRUE, scale. = TRUE, rank. = rank_val)
  k_knn <- min(10, ncol(toc) - 1)
  knn <- RANN::nn2(pca$x, k = k_knn)
  adj <- matrix(0, ncol(toc), ncol(toc))
  for (i in seq_len(ncol(toc))) {
    adj[i, knn$nn.idx[i, ]] <- 1
  }
  adj <- (adj + t(adj)) / 2
  adj[adj < 1] <- 0
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
  clusters <- igraph::cluster_louvain(g)$membership
  names(clusters) <- colnames(toc)
  cat(sprintf("[SoupX] quick clustering done: %d clusters\n", length(unique(clusters))))
}, error = function(e) {
  cat(sprintf("[SoupX] quick clustering failed (%s), using uniform clusters\n", e$message))
  clusters <<- NULL
})

# ===========================================================================
# 段 H：SoupChannel 构建
# ===========================================================================
# 第一参数 tod = 完整 raw droplet（含空液滴），第二参数 toc = filtered cells。
# 空液滴是环境 RNA 背景估计的唯一数据来源，绝不能裁剪。
cat("[SoupX] creating SoupChannel...\n")
sc <- SoupChannel(tod, toc)

if (!is.null(clusters) && length(clusters) > 0) {
  sc <- setClusters(sc, clusters)
  cat(sprintf("[SoupX] setClusters: %d clusters\n", length(unique(clusters))))
} else {
  cat("[SoupX] setClusters: skipped (clustering unavailable), autoEstCont will estimate without cluster info\n")
}

# ===========================================================================
# 段 I：autoEstCont 估计 + adjustCounts 校正
# ===========================================================================
# 三级渐进回退（合理估计鲁棒性策略，非"冒充成功"）：
#   1. 默认参数
#   2. 放宽 tfidfMin + soupQuantile
#   3. forceAccept=TRUE + priorRho=0.05
# 三级全失败 → fail() 硬退出，绝不返回未校正矩阵冒充成功。
# 为什么 fail 而非返回 filtered 冒充？
#   — returned filtered 在 I/O 形态上与 corrected 无法区分（都是 genes×cells）。
#     Python 侧读回后无自动检测手段，下游分析将建立在假校正数据上。
#   — 校正失败意味着污染估计完全不可靠，强行继续等于默认数据无污染，
#     这是比"不校正"更危险的隐性错误。

cat("[SoupX] running autoEstCont (estimating contamination fraction)...\n")
method_params <- "default"  # 记录实际生效的参数级别

sc <- tryCatch({
  autoEstCont(sc)
}, error = function(e1) {
  cat(sprintf("[SoupX] default autoEstCont failed: %s\n", e1$message))
  cat("[SoupX] retrying with tfidfMin=0.5, soupQuantile=0.5...\n")
  tryCatch({
    method_params <<- "relaxed"
    autoEstCont(sc, tfidfMin = 0.5, soupQuantile = 0.5)
  }, error = function(e2) {
    cat(sprintf("[SoupX] retry failed: %s\n", e2$message))
    cat("[SoupX] retrying with forceAccept=TRUE, priorRho=0.05...\n")
    tryCatch({
      method_params <<- "forceAccept"
      autoEstCont(sc, forceAccept = TRUE, priorRho = 0.05, doPlot = FALSE)
    }, error = function(e3) {
      fail(sprintf("autoEstCont failed after all fallback attempts: %s", e3$message))
    })
  })
})

# 提取污染比例估计
rho <- sc$metaData$rho
if (is.null(rho) || length(rho) == 0) {
  fail("autoEstCont succeeded but returned no rho estimates")
}

rho_global <- rho[1]  # per-gene rho 的第一值近似全局污染分数
cat(sprintf("[SoupX] estimated contamination (rho): median=%.4f, mean=%.4f, max=%.4f, global=%.4f\n",
            median(rho, na.rm = TRUE), mean(rho, na.rm = TRUE), max(rho, na.rm = TRUE), rho_global))

# 写 rho.txt（保留历史契约）
writeLines(as.character(rho_global), file.path(work_dir, "rho.txt"))
cat(sprintf("  rho (contamination fraction): %.4f\n", rho_global))

# adjustCounts 校正
cat("[SoupX] running adjustCounts (correcting counts)...\n")
corrected <- tryCatch({
  adjustCounts(sc)
}, error = function(e) {
  fail(sprintf("adjustCounts failed: %s", e$message))
})

cat(sprintf("[SoupX] corrected counts: %d genes x %d cells\n", nrow(corrected), ncol(corrected)))

# ===========================================================================
# 段 J：写出产物
# ===========================================================================
# 输出契约与 01_nancang notebook cell 49ea7e09 严格对齐：
#   - corrected_counts.mtx: 行=基因、列=cell（行列顺序 == 输入 filtered）
#     Python 侧读回 .T 得 cells×genes
#   - barcodes.tsv: 顺序 == 输入 filtered barcode（notebook 逐行比对）
#   - features.tsv: 三列无表头 tab，第 2 列 == 输入 var_names 顺序
#     （notebook 比对 iloc[:,1]）
#   - soupx_status.json: 供 P1-b 判定成功/失败并写契约

corrected_mtx_fp <- file.path(work_dir, "corrected_counts.mtx")
bc_fp_out        <- file.path(work_dir, "barcodes.tsv")
feat_fp_out      <- file.path(work_dir, "features.tsv")

Matrix::writeMM(corrected, corrected_mtx_fp)
writeLines(colnames(corrected), bc_fp_out)

# features.tsv: 行顺序跟随 corrected 矩阵（即 common_genes 过滤后的 filtered features 顺序）
feat_out <- filtered$features[match(rownames(corrected), filtered$features$feature_name), , drop = FALSE]
write.table(feat_out, feat_fp_out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat(sprintf("[SoupX] wrote: %s, %s, %s\n", corrected_mtx_fp, bc_fp_out, feat_fp_out))

# 写 soupx_status.json（成功路径）
n_genes <- nrow(corrected)
status_file <- file.path(work_dir, "soupx_status.json")
# 手工写 JSON 避免依赖 jsonlite 包
json_lines <- c(
  "{",
  '  "status": "success",',
  '  "reason": null,',
  sprintf('  "rho_global": %s,', as.character(rho_global)),
  sprintf('  "n_droplets": %d,', n_droplets),
  sprintf('  "n_cells_filtered": %d,', n_cells_filtered),
  sprintf('  "n_genes": %d,', n_genes),
  sprintf('  "n_common_genes": %d,', n_common_genes),
  sprintf('  "method_params": "%s"', method_params),
  "}"
)
writeLines(json_lines, status_file)
cat(sprintf("[SoupX] wrote: %s\n", status_file))

# ===========================================================================
# 段 K：摘要 + 正常退出
# ===========================================================================

total_umi_before <- sum(toc)
total_umi_after  <- sum(corrected)
umi_retained <- if (total_umi_before > 0) 100 * total_umi_after / total_umi_before else 0
cat(sprintf("[SoupX] summary: %d cells, %d genes, UMI retained: %.1f%% (before=%.0f, after=%.0f)\n",
            n_cells_filtered, n_genes, umi_retained, total_umi_before, total_umi_after))

# 显式 quit(status=0) 确保正常退出码
quit(save = "no", status = 0)
