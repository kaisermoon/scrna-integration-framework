#!/usr/bin/env Rscript --vanilla
# deseq2_contrast.R — Pseudobulk DESeq2 差异表达分析
#
# 用法（由 notebook subprocess 调用）:
#   Rscript --vanilla scripts/deseq2_contrast.R <work_dir> <factor_col> <numerator> <denominator>
#
# 参数:
#   work_dir    — 存放 counts.csv / metadata.csv 的工作目录（也在此产出结果）
#   factor_col  — metadata 中用作分组因子的列名（如 "disease"）
#   numerator   — 对比分子（如 "CAG"）
#   denominator — 对比分母（如 "normal"）
#
# 输出:
#   {work_dir}/deg_{numerator}_vs_{denominator}.csv  — DESeq2 results 全表

# ---- 入参解析 ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("usage: Rscript --vanilla deseq2_contrast.R <work_dir> <factor_col> <numerator> <denominator>")
}
work_dir    <- args[1]
factor_col  <- args[2]
numerator   <- args[3]
denominator <- args[4]

cat(sprintf("[DESeq2] work_dir=%s  factor=%s  contrast=%s_vs_%s\n",
            work_dir, factor_col, numerator, denominator))

# ---- 静默加载包 ----
suppressPackageStartupMessages({
  library(DESeq2)
})

# ---- 读入 pseudobulk 聚合后的 counts 与 metadata ----
counts_fp   <- file.path(work_dir, "counts.csv")
metadata_fp <- file.path(work_dir, "metadata.csv")

if (!file.exists(counts_fp))   stop("counts.csv not found: ", counts_fp)
if (!file.exists(metadata_fp)) stop("metadata.csv not found: ", metadata_fp)

counts   <- as.matrix(read.csv(counts_fp,   row.names = 1, check.names = FALSE))
metadata <- read.csv(metadata_fp, row.names = 1, check.names = FALSE)

cat(sprintf("[DESeq2] loaded: %d genes x %d pseudobulk samples\n", nrow(counts), ncol(counts)))

# ---- 取出参与对比的样本 ----
if (!(factor_col %in% colnames(metadata))) {
  stop(sprintf("factor_col '%s' not found in metadata. Columns: %s",
               factor_col, paste(colnames(metadata), collapse = ", ")))
}
keep <- metadata[[factor_col]] %in% c(numerator, denominator)
if (sum(keep) < 2) {
  stop(sprintf("only %d sample(s) match contrast groups '%s' / '%s'", sum(keep), numerator, denominator))
}
metadata <- metadata[keep, , drop = FALSE]
counts   <- counts[, rownames(metadata), drop = FALSE]

# 因子列为字符时强制转换；level 顺序确保"对比方向 = numerator vs denominator"
if (!is.factor(metadata[[factor_col]])) {
  metadata[[factor_col]] <- factor(metadata[[factor_col]], levels = c(denominator, numerator))
} else {
  metadata[[factor_col]] <- relevel(metadata[[factor_col]], ref = denominator)
}

cat(sprintf("[DESeq2] samples in contrast: %d (%s=%d, %s=%d)\n",
            nrow(metadata), numerator, sum(metadata[[factor_col]] == numerator),
            denominator, sum(metadata[[factor_col]] == denominator)))

# ---- DESeq2 ----
# pseudobulk counts 已经是整数和的形式，DESeq2 直接用。
# 注意: 不按常规加入 size factor 归一化——pseudobulk 的 library size 反映
# 每个 (sample x cell-type) 组合中的细胞数差异，DESeq2 内部 median-of-ratios
# 归一化已能处理。
dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData   = metadata,
  design    = as.formula(paste0("~", factor_col))
)
dds <- DESeq(dds, quiet = TRUE)

# 提取结果：numerator vs denominator
res <- results(dds, contrast = c(factor_col, numerator, denominator))
res <- as.data.frame(res[order(res$pvalue), , drop = FALSE])

# 补充 mean expression (baseMean 已在 results 中)
colnames(res)[colnames(res) == "baseMean"] <- "baseMean"

# ---- 写出结果 ----
out_fp <- file.path(work_dir, sprintf("deg_%s_vs_%s.csv", numerator, denominator))
write.csv(res, out_fp, row.names = TRUE)
cat(sprintf("[DESeq2] wrote %d genes -> %s\n", nrow(res), out_fp))

# ---- 快速摘要 ----
sig_alpha <- 0.05
n_sig <- sum(res$padj < sig_alpha, na.rm = TRUE)
n_up  <- sum(res$padj < sig_alpha & res$log2FoldChange > 0, na.rm = TRUE)
n_dn  <- sum(res$padj < sig_alpha & res$log2FoldChange < 0, na.rm = TRUE)
cat(sprintf("[DESeq2] sig genes (padj<%.2f): %d (up=%d, down=%d)\n", sig_alpha, n_sig, n_up, n_dn))
