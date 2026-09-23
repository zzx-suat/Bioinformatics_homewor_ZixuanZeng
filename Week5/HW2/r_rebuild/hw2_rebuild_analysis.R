# =============================================================================
# Week 5 Homework 2 — 平台导出缺失产物的 R 重建
# Zixuan Zeng / SUAT24000114 · 2026-09-23
#
# 背景：EMP-Web 的 Sync 把 diff_analysis / dimension / enrichment 三个结果文件
#       导成了空文件（4 字节，内容为 ""）。本脚本用同一批源数据在 R 里重建
#       这些产物。
#
# 重要声明：
#   本脚本的输出是「用 R 独立重跑得到的结果」，不是 EasyMultiProfiler-Web
#   平台的输出。两者用的是同一份输入和同一个对比方向，但统计实现、默认参数
#   与版本不同，数值不保证完全一致。报告中引用时必须标明来源为本次 R 重建。
#
# 输入（均为原始文件，未改动）：
#   EasyMultiProfiler-Web/tests/RNAseq_output.csv    原始整数 counts, 24393 x 24
#   EasyMultiProfiler-Web/tests/RNAseq_mapping.csv   样本分组
#   week5/RNAseq_output_assay.csv                    仅用于取 EMP 保留的基因名
#
# 输出：week5/hw2_rebuild/
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(ggrepel)
})

TESTS <- "D:/BioLession in total/Bioinformatics/week1/EasyMultiProfiler-Web/tests"
W5    <- "D:/BioLession in total/Bioinformatics/week5"
OUT   <- file.path(W5, "hw2_rebuild")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# --- 1. 读入 -----------------------------------------------------------------
counts <- read.csv(file.path(TESTS, "RNAseq_output.csv"),
                   row.names = 1, check.names = FALSE)
meta   <- read.csv(file.path(TESTS, "RNAseq_mapping.csv"),
                   fileEncoding = "UTF-8-BOM", check.names = FALSE)
rownames(meta) <- meta$SampleID

cm <- as.matrix(counts)

# --- 2. 输入断言：出问题要大声失败，不要静默继续 -----------------------------
stopifnot(all(cm == round(cm)), all(cm >= 0))                  # 非负整数
stopifnot(identical(colnames(cm), rownames(meta)))             # 列名与元数据严格对齐
stopifnot(!anyDuplicated(colnames(cm)), !anyDuplicated(rownames(cm)))
cat(sprintf("输入: %d 基因 x %d 样本, 全为非负整数, 样本顺序一致\n",
            nrow(cm), ncol(cm)))

# --- 3. 特征集：直接采用 EMP 实际保留的基因，避免猜测其过滤规则 --------------
assay_emp <- read.csv(file.path(W5, "RNAseq_output_assay.csv"),
                      row.names = 1, check.names = FALSE)
keep_genes <- rownames(assay_emp)
stopifnot(all(keep_genes %in% rownames(cm)))
cm <- cm[keep_genes, , drop = FALSE]
cat(sprintf("特征集: 采用 EMP 保留的 %d 个基因（从平台导出的 assay 读取，非本脚本重新过滤）\n",
            nrow(cm)))

meta$Group <- factor(meta$Group)
cat("\n分组结构:\n"); print(table(meta$Group))

# =============================================================================
# 4. PCA —— 全部 24 个样本，展示整个析因设计的结构
# =============================================================================
dds_all <- DESeqDataSetFromMatrix(cm, meta, design = ~ Group)
vsd <- vst(dds_all, blind = TRUE)

pca <- prcomp(t(assay(vsd)), scale. = FALSE)
pv  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pdat <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                   Group = meta$Group, Sample = rownames(meta))
pdat$Drug  <- sub("\\+LIPUS$", "", as.character(pdat$Group))
pdat$LIPUS <- ifelse(grepl("\\+LIPUS$", pdat$Group), "LIPUS", "no LIPUS")

p_pca <- ggplot(pdat, aes(PC1, PC2, colour = Drug, shape = LIPUS)) +
  geom_point(size = 3.6, alpha = .9) +
  scale_colour_manual(values = c(DMSO = "#6B6478", T4400 = "#B4341F",
                                 T3976 = "#2F6F4E")) +
  labs(title = "PCA of VST-transformed counts (all 24 samples)",
       subtitle = sprintf("EMP-retained feature set, %d genes | rebuilt in R, not the platform's own output",
                          nrow(cm)),
       x = sprintf("PC1 (%.1f%%)", pv[1]), y = sprintf("PC2 (%.1f%%)", pv[2])) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right")
ggsave(file.path(OUT, "hw2_pca_all24.png"), p_pca,
       width = 8.2, height = 5.4, dpi = 200, bg = "white")
write.csv(cbind(pdat, pca$x[, 1:4]), file.path(OUT, "hw2_pca_coordinates.csv"),
          row.names = FALSE)
cat(sprintf("\nPCA: PC1 %.1f%%, PC2 %.1f%%\n", pv[1], pv[2]))

# =============================================================================
# 5. 差异表达 —— T4400 vs DMSO，参考水平 DMSO（与平台上更正后的方向一致）
# =============================================================================
sel  <- meta$Group %in% c("DMSO", "T4400")
meta_sub <- droplevels(meta[sel, , drop = FALSE])
meta_sub$Group <- relevel(factor(meta_sub$Group), ref = "DMSO")
cm_sub <- cm[, rownames(meta_sub), drop = FALSE]

cat(sprintf("\n差异分析子集: %d 样本 (%s)\n", ncol(cm_sub),
            paste(names(table(meta_sub$Group)), table(meta_sub$Group),
                  sep = "=", collapse = ", ")))

dds <- DESeqDataSetFromMatrix(cm_sub, meta_sub, design = ~ Group)
mm  <- model.matrix(design(dds), colData(dds))
stopifnot(qr(mm)$rank == ncol(mm))                     # 满秩检查
cat(sprintf("模型矩阵满秩: %d/%d\n", qr(mm)$rank, ncol(mm)))

dds <- DESeq(dds, quiet = TRUE)
cat("resultsNames: ", paste(resultsNames(dds), collapse = " | "), "\n")

coef_name <- "Group_T4400_vs_DMSO"
stopifnot(coef_name %in% resultsNames(dds))            # 不假设系数名，显式核对
res <- results(dds, name = coef_name)
res_shr <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)

PADJ <- 0.05; LFC <- 1
rdf <- as.data.frame(res_shr)
rdf$feature <- rownames(rdf)
rdf <- rdf[order(rdf$padj, na.last = TRUE), ]
rdf$significant <- with(rdf, !is.na(padj) & padj < PADJ & abs(log2FoldChange) >= LFC)
rdf$direction <- with(rdf, ifelse(!significant, "ns",
                           ifelse(log2FoldChange > 0, "up in T4400", "down in T4400")))

n_up   <- sum(rdf$direction == "up in T4400")
n_down <- sum(rdf$direction == "down in T4400")
n_na   <- sum(is.na(rdf$padj))
n_padj_only <- sum(!is.na(rdf$padj) & rdf$padj < PADJ & abs(rdf$log2FoldChange) < LFC)

cat(sprintf("\n阈值 padj < %.2f 且 |log2FC| >= %g\n", PADJ, LFC))
cat(sprintf("  上调 (T4400 更高): %d\n", n_up))
cat(sprintf("  下调 (T4400 更低): %d\n", n_down))
cat(sprintf("  合计显著          : %d\n", n_up + n_down))
cat(sprintf("  padj 为 NA        : %d\n", n_na))
cat(sprintf("  仅过 padj 未过效应量: %d\n", n_padj_only))

write.csv(rdf[, c("feature", "baseMean", "log2FoldChange", "lfcSE",
                  "pvalue", "padj", "significant", "direction")],
          file.path(OUT, "hw2_deseq2_results_T4400_vs_DMSO.csv"), row.names = FALSE)

# --- 方向核验：收缩后的符号 vs 原始 counts 组均值 ----------------------------
top <- head(rdf[rdf$significant, ], 20)
if (nrow(top) > 0) {
  g_dmso  <- rowMeans(cm_sub[top$feature, meta_sub$Group == "DMSO",  drop = FALSE])
  g_t4400 <- rowMeans(cm_sub[top$feature, meta_sub$Group == "T4400", drop = FALSE])
  chk <- data.frame(feature = top$feature,
                    log2FC_shrunk = top$log2FoldChange,
                    mean_DMSO = g_dmso, mean_T4400 = g_t4400,
                    raw_higher_in = ifelse(g_t4400 > g_dmso, "T4400", "DMSO"))
  chk$agrees <- (chk$log2FC_shrunk > 0) == (chk$raw_higher_in == "T4400")
  write.csv(chk, file.path(OUT, "hw2_direction_check.csv"), row.names = FALSE)
  cat(sprintf("方向核验: top%d 中 %d 个与原始组均值一致\n",
              nrow(chk), sum(chk$agrees)))
}

# =============================================================================
# 6. 火山图
# =============================================================================
vdf <- rdf[!is.na(rdf$padj), ]
lab <- head(vdf[vdf$significant, ], 15)

p_vol <- ggplot(vdf, aes(log2FoldChange, -log10(padj), colour = direction)) +
  geom_point(alpha = .55, size = 1.5) +
  geom_hline(yintercept = -log10(PADJ), linetype = "22", colour = "#6B6478") +
  geom_vline(xintercept = c(-LFC, LFC), linetype = "22", colour = "#6B6478") +
  ggrepel::geom_text_repel(data = lab, aes(label = feature), size = 2.9,
                           max.overlaps = 20, seed = 1, show.legend = FALSE) +
  scale_colour_manual(values = c("ns" = "#C9C5D1",
                                 "up in T4400" = "#B4341F",
                                 "down in T4400" = "#2F6F4E"), name = NULL) +
  labs(title = "T4400 versus DMSO — apeglm-shrunken log2 fold change",
       subtitle = sprintf("padj < %.2f and |log2FC| >= %g | up %d, down %d | rebuilt in R",
                          PADJ, LFC, n_up, n_down),
       x = "log2 fold change (T4400 / DMSO)", y = "-log10 adjusted p") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
ggsave(file.path(OUT, "hw2_volcano_T4400_vs_DMSO.png"), p_vol,
       width = 8.2, height = 6.2, dpi = 200, bg = "white")

# =============================================================================
# 7. 存档
# =============================================================================
saveRDS(dds, file.path(OUT, "hw2_dds_T4400_vs_DMSO.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "hw2_session_info.txt"))

summary_tbl <- data.frame(
  item = c("source_counts", "source_mapping", "feature_set", "n_features",
           "contrast", "reference_level", "shrinkage", "padj_threshold",
           "lfc_threshold", "n_up", "n_down", "n_significant", "n_padj_NA",
           "n_padj_only", "PC1_pct", "PC2_pct"),
  value = c("EasyMultiProfiler-Web/tests/RNAseq_output.csv",
            "EasyMultiProfiler-Web/tests/RNAseq_mapping.csv",
            "EMP-retained genes read from RNAseq_output_assay.csv",
            nrow(cm), "T4400 vs DMSO", "DMSO", "apeglm", PADJ, LFC,
            n_up, n_down, n_up + n_down, n_na, n_padj_only, pv[1], pv[2])
)
write.csv(summary_tbl, file.path(OUT, "hw2_rebuild_summary.csv"), row.names = FALSE)

cat("\n完成。产物写入 ", OUT, "\n", sep = "")
