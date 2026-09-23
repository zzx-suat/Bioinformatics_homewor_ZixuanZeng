# ================================================================
# Week 5 Homework 1 — Bulk RNA-seq differential expression (DESeq2)
# Course : Bioinformatics: From Multi-Omics Data to Discovery
# Student: Zixuan Zeng (曾梓轩)
# Date   : 2026-09-22
#
# 运行方式：把本文件与两个输入 CSV 放在同一个文件夹，然后
#   - RStudio: 直接 Source
#   - 终端   : Rscript week5_deseq2_analysis.R
# 本脚本不依赖 rstudioapi（课程 starter 里的 getActiveDocumentContext()
# 在 RStudio 之外必定报错），路径自动探测。
# ================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(ggplot2)
  library(ggrepel)
})

set.seed(20260922)   # 可复现性：本脚本无随机步骤，仍固定种子以备后续扩展

# ----------------------------------------------------------------
# 1. 输入路径（自动探测，不覆盖原始输入文件）
# ----------------------------------------------------------------
find_inputs <- function() {
  cands <- c(".", "..", "Homework/for_student", "../Homework/for_student",
             "../../Homework/for_student")
  for (d in cands) {
    f <- file.path(d, "Week5_Homework_Count_Matrix.csv")
    if (file.exists(f)) return(normalizePath(d))
  }
  stop("找不到 Week5_Homework_Count_Matrix.csv，请把脚本放到与输入文件同级的目录。")
}
in_dir  <- find_inputs()
out_dir <- getwd()
count_file    <- file.path(in_dir, "Week5_Homework_Count_Matrix.csv")
metadata_file <- file.path(in_dir, "Week5_Homework_Sample_Metadata.csv")
message("输入目录: ", in_dir)

# ----------------------------------------------------------------
# 2. 导入
# ----------------------------------------------------------------
counts  <- read.csv(count_file,    row.names = 1, check.names = FALSE)
coldata <- read.csv(metadata_file, row.names = 1, check.names = FALSE)

# ----------------------------------------------------------------
# 3. 强制验证（对应 Verification Checklist 的 Inputs 段）
# ----------------------------------------------------------------
cmat <- as.matrix(counts)
stopifnot(ncol(counts) == nrow(coldata))
stopifnot(identical(colnames(counts), rownames(coldata)))   # 列名与行名同序
stopifnot(!any(duplicated(colnames(counts))))               # 无重复样本 ID
stopifnot(all(cmat >= 0))                                   # 非负
stopifnot(all(cmat == round(cmat)))                         # 整数

coldata$condition <- relevel(factor(coldata$condition), ref = "control")
coldata$batch     <- factor(coldata$batch)
stopifnot(levels(coldata$condition)[1] == "control")        # control 为参考水平

cat("\n== 设计表 ==\n");        print(table(coldata$batch, coldata$condition))
cat("\n== 文库大小 ==\n");      print(summary(colSums(counts)))
cat("\n== 因子水平 ==\n")
cat("condition:", paste(levels(coldata$condition), collapse = " | "), "\n")
cat("batch    :", paste(levels(coldata$batch),     collapse = " | "), "\n")

# ----------------------------------------------------------------
# 4. 构建 DESeq2 对象
# ----------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = coldata,
                              design    = ~ batch + condition)

# --- 学生判断 J1：为什么把 batch 写进设计公式 -------------------
# （原始判断见 学生独立判断记录_J1J2J3.md，此处为修订版）
# 样本分 A/B/C 三批完成，批次间存在试剂、仪器、操作条件造成的系统性差异。
# 这部分差异真实存在但不是本研究关心的对象；若不纳入模型，它会进入残差、
# 抬高离散度估计，从而降低检出真实处理效应的能力。
# 本设计每批内 control 与 treated 各 2 个，batch 与 condition 完全平衡、
# 不混杂，因此两者的效应可以分别估计，比较 condition 时能扣除批次项。
# 若某一批只含单一条件，两个效应将无法区分，模型矩阵也不满秩。
stopifnot(qr(model.matrix(design(dds), colData(dds)))$rank ==
            ncol(model.matrix(design(dds), colData(dds))))  # 满秩检查

# ----------------------------------------------------------------
# 5. 预过滤
# ----------------------------------------------------------------
# --- 学生判断 J2：过滤规则与阈值依据 ----------------------------
# （原始判断见 学生独立判断记录_J1J2J3.md，此处为修订版）
# 过滤对象是基因而非样本：仅保留在至少 3 个样本中达到 10 counts 的基因，
# 12 个样本全部保留。目的是移除近乎不表达的基因以减少检验个数，减轻 BH
# 校正的惩罚。阈值不得高于所关心的最小组的样本数：control 与 treated 各
# 6 个，若设为 7 则仅在单一条件下表达的基因会被误删；最小的 batch×condition
# 单元为 2 个。3 位于 2 与 6 之间，既滤掉偶然出现的基因，又不误删真实差异基因。
n_before <- nrow(dds)
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
cat(sprintf("\n== 过滤 ==\n规则: >=10 counts 于 >=3 个样本\n保留 %d / %d 个基因（滤除 %d）\n",
            nrow(dds), n_before, n_before - nrow(dds)))

# ----------------------------------------------------------------
# 6. 拟合与系数检查
# ----------------------------------------------------------------
dds <- DESeq(dds)
cat("\n== resultsNames(dds) ==\n"); print(resultsNames(dds))

target_coef <- "condition_treated_vs_control"
stopifnot(target_coef %in% resultsNames(dds))   # 不假设系数名，先核对再使用

# ----------------------------------------------------------------
# 7. 提取对比 + apeglm 收缩
# ----------------------------------------------------------------
res    <- results(dds, name = target_coef, alpha = 0.05)
resLFC <- lfcShrink(dds, coef = target_coef, type = "apeglm")
cat("\n== results() 元信息 ==\n"); print(mcols(res)$description)

# 收缩不改变 padj —— 明确验证这一点
stopifnot(identical(is.na(res$padj), is.na(resLFC$padj)))
stopifnot(all.equal(res$padj[!is.na(res$padj)], resLFC$padj[!is.na(resLFC$padj)]))

PADJ_CUT <- 0.05
LFC_CUT  <- 1
df <- as.data.frame(resLFC)
df$gene_id <- rownames(df)
df$significant <- !is.na(df$padj) & df$padj < PADJ_CUT & abs(df$log2FoldChange) >= LFC_CUT
df$direction <- ifelse(!df$significant, "not significant",
                       ifelse(df$log2FoldChange > 0, "up in treated", "down in treated"))

n_na   <- sum(is.na(df$padj))
n_up   <- sum(df$direction == "up in treated")
n_down <- sum(df$direction == "down in treated")
cat(sprintf("\n== 显著性（padj < %.2f 且 |log2FC| >= %g）==\n上调 %d ｜ 下调 %d ｜ 合计 %d\npadj 为 NA 的基因: %d\n",
            PADJ_CUT, LFC_CUT, n_up, n_down, n_up + n_down, n_na))

# ----------------------------------------------------------------
# 8a. PCA（VST 转换）
# ----------------------------------------------------------------
vsd <- if (nrow(dds) >= 1000) vst(dds, blind = TRUE) else
  varianceStabilizingTransformation(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = c("condition", "batch"), returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))

p_pca <- ggplot(pca, aes(PC1, PC2, colour = condition, shape = batch)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(aes(label = name), size = 3, show.legend = FALSE, max.overlaps = 20) +
  labs(title = "PCA of VST-transformed counts (all retained genes)",
       subtitle = "colour = condition, shape = batch",
       x = paste0("PC1: ", pv[1], "% variance"),
       y = paste0("PC2: ", pv[2], "% variance")) +
  theme_bw(base_size = 12)
ggsave(file.path(out_dir, "week5_pca.png"), plot = p_pca,
       width = 7, height = 5.2, dpi = 300)

# ----------------------------------------------------------------
# 8b. 火山图
# ----------------------------------------------------------------
# --- 学生判断 J3：火山图三个量的选择 ----------------------------
# x 轴 = apeglm 收缩后的 log2FoldChange（收缩抑制低表达基因的虚高效应量，
#        且不改变 padj，已在第 7 节用 stopifnot 验证）
# y 轴 = -log10(padj)，使用 BH 校正后的 p 值而非 raw p value
# 阈值线 = 水平 -log10(0.05) = 1.3010；垂直 log2FC = ±1
#          "显著"要求两个条件同时满足
vol <- df[!is.na(df$padj), ]
vol$negLog10Padj <- -log10(vol$padj)
lab <- head(vol[order(vol$padj), ][vol[order(vol$padj), ]$significant, ], 10)

p_vol <- ggplot(vol, aes(log2FoldChange, negLog10Padj, colour = direction)) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(PADJ_CUT),     linetype = "dashed", colour = "grey40") +
  geom_text_repel(data = lab, aes(label = gene_id), size = 3,
                  show.legend = FALSE, max.overlaps = 20) +
  scale_colour_manual(values = c("up in treated" = "#C0392B",
                                 "down in treated" = "#2471A3",
                                 "not significant" = "grey75")) +
  annotate("text", x = Inf, y = -log10(PADJ_CUT), hjust = 1.05, vjust = -0.6,
           size = 3, colour = "grey30",
           label = sprintf("padj = %.2f  (y = %.3f)", PADJ_CUT, -log10(PADJ_CUT))) +
  labs(title = "Volcano: treated versus control (reference = control)",
       subtitle = sprintf("apeglm-shrunken LFC; significance = padj < %.2f AND |log2FC| >= %g; %d up, %d down",
                          PADJ_CUT, LFC_CUT, n_up, n_down),
       x = "Shrunken log2 fold change (treated / control)",
       y = "-log10 adjusted p value (BH)",
       colour = NULL) +
  theme_bw(base_size = 12)
ggsave(file.path(out_dir, "week5_de_plot.png"), plot = p_vol,
       width = 7.2, height = 5.4, dpi = 300)

# ----------------------------------------------------------------
# 9. 导出结果（含非显著基因的完整表）
# ----------------------------------------------------------------
out <- df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE",
              "pvalue", "padj", "significant", "direction")]
out <- out[order(out$padj, na.last = TRUE), ]
write.csv(out, file.path(out_dir, "week5_deseq2_results.csv"), row.names = FALSE)

top20 <- head(out[out$significant, ], 20)
write.csv(top20, file.path(out_dir, "week5_top20_features.csv"), row.names = FALSE)

# 方向核验所需的原始 counts（供人工比对，对应 checklist 的 direction 验证）
ctrl_s <- rownames(coldata)[coldata$condition == "control"]
trt_s  <- rownames(coldata)[coldata$condition == "treated"]
chk <- data.frame(
  gene_id          = top20$gene_id,
  shrunken_log2FC  = round(top20$log2FoldChange, 3),
  padj             = signif(top20$padj, 3),
  mean_raw_control = round(rowMeans(cmat[top20$gene_id, ctrl_s, drop = FALSE]), 1),
  mean_raw_treated = round(rowMeans(cmat[top20$gene_id, trt_s,  drop = FALSE]), 1)
)
chk$sign_consistent <- (chk$shrunken_log2FC > 0) ==
  (chk$mean_raw_treated > chk$mean_raw_control)
write.csv(chk, file.path(out_dir, "week5_direction_check.csv"), row.names = FALSE)
cat("\n== 方向一致性（收缩 LFC 符号 vs 原始 counts 组均值）==\n")
cat(sprintf("top20 中一致: %d / %d\n", sum(chk$sign_consistent), nrow(chk)))

# ----------------------------------------------------------------
# 10. 设计与检查审计表
# ----------------------------------------------------------------
audit <- data.frame(
  item = c("RNA class / assay", "value type entering DESeq2", "biological unit",
           "n samples", "conditions", "batches", "design formula",
           "reference level", "contrast extracted", "shrinkage",
           "filter rule", "genes before filter", "genes after filter",
           "colnames == rownames(metadata)", "duplicated sample IDs",
           "all counts non-negative integers", "model matrix full rank",
           "library size min/max", "padj cutoff", "abs(log2FC) cutoff",
           "significant up", "significant down", "padj = NA"),
  value = c("bulk gene-level RNA-seq (teaching synthetic data)",
            "raw integer counts (no TPM/CPM/z-score)",
            "one library per sample; 2 replicates per batch x condition",
            ncol(counts), paste(levels(coldata$condition), collapse = ", "),
            paste(levels(coldata$batch), collapse = ", "),
            "~ batch + condition", "control", target_coef, "apeglm",
            ">=10 counts in >=3 samples", n_before, nrow(dds),
            identical(colnames(counts), rownames(coldata)),
            any(duplicated(colnames(counts))),
            all(cmat >= 0) && all(cmat == round(cmat)), TRUE,
            paste(range(colSums(counts)), collapse = " / "),
            PADJ_CUT, LFC_CUT, n_up, n_down, n_na),
  stringsAsFactors = FALSE
)
write.table(audit, file.path(out_dir, "design_and_checks.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ----------------------------------------------------------------
# 11. 保存对象与环境
# ----------------------------------------------------------------
saveRDS(dds, file.path(out_dir, "week5_deseq2_object.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))

cat("\n== 完成 ==\n")
cat(paste(" -", c("week5_pca.png", "week5_de_plot.png", "week5_deseq2_results.csv",
                  "week5_top20_features.csv", "week5_direction_check.csv",
                  "design_and_checks.tsv", "week5_deseq2_object.rds",
                  "session_info.txt"), collapse = "\n"), "\n")
