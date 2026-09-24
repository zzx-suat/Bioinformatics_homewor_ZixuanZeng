# =============================================================================
# Week 5 Homework 2 — 完整析因设计分析
#
# 作业原文要求 "perform a complete analysis of the selected dataset"。
# 平台上只跑了 T4400 vs DMSO 一个两两对比，而数据本身是 3 药 × 2 超声的
# 完整析因设计（6 组 × 4 重复）。本脚本把设计跑全。
#
# 三个层次的问题：
#   1. 药物主效应      T4400 vs DMSO, T3976 vs DMSO（无超声条件下）
#   2. 超声主效应      DMSO+LIPUS vs DMSO（在无药条件下，超声单独有没有作用）
#   3. 药物 × 超声交互  超声是否改变了药物的效应（平台界面无此选项）
#
# 输出：week5/hw2_rebuild/full_factorial/
# 声明：R 独立分析，不是 EasyMultiProfiler-Web 平台输出。
# =============================================================================

suppressPackageStartupMessages({ library(DESeq2); library(ggplot2) })

TESTS <- "D:/BioLession in total/Bioinformatics/week1/EasyMultiProfiler-Web/tests"
W5    <- "D:/BioLession in total/Bioinformatics/week5"
OUT   <- file.path(W5, "hw2_rebuild", "full_factorial")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

PADJ <- 0.05; LFC <- 1

counts <- read.csv(file.path(TESTS, "RNAseq_output.csv"), row.names = 1, check.names = FALSE)
meta   <- read.csv(file.path(TESTS, "RNAseq_mapping.csv"),
                   fileEncoding = "UTF-8-BOM", check.names = FALSE)
rownames(meta) <- meta$SampleID
cm <- as.matrix(counts)

stopifnot(all(cm == round(cm)), all(cm >= 0))
stopifnot(identical(colnames(cm), rownames(meta)))

# 沿用 EMP 保留的特征集，与平台及前一份重建保持一致
keep <- rownames(read.csv(file.path(W5, "RNAseq_output_assay.csv"),
                          row.names = 1, check.names = FALSE))
cm <- cm[keep, , drop = FALSE]
cat(sprintf("特征集 %d 个基因（沿用平台保留的集合）\n\n", nrow(cm)))

# --- 把 Group 拆成 Drug 与 LIPUS 两个因子 ------------------------------------
meta$Drug  <- factor(sub("\\+LIPUS$", "", meta$Group), levels = c("DMSO", "T4400", "T3976"))
meta$LIPUS <- factor(ifelse(grepl("\\+LIPUS$", meta$Group), "yes", "no"), levels = c("no", "yes"))
cat("设计结构：\n"); print(table(meta$Drug, meta$LIPUS)); cat("\n")

# =============================================================================
# 1. 成对对比 —— 每个对比单独建模，避免跨条件借用离散度估计引起的争议
# =============================================================================
pairwise <- list(
  list(name = "T4400_vs_DMSO",            a = "DMSO",        b = "T4400",       q = "药物 T4400 的效应（无超声）"),
  list(name = "T3976_vs_DMSO",            a = "DMSO",        b = "T3976",       q = "药物 T3976 的效应（无超声）"),
  list(name = "DMSOLIPUS_vs_DMSO",        a = "DMSO",        b = "DMSO+LIPUS",  q = "超声单独的效应（无药）"),
  list(name = "T4400LIPUS_vs_T4400",      a = "T4400",       b = "T4400+LIPUS", q = "超声叠加在 T4400 之上的效应"),
  list(name = "T3976LIPUS_vs_T3976",      a = "T3976",       b = "T3976+LIPUS", q = "超声叠加在 T3976 之上的效应")
)

summ <- list()
for (p in pairwise) {
  sel <- meta$Group %in% c(p$a, p$b)
  md  <- droplevels(meta[sel, , drop = FALSE])
  md$Grp <- relevel(factor(make.names(md$Group)), ref = make.names(p$a))
  dsub <- cm[, rownames(md), drop = FALSE]

  dds <- DESeqDataSetFromMatrix(dsub, md, design = ~ Grp)
  mm  <- model.matrix(design(dds), colData(dds))
  stopifnot(qr(mm)$rank == ncol(mm))
  dds <- DESeq(dds, quiet = TRUE)

  cf <- setdiff(resultsNames(dds), "Intercept")
  stopifnot(length(cf) == 1L)
  r <- lfcShrink(dds, coef = cf, type = "apeglm", quiet = TRUE)
  rd <- as.data.frame(r); rd$feature <- rownames(rd)

  sig  <- !is.na(rd$padj) & rd$padj < PADJ & abs(rd$log2FoldChange) >= LFC
  up   <- sum(sig & rd$log2FoldChange > 0)
  down <- sum(sig & rd$log2FoldChange < 0)

  rd <- rd[order(rd$padj, na.last = TRUE), ]
  write.csv(rd[, c("feature","baseMean","log2FoldChange","lfcSE","pvalue","padj")],
            file.path(OUT, paste0("de_", p$name, ".csv")), row.names = FALSE)

  summ[[length(summ) + 1L]] <- data.frame(
    contrast = p$name, question = p$q, n_samples = ncol(dsub),
    coef = cf, up = up, down = down, total = up + down,
    padj_NA = sum(is.na(rd$padj))
  )
  cat(sprintf("  %-22s %-28s 上调 %4d  下调 %4d  合计 %4d\n",
              p$name, p$q, up, down, up + down))
}
summary_tbl <- do.call(rbind, summ)
write.csv(summary_tbl, file.path(OUT, "contrast_summary.csv"), row.names = FALSE)

# =============================================================================
# 2. 交互项 —— 超声是否改变了药物的效应（平台界面做不了）
# =============================================================================
cat("\n交互模型 ~ Drug * LIPUS （全部 24 样本）\n")
dds_i <- DESeqDataSetFromMatrix(cm, meta, design = ~ Drug * LIPUS)
mm_i  <- model.matrix(design(dds_i), colData(dds_i))
stopifnot(qr(mm_i)$rank == ncol(mm_i))
cat(sprintf("  模型矩阵满秩 %d/%d\n", qr(mm_i)$rank, ncol(mm_i)))
dds_i <- DESeq(dds_i, quiet = TRUE)
cat("  resultsNames: ", paste(resultsNames(dds_i), collapse = " | "), "\n")

inter <- grep("^Drug.*\\.LIPUS", resultsNames(dds_i), value = TRUE)
int_rows <- list()
for (tn in inter) {
  ri <- results(dds_i, name = tn)
  n_sig <- sum(!is.na(ri$padj) & ri$padj < PADJ)
  rd <- as.data.frame(ri); rd$feature <- rownames(rd)
  rd <- rd[order(rd$padj, na.last = TRUE), ]
  write.csv(rd[, c("feature","baseMean","log2FoldChange","lfcSE","pvalue","padj")],
            file.path(OUT, paste0("interaction_", tn, ".csv")), row.names = FALSE)
  int_rows[[length(int_rows) + 1L]] <- data.frame(term = tn, n_padj_lt_0.05 = n_sig)
  cat(sprintf("  %-34s padj<0.05 的基因数: %d\n", tn, n_sig))
}
write.csv(do.call(rbind, int_rows), file.path(OUT, "interaction_summary.csv"),
          row.names = FALSE)

# =============================================================================
# 3. 汇总图
# =============================================================================
pd <- summary_tbl
pd$contrast <- factor(pd$contrast, levels = rev(pd$contrast))
pl <- data.frame(contrast = rep(pd$contrast, 2),
                 n = c(pd$up, -pd$down),
                 dir = rep(c("上调", "下调"), each = nrow(pd)))

p <- ggplot(pl, aes(contrast, n, fill = dir)) +
  geom_col(width = .62) +
  geom_hline(yintercept = 0, colour = "#4A4A55") +
  geom_text(aes(label = abs(n)), hjust = ifelse(pl$n >= 0, -0.25, 1.25), size = 3.2) +
  coord_flip() +
  scale_fill_manual(values = c("上调" = "#B4341F", "下调" = "#2F6F4E"), name = NULL) +
  labs(title = "完整析因设计下的五个成对对比",
       subtitle = sprintf("padj < %.2f 且 |log2FC| >= %g | %d 个基因 | R 独立分析，非平台输出",
                          PADJ, LFC, nrow(cm)),
       x = NULL, y = "显著基因数（左=下调，右=上调）") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
ggsave(file.path(OUT, "contrast_summary.png"), p, width = 8.6, height = 4.8,
       dpi = 200, bg = "white")

writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info.txt"))
cat("\n完成。产物写入 ", OUT, "\n", sep = "")
