# =============================================================================
# Week 6 — 04 差异丰度（受试者内配对）
#
# 方法选择说明：
#   阅读材料列了 ANCOM-BC2 / MaAsLin3 / ALDEx2 / LinDA 等。本机未安装这些
#   Bioconductor 包，且 Week 5 已验证在本环境大规模安装 Bioconductor 包不可靠。
#   故采用阅读材料同样列出的 "Wilcoxon / simple regression —— 透明的探索性分析"，
#   但配合两项关键处理，使其适配本数据的设计：
#     1. 在 CLR 尺度上做（处理组成性），不是相对丰度
#     2. 受试者内配对检验（处理非独立性），不是两独立样本检验
#   局限：Wilcoxon 本身不解决组成性，只是 CLR 先解决了；也不提供协变量调整。
#   这一点在报告局限一节明确写出。
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(ggrepel) })
set.seed(1)

PROJ <- "D:/BioLession in total/Bioinformatics/week6/week6_project"
d    <- readRDS(file.path(PROJ, "results", "01_clean_data.rds"))
o    <- readRDS(file.path(PROJ, "results", "03_clr_objects.rds"))
meta <- d$meta; CLR <- o$CLR; delta <- o$delta; sub_meta <- o$sub_meta

LOG <- file.path(PROJ, "results", "04_diff_log.txt")
con <- file(LOG, open = "wt", encoding = "UTF-8")
say <- function(...) { m <- paste0(...); cat(m, "\n"); writeLines(m, con) }

say("=========== Week 6 差异丰度 ===========")
say("时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("尺度: CLR；检验: 受试者内配对 Wilcoxon signed-rank；多重校正: BH")
say("特征数: ", ncol(CLR))
say("")

short <- function(x) {
  # 把 k__Bacteria;p__...;g__Xxx;s__yyy 压成最末两级
  parts <- strsplit(x, ";")[[1]]
  parts <- parts[parts != "__" & parts != ""]
  if (!length(parts)) return(x)
  paste(tail(parts, 2), collapse = ";")
}

# --- 1. 每个疾病内：before vs after（配对） ---------------------------------
res_all <- list()
for (dis in c("IBS", "UC")) {
  subs <- sub_meta$Subject[sub_meta$Disease == dis]
  rows_b <- sapply(subs, function(s) rownames(meta)[meta$Subject == s & meta$Timepoint == "before"])
  rows_a <- sapply(subs, function(s) rownames(meta)[meta$Subject == s & meta$Timepoint == "after"])
  B <- CLR[rows_b, , drop = FALSE]; A <- CLR[rows_a, , drop = FALSE]
  stopifnot(nrow(B) == nrow(A), length(subs) == nrow(B))

  out <- data.frame(taxon = colnames(CLR), stringsAsFactors = FALSE)
  out$n_pairs <- nrow(B)
  out$median_before <- apply(B, 2, median)
  out$median_after  <- apply(A, 2, median)
  out$median_delta  <- apply(A - B, 2, median)
  out$p <- sapply(seq_len(ncol(CLR)), function(j)
    tryCatch(wilcox.test(A[, j], B[, j], paired = TRUE, exact = FALSE)$p.value,
             error = function(e) NA_real_))
  out$p_adj <- p.adjust(out$p, method = "BH")
  out$disease <- dis
  out <- out[order(out$p), ]
  res_all[[dis]] <- out

  sig <- sum(out$p_adj < 0.05, na.rm = TRUE)
  say(sprintf("---- %s：before vs after（%d 对）----", dis, nrow(B)))
  say(sprintf("  BH 校正后 p.adj < 0.05 的 taxa: %d / %d", sig, nrow(out)))
  say(sprintf("  未校正 p < 0.05 的 taxa      : %d（按随机预期约 %.0f 个）",
              sum(out$p < 0.05, na.rm = TRUE), 0.05 * nrow(out)))
  top <- head(out, 8)
  for (i in seq_len(nrow(top))) {
    say(sprintf("   %-52s delta=%+6.3f  p=%.4f  p.adj=%.3f",
                substr(short(top$taxon[i]), 1, 52), top$median_delta[i],
                top$p[i], top$p_adj[i]))
  }
  say("")
}
diff_tbl <- bind_rows(res_all)
write.csv(diff_tbl, file.path(PROJ, "results", "04_paired_taxa.csv"), row.names = FALSE)

# --- 2. 变化量在应答组间是否不同（组间独立） --------------------------------
say("---- 受试者内变化量 delta 在 poor vs great 之间的差异 ----")
say("  每受试者一个 delta，组间独立，用 Mann-Whitney。")
resp_all <- list()
for (dis in c("IBS", "UC")) {
  sel <- sub_meta$Disease == dis
  dd <- delta[sel, , drop = FALSE]; mm <- droplevels(sub_meta[sel, ])
  if (n_distinct(mm$Response) < 2) next
  g1 <- dd[mm$Response == "great", , drop = FALSE]
  g0 <- dd[mm$Response == "poor",  , drop = FALSE]
  out <- data.frame(taxon = colnames(dd), stringsAsFactors = FALSE)
  out$n_great <- nrow(g1); out$n_poor <- nrow(g0)
  out$median_delta_great <- apply(g1, 2, median)
  out$median_delta_poor  <- apply(g0, 2, median)
  out$p <- sapply(seq_len(ncol(dd)), function(j)
    tryCatch(wilcox.test(g1[, j], g0[, j], exact = FALSE)$p.value,
             error = function(e) NA_real_))
  out$p_adj <- p.adjust(out$p, method = "BH")
  out$disease <- dis
  out <- out[order(out$p), ]
  resp_all[[dis]] <- out
  say(sprintf("  %s: BH 后显著 %d 个；未校正 p<0.05 的 %d 个（随机预期约 %.0f）",
              dis, sum(out$p_adj < 0.05, na.rm = TRUE),
              sum(out$p < 0.05, na.rm = TRUE), 0.05 * nrow(out)))
  top <- head(out, 5)
  for (i in seq_len(nrow(top))) {
    say(sprintf("   %-52s great%+6.3f poor%+6.3f  p=%.4f p.adj=%.3f",
                substr(short(top$taxon[i]), 1, 52),
                top$median_delta_great[i], top$median_delta_poor[i],
                top$p[i], top$p_adj[i]))
  }
}
resp_tbl <- bind_rows(resp_all)
write.csv(resp_tbl, file.path(PROJ, "results", "04_response_taxa.csv"), row.names = FALSE)
say("")

# --- 图：火山式（delta vs -log10 p），两个疾病并排 --------------------------
pd <- diff_tbl %>% filter(!is.na(p)) %>%
  mutate(sig = ifelse(p_adj < 0.05, "BH<0.05",
               ifelse(p < 0.05, "p<0.05 未过校正", "ns")),
         lab = sapply(taxon, short))
lab <- pd %>% group_by(disease) %>% slice_min(p, n = 6) %>% ungroup()

p1 <- ggplot(pd, aes(median_delta, -log10(p), colour = sig)) +
  geom_hline(yintercept = -log10(0.05), linetype = "22", colour = "#6B6478") +
  geom_vline(xintercept = 0, colour = "#9A9AA6") +
  geom_point(alpha = .75, size = 1.9) +
  geom_text_repel(data = lab, aes(label = substr(lab, 1, 30)), size = 2.5,
                  max.overlaps = 18, seed = 1, show.legend = FALSE) +
  facet_wrap(~ disease) +
  scale_colour_manual(values = c("BH<0.05" = "#B4341F",
                                 "p<0.05 未过校正" = "#C58A2E", "ns" = "#C9C5D1"),
                      name = NULL) +
  labs(title = "受试者内配对差异（CLR 尺度）",
       subtitle = "横轴为每个 taxon 的 after − before 中位变化；配对 Wilcoxon + BH 校正",
       x = "median Δ CLR (after − before)", y = "-log10 p（未校正）") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

ggsave(file.path(PROJ, "figures", "04_paired_diff.png"), p1,
       width = 10, height = 5.6, dpi = 200, bg = "white")

say("已写出 results/04_*.csv, figures/04_paired_diff.png")
close(con)
cat("\n04 完成\n")
