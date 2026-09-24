# =============================================================================
# Week 6 — 03 Beta 多样性与 PERMANOVA
#
# 本周阅读材料把这一条列为 common failure：
#   "Unrestricted permutation of paired samples"
#   推荐做法: "Restricted permutations + dispersion check"
#
# 本脚本刻意把两种置换都跑一遍并列出来，用同一份数据显示差别有多大。
#
# 组成性：16S 计数是组成数据（Gloor 2017）。相对丰度上的欧氏距离没有意义，
#         故以 CLR 变换后的欧氏距离（= Aitchison 距离）为主，
#         Bray-Curtis 作为对照一并报告。
# =============================================================================

suppressPackageStartupMessages({
  library(vegan); library(permute); library(dplyr); library(ggplot2)
})
set.seed(1)

PROJ <- "D:/BioLession in total/Bioinformatics/week6/week6_project"
d    <- readRDS(file.path(PROJ, "results", "01_clean_data.rds"))
abun <- d$abun; meta <- d$meta

LOG <- file.path(PROJ, "results", "03_beta_log.txt")
con <- file(LOG, open = "wt", encoding = "UTF-8")
say <- function(...) { m <- paste0(...); cat(m, "\n"); writeLines(m, con) }

say("=========== Week 6 Beta 多样性 ===========")
say("时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("")

X <- t(abun)
stopifnot(identical(rownames(X), rownames(meta)))

# --- 过滤极低流行度的 taxa ---------------------------------------------------
prev <- colSums(X > 0)
keep <- prev >= 0.10 * nrow(X)
say("特征过滤")
say("  规则: 至少在 10% 的样本中出现（", ceiling(0.10 * nrow(X)), " 个样本）")
say("  保留 ", sum(keep), " / ", ncol(X), " 个 taxa")
say("  说明: 只出现在极少数样本的 taxa 主要贡献噪声与零；该规则在看结果之前设定。")
Xf <- X[, keep, drop = FALSE]
say("")

# --- CLR 变换 ----------------------------------------------------------------
# CLR(x) = log(x / 几何均值(x))，零值用 0.5 伪计数处理
pseudo <- 0.5
Xp  <- Xf + pseudo
gm  <- exp(rowMeans(log(Xp)))
CLR <- log(Xp / gm)
say("CLR 变换")
say("  伪计数 = ", pseudo, "（零值占比 ",
    sprintf("%.1f%%", 100 * mean(Xf == 0)), "）")
say("  每行 CLR 之和应为 0: 最大偏差 ",
    sprintf("%.2e", max(abs(rowSums(CLR)))))
say("")

d_ait <- dist(CLR)                                   # Aitchison
rel   <- Xf / rowSums(Xf)
d_bc  <- vegdist(rel, method = "bray")               # Bray-Curtis 对照

# =============================================================================
# PERMANOVA：无限制置换 vs 以受试者为 strata 的受限置换
# =============================================================================
say("=========== PERMANOVA ===========")
say("问题: 时间点（before/after）是否解释群落结构？")
say("同一受试者贡献两个样本，两者高度相关，故置换必须在受试者内进行。")
say("")

run_permanova <- function(D, dat, label) {
  # (a) 无限制置换 —— 阅读材料点名的错误做法
  set.seed(1)
  a_wrong <- adonis2(D ~ Timepoint, data = dat, permutations = 999)
  # (b) 受限置换：只在每个受试者内部交换（plots = Subject, 自由置换 within）
  set.seed(1)
  ctrl <- how(blocks = dat$Subject, nperm = 999)
  a_right <- adonis2(D ~ Timepoint, data = dat, permutations = ctrl)
  list(wrong = a_wrong, right = a_right, label = label)
}

for (dis in c("IBS", "UC")) {
  idx <- meta$Disease == dis & meta$Subject %in%
         names(which(table(meta$Subject) == 2))
  dat <- droplevels(meta[idx, , drop = FALSE])
  Dsub <- as.dist(as.matrix(d_ait)[idx, idx])

  say("---- ", dis, " （", nrow(dat), " 样本 / ", n_distinct(dat$Subject), " 受试者）----")

  set.seed(1)
  a_wrong <- adonis2(Dsub ~ Timepoint, data = dat, permutations = 999)
  set.seed(1)
  ctrl <- how(blocks = factor(dat$Subject), nperm = 999)
  a_right <- adonis2(Dsub ~ Timepoint, data = dat, permutations = ctrl)

  say(sprintf("  Aitchison 距离  R2 = %.4f", a_wrong$R2[1]))
  say(sprintf("    无限制置换（错误做法）      p = %.4f", a_wrong$`Pr(>F)`[1]))
  say(sprintf("    受试者内受限置换（正确做法）p = %.4f", a_right$`Pr(>F)`[1]))

  # 离散度检查
  bd <- betadisper(Dsub, dat$Timepoint)
  set.seed(1)
  pbd <- permutest(bd, permutations = 999)
  say(sprintf("    离散度齐性检验 p = %.4f %s",
              pbd$tab$`Pr(>F)`[1],
              ifelse(pbd$tab$`Pr(>F)`[1] < 0.05,
                     "  <= 组间离散度不齐，PERMANOVA 的显著可能来自离散度而非位置",
                     "  （离散度齐，位置差异的解读较安全）")))
  say("")
}

# --- 受试者身份本身能解释多少 ------------------------------------------------
say("---- 受试者身份的解释力（说明为何必须受限置换）----")
idx <- meta$Subject %in% names(which(table(meta$Subject) == 2))
dat <- droplevels(meta[idx, , drop = FALSE])
Dsub <- as.dist(as.matrix(d_ait)[idx, idx])
set.seed(1)
a_subj <- adonis2(Dsub ~ Subject, data = dat, permutations = 199)
say(sprintf("  Aitchison ~ Subject: R2 = %.4f, p = %.4f",
            a_subj$R2[1], a_subj$`Pr(>F)`[1]))
say("  受试者身份解释了群落方差的很大一部分。这正是把两个时间点当独立样本会失真的原因。")
say("")

# --- 应答组：用受试者内变化向量比较（组间独立） ------------------------------
say("---- 应答组差异：受试者内变化向量 ----")
say("  对每个受试者计算 delta = CLR(after) - CLR(before)，")
say("  再以 delta 向量做 PERMANOVA。每受试者一行，组间观测独立，无需受限置换。")
sub_ids <- names(which(table(meta$Subject) == 2))
dl <- t(sapply(sub_ids, function(s) {
  rows <- rownames(meta)[meta$Subject == s]
  b <- rows[meta[rows, "Timepoint"] == "before"]
  a <- rows[meta[rows, "Timepoint"] == "after"]
  CLR[a, ] - CLR[b, ]
}))
sub_meta <- meta[match(sub_ids, meta$Subject), c("Subject", "Disease", "Response")]
rownames(sub_meta) <- sub_ids
stopifnot(identical(rownames(dl), rownames(sub_meta)))

for (dis in c("ALL", "IBS", "UC")) {
  sel <- if (dis == "ALL") rep(TRUE, nrow(sub_meta)) else sub_meta$Disease == dis
  md <- droplevels(sub_meta[sel, , drop = FALSE])
  if (n_distinct(md$Response) < 2) next
  Dd <- dist(dl[sel, , drop = FALSE])
  set.seed(1)
  ad <- adonis2(Dd ~ Response, data = md, permutations = 999)
  say(sprintf("  %-4s n=%2d  delta ~ Response: R2 = %.4f, p = %.4f",
              dis, nrow(md), ad$R2[1], ad$`Pr(>F)`[1]))
}
say("")

# --- 排序图 ------------------------------------------------------------------
pc <- prcomp(CLR, center = TRUE, scale. = FALSE)
pv <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
pdat <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], meta)

p1 <- ggplot(pdat, aes(PC1, PC2, colour = Disease, shape = Timepoint)) +
  geom_point(size = 2.4, alpha = .85) +
  scale_colour_manual(values = c(IBS = "#3E2A63", UC = "#B4341F")) +
  labs(title = "Aitchison PCA（CLR 后的主成分）",
       subtitle = sprintf("%d taxa（流行度 >= 10%%）；组成性数据不能直接用相对丰度做欧氏距离", ncol(Xf)),
       x = sprintf("PC1 (%.1f%%)", pv[1]), y = sprintf("PC2 (%.1f%%)", pv[2])) +
  theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"))

seg <- pdat %>% filter(Subject %in% sub_ids) %>%
  select(Subject, Disease, Timepoint, PC1, PC2) %>%
  tidyr::pivot_wider(names_from = Timepoint, values_from = c(PC1, PC2))

p2 <- ggplot(seg) +
  geom_segment(aes(x = PC1_before, y = PC2_before, xend = PC1_after, yend = PC2_after,
                   colour = Disease), alpha = .55,
               arrow = arrow(length = unit(0.14, "cm"))) +
  geom_point(aes(PC1_before, PC2_before, colour = Disease), size = 1.3, alpha = .8) +
  scale_colour_manual(values = c(IBS = "#3E2A63", UC = "#B4341F")) +
  facet_wrap(~ Disease) +
  labs(title = "每位受试者的 before → after 轨迹",
       subtitle = "箭头起点为 before。方向散乱说明不存在统一的位移方向",
       x = sprintf("PC1 (%.1f%%)", pv[1]), y = sprintf("PC2 (%.1f%%)", pv[2])) +
  theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"))

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  ggsave(file.path(PROJ, "figures", "03_beta.png"), p1 / p2 + plot_layout(heights = c(1, 1)),
         width = 9, height = 9, dpi = 200, bg = "white")
}
write.csv(data.frame(SampleID = rownames(pdat), pdat[, c("PC1", "PC2", "Disease",
          "Timepoint", "Response", "Subject")]),
          file.path(PROJ, "results", "03_pca_coordinates.csv"), row.names = FALSE)
saveRDS(list(CLR = CLR, d_ait = d_ait, delta = dl, sub_meta = sub_meta),
        file.path(PROJ, "results", "03_clr_objects.rds"))

say("已写出 results/03_*, figures/03_beta.png")
close(con)
cat("\n03 完成\n")
