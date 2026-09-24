# =============================================================================
# Week 6 — 02 Alpha 多样性
#
# 设计要点：同一受试者的 before/after 不是独立观测。
#   · 组间比较一律用受试者内配对（Wilcoxon signed-rank），不用两独立样本检验
#   · 应答组比较用「受试者内变化量」，每个受试者贡献一个数，因此彼此独立
#   · 另用 lme4 混合模型（Subject 随机截距）做一次对照，p 值走似然比检验
#
# 深度问题：样本深度相差 12.5 倍。丰富度类指标对深度极敏感，故先抽平到最小深度；
#           Shannon 对深度不敏感，两种都报，并说明差异。
# =============================================================================

suppressPackageStartupMessages({
  library(vegan); library(dplyr); library(ggplot2); library(lme4)
})
set.seed(1)

PROJ <- "D:/BioLession in total/Bioinformatics/week6/week6_project"
d    <- readRDS(file.path(PROJ, "results", "01_clean_data.rds"))
abun <- d$abun; meta <- d$meta

LOG <- file.path(PROJ, "results", "02_alpha_log.txt")
con <- file(LOG, open = "wt", encoding = "UTF-8")
say <- function(...) { m <- paste0(...); cat(m, "\n"); writeLines(m, con) }

say("=========== Week 6 Alpha 多样性 ===========")
say("时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("")

X <- t(abun)                                   # 样本 x taxa
stopifnot(identical(rownames(X), rownames(meta)))
depth <- rowSums(X)

# --- 抽平 --------------------------------------------------------------------
dmin <- min(depth)
say("抽平设置")
say("  最小深度 = ", dmin, "，全部样本抽平至该深度（无样本被丢弃）")
Xr <- rrarefy(X, dmin)
say("  抽平后每样本深度是否一致: ", length(unique(rowSums(Xr))) == 1)
say("")

alpha <- data.frame(
  SampleID  = rownames(X),
  depth     = as.integer(depth),
  Shannon   = diversity(X,  index = "shannon"),      # 原始计数
  Shannon_r = diversity(Xr, index = "shannon"),      # 抽平后
  Observed  = rowSums(X  > 0),
  Observed_r= rowSums(Xr > 0),
  Simpson_r = diversity(Xr, index = "simpson"),
  stringsAsFactors = FALSE
) %>% left_join(meta, by = "SampleID")

say("深度与 alpha 的相关（说明为何要抽平）")
say(sprintf("  Observed（原始） 对 depth 的 Spearman rho = %.3f",
            cor(alpha$Observed, alpha$depth, method = "spearman")))
say(sprintf("  Observed（抽平） 对 depth 的 Spearman rho = %.3f",
            cor(alpha$Observed_r, alpha$depth, method = "spearman")))
say(sprintf("  Shannon （原始） 对 depth 的 Spearman rho = %.3f",
            cor(alpha$Shannon, alpha$depth, method = "spearman")))
say(sprintf("  Shannon （抽平） 对 depth 的 Spearman rho = %.3f",
            cor(alpha$Shannon_r, alpha$depth, method = "spearman")))
say("  => 抽平后若相关大幅下降，说明原始丰富度差异有相当部分只是深度差异。")
say("")

write.csv(alpha, file.path(PROJ, "results", "02_alpha_table.csv"), row.names = FALSE)

# --- 受试者内配对：before -> after -------------------------------------------
say("检验 1 — 受试者内 before vs after（配对 Wilcoxon signed-rank）")
paired <- alpha %>% group_by(Subject) %>% filter(n() == 2) %>% ungroup()
say("  纳入受试者: ", n_distinct(paired$Subject), " 对（仅单次采样者已排除）")

pair_res <- list()
for (metric in c("Shannon_r", "Observed_r", "Simpson_r")) {
  for (dis in c("ALL", "IBS", "UC")) {
    sub <- if (dis == "ALL") paired else paired %>% filter(Disease == dis)
    w <- sub %>% select(Subject, Timepoint, val = all_of(metric)) %>%
      tidyr::pivot_wider(names_from = Timepoint, values_from = val) %>%
      filter(!is.na(before), !is.na(after))
    if (nrow(w) < 5) next
    tt <- wilcox.test(w$after, w$before, paired = TRUE, exact = FALSE)
    pair_res[[length(pair_res) + 1L]] <- data.frame(
      metric = metric, disease = dis, n_pairs = nrow(w),
      median_before = median(w$before), median_after = median(w$after),
      median_change = median(w$after - w$before),
      p = tt$p.value)
  }
}
pair_tbl <- bind_rows(pair_res)
pair_tbl$p_adj <- p.adjust(pair_tbl$p, method = "BH")
for (i in seq_len(nrow(pair_tbl))) {
  r <- pair_tbl[i, ]
  say(sprintf("  %-11s %-4s n=%2d  中位 %6.3f -> %6.3f  变化 %+6.3f  p=%.4f  p.adj=%.4f",
              r$metric, r$disease, r$n_pairs, r$median_before, r$median_after,
              r$median_change, r$p, r$p_adj))
}
write.csv(pair_tbl, file.path(PROJ, "results", "02_paired_before_after.csv"), row.names = FALSE)
say("")

# --- 受试者内变化量：poor vs great（组间独立） -------------------------------
say("检验 2 — 受试者内变化量在应答组间的差异（Mann-Whitney）")
say("  每个受试者贡献一个 delta = after - before，故组间观测互相独立。")

chg_res <- list()
for (metric in c("Shannon_r", "Observed_r", "Simpson_r")) {
  for (dis in c("ALL", "IBS", "UC")) {
    sub <- if (dis == "ALL") paired else paired %>% filter(Disease == dis)
    w <- sub %>% select(Subject, Response, Timepoint, val = all_of(metric)) %>%
      tidyr::pivot_wider(names_from = Timepoint, values_from = val) %>%
      filter(!is.na(before), !is.na(after)) %>%
      mutate(delta = after - before)
    if (n_distinct(w$Response) < 2) next
    g <- split(w$delta, w$Response)
    if (min(lengths(g)) < 5) next
    tt <- wilcox.test(g$great, g$poor, exact = FALSE)
    chg_res[[length(chg_res) + 1L]] <- data.frame(
      metric = metric, disease = dis,
      n_poor = length(g$poor), n_great = length(g$great),
      median_delta_poor = median(g$poor), median_delta_great = median(g$great),
      p = tt$p.value)
  }
}
chg_tbl <- bind_rows(chg_res)
chg_tbl$p_adj <- p.adjust(chg_tbl$p, method = "BH")
for (i in seq_len(nrow(chg_tbl))) {
  r <- chg_tbl[i, ]
  say(sprintf("  %-11s %-4s poor n=%2d delta中位 %+6.3f | great n=%2d delta中位 %+6.3f | p=%.4f p.adj=%.4f",
              r$metric, r$disease, r$n_poor, r$median_delta_poor,
              r$n_great, r$median_delta_great, r$p, r$p_adj))
}
write.csv(chg_tbl, file.path(PROJ, "results", "02_change_by_response.csv"), row.names = FALSE)
say("")

# --- 混合模型对照 ------------------------------------------------------------
say("检验 3 — 线性混合模型对照（Subject 随机截距，似然比检验）")
say("  模型: Shannon_r ~ Timepoint * Response + (1 | Subject)")
for (dis in c("IBS", "UC")) {
  sub <- paired %>% filter(Disease == dis)
  full <- lmer(Shannon_r ~ Timepoint * Response + (1 | Subject), data = sub, REML = FALSE)
  noint <- lmer(Shannon_r ~ Timepoint + Response + (1 | Subject), data = sub, REML = FALSE)
  notp  <- lmer(Shannon_r ~ Response + (1 | Subject), data = sub, REML = FALSE)
  a_int <- anova(full, noint); a_tp <- anova(noint, notp)
  vc <- as.data.frame(VarCorr(full))
  icc <- vc$vcov[vc$grp == "Subject"] / sum(vc$vcov)
  say(sprintf("  %s: 交互项 LRT p=%.4f | Timepoint 主效应 LRT p=%.4f | 受试者 ICC=%.3f",
              dis, a_int$`Pr(>Chisq)`[2], a_tp$`Pr(>Chisq)`[2], icc))
}
say("  ICC 越高，说明受试者身份解释的方差越大，把样本当独立处理的偏差就越严重。")
say("")

# --- 图 ----------------------------------------------------------------------
pd <- paired %>%
  select(Subject, Disease, Response, Timepoint, Shannon_r)

p1 <- ggplot(pd, aes(Timepoint, Shannon_r, group = Subject)) +
  geom_line(alpha = .35, colour = "#8A8596") +
  geom_point(aes(colour = Response), size = 1.8, alpha = .85) +
  facet_grid(Disease ~ Response) +
  scale_colour_manual(values = c(poor = "#B4341F", great = "#2F6F4E"), guide = "none") +
  labs(title = "Shannon 多样性：每条线是一个受试者",
       subtitle = "抽平至最小深度；配对结构直接画出来，避免把两个时间点当独立组",
       x = NULL, y = "Shannon (rarefied)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

delta <- pd %>% tidyr::pivot_wider(names_from = Timepoint, values_from = Shannon_r) %>%
  filter(!is.na(before), !is.na(after)) %>% mutate(delta = after - before)

p2 <- ggplot(delta, aes(Response, delta, fill = Response)) +
  geom_hline(yintercept = 0, linetype = "22", colour = "#6B6478") +
  geom_boxplot(width = .55, outlier.shape = NA, alpha = .75) +
  geom_jitter(width = .12, size = 1.6, alpha = .7) +
  facet_wrap(~ Disease) +
  scale_fill_manual(values = c(poor = "#B4341F", great = "#2F6F4E"), guide = "none") +
  labs(title = "受试者内变化量 Δ Shannon（after − before）",
       subtitle = "每点一个受试者，组间独立；这是可以做两组比较的量",
       x = NULL, y = "Δ Shannon") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  ggsave(file.path(PROJ, "figures", "02_alpha.png"), p1 / p2 + plot_layout(heights = c(1.15, 1)),
         width = 9, height = 9.5, dpi = 200, bg = "white")
} else {
  ggsave(file.path(PROJ, "figures", "02_alpha_lines.png"), p1, width = 9, height = 5, dpi = 200, bg = "white")
  ggsave(file.path(PROJ, "figures", "02_alpha_delta.png"), p2, width = 9, height = 4.5, dpi = 200, bg = "white")
}

say("已写出 results/02_*.csv, figures/02_alpha.png")
close(con)
cat("\n02 完成\n")
