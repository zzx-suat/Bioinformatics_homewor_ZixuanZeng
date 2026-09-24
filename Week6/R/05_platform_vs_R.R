# =============================================================================
# Week 6 — 05 EMP 平台结果与 R 结果对照
#
# 平台产出（results/EMPresult/）：
#   diff-result.csv      IBS_before vs IBS_after，分析页，Wilcoxon 秩和（不配对）
#   UC_diff-result.csv   UC_before  vs UC_after，取自 EMP 同步到 GitHub 的文件
#   20260925-004748_RunAll_全部132样本.zip  全样本 Run All（alpha / beta / 热图）
#
# 对照三件事：
#   1. 差异丰度的方向是否一致（平台不配对 vs R 配对）
#   2. 两种检验的灵敏度（配对是否更灵敏）
#   3. 平台算的 Shannon，用配对检验能否复现 R 看到的 UC 多样性下降
#
# 平台 log2FC 的方向：已用原始数据核验，正值 = 参照组（*_before）更高。
# 本脚本不依赖 log2FC 的正负号，而用平台的 sign_group 列判断哪组更高。
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
set.seed(1)

PROJ <- "D:/BioLession in total/Bioinformatics/week6/week6_project"
EMP  <- file.path(PROJ, "results", "EMPresult")
stopifnot(dir.exists(EMP))

LOG <- file.path(PROJ, "results", "05_platform_vs_R_log.txt")
con <- file(LOG, open = "wt", encoding = "UTF-8")
say <- function(...) { m <- paste0(...); cat(m, "\n"); writeLines(m, con) }

say("=========== Week 6 平台结果与 R 结果对照 ===========")
say("时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("")

meta <- readRDS(file.path(PROJ, "results", "01_clean_data.rds"))$meta
r_dif <- read.csv(file.path(PROJ, "results", "04_paired_taxa.csv"), check.names = FALSE)

short <- function(x) {
  p <- strsplit(x, ";")[[1]]; p <- p[!p %in% c("__", "")]
  if (!length(p)) return(x); paste(tail(p, 2), collapse = ";")
}

# --- 1 & 2. 差异丰度对照 -----------------------------------------------------
plat_files <- c(IBS = "diff-result.csv", UC = "UC_diff-result.csv")
conc_rows <- list(); sens_rows <- list(); plot_rows <- list()

for (dis in names(plat_files)) {
  pf <- read.csv(file.path(EMP, plat_files[[dis]]), check.names = FALSE)
  stopifnot(all(c("feature", "pvalue", "fdr", "sign_group") %in% names(pf)))
  before_lab <- paste0(dis, "_before"); after_lab <- paste0(dis, "_after")
  stopifnot(all(pf$sign_group %in% c(before_lab, after_lab)))

  rr <- r_dif %>% filter(disease == dis)
  m  <- inner_join(pf, rr, by = c("feature" = "taxon"), suffix = c(".plat", ".R"))

  say("---- ", dis, " ----")
  say(sprintf("  平台特征 %d | R 特征 %d | 两边共有 %d", nrow(pf), nrow(rr), nrow(m)))

  # 方向：平台看 sign_group；R 看受试者内中位变化（after - before）的正负
  m <- m %>% mutate(
    plat_higher = ifelse(sign_group == after_lab, "after", "before"),
    R_higher    = ifelse(median_delta > 0, "after",
                  ifelse(median_delta < 0, "before", NA)))
  mm <- m %>% filter(!is.na(R_higher))
  agree_all <- mean(mm$plat_higher == mm$R_higher)

  nom <- mm %>% filter(pvalue < 0.05 | p < 0.05)
  agree_nom <- if (nrow(nom)) mean(nom$plat_higher == nom$R_higher) else NA

  say(sprintf("  方向一致率（R 中位变化非零的 %d 个共有特征）: %.1f%%",
              nrow(mm), 100 * agree_all))
  say(sprintf("  方向一致率（任一方法未校正 p<0.05 的 %d 个特征）: %.1f%%",
              nrow(nom), 100 * agree_nom))

  sens <- data.frame(
    disease = dis,
    method = c("平台：Wilcoxon 秩和（不配对）", "R：Wilcoxon 符号秩（配对）"),
    n_features = c(nrow(pf), nrow(rr)),
    min_p   = c(min(pf$pvalue, na.rm = TRUE), min(rr$p, na.rm = TRUE)),
    min_fdr = c(min(pf$fdr,    na.rm = TRUE), min(rr$p_adj, na.rm = TRUE)),
    n_p_lt_0.05   = c(sum(pf$pvalue < 0.05, na.rm = TRUE), sum(rr$p < 0.05, na.rm = TRUE)),
    n_fdr_lt_0.05 = c(sum(pf$fdr    < 0.05, na.rm = TRUE), sum(rr$p_adj < 0.05, na.rm = TRUE)))
  for (i in seq_len(nrow(sens))) {
    s <- sens[i, ]
    say(sprintf("  %-28s 最小p=%.4f 最小FDR=%.3f  p<0.05:%2d  FDR<0.05:%d",
                s$method, s$min_p, s$min_fdr, s$n_p_lt_0.05, s$n_fdr_lt_0.05))
  }

  top <- nom %>% arrange(pmin(pvalue, p)) %>% head(8)
  say("  主要特征对照（平台 sign_group vs R 中位变化）:")
  for (i in seq_len(nrow(top))) {
    t <- top[i, ]
    say(sprintf("    %-40s 平台:%-6s p=%.4f | R:%-6s Δ=%+.3f p=%.4f  %s",
                substr(short(t$feature), 1, 40), t$plat_higher, t$pvalue,
                t$R_higher, t$median_delta, t$p,
                ifelse(t$plat_higher == t$R_higher, "一致", "不一致")))
  }
  say("")

  conc_rows[[dis]] <- data.frame(disease = dis, n_common = nrow(m),
                                 n_direction_scored = nrow(mm),
                                 agree_all_pct = round(100 * agree_all, 1),
                                 n_nominal = nrow(nom),
                                 agree_nominal_pct = round(100 * agree_nom, 1))
  sens_rows[[dis]] <- sens
  plot_rows[[dis]] <- mm %>% transmute(disease = dis, feature,
                                       plat_signed = ifelse(plat_higher == "after", 1, -1) * -log10(pvalue),
                                       R_signed    = sign(median_delta) * -log10(p))
}

write.csv(bind_rows(conc_rows), file.path(PROJ, "results", "05_direction_concordance.csv"), row.names = FALSE)
write.csv(bind_rows(sens_rows), file.path(PROJ, "results", "05_sensitivity.csv"), row.names = FALSE)

# --- 3. 平台 Shannon 的配对复核 ----------------------------------------------
say("---- 平台计算的 Shannon（全样本 Run All），用配对检验复核 ----")
tmp <- tempfile(); dir.create(tmp)
utils::unzip(file.path(EMP, "20260925-004748_RunAll_全部132样本.zip"),
             files = "tables/02_alpha_indices.csv", exdir = tmp)
pa <- read.csv(file.path(tmp, "tables", "02_alpha_indices.csv"), check.names = FALSE)
names(pa)[1] <- "SampleID"
say("  平台 alpha 表样本数: ", nrow(pa), "（含 2 个无元数据样本，下方配对时自然排除）")

pa <- inner_join(pa, meta[, c("SampleID", "Subject", "Disease", "Timepoint")], by = "SampleID")
sh_rows <- list()
for (dis in c("IBS", "UC")) {
  w <- pa %>% filter(Disease == dis) %>% select(Subject, Timepoint, shannon) %>%
    tidyr::pivot_wider(names_from = Timepoint, values_from = shannon) %>%
    filter(!is.na(before), !is.na(after))
  tt <- wilcox.test(w$after, w$before, paired = TRUE, exact = FALSE)
  tu <- wilcox.test(w$after, w$before, paired = FALSE, exact = FALSE)
  say(sprintf("  %-3s n=%2d 对  中位 %.3f -> %.3f  配对 p=%.4f | 若当作独立两组 p=%.4f",
              dis, nrow(w), median(w$before), median(w$after), tt$p.value, tu$p.value))
  sh_rows[[dis]] <- data.frame(disease = dis, n_pairs = nrow(w),
                               median_before = median(w$before), median_after = median(w$after),
                               p_paired = tt$p.value, p_unpaired = tu$p.value)
}
write.csv(bind_rows(sh_rows), file.path(PROJ, "results", "05_platform_shannon_paired.csv"), row.names = FALSE)
say("  说明：平台的 alpha 在属水平聚合后、未抽平的计数上计算，R 的是物种水平抽平后计算，")
say("        两者数值不直接可比，这里只比较「配对与否」对结论的影响。")
say("")

# --- 图 ---------------------------------------------------------------------
pd <- bind_rows(plot_rows)
p <- ggplot(pd, aes(R_signed, plat_signed)) +
  geom_hline(yintercept = 0, colour = "#9A9AA6") +
  geom_vline(xintercept = 0, colour = "#9A9AA6") +
  geom_point(alpha = .7, size = 1.9, colour = "#3E2A63") +
  facet_wrap(~ disease) +
  labs(title = "平台差异分析与 R 配对分析的方向对照",
       subtitle = "每点一个 taxon；坐标 = 方向 × −log10(p)，右上/左下象限为两种方法方向一致",
       x = "R（配对 Wilcoxon）：方向 × −log10 p", y = "平台（不配对 Wilcoxon）：方向 × −log10 p") +
  theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(PROJ, "figures", "05_platform_vs_R.png"), p, width = 9.5, height = 5, dpi = 200, bg = "white")

say("已写出 results/05_*.csv, figures/05_platform_vs_R.png")
close(con)
cat("\n05 完成\n")
