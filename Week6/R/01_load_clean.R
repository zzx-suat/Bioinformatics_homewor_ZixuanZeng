# =============================================================================
# Week 6 — 01 载入与清洗
# 16S 微生物组（EMP tests: 16S_level-7.csv + 16S_mapping.csv）
#
# 这一步只做一件事：把数据变成可信的分析对象，并把每一处不一致都记录下来。
# 任何静默的修改都不允许——删了什么、改了什么，全部写进 results/01_cleaning_log.txt
# =============================================================================

suppressPackageStartupMessages({ library(dplyr) })

TESTS <- "D:/BioLession in total/Bioinformatics/week1/EasyMultiProfiler-Web/tests"
PROJ  <- "D:/BioLession in total/Bioinformatics/week6/week6_project"
dir.create(file.path(PROJ, "results"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PROJ, "figures"), showWarnings = FALSE, recursive = TRUE)

LOG <- file.path(PROJ, "results", "01_cleaning_log.txt")
con <- file(LOG, open = "wt", encoding = "UTF-8")
say <- function(...) { msg <- paste0(...); cat(msg, "\n"); writeLines(msg, con) }

say("=========== Week 6 数据清洗日志 ===========")
say("时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("")

# --- 读入 --------------------------------------------------------------------
abun_raw <- read.csv(file.path(TESTS, "16S_level-7.csv"),
                     row.names = 1, check.names = FALSE, fileEncoding = "UTF-8-BOM")
meta_raw <- read.csv(file.path(TESTS, "16S_mapping.csv"),
                     check.names = FALSE, fileEncoding = "UTF-8-BOM",
                     stringsAsFactors = FALSE)

say("原始输入")
say("  丰度表 : ", nrow(abun_raw), " taxa x ", ncol(abun_raw), " 样本")
say("  元数据 : ", nrow(meta_raw), " 行, 字段 = ", paste(names(meta_raw), collapse = ", "))
say("")

# --- 问题 1：样本集合不一致 --------------------------------------------------
a_ids <- colnames(abun_raw); m_ids <- meta_raw$SampleID
only_a <- setdiff(a_ids, m_ids); only_m <- setdiff(m_ids, a_ids)

say("问题 1 — 样本集合")
say("  两边都有        : ", length(intersect(a_ids, m_ids)))
say("  仅在丰度表 (孤儿): ", length(only_a),
    if (length(only_a)) paste0("  -> ", paste(only_a, collapse = ", ")) else "")
say("  仅在元数据      : ", length(only_m))
say("  顺序是否一致    : ", identical(a_ids, m_ids))
say("  处置：丢弃无元数据的孤儿样本（无分组信息则无法进入任何比较），")
say("        并把丰度表按元数据顺序重排。两项均记录于此，不静默进行。")
say("")

# --- 问题 2：命名不一致 ------------------------------------------------------
suffix <- sub("^.*_", "", m_ids)
say("问题 2 — 样本编号后缀不统一")
say("  后缀分布: ", paste(sprintf("_%s=%d", names(table(suffix)), table(suffix)),
                          collapse = ", "))
say("  '_2' 与 '_02' 指同一时间点，但字符串不同；解析受试者/时间点时按数值处理。")
say("")

# --- 对齐 --------------------------------------------------------------------
keep <- intersect(m_ids, a_ids)
meta <- meta_raw[match(keep, meta_raw$SampleID), , drop = FALSE]
abun <- as.matrix(abun_raw[, keep, drop = FALSE])
rownames(meta) <- meta$SampleID
stopifnot(identical(colnames(abun), rownames(meta)))
say("对齐后: ", nrow(abun), " taxa x ", ncol(abun), " 样本；列名与元数据行名严格一致")
say("")

# --- 解析受试者与时间点 ------------------------------------------------------
meta$Subject   <- sub("_[^_]+$", "", meta$SampleID)
meta$Visit     <- as.integer(sub("^.*_", "", meta$SampleID))
meta$Disease   <- sub("_.*$", "", meta$Group)                 # IBS / UC
meta$Timepoint <- factor(sub("^.*_", "", meta$Group), levels = c("before", "after"))
meta$Response  <- factor(sub("^.*_", "", meta$Group_sub), levels = c("poor", "great"))

say("解析出的设计变量")
say("  Disease   : ", paste(sprintf("%s=%d", names(table(meta$Disease)),
                                    table(meta$Disease)), collapse = ", "))
say("  Timepoint : ", paste(sprintf("%s=%d", names(table(meta$Timepoint)),
                                    table(meta$Timepoint)), collapse = ", "))
say("  Response  : ", paste(sprintf("%s=%d", names(table(meta$Response)),
                                    table(meta$Response)), collapse = ", "))
say("  受试者数  : ", length(unique(meta$Subject)))
say("")

# --- 问题 3：配对结构（这是整份作业的关键） ----------------------------------
tab <- table(meta$Subject)
paired_subj <- names(tab)[tab == 2]
single_subj <- names(tab)[tab == 1]

say("问题 3 — 配对结构（决定统计方法）")
say("  有 before+after 两次采样的受试者: ", length(paired_subj))
say("  只有单次采样的受试者            : ", length(single_subj),
    if (length(single_subj)) paste0("  -> ", paste(head(single_subj, 6), collapse = ", ")) else "")

chk <- meta %>% filter(Subject %in% paired_subj) %>%
  group_by(Subject) %>%
  summarise(tp = paste(sort(as.character(Timepoint)), collapse = "+"),
            resp = n_distinct(Response), dis = n_distinct(Disease), .groups = "drop")
say("  配对内 timepoint 组合: ",
    paste(sprintf("%s=%d", names(table(chk$tp)), table(chk$tp)), collapse = ", "))
say("  配对内 Response 前后不一致的受试者数: ", sum(chk$resp > 1))
say("  配对内 Disease  前后不一致的受试者数: ", sum(chk$dis  > 1))
say("")
say("  => 同一受试者的两个样本不是独立观测。")
say("     任何把 before/after 当作两个独立组的检验都是伪重复。")
say("     后续采用：受试者内对比 + 以 Subject 为 strata 的受限置换。")
say("")

# --- 测序深度与稀疏度 --------------------------------------------------------
depth <- colSums(abun)
say("测序深度（每样本总计数）")
say(sprintf("  min %.0f | Q1 %.0f | 中位 %.0f | Q3 %.0f | max %.0f",
            min(depth), quantile(depth, .25), median(depth),
            quantile(depth, .75), max(depth)))
say("  最深/最浅 = ", sprintf("%.1f 倍", max(depth) / min(depth)))
low <- names(depth)[depth < 1000]
say("  深度 < 1000 的样本: ", length(low),
    if (length(low)) paste0(" -> ", paste(head(low, 8), collapse = ", ")) else "")
say("")
say("稀疏度")
say(sprintf("  零值占比 %.1f%%", 100 * mean(abun == 0)))
prev <- rowSums(abun > 0)
say(sprintf("  仅出现在 1 个样本中的 taxa: %d / %d", sum(prev == 1), length(prev)))
say(sprintf("  出现在 >= 10%% 样本中的 taxa: %d", sum(prev >= 0.1 * ncol(abun))))
say("")

# --- 存档 --------------------------------------------------------------------
saveRDS(list(abun = abun, meta = meta), file.path(PROJ, "results", "01_clean_data.rds"))
write.csv(data.frame(SampleID = names(depth), depth = as.integer(depth),
                     row.names = NULL),
          file.path(PROJ, "results", "01_read_depth.csv"), row.names = FALSE)
write.csv(meta, file.path(PROJ, "results", "01_sample_metadata_parsed.csv"), row.names = FALSE)

say("已写出:")
say("  results/01_clean_data.rds")
say("  results/01_read_depth.csv")
say("  results/01_sample_metadata_parsed.csv")
say("  results/01_cleaning_log.txt  (本文件)")
close(con)
cat("\n01 完成\n")
