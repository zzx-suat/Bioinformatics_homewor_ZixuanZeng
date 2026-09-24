# =============================================================================
# Week 6 — 一键运行全部分析
# 依次执行 01 清洗 -> 02 alpha -> 03 beta/PERMANOVA -> 04 差异丰度 -> 05 平台对照
# 用法: Rscript R/run_all.R   （工作目录为 week6_project）
# =============================================================================
t0 <- Sys.time()
cat("=========== Week 6 全流程 ===========\n")
cat("开始时间 :", format(t0, "%Y-%m-%d %H:%M:%S"), "\n")
cat("R 版本   :", R.version.string, "\n")
cat("工作目录 :", getwd(), "\n\n")

steps <- c("R/01_load_clean.R", "R/02_alpha_diversity.R",
           "R/03_beta_diversity.R", "R/04_differential.R")
# 05 需要 EMP 平台导出的结果（results/EMPresult/），存在时才运行
if (dir.exists("results/EMPresult")) steps <- c(steps, "R/05_platform_vs_R.R")

for (s in steps) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat(">>> ", s, "\n", sep = "")
  cat(strrep("=", 70), "\n", sep = "")
  t1 <- Sys.time()
  source(s, echo = FALSE, local = new.env())
  cat(sprintf("[%s 用时 %.1f 秒]\n", basename(s),
              as.numeric(difftime(Sys.time(), t1, units = "secs"))))
}

cat("\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("全流程完成，总用时 %.1f 秒\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
writeLines(capture.output(sessionInfo()), "results/session_info.txt")
cat("sessionInfo 已写入 results/session_info.txt\n")
