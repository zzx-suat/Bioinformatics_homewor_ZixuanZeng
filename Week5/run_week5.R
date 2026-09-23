# ================================================================
# Week 5 一键运行引导脚本（由 一键运行.bat 调用，不要单独双击）
# 作用：检查并安装 DESeq2 / apeglm / ggrepel，然后运行分析脚本
# ================================================================
cat("\n========== Week 5 作业自动运行 ==========\n")
cat("R 版本 :", R.version.string, "\n")
cat("工作目录:", getwd(), "\n")
cat("时间   :", format(Sys.time()), "\n\n")

pick_mirror <- function(cands) {
  for (u in cands) {
    ok <- tryCatch({ con <- url(paste0(u, "/src/contrib/PACKAGES.gz"), "rb")
                     on.exit(close(con)); length(readBin(con, "raw", 50)) > 0 },
                   error = function(e) FALSE)
    if (ok) { cat("CRAN 镜像:", u, "\n"); return(u) }
  }
  cat("所有 CRAN 镜像都不通，仍尝试默认源\n"); cands[[1]]
}
cran <- pick_mirror(c("https://mirrors.tuna.tsinghua.edu.cn/CRAN",
                      "https://mirrors.ustc.edu.cn/CRAN",
                      "https://cloud.r-project.org"))
options(repos = c(CRAN = cran),
        BioC_mirror = if (grepl("tuna", cran)) "https://mirrors.tuna.tsinghua.edu.cn/bioconductor"
                      else "https://bioconductor.org")

inst <- rownames(installed.packages())
cran_need <- setdiff(c("ggplot2", "ggrepel", "BiocManager"), inst)
if (length(cran_need)) {
  cat("安装 CRAN 包:", paste(cran_need, collapse = ", "), "\n")
  install.packages(cran_need)
}
bioc_need <- setdiff(c("DESeq2", "apeglm"), rownames(installed.packages()))
if (length(bioc_need)) {
  cat("安装 Bioconductor 包:", paste(bioc_need, collapse = ", "), "（首次可能要几分钟）\n")
  BiocManager::install(bioc_need, ask = FALSE, update = FALSE)
} else {
  cat("DESeq2 / apeglm 已安装\n")
}

miss <- c("DESeq2", "apeglm", "ggplot2", "ggrepel")[
  !sapply(c("DESeq2", "apeglm", "ggplot2", "ggrepel"), requireNamespace, quietly = TRUE)]
if (length(miss)) {
  cat("\n!!! 以下包仍然缺失，无法继续:", paste(miss, collapse = ", "), "\n")
  cat("请把本文件 run_log.txt 全文发给助手。\n")
  quit(status = 1)
}

cat("\n---------- 开始分析 ----------\n")
source("week5_deseq2_analysis.R", echo = FALSE)
cat("\n---------- 分析结束 ----------\n")
cat("完成时间:", format(Sys.time()), "\n")
