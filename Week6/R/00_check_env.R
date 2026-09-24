# Week 6 环境检查：微生物组分析需要哪些包，装了哪些
pkgs <- c(
  "vegan",        # PERMANOVA (adonis2)、alpha/beta 多样性、受限置换
  "permute",      # 受限置换设计 how()
  "lme4",         # 线性混合模型（受试者随机效应）
  "lmerTest",     # 混合模型的 p 值
  "ggplot2", "ggrepel", "dplyr", "tidyr", "patchwork",
  "compositions", # CLR
  "phyloseq",     # 微生物组数据结构
  "ANCOMBC",      # 差异丰度（支持重复测量）
  "Maaslin2"      # 差异丰度（纵向设计）
)
for (p in pkgs) {
  cat(sprintf("  %-14s %s\n", p, requireNamespace(p, quietly = TRUE)))
}
cat("\nR: ", R.version.string, "\n", sep = "")
cat("Bioconductor: ",
    tryCatch(as.character(BiocManager::version()), error = function(e) "BiocManager 未装"),
    "\n", sep = "")
