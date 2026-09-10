## 00_setup.R -- one-time package installation + project folders
## Week 3 homework | Zixuan Zeng (SUAT24000114)
## Run this ONCE. Safe to re-run: it skips anything already installed.

pkgs_cran <- c("ggplot2", "reshape2", "pheatmap")
pkgs_bioc <- c("GEOquery", "Biobase", "SummarizedExperiment")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

need_cran <- pkgs_cran[!vapply(pkgs_cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(need_cran))
  install.packages(need_cran, repos = "https://cloud.r-project.org")

need_bioc <- pkgs_bioc[!vapply(pkgs_bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(need_bioc))
  BiocManager::install(need_bioc, ask = FALSE, update = FALSE)

## project folders
for (d in c("data_raw", "data_clean", "metadata", "results", "R"))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

cat("\n--- package status ---\n")
for (p in c(pkgs_cran, pkgs_bioc))
  cat(sprintf("%-24s %s\n", p,
      if (requireNamespace(p, quietly = TRUE))
        paste("OK", as.character(packageVersion(p))) else "MISSING"))
cat("\n00_setup.R done\n")
