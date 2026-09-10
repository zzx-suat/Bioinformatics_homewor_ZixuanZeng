## 03_qc.R -- library size, distribution, PCA by group and by batch/covariate
## Week 3 homework | Zixuan Zeng (SUAT24000114)
## Normalization here is for VISUALIZATION ONLY. Statistics stay on raw counts.

suppressPackageStartupMessages({ library(ggplot2) })
obj    <- readRDS(file.path("data_clean", "counts_filtered.rds"))
counts <- obj$counts
dict   <- obj$dict
dir.create("results", showWarnings = FALSE)

cat("\n================ LIBRARY SIZE ================\n")
lib <- data.frame(sample = colnames(counts), libsize = colSums(counts),
                  group = dict$diagnosis, stringsAsFactors = FALSE)
print(summary(lib$libsize))
cat("min / median / max library size:",
    format(c(min(lib$libsize), median(lib$libsize), max(lib$libsize)), big.mark = ","), "\n")
cat("fold range (max/min):", round(max(lib$libsize) / min(lib$libsize), 1), "\n")
low <- lib[lib$libsize < 0.5 * median(lib$libsize), ]
cat("samples below 50% of median depth (flag, do NOT delete):", nrow(low), "\n")
if (nrow(low)) print(low[order(low$libsize), ])
write.csv(lib, "results/library_sizes.csv", row.names = FALSE)

p1 <- ggplot(lib, aes(x = reorder(sample, libsize), y = libsize / 1e6, fill = group)) +
  geom_col() + coord_flip() +
  labs(title = "GSE111889 library size per sample", x = NULL, y = "million counts") +
  theme_bw(base_size = 7) + theme(legend.position = "top")
ggsave("results/qc_library_size.png", p1, width = 6, height = 10, dpi = 150)

cat("\n================ ZERO-DEPTH SAMPLES: DOCUMENT BEFORE EXCLUDING ================\n")
## A sample with 0 total counts cannot be CPM-normalised (0/0 -> NaN) and would
## crash PCA. Per the lecture: never silently drop a sample. Every exclusion is
## logged, such samples stay in the library-size figure above, and they are
## excluded ONLY from the normalisation-dependent plots below.
zero_s <- colnames(counts)[colSums(counts) == 0]
cat("samples with ZERO total counts after filtering:", length(zero_s), "\n")
excl <- data.frame(sample    = zero_s,
                   gsm       = dict$gsm[match(zero_s, dict$sample_col)],
                   diagnosis = dict$diagnosis[match(zero_s, dict$sample_col)],
                   reason    = rep("zero total counts; cannot be CPM-normalised",
                                   length(zero_s)),
                   stringsAsFactors = FALSE)
if (nrow(excl)) print(excl)
write.csv(excl, "results/excluded_samples.csv", row.names = FALSE)
cat("-> results/excluded_samples.csv written (declare these in the report)\n")
cat("very low depth (<50% of median) but RETAINED:",
    paste(setdiff(low$sample, zero_s), collapse = ", "), "\n")

keep_n <- colSums(counts) > 0
counts <- counts[, keep_n, drop = FALSE]
dict   <- dict[keep_n, , drop = FALSE]
stopifnot(identical(colnames(counts), dict$sample_col))
cat("samples carried into normalisation/PCA:", ncol(counts), "of",
    length(keep_n), "\n")

cat("\n================ DISTRIBUTION (log2 CPM density) ================\n")
cpm    <- t(t(counts) / colSums(counts)) * 1e6
logcpm <- log2(cpm + 1)
df <- data.frame(value = as.vector(logcpm),
                 sample = rep(colnames(logcpm), each = nrow(logcpm)),
                 group  = rep(dict$diagnosis, each = nrow(logcpm)))
p2 <- ggplot(df, aes(value, group = sample, colour = group)) +
  geom_density(alpha = .3, linewidth = .2) +
  labs(title = "GSE111889 log2(CPM+1) density, filtered genes", x = "log2(CPM+1)") +
  theme_bw(base_size = 9)
ggsave("results/qc_density.png", p2, width = 7, height = 4.5, dpi = 150)
cat("density plot written. Any single curve far from the bulk = candidate outlier.\n")

cat("\n================ PCA (on log2 CPM, visualization only) ================\n")
v   <- apply(logcpm, 1, var)
top <- head(order(v, decreasing = TRUE), 2000)
pca <- prcomp(t(logcpm[top, ]), center = TRUE, scale. = FALSE)
pv  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
cat("variance explained PC1-PC4:", paste0(pv[1:4], "%", collapse = "  "), "\n")

pdf_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                     sample = rownames(pca$x),
                     group = dict$diagnosis, site = dict$biopsy_site,
                     subject = dict$subject_id, libsize = colSums(counts))

mk <- function(colvar, ttl) {
  ggplot(pdf_df, aes(PC1, PC2, colour = .data[[colvar]])) +
    geom_point(size = 1.8, alpha = .85) +
    labs(title = ttl, x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)"),
         colour = colvar) +
    theme_bw(base_size = 9)
}
ggsave("results/pca_by_group.png",   mk("group", "PCA coloured by diagnosis (the comparison of interest)"),
       width = 6.5, height = 4.5, dpi = 150)
ggsave("results/pca_by_site.png",    mk("site",  "PCA coloured by biopsy site (technical/anatomical covariate)"),
       width = 6.5, height = 4.5, dpi = 150)
p5 <- ggplot(pdf_df, aes(PC1, PC2, colour = log10(libsize))) + geom_point(size = 1.8) +
  scale_colour_viridis_c() +
  labs(title = "PCA coloured by sequencing depth", x = paste0("PC1 (", pv[1], "%)"),
       y = paste0("PC2 (", pv[2], "%)")) + theme_bw(base_size = 9)
ggsave("results/pca_by_depth.png", p5, width = 6.5, height = 4.5, dpi = 150)

cat("\n-- association of PC1/PC2 with candidate drivers (ANOVA p-values) --\n")
for (vn in c("group", "site")) {
  x <- pdf_df[[vn]]
  if (length(unique(na.omit(x))) > 1) {
    for (pc in c("PC1", "PC2")) {
      p <- tryCatch(summary(aov(pdf_df[[pc]] ~ factor(x)))[[1]][["Pr(>F)"]][1],
                    error = function(e) NA)
      cat(sprintf("  %s ~ %-6s : p = %s\n", pc, vn, format.pval(p, digits = 3)))
    }
  }
}
cat("  PC1 ~ log10(libsize) : r =",
    round(cor(pdf_df$PC1, log10(pdf_df$libsize)), 3), "\n")

write.csv(pdf_df, "results/pca_coordinates.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
cat("\n03_qc.R done. Figures in results/ ; sessionInfo.txt written\n")
