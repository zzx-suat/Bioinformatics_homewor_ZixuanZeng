## 01_fetch.R -- fetch GSE111889 (Task A primary dataset) from GEO
## Week 3 homework | Zixuan Zeng (SUAT24000114)
## Principle from the lecture: LIST supplementary files first, decide, THEN download.

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})
source(file.path("R", "utils_download.R"))   # retry helpers: this network's DNS is flaky
options(timeout = 3600)          # 10.6 MB over a slow link needs more than the 60 s default
GSE <- "GSE111889"
for (d in c("data_raw", "data_clean", "metadata", "results")) dir.create(d, showWarnings = FALSE)

cat("\n================ STEP 1: list supplementary files (no download) ================\n")
supp <- robust_supp_list(GSE)
print(supp[, intersect(c("fname", "size", "url"), colnames(supp))])
write.csv(supp, file.path("metadata", paste0(GSE, "_supp_manifest.csv")), row.names = FALSE)

cat("\n================ STEP 2: download the series matrix (metadata layer) ==========\n")
eset <- load_series_matrix(GSE, "data_raw")
cat("class:", class(eset), "\n")
cat("dim (features x samples):", paste(dim(eset), collapse = " x "), "\n")
cat("platform:", unique(as.character(eset$platform_id)), "\n")

pheno <- as.data.frame(pData(eset))
cat("n samples in series matrix:", nrow(pheno), "\n")
cat("GSM range:", pheno$geo_accession[1], "...", pheno$geo_accession[nrow(pheno)], "\n")
cat("\n-- phenoData columns --\n"); print(colnames(pheno))
cat("\n-- first 3 sample titles --\n"); print(head(pheno$title, 3))
cat("\n-- characteristics columns, first sample --\n")
ch <- grep("^characteristics_ch1", colnames(pheno), value = TRUE)
for (c1 in ch) cat(sprintf("  %-24s : %s\n", c1, as.character(pheno[1, c1])))
cat("\n-- relation fields (BioProject / SRA links) --\n")
rel <- grep("^relation", colnames(pheno), value = TRUE)
for (r in rel) cat(sprintf("  %-24s : %s\n", r, as.character(pheno[1, r])))
cat("\n-- supplementary_file_1 of first 2 samples (raw-read layer) --\n")
sfc <- grep("^supplementary_file", colnames(pheno), value = TRUE)
if (length(sfc)) print(head(pheno[, sfc[1], drop = FALSE], 2))

cat("\n-- does the series matrix itself carry an expression table? --\n")
em <- exprs(eset)
cat("exprs() dim:", paste(dim(em), collapse = " x "),
    "  (0 rows means: NO data table in the series matrix -> counts must come from the",
    "supplementary file. This is itself an answer for Table 0 Step 4.)\n")

write.csv(pheno, file.path("metadata", paste0(GSE, "_pheno_raw.csv")), row.names = FALSE)

cat("\n================ STEP 3: download the supplementary count matrix ===============\n")
i <- grep("host_tx_counts", supp$fname)
if (!length(i)) i <- grep("counts", supp$fname)[1]
stopifnot(length(i) >= 1)
i <- i[1]
target <- file.path("data_raw", supp$fname[i])
cat("chosen file :", supp$fname[i], "\n")
cat("reported size:", as.character(supp$size[i]), "bytes\n")
cat("url          :", as.character(supp$url[i]), "\n")
robust_download(as.character(supp$url[i]), target)
cat("downloaded bytes:", file.size(target), "\n")
cat("md5             :", unname(tools::md5sum(target)), "\n")
cat("download date   :", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

saveRDS(list(eset = eset, pheno = pheno, supp = supp, counts_file = target,
             md5 = unname(tools::md5sum(target)),
             fetched_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
        file.path("data_clean", "01_fetch_objects.rds"))

cat("\n01_fetch.R done. Provenance saved to data_clean/01_fetch_objects.rds\n")
