## 02_inspect_clean.R -- integrity checks, metadata alignment, low-count filter
## Week 3 homework | Zixuan Zeng (SUAT24000114)

suppressPackageStartupMessages({ library(Biobase) })
obj   <- readRDS(file.path("data_clean", "01_fetch_objects.rds"))
pheno <- obj$pheno
cfile <- obj$counts_file

cat("\n================ READ THE SUPPLEMENTARY COUNT MATRIX ================\n")
cat("file:", cfile, " md5:", obj$md5, "\n")
raw <- read.delim(gzfile(cfile), header = TRUE, check.names = FALSE,
                  row.names = 1, stringsAsFactors = FALSE)
expr <- as.matrix(raw)
cat("dim (genes x samples):", paste(dim(expr), collapse = " x "), "\n")
cat("gene ID examples :", paste(head(rownames(expr), 4), collapse = ", "), "\n")
cat("column examples  :", paste(head(colnames(expr), 4), collapse = ", "), "\n")
cat("storage mode     :", storage.mode(expr), "\n")

cat("\n================ 5 INTEGRITY CHECKS (stopifnot = failures are LOUD) =========\n")
## 1. duplicated sample columns
cat("1. duplicated sample columns :", anyDuplicated(colnames(expr)), "(0 = none)\n")
stopifnot(!anyDuplicated(colnames(expr)))
## 2. duplicated / missing gene IDs
cat("2. duplicated gene IDs       :", anyDuplicated(rownames(expr)), "(0 = none)\n")
cat("   NA values in matrix       :", sum(is.na(expr)), "\n")
stopifnot(!anyDuplicated(rownames(expr)))
## 3. ARE THE VALUES COUNTS? (integer test -- Table 1 Q1 logic applied to Task A too)
is_int <- all(expr == floor(expr), na.rm = TRUE)
cat("3. all values integers?      :", is_int, "\n")
cat("   min / max                 :", min(expr, na.rm = TRUE), "/", max(expr, na.rm = TRUE), "\n")
cat("   any negative?             :", any(expr < 0, na.rm = TRUE), "\n")
print(summary(as.vector(expr[, 1:min(5, ncol(expr))])))
if (!is_int) cat("   !! NOT integers -> these are normalized values, NOT raw counts.",
                 "Report this honestly; do NOT round them.\n")
## 4. dimensions vs metadata
cat("4. ncol(expr) =", ncol(expr), " vs nrow(pheno) =", nrow(pheno), "\n")

cat("\n================ MAP COUNT COLUMNS TO GSM METADATA ================\n")
cn  <- colnames(expr)
key <- list(title = as.character(pheno$title),
            geo   = as.character(pheno$geo_accession),
            desc  = if ("description" %in% colnames(pheno)) as.character(pheno$description) else NA)
hits <- vapply(key, function(k) sum(cn %in% k), integer(1))
print(hits)
use <- names(which.max(hits))
cat("best matching metadata column:", use, "-> matched", max(hits), "of", ncol(expr), "columns\n")
if (max(hits) == 0) {
  cat("!! No exact match. Showing both sides so the mapping can be built by hand:\n")
  print(head(cn, 10)); print(head(key$title, 10))
  stop("Stopping: sample-ID mapping must be resolved before any analysis (lecture pitfall #6).")
}
mkey <- key[[use]]
ord  <- match(cn, mkey)
meta <- pheno[ord, , drop = FALSE]
keep_s <- !is.na(ord)
cat("samples with metadata:", sum(keep_s), " | without:", sum(!keep_s), "\n")
expr <- expr[, keep_s, drop = FALSE]
meta <- meta[keep_s, , drop = FALSE]
rownames(meta) <- colnames(expr)

## 5. order must be identical -- never assume
cat("5. identical(colnames(expr), rownames(meta)) :",
    identical(colnames(expr), rownames(meta)), "\n")
stopifnot(identical(colnames(expr), rownames(meta)))

cat("\n================ BUILD THE SAMPLE DICTIONARY (design reconstruction) ========\n")
ch <- grep("^characteristics_ch1", colnames(meta), value = TRUE)
getch <- function(tag) {
  out <- rep(NA_character_, nrow(meta))
  for (c1 in ch) {
    v <- as.character(meta[[c1]])
    j <- grepl(tag, v, ignore.case = TRUE) & is.na(out)
    out[j] <- trimws(sub("^[^:]*:\\s*", "", v[j]))
  }
  out
}
dict <- data.frame(
  sample_col    = colnames(expr),
  gsm           = as.character(meta$geo_accession),
  title         = as.character(meta$title),
  diagnosis     = getch("diagnosis|disease"),
  subject_id    = getch("subject|participant|patient|donor"),
  biopsy_site   = getch("biopsy|location|site|tissue"),
  visit_week    = getch("week|visit|time"),
  sex           = getch("sex|gender"),
  age           = getch("age"),
  stringsAsFactors = FALSE
)
cat("\n-- group sizes (diagnosis) --\n");   print(table(dict$diagnosis, useNA = "ifany"))
cat("\n-- biopsy site --\n");               print(table(dict$biopsy_site, useNA = "ifany"))
cat("\n-- unique subjects --\n");           cat(length(unique(na.omit(dict$subject_id))), "\n")
cat("\n-- samples per subject (repeated-measures / paired structure) --\n")
print(table(table(dict$subject_id[!is.na(dict$subject_id)])))
cat("\n-- head of dictionary --\n");        print(head(dict, 5))
write.csv(dict, file.path("metadata", "sample_dictionary.csv"), row.names = FALSE)

cat("\n================ LOW-COUNT FILTER ================\n")
min_count <- 10
min_samp  <- max(3, min(table(dict$diagnosis, useNA = "no")))
cat("threshold: >=", min_count, "counts in >=", min_samp,
    "samples (min_samp = smallest group size, floor 3)\n")
keep   <- rowSums(expr >= min_count) >= min_samp
expr_f <- expr[keep, , drop = FALSE]
cat(sprintf("genes retained: %d of %d (%.1f%%); removed %d\n",
            nrow(expr_f), nrow(expr), 100 * nrow(expr_f) / nrow(expr),
            nrow(expr) - nrow(expr_f)))
cat("all-zero genes before filtering:", sum(rowSums(expr) == 0), "\n")

saveRDS(list(counts = expr_f, meta = meta, dict = dict,
             is_integer = is_int, min_count = min_count, min_samp = min_samp,
             n_genes_before = nrow(expr), n_genes_after = nrow(expr_f),
             md5 = obj$md5, fetched_at = obj$fetched_at),
        file.path("data_clean", "counts_filtered.rds"))
cat("\n02_inspect_clean.R done -> data_clean/counts_filtered.rds\n")
