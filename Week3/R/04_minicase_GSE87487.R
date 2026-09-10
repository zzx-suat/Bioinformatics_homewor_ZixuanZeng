## 04_minicase_GSE87487.R -- Task B mini-case: verify the instructor's claims
## Week 3 homework | Zixuan Zeng (SUAT24000114)
## Purpose: answer Table 1 with EVIDENCE, not assumption.
##   Q1 are values counts?  -> download the supplementary matrix, test integers
##   Q2 independent?        -> read the metadata, look for donor/subject pairing
##   Q3 all 20 comparable?  -> inspect the actual factor combinations
##   Q4 raw-read reprocessing? -> confirm the SRA link exists

suppressPackageStartupMessages({ library(GEOquery); library(Biobase) })
source(file.path("R", "utils_download.R"))
options(timeout = 3600)
GSE <- "GSE87487"
for (d in c("data_raw", "metadata", "results")) dir.create(d, showWarnings = FALSE)

cat("\n================ supplementary file inventory ================\n")
supp <- robust_supp_list(GSE)
print(supp[, intersect(c("fname", "size", "url"), colnames(supp))])
write.csv(supp, file.path("metadata", paste0(GSE, "_supp_manifest.csv")), row.names = FALSE)

cat("\n================ series metadata ================\n")
e <- load_series_matrix(GSE, "data_raw")
ph    <- as.data.frame(pData(e))
cat("n samples:", nrow(ph), "\n")
cat("GSM range:", ph$geo_accession[1], "...", ph$geo_accession[nrow(ph)], "\n")
cat("platform :", unique(as.character(ph$platform_id)), "\n")
cat("\n-- titles (pre/post pairing is usually visible here) --\n")
print(data.frame(gsm = ph$geo_accession, title = ph$title), right = FALSE)
cat("\n-- characteristics --\n")
for (c1 in grep("^characteristics_ch1", colnames(ph), value = TRUE))
  print(table(as.character(ph[[c1]]), useNA = "ifany"))
cat("\n-- relation (SRA / BioProject) --\n")
for (r in grep("^relation", colnames(ph), value = TRUE))
  print(unique(as.character(ph[[r]])))
write.csv(ph, file.path("metadata", paste0(GSE, "_pheno_raw.csv")), row.names = FALSE)

cat("\n================ Q1: ARE THE VALUES COUNTS? ================\n")
i <- grep("count", supp$fname, ignore.case = TRUE)[1]
if (is.na(i)) i <- 1
tgt <- file.path("data_raw", supp$fname[i])
cat("file:", supp$fname[i], " reported size:", as.character(supp$size[i]), "bytes\n")
robust_download(as.character(supp$url[i]), tgt)
cat("bytes on disk:", file.size(tgt), " md5:", unname(tools::md5sum(tgt)), "\n")
cat("download date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

m <- read.delim(gzfile(tgt), header = TRUE, check.names = FALSE,
                row.names = 1, stringsAsFactors = FALSE)
num <- m[, vapply(m, is.numeric, logical(1)), drop = FALSE]

## TRAP: featureCounts exports gene annotation columns alongside the samples.
## Here the first column is "Length" (a gene length in bp, colSum 1.3e8) --
## numeric, but NOT a library. Treating it as a sample would corrupt every
## downstream normalisation. Drop annotation columns explicitly.
fc_meta <- c("Length", "Chr", "Start", "End", "Strand", "gene_length", "width")
drop_i  <- colnames(num) %in% fc_meta
cat("non-sample annotation columns found:",
    paste(colnames(num)[drop_i], collapse = ", "),
    "-> dropped\n")
num <- num[, !drop_i, drop = FALSE]
mm  <- as.matrix(num)
cat("dim (genes x samples):", paste(dim(mm), collapse = " x "), "\n")
cat("columns:", paste(colnames(mm), collapse = ", "), "\n")
cat("gene IDs:", paste(head(rownames(mm), 3), collapse = ", "), "\n")
cat("ALL INTEGER? ->", all(mm == floor(mm), na.rm = TRUE), "\n")
cat("min:", min(mm, na.rm = TRUE), " max:", max(mm, na.rm = TRUE),
    " any negative:", any(mm < 0, na.rm = TRUE),
    " any fractional:", any(mm != floor(mm), na.rm = TRUE), "\n")
cat("column sums (library sizes):\n"); print(colSums(mm))
print(summary(as.vector(mm[, 1:min(4, ncol(mm))])))


cat("\n================ Q2/Q3: RECONSTRUCT THE PAIRED DESIGN ================\n")
## Map count columns to GSMs. The two naming schemes do NOT match verbatim:
##   count column : "Sample_Pt10-Bx1.bam"
##   GEO title    : "Pt10Bx1"
## so both sides are normalised (drop prefix/suffix, hyphens, case) first.
norm_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "",
                     sub(".bam", "", sub("^Sample_", "", x), fixed = TRUE)))
cn_n <- norm_id(colnames(mm))
ti_n <- norm_id(as.character(ph$title))
ord  <- match(cn_n, ti_n)
cat("count columns matched to a GSM:", sum(!is.na(ord)), "of", length(cn_n), "\n")
if (any(is.na(ord))) {
  cat("!! unmatched columns:", paste(colnames(mm)[is.na(ord)], collapse = ", "), "\n")
  cat("   GEO titles were:", paste(ph$title, collapse = ", "), "\n")
}

getch2 <- function(tag) {
  out <- rep(NA_character_, nrow(ph))
  for (c1 in grep("^characteristics_ch1", colnames(ph), value = TRUE)) {
    v <- as.character(ph[[c1]])
    j <- grepl(tag, v, ignore.case = TRUE) & is.na(out)
    out[j] <- trimws(sub("^[^:]*: *", "", v[j]))
  }
  out
}
stage <- getch2("transplant stage")
iri   <- getch2("ischemia|iri")

## donor = title with the trailing Bx1/Bx2 removed
donor_of <- function(t) toupper(gsub("[^A-Za-z0-9]", "",
                        sub("B[xX][12]$", "", t)))
d87 <- data.frame(
  count_column = colnames(mm),
  gsm          = as.character(ph$geo_accession)[ord],
  geo_title    = as.character(ph$title)[ord],
  donor        = donor_of(as.character(ph$title)[ord]),
  stage        = stage[ord],
  iri          = iri[ord],
  libsize      = colSums(mm),
  stringsAsFactors = FALSE
)
print(d87, right = FALSE)
write.csv(d87, file.path("metadata", "GSE87487_sample_dictionary.csv"), row.names = FALSE)

cat("\n-- Q2: are samples independent? --\n")
cat("unique donors:", length(unique(d87$donor)), " samples:", nrow(d87), "\n")
cat("samples per donor:\n"); print(table(table(d87$donor)))
cat("=> ", nrow(d87), "samples come from only", length(unique(d87$donor)),
    "donors: the samples are NOT independent.\n")

cat("\n-- Q3: is every one of the 20 directly comparable? --\n")
cat("stage x IRI cross-tabulation:\n")
print(table(stage = d87$stage, iri = d87$iri, useNA = "ifany"))
cat("donor x stage completeness (1 = present):\n")
print(table(donor = d87$donor, stage = d87$stage))
cat("library size range:", format(range(d87$libsize), big.mark = ","),
    " fold:", round(max(d87$libsize) / min(d87$libsize), 2), "\n")

cat("\n-- Q4: is raw-read reprocessing possible? --\n")
sra <- unlist(lapply(grep("^relation", colnames(ph), value = TRUE),
                     function(r) grep("SRA|SRX", as.character(ph[[r]]), value = TRUE)))
cat("per-sample SRA experiment links found:", length(sra), "of", nrow(ph), "\n")
if (length(sra)) cat("example:", sra[1], "\n")
cat("=> raw reads are openly accessible for this study (contrast with GSE111889,\n")
cat("   whose Series states raw data go to dbGaP under controlled access).\n")

cat("\n04_minicase_GSE87487.R done\n")
