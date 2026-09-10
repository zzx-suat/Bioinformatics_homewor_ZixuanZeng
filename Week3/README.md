# Week 3 — Data Fetching and Cleaning from Public Repositories

**Student:** Zixuan Zeng (曾梓轩) · SUAT24000114 · Bioinformatics: From Multi-Omics Data to Discovery (Dr. Liwei Xie, SUAT, Fall 2026)

## 1. The data-driven question (Week 2 → Week 3)

Week 2 Part B, Paper 1 was Priya *et al.*, *Nature Microbiology* 2022, "Identification of
shared and disease-specific host gene–microbiome associations across human diseases using
multi-omic integration" (doi:10.1038/s41564-022-01121-z). Its Data availability statement
names the host RNA-seq accession for the IBD cohort: **GSE111889**. That is the GSE carried
forward into Week 3.

**PECO question:** In colonic/ileal mucosal biopsies from patients enrolled in the HMP2/IBDMDB
cohort (**P**), is Crohn's disease or ulcerative colitis status (**E**), compared with non-IBD
controls (**C**), associated with differences in host mucosal gene expression measured as bulk
RNA-seq counts (**O**)?

## 2. Accessions

| Layer | Accession | Note |
|---|---|---|
| GEO Series (host RNA-seq) | GSE111889 | 251 samples, GSM3043377–GSM3043627 |
| Platform | GPL11154 | Illumina HiSeq 2000, *Homo sapiens* |
| Publication of the series | PMID 31142855 — Lloyd-Price *et al.*, *Nature* 2019, doi:10.1038/s41586-019-1237-2 | HMP2/IBDMDB multi-omics |
| Re-use publication | Priya *et al.*, *Nat Microbiol* 2022, doi:10.1038/s41564-022-01121-z | Week 2 paper |
| BioProject of the RNA-seq series | PRJNA438663 | raw-read layer |
| 16S BioProject of the same cohort | PRJNA398089 | named in Priya *et al.*; **not** the RNA-seq project |
| ArrayExpress cross-reference | E-GEOD-111889 | |
| Mini-case (Task B) | GSE87487, SRA SRP090633 | 20 liver allograft biopsies |

## 3. Folder meanings

```
week3_project/
├── README.md            this file: question, accessions, run order, limitations
├── R/                   numbered scripts; run in order
│   ├── 00_setup.R       install GEOquery etc. + create folders (run once)
│   ├── utils_download.R retry helpers — this network resolves ftp.ncbi.nlm.nih.gov
│   │                    only intermittently, so every fetch is retried with backoff
│   │                    and getGEO() is fed a local file instead of an accession
│   ├── 01_fetch.R       list supplementary files FIRST, then download; records md5 + date
│   ├── 02_inspect_clean.R  5 integrity checks, metadata alignment, low-count filter
│   ├── 03_qc.R          library size, log2-CPM density, PCA by group / site / depth
│   ├── 04_minicase_GSE87487.R  verifies the Table 1 claims with real data
│   └── run_all.R        runs all of the above and writes results/run_log.txt
├── data_raw/            untouched downloads — treat as read-only after fetching
├── data_clean/          derived objects (.rds) only; regenerable from data_raw
├── metadata/            supp manifests, raw phenoData, sample_dictionary.csv
├── results/             QC figures (.png), tables (.csv), run_log.txt
└── sessionInfo.txt      exact package versions of the run
```

## 4. Run order

```r
setwd("D:/BioLession in total/Bioinformatics/week3/week3_project")
source("R/00_setup.R")     # once, installs packages
source("R/run_all.R")      # the whole pipeline, logged to results/run_log.txt
```

## 5. Rules followed (from the lecture)

- Supplementary files are **listed before downloading**; nothing is blind-downloaded.
- Provenance recorded for every download: filename, byte size, URL, md5, date.
- Sample order is never assumed — `stopifnot(identical(colnames(counts), rownames(meta)))`.
- Counts are verified as integers before being treated as counts. Normalized values are never
  rounded into counts.
- CPM/log transforms are used for **visualization only**; statistics would stay on raw counts.
- Low-count filter is stated explicitly: ≥10 counts in ≥ (smallest group size, floor 3) samples,
  with gene numbers before and after logged.
- No sample is deleted to improve a PCA. Outliers are flagged and documented, not removed.

## 6. Known limitations

- The GSE111889 series matrix carries metadata but no expression table, so counts must come
  from the supplementary file `GSE111889_host_tx_counts.tsv.gz` (10.6 MB).
- The cohort is longitudinal with repeated biopsies per subject, so samples are **not**
  independent; any downstream model needs subject as a blocking factor.
- Sequencing-run/batch labels are not exposed in the GEO metadata; biopsy site and library
  depth are used as the observable technical covariates instead.
- Fetch + clean + QC only. Differential expression is out of scope for Week 3.
