## run_all.R -- run the whole pipeline and capture ONE log file to hand back
## Week 3 homework | Zixuan Zeng (SUAT24000114)
## In RStudio: Session > Set Working Directory > To Source File Location's PARENT
## (working directory must be week3_project), then: source("R/run_all.R")

stopifnot(dir.exists("R"))
dir.create("results", showWarnings = FALSE)
log_file <- file.path("results", "run_log.txt")
con <- file(log_file, open = "wt")
sink(con, split = TRUE); sink(con, type = "message")

cat("### Week 3 pipeline run\n")
cat("started :", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("R       :", R.version.string, "\n")
cat("wd      :", getwd(), "\n")

for (s in c("R/01_fetch.R", "R/02_inspect_clean.R", "R/03_qc.R", "R/04_minicase_GSE87487.R")) {
  cat("\n\n##############################################################\n")
  cat("### SOURCING", s, "\n")
  cat("##############################################################\n")
  res <- try(source(s, echo = FALSE), silent = TRUE)
  if (inherits(res, "try-error")) cat("\n!!! ERROR in", s, ":\n", conditionMessage(attr(res, "condition")), "\n")
}

cat("\n\n### sessionInfo\n")
print(sessionInfo())
cat("\nfinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
sink(type = "message"); sink(); close(con)
cat("\nDONE. Send back:", log_file, "and the PNGs in results/\n")
