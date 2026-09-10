## utils_download.R -- resilient download helper
## Week 3 homework | Zixuan Zeng (SUAT24000114)
## This network resolves ftp.ncbi.nlm.nih.gov only intermittently
## ("Could not resolve hostname"). Every fetch is therefore retried with
## backoff instead of failing the whole pipeline on one bad DNS lookup.

robust_download <- function(url, dest, tries = 10, wait = 6, quiet = TRUE) {
  options(timeout = 3600)
  if (file.exists(dest) && file.size(dest) > 0) {
    cat("  [cache] already on disk:", basename(dest),
        file.size(dest), "bytes\n")
    return(invisible(TRUE))
  }
  for (k in seq_len(tries)) {
    r <- try(suppressWarnings(
           download.file(url, dest, mode = "wb", quiet = quiet)),
         silent = TRUE)
    ok <- !inherits(r, "try-error") && file.exists(dest) && file.size(dest) > 0
    if (ok) {
      cat("  [ok] attempt", k, "->", basename(dest),
          file.size(dest), "bytes\n")
      return(invisible(TRUE))
    }
    if (file.exists(dest)) unlink(dest)
    cat("  [retry", k, "of", tries, "] DNS/connection failed; waiting",
        wait, "s\n")
    Sys.sleep(wait)
  }
  stop("robust_download failed after ", tries, " attempts: ", url)
}

## GEOquery reuses a series matrix already present in destdir, so pre-fetching
## it with the retry helper makes getGEO() itself offline-safe.
prefetch_series_matrix <- function(gse, destdir = "data_raw", ...) {
  stub <- paste0(substr(gse, 1, nchar(gse) - 3), "nnn")
  url  <- sprintf(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s_series_matrix.txt.gz",
    stub, gse, gse)
  dest <- file.path(destdir, paste0(gse, "_series_matrix.txt.gz"))
  dir.create(destdir, showWarnings = FALSE, recursive = TRUE)
  cat("prefetch series matrix:", url, "\n")
  robust_download(url, dest, ...)
  invisible(dest)
}

## The supplementary-file listing scrapes the same flaky host, so retry it too.
robust_supp_list <- function(gse, tries = 10, wait = 6) {
  for (k in seq_len(tries)) {
    s <- try(GEOquery::getGEOSuppFiles(gse, fetch_files = FALSE), silent = TRUE)
    if (!inherits(s, "try-error") && !is.null(s) && NROW(s) > 0) {
      cat("  [ok] supp listing on attempt", k, "->", NROW(s), "file(s)\n")
      return(s)
    }
    cat("  [retry", k, "of", tries, "] supp listing empty/failed; waiting",
        wait, "s\n")
    Sys.sleep(wait)
  }
  stop("robust_supp_list failed for ", gse)
}

## getGEO(<accession>) always lists the remote matrix directory first, so it
## fails on a bad DNS lookup even when the file is already cached. Fetching the
## file ourselves (with retries) and handing GEOquery the local path via
## getGEO(filename=) removes that network dependency entirely.
load_series_matrix <- function(gse, destdir = "data_raw", ...) {
  dest <- prefetch_series_matrix(gse, destdir, ...)
  cat("parsing local series matrix:", dest, "\n")
  e <- GEOquery::getGEO(filename = dest, getGPL = FALSE)
  if (is.list(e) && !methods::is(e, "ExpressionSet")) e <- e[[1]]
  cat("  parsed:", class(e), " dim:", paste(dim(e), collapse = " x "), "\n")
  e
}
