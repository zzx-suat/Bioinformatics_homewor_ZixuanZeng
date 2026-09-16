# =============================================================================
# Week 4 — Q4: AI-Assisted Variant Prioritization
# Bioinformatics: From Multi-Omics Data to Discovery (SUAT, Fall 2026)
# Zixuan Zeng / SUAT24000114
#
# Data : data/variants_q4.tsv  (synthetic teaching table — NOT real patient data)
# Out  : results/q4_*.tsv, results/q4_filter_audit.tsv, figures/Q4_prioritization.png
#
# Design principle (Week 4 reading, s.6):
#   "A PASS flag indicates that a record passed the caller's filters; it does not
#    establish rarity, pathogenicity, regulatory function, or causality."
#   => Technical callability is a GATE. Clinical annotation is a RANK, never a gate.
#      Ranking on CLINVAR_SIG first would promote two uncallable records.
# =============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(tidyr); library(stringr)
})

set.seed(1)

# Resolve the project root whether this is run from the project dir or from R/
proj <- if (file.exists(file.path("data", "variants_q4.tsv"))) "." else ".."
stopifnot(dir.exists(file.path(proj, "data")))

variants_path <- file.path(proj, "data", "variants_q4.tsv")
stopifnot(file.exists(variants_path))

v <- read_tsv(variants_path, comment = "#", show_col_types = FALSE) %>%
  mutate(
    variant_id = sprintf("%s:%.0f %s>%s", CHROM, as.numeric(POS), REF, ALT),
    GENE        = na_if(GENE, "."),
    CLINVAR_SIG = na_if(CLINVAR_SIG, ".")
  )

# --- schema assertions: fail loudly rather than silently mis-filter -----------
required_cols <- c("CHROM","POS","REF","ALT","FILTER","DP","GQ","AF",
                   "GENE","CONSEQUENCE","CLINVAR_SIG","CLINVAR_ID","NOTE")
stopifnot(all(required_cols %in% names(v)))
stopifnot(nrow(v) == 12L, !anyDuplicated(v$variant_id))
stopifnot(all(v$DP >= 0), all(v$GQ >= 0), all(v$AF >= 0 & v$AF <= 1))

message(sprintf("Loaded %d variants, %d columns", nrow(v), ncol(v)))

# =============================================================================
# Thresholds — declared BEFORE looking at which variants they remove.
# =============================================================================
MIN_DP  <- 20      # depth floor: below this, allele balance is not assessable
MIN_GQ  <- 30      # genotype quality floor: GQ30 ~ P(wrong genotype) <= 0.001
MAX_AF  <- 1e-3    # rarity ceiling under a severe, penetrant Mendelian model
IMPACT  <- c("stop_gained", "frameshift_variant",
             "splice_acceptor_variant", "splice_donor_variant",
             "missense_variant")

# Clinical evidence is a RANK, not a gate (see header note).
clinvar_rank <- c(
  "Pathogenic" = 6, "Likely_pathogenic" = 5,
  "Conflicting_interpretations_of_pathogenicity" = 4,
  "Uncertain_significance" = 3, "Likely_benign" = 2, "Benign" = 1
)

# =============================================================================
# Tiered filtering. Row counts recorded at EVERY step.
# =============================================================================
audit <- list()
keep_track <- function(df, stage, rule) {
  audit[[length(audit) + 1L]] <<- tibble(
    stage = stage, rule = rule, n_remaining = nrow(df)
  )
  df
}

t0 <- v %>% keep_track("0. All records", "—")

# --- Tier 1: technical callability ------------------------------------------
t1a <- t0  %>% filter(FILTER == "PASS") %>%
  keep_track("1a. FILTER", 'FILTER == "PASS"')
t1b <- t1a %>% filter(DP >= MIN_DP) %>%
  keep_track("1b. Depth", sprintf("DP >= %d", MIN_DP))
t1  <- t1b %>% filter(GQ >= MIN_GQ) %>%
  keep_track("1c. Genotype quality", sprintf("GQ >= %d", MIN_GQ))

# --- Tier 2: population frequency -------------------------------------------
t2 <- t1 %>% filter(AF <= MAX_AF) %>%
  keep_track("2. Population AF", sprintf("AF <= %g", MAX_AF))

# --- Tier 3: functional consequence -----------------------------------------
t3 <- t2 %>% filter(CONSEQUENCE %in% IMPACT) %>%
  keep_track("3. Consequence", "protein/splice-impacting classes")

filter_audit <- bind_rows(audit) %>%
  mutate(n_removed = lag(n_remaining) - n_remaining)

# --- Tier 4: rank survivors on clinical + mechanistic evidence ---------------
shortlist <- t3 %>%
  mutate(
    clin_score  = coalesce(clinvar_rank[CLINVAR_SIG], 0),
    # splice/LoF classes carry more mechanistic prior than missense
    mech_score  = case_when(
      CONSEQUENCE %in% c("splice_acceptor_variant","splice_donor_variant",
                         "stop_gained","frameshift_variant") ~ 2,
      CONSEQUENCE == "missense_variant" ~ 1,
      TRUE ~ 0
    ),
    rarity_score = case_when(AF <= 1e-5 ~ 2, AF <= 1e-4 ~ 1, TRUE ~ 0),
    total_score  = clin_score + mech_score + rarity_score
  ) %>%
  arrange(desc(total_score), desc(GQ), desc(DP)) %>%
  mutate(rank = row_number())

# --- what was dropped, and exactly why (this is the graded part) -------------
reason_for_drop <- function(row) {
  r <- c()
  if (row$FILTER != "PASS")          r <- c(r, sprintf("FILTER=%s", row$FILTER))
  if (row$DP < MIN_DP)               r <- c(r, sprintf("DP=%g < %d", row$DP, MIN_DP))
  if (row$GQ < MIN_GQ)               r <- c(r, sprintf("GQ=%g < %d", row$GQ, MIN_GQ))
  if (row$AF > MAX_AF)               r <- c(r, sprintf("AF=%g > %g", row$AF, MAX_AF))
  if (!row$CONSEQUENCE %in% IMPACT)  r <- c(r, sprintf("consequence=%s", row$CONSEQUENCE))
  if (!length(r)) "kept" else paste(r, collapse = "; ")
}

excluded <- v %>%
  filter(!variant_id %in% shortlist$variant_id) %>%
  rowwise() %>%
  mutate(drop_reason = reason_for_drop(pick(everything()))) %>%
  ungroup() %>%
  mutate(
    # flag the records whose ClinVar label conflicts with their technical quality
    trap = CLINVAR_SIG %in% c("Pathogenic","Likely_pathogenic") &
           str_detect(drop_reason, "FILTER=|DP=|GQ=")
  ) %>%
  arrange(desc(trap), CHROM)

# =============================================================================
# Outputs
# =============================================================================
dir.create(file.path(proj, "results"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(proj, "figures"), showWarnings = FALSE, recursive = TRUE)

write_tsv(filter_audit, file.path(proj, "results", "q4_filter_audit.tsv"))
write_tsv(shortlist %>% select(rank, variant_id, GENE, CONSEQUENCE, FILTER,
                               DP, GQ, AF, CLINVAR_SIG, CLINVAR_ID, total_score),
          file.path(proj, "results", "q4_shortlist.tsv"))
write_tsv(excluded %>% select(variant_id, GENE, CONSEQUENCE, CLINVAR_SIG,
                              FILTER, DP, GQ, AF, drop_reason, trap),
          file.path(proj, "results", "q4_excluded.tsv"))

cat("\n================ FILTER CASCADE ================\n")
print(as.data.frame(filter_audit), row.names = FALSE)
cat("\n================ SHORTLIST (ranked) ================\n")
print(as.data.frame(shortlist %>% select(rank, variant_id, GENE, CONSEQUENCE,
                                         DP, GQ, AF, CLINVAR_SIG, total_score)),
      row.names = FALSE)
cat("\n======= EXCLUDED — note the ClinVar-Pathogenic traps =======\n")
print(as.data.frame(excluded %>% select(variant_id, GENE, CLINVAR_SIG,
                                        drop_reason, trap)), row.names = FALSE)

# =============================================================================
# Figure: filter cascade + why the ClinVar-first shortcut fails
# =============================================================================
pal <- c(kept = "#2F6F4E", dropped = "#9A9AA6", trap = "#B4341F")

casc <- filter_audit %>%
  mutate(stage = factor(stage, levels = stage),
         lbl   = ifelse(is.na(n_removed) | n_removed == 0, "",
                        sprintf("-%d", n_removed)))

p_casc <- ggplot(casc, aes(stage, n_remaining)) +
  geom_col(fill = "#3E2A63", width = .62) +
  geom_text(aes(label = n_remaining), vjust = -0.55, size = 4.2,
            fontface = "bold", colour = "#241A38") +
  geom_text(aes(label = lbl), y = 0.6, colour = "#B4341F",
            size = 3.6, fontface = "bold") +
  scale_y_continuous(limits = c(0, 13.5), expand = c(0, 0)) +
  labs(title = "Q4 filter cascade — technical callability gates first",
       subtitle = sprintf("12 synthetic variants -> %d prioritized | DP>=%d, GQ>=%d, AF<=%g",
                          nrow(shortlist), MIN_DP, MIN_GQ, MAX_AF),
       x = NULL, y = "variants remaining") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor  = element_blank(),
        axis.text.x = element_text(angle = 18, hjust = 1),
        plot.title = element_text(face = "bold"))

trapdat <- v %>%
  mutate(
    status = case_when(
      variant_id %in% shortlist$variant_id ~ "kept",
      CLINVAR_SIG %in% c("Pathogenic","Likely_pathogenic") ~ "trap",
      TRUE ~ "dropped"),
    label = ifelse(is.na(GENE), "(no gene)", GENE)
  )

p_trap <- ggplot(trapdat, aes(DP, GQ, colour = status)) +
  annotate("rect", xmin = MIN_DP, xmax = Inf, ymin = MIN_GQ, ymax = Inf,
           fill = "#2F6F4E", alpha = .07) +
  geom_vline(xintercept = MIN_DP, linetype = "22", colour = "#6B6478") +
  geom_hline(yintercept = MIN_GQ, linetype = "22", colour = "#6B6478") +
  geom_point(size = 3.1) +
  ggrepel::geom_text_repel(aes(label = label), size = 3.1, seed = 1,
                           max.overlaps = 20, show.legend = FALSE) +
  scale_colour_manual(values = pal, name = NULL,
                      labels = c(kept = "prioritized",
                                 dropped = "excluded",
                                 trap = "ClinVar Pathogenic but NOT callable")) +
  labs(title = "Why ClinVar cannot be the first filter",
       subtitle = "MSH2 and MECP2 are ClinVar Pathogenic yet fail depth / genotype quality",
       x = "DP (read depth)", y = "GQ (genotype quality)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  combo <- p_casc / p_trap + plot_layout(heights = c(1, 1.15))
  ggsave(file.path(proj, "figures", "Q4_prioritization.png"), combo,
         width = 9, height = 10, dpi = 200, bg = "white")
} else {
  ggsave(file.path(proj, "figures", "Q4_cascade.png"), p_casc,
         width = 9, height = 4.6, dpi = 200, bg = "white")
  ggsave(file.path(proj, "figures", "Q4_traps.png"), p_trap,
         width = 9, height = 5.4, dpi = 200, bg = "white")
}

writeLines(capture.output(sessionInfo()),
           file.path(proj, "results", "q4_sessionInfo.txt"))
message("Q4 done.")
