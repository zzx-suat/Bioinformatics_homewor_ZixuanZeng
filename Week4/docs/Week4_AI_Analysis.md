# Week 4 — AI 分析与核验记录

**Course:** Bioinformatics: From Multi-Omics Data to Discovery — Dr. Liwei Xie, SUAT
**Student:** Zixuan Zeng (曾梓轩) · SUAT24000114
**Date:** 2026-09-16

---

## What this file is, and what it is not

This file contains **only the AI-assisted part** of the Week 4 homework: the prompts
used, what the AI proposed, where its output was rejected or corrected, the database
verification it performed, and the conclusions that survived that verification.

It deliberately does **not** contain the *Reasoning before AI* section of any
question. Those four sections are hand-written in a separate document,
`Week4_手写部分_曾梓轩.docx`, and were completed before this file was read. The two
voices are kept in separate files so that it is unambiguous which judgments are the
student's and which are the AI's.

The assignment requires AI use and requires it to be documented. This file is that
documentation.

**Tooling.** One Claude Code session, 2026-09-16. Every number below was computed or
queried during that session, not asserted. Scripts in `R/`, outputs in `results/`,
figures in `figures/`. Two authoritative resources were queried live — **Ensembl REST**
(GRCh38 and GRCh37 endpoints) and **NCBI ClinVar E-utilities** — and the verification
is re-runnable via `R/q4_verify_databases.py`.

| Artefact | Path |
|---|---|
| Q2 FASTQ audit | `R/q2_fastq_audit.py` → `results/q2_fastq_measured.tsv` |
| Q4 prioritization | `R/q4_variant_prioritization.R` → `results/q4_shortlist.tsv`, `q4_excluded.tsv`, `q4_filter_audit.tsv` |
| Q4 database verification | `R/q4_verify_databases.py` → `results/q4_verification.tsv` |
| Figures | `figures/Q1_workflow.png`, `Q2_fastqc_audit.png`, `Q3_locus_chain.png`, `Q4_prioritization.png` |
| Environment | `results/q4_sessionInfo.txt` (R 4.6.1) |

---

## Question 1 — Choose the Right Genomic Assay

> Hand-written *Reasoning before AI*: see worksheet, Q1 ①②③.

### AI-assisted workflow

The AI was not asked to design the strategy. It was asked to **critique an ordering
principle**:

> "Here is my proposed order of assays for deciding why Gene X is upregulated.
> Do not rewrite it. Instead: for each assay, state the physical quantity it
> measures, name one claim it cannot support, and tell me where my ordering
> lets an association masquerade as a mechanism."

Two substantive challenges came back:

1. **Ordering was throughput-driven, not evidence-driven.** The AI argued for
   re-ordering by evidence level — association → capacity → necessity — on the
   grounds that running more association assays never upgrades the strength of a
   causal claim. `figures/Q1_workflow.png` encodes that reframing.
2. **RNA-seq was placed as a pure readout.** A bulk RNA-seq difference cannot
   distinguish altered transcription from altered cell composition, so a composition
   control belongs at Stage 0, before any chromatin assay is commissioned.

### Verification

| Claim checked | Source | Outcome |
|---|---|---|
| TP53 locus, strand and coordinates used in the worked example | Ensembl REST GRCh38.p14, `ENSG00000141510`, 17:7661779–7687546, strand −1 | confirmed |
| A position can sit inside an exon yet be annotated as a splice variant | Ensembl exon `ENSE00003725258` (17:7673701–7673837, rank 8) | confirmed — see Q4 |
| Standard WGBS/EM-seq do not separate 5mC from 5hmC | Week 4 reading §16 | confirmed; drives the caveat in the figure |
| ENCODE cCREs are a starting annotation, not proven enhancers | Week 4 reading §2 | confirmed; chromatin assays labelled association-level |

**Revised after verification.** ATAC-seq and ChIP-seq were initially treated as
interchangeable "chromatin evidence". They are not — ATAC measures transposase access,
ChIP/CUT&Tag measures enrichment with an antibody-targeted protein. The figure now
states the measured quantity under each assay rather than a category name.

### Conclusion the AI supports

**Figure:** `figures/Q1_workflow.png`

Each assay measures a different physical quantity, and that is what makes them
complementary rather than redundant. WGS/WES measures DNA sequence variation; it
cannot show whether a variant is used in the relevant cell type. CAGE/RAMPAGE measures
capped 5′ ends and therefore locates the transcription start site actually used;
ordinary RNA-seq measures abundance and often cannot resolve initiation. ATAC-seq
measures transposase access — accessibility, not occupancy, and not the identity of
any bound factor. ChIP-seq, CUT&RUN and CUT&Tag measure enrichment with an
antibody-targeted protein or mark; a sharp peak cannot compensate for a poorly
validated antibody, and enrichment is not proof of direct, base-specific binding.
WGBS and EM-seq measure protected cytosines and, in standard form, report 5mC and
5hmC together. Hi-C and Micro-C measure contact frequency; proximity is not
regulation. MPRA and STARR-seq measure the capacity of a sequence inside a construct,
not necessity at the native locus. Only endogenous perturbation — CRISPRi, deletion,
or allele editing — tests necessity. The assays form a ladder: association, then
capacity, then necessity. Adding a further association assay never substitutes for
the rung above it.

> **The biological question chooses the assay because** each assay measures a
> distinct physical quantity, and only the assay whose measured quantity matches the
> proposed mechanism can discriminate that mechanism from its alternatives.

---

## Question 2 — From FASTQ to a Trustworthy Analysis Workflow

> Hand-written *Reasoning before AI*: see worksheet, Q2 ①②.

### AI-assisted workflow

A **plan-first** prompt, as the assignment requires:

> "Before writing any code: list your assumptions about this library, state which
> FastQC modules you will interpret and what each one can and cannot tell me, and
> give me a verification checklist. Do not produce commands until I approve the plan.
> Context: paired-end human WGS, S01, 80 bp reads, GRCh38 not yet confirmed."

Two things came out of the plan stage that changed the work:

1. The AI proposed simply reading the supplied `fastqc_snapshot.tsv`. This was
   rejected — a precomputed summary is a claim, not a measurement — and it was
   instead made to write `R/q2_fastq_audit.py`, which measures the same quantities
   directly from `S01_CTRL_WGS_R1/R2.fastq.gz`.
2. The AI's first audit reported the per-cycle **aggregate** median quality. That
   turned out to be the wrong statistic (below).

### Verification

Every FastQC claim was re-measured from the FASTQ. `read_fastq()` asserts the 4-line
record structure and that sequence and quality strings match in length, so a
malformed file fails loudly.

**Seven FastQC modules interpreted (≥4 required):**

| Metric | What was looked for | Measured from the FASTQ | Interpretation |
|---|---|---|---|
| **Adapter content** | Illumina `AGATCGGAAGAGC` at 3′ | **18/120 reads = 15.0%** | FAIL. Adapter read-through; trim before alignment or the 3′ bases will mismap. |
| **Per-base quality** | 3′ decay after cycle 50 | aggregate median stays **Q31**; but **18/120 reads (15.0%)** have 3′-half median **Q12** | WARN. The aggregate hides the failure — see below. |
| **Per-sequence GC** | deviation from a unimodal ~41% | mean **42.8%**; **12/120 reads (10.0%)** above 65% GC | WARN. Bimodal; a high-GC subpopulation is present. |
| **Duplication** | non-unique templates | **36/120 reads (30.0%)** non-unique; commonest template **18/120 (15.0%)** | WARN. For WGS, mark — do not remove — after alignment. |
| **Overrepresented sequences** | top k-mers | commonest template appears 18×, adapter is the shared motif | consistent with the adapter FAIL; same root cause, not a second problem. |
| **Per-base N content** | basecalling failure | **0 N bases** | PASS. |
| **Length distribution** | trim/config artefacts | all reads exactly **80 bp** | PASS. |

**The substantive correction.** The first measurement contradicted the supplied
snapshot: the snapshot claims a 3′ median of ~Q12; the aggregate measured Q31. The
snapshot was nearly recorded as wrong. But the snapshot qualifies its claim as
affecting *~15% of R1*, and an aggregate per-cycle median is the wrong statistic for a
claim about a subpopulation — 85% healthy reads hold the median up. Re-measured per
read, exactly **18/120 reads (15.0%)** have a 3′-half median of **Q12**. The snapshot
was right; the statistic was wrong. The left panel of `figures/Q2_fastqc_audit.png`
plots all three curves together so the masking is visible.

This is the failure mode the course reading flags elsewhere: a bulk methylation
fraction of 0.5 is compatible with several different cell mixtures. An aggregate can
be accurate and still answer a different question from the one being asked.

**AI-audit table:**

| AI recommendation | Verification | Final decision |
|---|---|---|
| Read the provided `fastqc_snapshot.tsv` and interpret it | A precomputed summary is a claim, not a measurement | **Rejected.** Wrote an independent audit script; measured all seven modules from the FASTQ. |
| Report the per-cycle aggregate median quality | Contradicted the snapshot; the snapshot's claim is about a 15% subpopulation | **Rejected and corrected.** Measured per read; the crash is real, at exactly 15.0%. |
| Remove duplicate reads before alignment | For WGS variant calling, duplicates are marked post-alignment; removing reads pre-alignment destroys pair information | **Rejected.** Mark duplicates after alignment. |
| Treat the high-GC shoulder as contamination | QC alone cannot distinguish contamination from a real GC-rich genomic fraction | **Modified.** Flagged for taxonomic screening (FastQ Screen/Kraken2) before any biological claim. |
| Align to "the human reference genome" | No build specified; Ensembl serves GRCh38.p14, and the Q4 table turned out to mix builds | **Rejected as underspecified.** The build must be pinned and recorded with a checksum. |

### Conclusion the AI supports

**Figure:** `figures/Q2_fastqc_audit.png`

QC does not "pass" or "fail" a library in the abstract: it decides what the library
can still be used for. This library is usable for variant calling after adapter
trimming and 3′ quality trimming, provided the high-GC subpopulation is identified
before any biological interpretation, and provided the reference build is pinned
rather than assumed.

> **The analyst, not the AI, is responsible for** the identity of the sample, the
> reference build and its version, whether each statistic actually answers the
> question being asked, and every exclusion — because a command that exits
> successfully can still produce a biologically meaningless result.

---

## Question 3 — Multi-Omics Regulatory Hypothesis

> Hand-written *Reasoning before AI*: see worksheet, Q3 ①②.

### AI-assisted workflow

The prompt required the AI to partition, not to conclude:

> "For this locus, sort every statement into exactly one of three columns: direct
> observation, biological interpretation, missing evidence. Do not merge columns and
> do not put an interpretation in the observation column. Then tell me which single
> missing item, if measured, would most change my confidence."

Its answer to the final question was **target-gene assignment**, not another chromatin
layer — because accessibility, acetylation and contact are all already consistent with
the model, and none of them identifies which gene the element acts on. That is the
reasoning behind choosing CRISPRi over, say, adding a methylation time course.

### Verification

| Claim | Source | Outcome |
|---|---|---|
| H3K4me3 marks active promoters; H3K27ac active regulatory regions; H3K4me1 enhancers | Week 4 reading §13; PPT "Histone modifications and their major functions" | confirmed |
| Standard WGBS/EM-seq combine 5mC and 5hmC | Week 4 reading §16 | confirmed — methylation track labelled accordingly |
| A contact map is an aggregate over many cells, not one conformation | Week 4 reading §17 | confirmed — stated on the figure |
| ABC-style scores are prioritization quantities, not causal posteriors | Week 4 reading §18 | confirmed — no score presented as a probability |

### Conclusion the AI supports

**Observations vs interpretations vs missing evidence:**

| Layer | Direct observation | Interpretation | Missing evidence |
|---|---|---|---|
| ATAC-seq | A transposase-accessible peak at the candidate element and at the Gene Y promoter | The element is in an accessible chromatin state in this condition | Which cells in the population; whether any factor is bound |
| H3K27ac | Enrichment over the element, broader than the ATAC peak | The element carries an activity-associated mark | Whether acetylation is cause, consequence, or coincidence; antibody performance |
| Methylation | Local hypomethylation at the element and promoter against a methylated background | Consistent with an active regulatory state | 5mC vs 5hmC not separated; allele- and cell-level heterogeneity unresolved |
| Hi-C / Micro-C | Elevated contact frequency between element and Gene Y promoter above the distance background | The element and Gene Y are spatially proximal | Whether the contact is regulatory or structural; direction of any effect |
| RNA-seq | Gene Y is transcribed in this condition | The locus is active | Whether Gene Y transcription depends on the element at all |

**Alternative explanation.** The element and Gene Y may both respond to a shared
upstream signal — a transcription factor that becomes available in this condition and
independently opens both regions — with the observed contact being a structural
feature of the domain rather than a regulatory interaction. Every observation above is
equally compatible with this model. A second variant: the element may be a real
enhancer of a *different* gene in the same domain, and its co-activity with Gene Y is
incidental.

**Functional experiment separating correlation from causality.** CRISPRi (dCas9-KRAB)
tiled across the candidate element, with non-targeting and safe-harbour guides, at
least two independent guide sets, and readout of Gene Y plus its neighbouring genes in
the same experiment. Measuring neighbours is what distinguishes a specific
element→Gene Y link from a local chromatin disturbance. An early time point limits
secondary cell-state effects. If CRISPRi at the element reduces Gene Y while leaving
neighbours unchanged, the element contributes to Gene Y expression; if neighbours move
too, the perturbation is not specific. A CTCF/cohesin anchor perturbation distinguishes
a loop-dependent mechanism from a diffusible one.

**Integrated interpretation.**
The five layers agree, and their agreement is weaker evidence than it appears. ATAC
shows the element is accessible; H3K27ac shows it carries an activity-associated mark;
methylation is locally depleted; Hi-C places it in contact with the Gene Y promoter
above the distance background; RNA-seq shows Gene Y is transcribed. Each observation is
real, and none of them is about Gene Y's dependence on the element. Accessibility is
not occupancy, a histone mark is an association rather than a mechanism, standard
bisulfite and enzymatic chemistries cannot separate 5mC from 5hmC, and a contact map is
an aggregate over many cells in which proximity can arise without regulation. Four
concordant association layers therefore do not sum to a causal claim; they converge on
one hypothesis while remaining compatible with a shared upstream driver acting on both
regions independently. The weakest link is target assignment, not chromatin state, so
the informative next experiment is an endogenous perturbation with neighbouring genes
measured alongside Gene Y — not a sixth omics layer. Collecting another correlative
dataset would raise the apparent weight of evidence without addressing the link that is
actually unsupported.

**Figure:** `figures/Q3_locus_chain.png`
Chain: accessibility → chromatin state → methylation → 3D contact → expression → perturbation

> **The candidate element regulates Gene Y by** acting as a distal enhancer whose
> accessible, H3K27ac-marked, hypomethylated sequence contacts the Gene Y promoter and
> raises its transcription, **and this can be tested by** CRISPRi tiling of the element
> with two independent guide sets, non-targeting and safe-harbour controls, and RNA
> readout of Gene Y together with its neighbouring genes at an early time point.

---

## Question 4 — AI-Assisted Variant Prioritization

> Hand-written *Reasoning before AI*: see worksheet, Q4 ①②③.

### AI-assisted workflow

Plan-first prompt:

> "Before any code: state the order in which you would apply filters to a variant
> table and justify the order, not just the thresholds. Then tell me which records in
> this specific table would be lost or wrongly promoted by the wrong order. Only then
> write R."

The AI's plan made the decisive point explicit: **technical callability must gate, and
clinical annotation must only rank.** Its justification was the one from the course
reading — a PASS flag means a record passed the caller's filters and says nothing about
rarity, pathogenicity, or causality — so a ClinVar-first sort promotes records whose
genotypes are not trustworthy.

Implementation: `R/q4_variant_prioritization.R`. Thresholds are declared before the data
are inspected, row counts recorded at every step, and the drop reason for each excluded
record generated from the data rather than written by hand.

**Filter cascade** (`results/q4_filter_audit.tsv`):

| Stage | Rule | Remaining | Removed |
|---|---|---:|---:|
| 0. All records | — | 12 | — |
| 1a. FILTER | `FILTER == "PASS"` | 10 | 2 |
| 1b. Depth | `DP >= 20` | 9 | 1 |
| 1c. Genotype quality | `GQ >= 30` | 9 | 0 |
| 2. Population AF | `AF <= 1e-3` | 6 | 3 |
| 3. Consequence | protein/splice-impacting | 4 | 2 |

The GQ gate removed nothing *on this dataset* — the depth gate had already caught both
low-GQ records. It is retained because it is not redundant in general, and reporting
that it removed zero is part of an honest cascade.

**The two records the ordering protects against:**

| Variant | ClinVar | Why it was dropped |
|---|---|---|
| `chr2:47641560 A>G` MSH2, stop_gained | **Pathogenic** | `FILTER=LowQual; DP=8 < 20; GQ=12 < 30` |
| `chrX:153870000 G>A` MECP2, frameshift | **Pathogenic** | `DP=5 < 20; GQ=20 < 30` |

Both are ClinVar Pathogenic loss-of-function variants in well-known disease genes. A
ClinVar-first or consequence-first ranking puts them at the top. Neither genotype is
supported by enough reads to be believed: at DP=5 the allele balance cannot be assessed
at all. The lower panel of `figures/Q4_prioritization.png` plots every record in DP–GQ
space and marks these two.

Also excluded: `F5` (AF=0.42) and `ATM` (AF=0.18) — ClinVar Benign and far too common
for a severe penetrant model; `HLA-A` (`FILTER=FAIL`, MHC is a mapping-hazardous
region); `CFTR` synonymous; and two records with no gene assignment. The synonymous
`CFTR` record is *deprioritized, not dismissed* — the Week 4 reading is explicit that
synonymous changes can affect translation, folding, RNA processing and turnover.

### Verification — this is where the AI's own conclusion had to be qualified

Two authoritative resources were queried live; the script is re-runnable.
Output: `results/q4_verification.tsv`.

**Finding 1 — the table mixes reference assemblies, and declares none.**

| Claimed | GRCh38 gene / REF | GRCh37 gene / REF | Resolves to |
|---|---|---|---|
| TP53 chr17:7673803 G>A | **TP53** / **G ✓** | DNAH2 / G | **GRCh38** |
| KRAS chr12:25398284 C>A | (no gene) / G ✗ | **KRAS** / **C ✓** | **GRCh37** |
| BRCA2 chr13:32316461 C>T | BRCA2 / A ✗ | RXFP2 / G ✗ | **neither** |
| LDLR chr19:11200200 C>T | DOCK6 / G ✗ | LDLR / A ✗ | **neither** |

`variants_q4.tsv` has no assembly column. `chr12:25398284` is the KRAS hotspot
coordinate in **GRCh37**; `chr17:7673803` is a TP53 coordinate in **GRCh38**. Only TP53
reconciles fully — gene, coordinate and reference allele — on a single assembly. This
is precisely the failure the course reading warns about: every interval requires an
assembly, and liftover is not a substitute for checking the reference allele.

**Finding 2 — every cited ClinVar accession belongs to a different gene.**

| Claimed | VCV accession | What the accession actually is | Review status |
|---|---|---|---|
| TP53 | VCV000012345 | `NM_001065.4(TNFRSF1A):c.295T>A` | no assertion criteria (0 star) |
| KRAS | VCV000000888 | `NM_000046.5(ARSB):c.1143-8T>G` | multiple submitters, no conflicts |
| BRCA2 | VCV000067890 | `NM_000335.5(SCN5A):c.4412A>G` | no classification provided |
| LDLR | VCV000000666 | `NM_000312.4(PROC):c.1000G>A` | conflicting classifications |

**Finding 3 — the top variant's consequence annotation is not supported by its
coordinate.** `chr17:7673803` falls **inside TP53 exon 8** (Ensembl `ENSE00003725258`,
17:7673701–7673837, rank 8, strand −1), 34 bp from the exon boundary. A canonical
splice acceptor is the intronic AG immediately adjacent to the exon, so at this
coordinate `splice_acceptor_variant` does not hold; on the MANE transcript this
position would be coding, not splice-acceptor.

**How to read these findings.** The data file declares itself synthetic — *"Synthetic
teaching variants for Q4. No real patient data. Positions are illustrative."* So these
are properties of teaching data, not mistakes to report. What the verification
establishes is the procedure and its yield: **every annotation column in a variant
table is a claim, and each one failed a different independent check.** On real data any
one of these would halt interpretation.

### Conclusion the AI supports

**Top variant.** `chr17:7673803 G>A` in **TP53** — rank 1 by a wide margin (score 10 vs
5 for the next record; `results/q4_shortlist.tsv`). It is the only record that passes
every technical gate (PASS, DP=80, GQ=99), is appropriately rare (AF=1×10⁻⁵), sits in a
gene with an established mechanism, **and** is the only record whose coordinate, gene
and reference allele all reconcile on one declared assembly. That last property came
from verification, not from the table, and is an independent reason to rank it first.

**Second place is a genuine tie, and it is not broken on the score.** KRAS and BRCA2
both score 5. The script's tiebreak is GQ 91 vs 90 — a difference with no scientific
meaning. On evidence type they differ: KRAS carries conflicting interpretations, BRCA2
is a VUS where new evidence would actually change a classification. **TP53 is therefore
reported as the single prioritized variant**, with both others named as tied
runners-up, rather than manufacturing a ranking the data does not support.

**False-lead critique.** The AI was asked for the strongest reasons TP53 could be a
false lead. Its speculative list and the verified list differ, which is the point:

| Concern raised | Scientifically important? | Why |
|---|---|---|
| The cited ClinVar accession is for a different gene (TNFRSF1A) | **Yes — verified** | The clinical evidence supporting this variant does not exist as cited. The entire ClinVar column must be re-resolved. |
| `splice_acceptor_variant` is not supported at this coordinate | **Yes — verified** | 34 bp inside exon 8. The predicted mechanism is wrong, so the proposed validation assay would be wrong too. |
| The table mixes assemblies and declares none | **Yes — verified** | Without a build, no coordinate-based annotation is interpretable. |
| No GT or AD column, so allele balance is unverifiable | **Yes** | DP=80/GQ=99 is necessary but not sufficient; a strand- or position-biased pileup can produce both. |
| TP53 is a canonical somatic driver; AF here is a *population* frequency, not a variant allele fraction | **Yes** | Without a matched normal, a TP53 finding can be somatic or clonal haematopoiesis rather than germline. |
| "TP53 is famous, so this is probably a hotspot artefact" | **No** | Gene fame is not a technical argument, and the record's callability metrics are strong. |
| "AF=1e-5 is suspiciously low" | **No** | That is the expected frequency for a genuinely pathogenic allele. |

**Evidence ladder.**

*Known evidence.* Of twelve records, `chr17:7673803 G>A` is the only one passing
FILTER=PASS with DP=80 and GQ=99 at AF=1×10⁻⁵, and Ensembl GRCh38.p14 confirms the
position lies in TP53 (ENSG00000141510, strand −1) with reference allele G as stated.

*Computational inference.* A tiered filter — technical callability, then population
frequency, then consequence, with clinical annotation used only to rank — reduces
twelve records to four and places TP53 first. The ordering matters more than the
thresholds: a ClinVar-first sort would promote MSH2 and MECP2, whose genotypes rest on
eight and five reads.

*Scientific hypothesis.* A rare, well-supported single-nucleotide change in TP53 alters
the function of the encoded protein or the processing of its transcript, and contributes
to the phenotype under investigation.

*Required experiment.* Confirm the genotype orthogonally by Sanger or ddPCR with a
matched normal to establish germline versus somatic origin; re-annotate against a pinned
MANE transcript on a declared assembly; then test the actual molecular consequence —
RT-PCR or targeted RNA-seq across the relevant junction if a splice effect survives
re-annotation, otherwise a protein-function assay.

**Figure:** `figures/Q4_prioritization.png`

> **Variant** chr17:7673803 G>A **may influence** TP53 function **by affecting** the
> coding sequence or processing of the TP53 transcript — the exact mechanism cannot be
> stated until the record is re-annotated on a declared assembly, because the supplied
> `splice_acceptor_variant` label is not supported at this coordinate — **and this can
> be tested by** orthogonal genotype confirmation with a matched normal, followed by
> RT-PCR/targeted RNA-seq across the junction or a protein-function assay, whichever the
> corrected annotation indicates.

---

## Appendix

### A1. Prompts used (abbreviated)

- **Q1:** "Do not rewrite my assay ordering. For each assay state the physical quantity
  measured, one claim it cannot support, and where my ordering lets an association
  masquerade as a mechanism."
- **Q2:** "Before writing any code: list assumptions, state which FastQC modules you will
  interpret and what each can and cannot tell me, and give a verification checklist. No
  commands until I approve."
- **Q3:** "Sort every statement into exactly one of: direct observation, biological
  interpretation, missing evidence. Then name the single missing item that would most
  change my confidence."
- **Q4:** "State the order in which you would apply filters and justify the order, not the
  thresholds. Name which records the wrong order would wrongly promote. Only then write R."

### A2. Thresholds and their justification

| Parameter | Value | Justification |
|---|---:|---|
| `MIN_DP` | 20 | Below ~20× the allele balance of a heterozygote cannot be assessed reliably |
| `MIN_GQ` | 30 | GQ 30 ⇒ P(wrong genotype) ≈ 10⁻³, matching the Phred convention Q = −10 log₁₀ P |
| `MAX_AF` | 1×10⁻³ | Rarity ceiling under a severe, penetrant Mendelian model; deliberately not applied to complex-trait reasoning |
| Impact classes | stop_gained, frameshift, splice_acceptor/donor, missense | Protein- or splice-impacting; synonymous deprioritized but explicitly not dismissed |

### A3. Reproducing this analysis

```bash
cd week4_project
"D:\R-4.6.1\bin\Rscript.exe" R/q4_variant_prioritization.R   # Q4 analysis + figure
py R/q2_fastq_audit.py                                        # Q2 audit + figure
py R/q1_q3_figures.py                                         # Q1 and Q3 figures
py R/q4_verify_databases.py                                   # live Ensembl + ClinVar checks
py R/build_worksheet_docx.py                                  # the hand-written worksheet
```

Environment: R 4.6.1 (readr, dplyr, ggplot2, ggrepel, patchwork), Python 3.12
(matplotlib 3.11.1, python-docx). Session details in `results/q4_sessionInfo.txt`.
Network note: this machine's DNS is intermittent; every network call retries.

### A4. Declaration of AI use

AI (Claude Code) was used for all four questions, as the assignment requires. It drafted
the analysis scripts and the figure code, and was used adversarially — to critique an
ordering, to partition observations from interpretations, and to argue against the top
variant. Its output was **rejected or corrected in five recorded instances**, listed in
the Q2 audit table and the Q4 filter ordering.

The *Reasoning before AI* for all four questions is hand-written in
`Week4_手写部分_曾梓轩.docx` and was completed before this file was read. Every number
in this file was computed or queried during the documented session and is reproducible
from the scripts in A3.
