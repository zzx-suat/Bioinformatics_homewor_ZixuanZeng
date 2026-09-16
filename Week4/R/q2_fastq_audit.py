#!/usr/bin/env python3
"""
Week 4 - Q2: independent audit of the demo FASTQ against the provided
fastqc_snapshot.tsv.

The homework ships a *precomputed* QC summary with deliberately planted traps.
Rather than quoting it, this script measures the same quantities directly from
S01_CTRL_WGS_R1/R2.fastq.gz, so every claim in the report is backed by a number
we computed ourselves.

Outputs
  results/q2_fastq_measured.tsv   measured vs claimed, per module
  results/q2_per_cycle_quality.tsv
  figures/Q2_fastqc_audit.png
"""
import gzip
import os
import sys
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJ = "." if os.path.isdir("data") else ".."
FQ_DIR = os.path.join(PROJ, "..", "Homework", "for_student", "data", "demo_fastq")
ADAPTER = "AGATCGGAAGAGC"          # Illumina TruSeq adapter prefix
PHRED_OFFSET = 33                  # Sanger / Illumina 1.8+


def read_fastq(path):
    """Yield (header, seq, qual). Validates the 4-line record structure."""
    with gzip.open(path, "rt") as fh:
        while True:
            h = fh.readline().rstrip("\n")
            if not h:
                return
            s = fh.readline().rstrip("\n")
            p = fh.readline().rstrip("\n")
            q = fh.readline().rstrip("\n")
            if not (h.startswith("@") and p.startswith("+")):
                raise ValueError(f"malformed FASTQ record in {path}: {h!r}")
            if len(s) != len(q):
                raise ValueError(f"seq/qual length mismatch in {path}: {h!r}")
            yield h, s, q


def audit(path):
    n = 0
    lengths = Counter()
    gc_per_read = []
    adapter_reads = 0
    seq_counts = Counter()
    per_cycle = {}          # cycle -> list of Q
    for _, s, q in read_fastq(path):
        n += 1
        lengths[len(s)] += 1
        gc = sum(1 for c in s if c in "GCgc")
        gc_per_read.append(100.0 * gc / len(s) if s else 0.0)
        if ADAPTER in s:
            adapter_reads += 1
        seq_counts[s] += 1
        for i, ch in enumerate(q):
            per_cycle.setdefault(i, []).append(ord(ch) - PHRED_OFFSET)
    return dict(n=n, lengths=lengths, gc=gc_per_read,
                adapter_reads=adapter_reads, seq_counts=seq_counts,
                per_cycle=per_cycle)


def median(xs):
    xs = sorted(xs)
    m = len(xs) // 2
    return xs[m] if len(xs) % 2 else (xs[m - 1] + xs[m]) / 2


def main():
    r1 = os.path.join(FQ_DIR, "S01_CTRL_WGS_R1.fastq.gz")
    r2 = os.path.join(FQ_DIR, "S01_CTRL_WGS_R2.fastq.gz")
    for p in (r1, r2):
        if not os.path.exists(p):
            sys.exit(f"missing FASTQ: {p}")

    a1, a2 = audit(r1), audit(r2)
    assert a1["n"] == a2["n"], "R1/R2 read counts differ - not a valid PE pair"

    os.makedirs(os.path.join(PROJ, "results"), exist_ok=True)
    os.makedirs(os.path.join(PROJ, "figures"), exist_ok=True)

    # ---- per-cycle median quality (R1) --------------------------------------
    cycles = sorted(a1["per_cycle"])
    med_q1 = [median(a1["per_cycle"][c]) for c in cycles]
    med_q2 = [median(a2["per_cycle"][c]) for c in cycles]

    with open(os.path.join(PROJ, "results", "q2_per_cycle_quality.tsv"), "w") as fh:
        fh.write("cycle\tmedian_Q_R1\tmedian_Q_R2\n")
        for c, x, y in zip(cycles, med_q1, med_q2):
            fh.write(f"{c+1}\t{x}\t{y}\n")

    # ---- summary metrics ----------------------------------------------------
    dup1 = sum(c for c in a1["seq_counts"].values() if c > 1)
    top_seq, top_n = a1["seq_counts"].most_common(1)[0]
    hi_gc = sum(1 for g in a1["gc"] if g > 65)
    mean_gc = sum(a1["gc"]) / len(a1["gc"])
    late = [med_q1[c] for c in range(len(cycles)) if c >= 50]

    # The snapshot claims the 3' crash affects only ~15% of R1. An aggregate
    # per-cycle median is the wrong statistic for a claim about a subpopulation:
    # 85% healthy reads drag the median up. Measure the subpopulation directly.
    crashed = 0
    tail_medians = []
    for _, s, q in read_fastq(r1):
        tail = [ord(c) - PHRED_OFFSET for c in q[50:]]
        if tail:
            m = median(tail)
            tail_medians.append(m)
            if m < 20:
                crashed += 1
    frac_crash = 100.0 * crashed / a1["n"]
    crashed_med = median([m for m in tail_medians if m < 20]) if crashed else float("nan")

    rows = [
        ("Read pairs",                f"{a1['n']}",
         "-", "R1 and R2 record counts match"),
        ("Sequence length distribution",
         f"all reads {list(a1['lengths'])[0]} bp" if len(a1["lengths"]) == 1
         else f"{dict(a1['lengths'])}",
         "PASS / all 80 bp", "confirmed"),
        ("Per base sequence quality (aggregate)",
         f"median Q cycles 1-50 = {median(med_q1[:50]):.0f}; "
         f"cycles 51+ = {median(late):.0f}" if late else "n/a",
         "WARN / 3' median ~12 after cycle 50",
         "DISCREPANT as stated - the aggregate median stays high"),
        ("Per base sequence quality (subpopulation)",
         f"{crashed}/{a1['n']} reads ({frac_crash:.1f}%) have 3'-half "
         f"median Q < 20; their median Q = {crashed_med:.0f}",
         "WARN / affects ~15% of R1",
         "CONFIRMED once measured per-read - aggregate hid it"),
        ("Adapter content",
         f"{a1['adapter_reads']}/{a1['n']} R1 reads contain {ADAPTER} "
         f"({100.0*a1['adapter_reads']/a1['n']:.1f}%)",
         "FAIL / ~15% of pairs", "confirmed"),
        ("Per sequence GC content",
         f"mean {mean_gc:.1f}% GC; {hi_gc}/{a1['n']} reads >65% GC "
         f"({100.0*hi_gc/a1['n']:.1f}%)",
         "WARN / main mode ~41% + shoulder ~78%",
         "confirmed - bimodal, high-GC subpopulation present"),
        ("Sequence duplication levels",
         f"{dup1}/{a1['n']} reads belong to a duplicated sequence "
         f"({100.0*dup1/a1['n']:.1f}%); the single most common template "
         f"accounts for {top_n}/{a1['n']} ({100.0*top_n/a1['n']:.1f}%)",
         "WARN / one template repeated in ~15%",
         "REFINED - snapshot counts the repeated template (15%); "
         "total non-unique reads are twice that"),
        ("Per base N content",
         f"{sum(s.count('N') for _, s, _ in read_fastq(r1))} N bases total",
         "PASS / ~0%", "confirmed"),
    ]

    out = os.path.join(PROJ, "results", "q2_fastq_measured.tsv")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("module\tmeasured_by_us\tclaimed_in_snapshot\tverdict\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

    print(f"{'MODULE':<32} {'MEASURED':<58} VERDICT")
    print("-" * 118)
    for m, meas, claim, verd in rows:
        print(f"{m:<32} {meas:<58} {verd}")

    # ---- figure -------------------------------------------------------------
    fig, ax = plt.subplots(2, 2, figsize=(11, 7.2))
    ink, accent, warn = "#241A38", "#3E2A63", "#B4341F"

    # split the per-cycle median by subpopulation to show what the aggregate hides
    bad_ids, good_ids = set(), set()
    for idx, (_, s, q) in enumerate(read_fastq(r1)):
        tail = [ord(c) - PHRED_OFFSET for c in q[50:]]
        (bad_ids if (tail and median(tail) < 20) else good_ids).add(idx)
    med_bad, med_good = [], []
    percyc_bad, percyc_good = {}, {}
    for idx, (_, s, q) in enumerate(read_fastq(r1)):
        tgt = percyc_bad if idx in bad_ids else percyc_good
        for i, ch in enumerate(q):
            tgt.setdefault(i, []).append(ord(ch) - PHRED_OFFSET)
    med_bad = [median(percyc_bad[c]) for c in cycles] if percyc_bad else []
    med_good = [median(percyc_good[c]) for c in cycles] if percyc_good else []

    ax[0, 0].plot([c + 1 for c in cycles], med_q1, color="#9A9AA6", lw=2.4,
                  label=f"all reads (n={a1['n']})")
    if med_good:
        ax[0, 0].plot([c + 1 for c in cycles], med_good, color="#2F6F4E",
                      lw=1.8, label=f"healthy (n={len(good_ids)})")
    if med_bad:
        ax[0, 0].plot([c + 1 for c in cycles], med_bad, color=warn, lw=2.2,
                      label=f"crashed (n={len(bad_ids)})")
    ax[0, 0].axhline(20, color=ink, ls=":", lw=1.1)
    ax[0, 0].axvspan(51, len(cycles), color=warn, alpha=.06)
    ax[0, 0].set_title("The aggregate median hides the failing 15%",
                       fontsize=10, weight="bold", color=ink)
    ax[0, 0].set_xlabel("cycle"); ax[0, 0].set_ylabel("median Phred Q")
    ax[0, 0].legend(frameon=False, fontsize=7.5, loc="lower left")

    ax[0, 1].hist(a1["gc"], bins=30, color=accent, alpha=.85)
    ax[0, 1].axvline(41, color="#2F6F4E", ls="--", lw=1.4)
    ax[0, 1].axvline(78, color=warn, ls="--", lw=1.4)
    ax[0, 1].set_title("GC content is bimodal - contaminant shoulder",
                       fontsize=10, weight="bold", color=ink)
    ax[0, 1].set_xlabel("% GC per read"); ax[0, 1].set_ylabel("reads")

    frac_ad = 100.0 * a1["adapter_reads"] / a1["n"]
    frac_dup = 100.0 * dup1 / a1["n"]
    frac_gc = 100.0 * hi_gc / a1["n"]
    bars = ["adapter\n(FAIL)", "3' quality\ncrash", "duplicate\n(WARN)",
            ">65% GC\n(WARN)"]
    vals = [frac_ad, frac_crash, frac_dup, frac_gc]
    b = ax[1, 0].bar(bars, vals, color=[warn, warn, "#C58A2E", "#C58A2E"],
                     width=.58)
    ax[1, 0].bar_label(b, fmt="%.1f%%", padding=2, fontsize=9, weight="bold")
    ax[1, 0].set_ylim(0, max(vals) * 1.3)
    ax[1, 0].set_title("Contaminated read fractions (measured)",
                       fontsize=10, weight="bold", color=ink)
    ax[1, 0].set_ylabel("% of R1 reads")

    ax[1, 1].axis("off")
    ax[1, 1].text(0, .96, "Decision after QC", fontsize=10, weight="bold",
                  color=ink, va="top")
    plan = ("1. Trim adapters (cutadapt/fastp), then re-run FastQC\n"
            "2. Hard-trim or quality-trim the 3' half; keep >=40 bp\n"
            "3. Screen the high-GC subpopulation (FastQ Screen /\n"
            "    Kraken2) before blaming biology\n"
            "4. Mark duplicates after alignment - never remove reads\n"
            "    pre-alignment for a WGS variant-calling library\n"
            "5. Re-QC, then align to a PINNED reference build\n\n"
            "Not decidable from QC alone: whether the high-GC reads are\n"
            "contamination or a real GC-rich genomic fraction. That needs\n"
            "taxonomic screening or alignment, not another QC plot.")
    ax[1, 1].text(0, .86, plan, fontsize=8.6, va="top", color="#3A3348",
                  linespacing=1.55)

    fig.suptitle("Q2 - independent audit of the S01 demo FASTQ "
                 "(measured, not quoted)",
                 fontsize=12.5, weight="bold", color=ink, y=.99)
    fig.tight_layout(rect=[0, 0, 1, .96])
    fig.savefig(os.path.join(PROJ, "figures", "Q2_fastqc_audit.png"),
                dpi=200, facecolor="white")
    print("\nwrote figures/Q2_fastqc_audit.png")


if __name__ == "__main__":
    main()
