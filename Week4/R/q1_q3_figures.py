#!/usr/bin/env python3
"""
Week 4 - figures for Q1 (assay strategy) and Q3 (integrated regulatory locus).

Both figures are SCHEMATIC. The Q3 track shapes are drawn from a fixed seed to
illustrate analytical logic; they are not experimental measurements. This mirrors
the course reading, which labels its own toy contact matrix the same way.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.path import Path
import matplotlib.patches as mpatches

PROJ = "." if os.path.isdir("data") else ".."
FIG = os.path.join(PROJ, "figures")
os.makedirs(FIG, exist_ok=True)

INK = "#241A38"
MUTED = "#6B6478"
OBS = "#3E2A63"     # observation-level assay
CAP = "#C58A2E"     # capacity-level assay (reporter)
NEC = "#2F6F4E"     # necessity-level assay (endogenous perturbation)
SEQ = "#2B5D8A"     # sequence-level assay
WARN = "#B4341F"


def box(ax, x, y, w, h, text, fc, tc="white", fs=8.2, weight="normal", alpha=1.0):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.02",
                       linewidth=0, facecolor=fc, alpha=alpha, zorder=2)
    ax.add_patch(p)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fs, color=tc, zorder=3, linespacing=1.35, weight=weight)


def arrow(ax, x1, y1, x2, y2, color=MUTED, style="-|>", lw=1.3, rad=0.0, ls="-"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                 mutation_scale=11, linewidth=lw, color=color,
                                 connectionstyle=f"arc3,rad={rad}",
                                 linestyle=ls, zorder=1))


# =============================================================================
# Q1 — assay strategy: the question drives a staged, branching design
# =============================================================================
def figure_q1():
    fig, ax = plt.subplots(figsize=(12.4, 8.4))
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

    ax.text(50, 97.6, "Q1 — Which mechanism explains Gene X upregulation?",
            ha="center", fontsize=14, weight="bold", color=INK)
    ax.text(50, 94.6,
            "Each stage answers one question and sets a decision rule for the next. "
            "Evidence level, not throughput, orders the assays.",
            ha="center", fontsize=9, color=MUTED)

    L = 21.0          # left edge of all content boxes
    RAIL = 98.0

    # --- stage rail: label sits in the left gutter, on the divider ----------
    stages = [
        (91.6, "STAGE 0", "Establish\nthe phenotype"),
        (76.6, "STAGE 1", "Sequence\nor state?"),
        (59.6, "STAGE 2", "Which\nregulatory layer?"),
        (42.6, "STAGE 3", "Which\ntarget gene?"),
        (26.0, "STAGE 4", "Capacity,\nthen necessity"),
    ]
    for yline, tag, label in stages:
        ax.plot([L - 1.5, RAIL], [yline, yline], color="#E3DFEA", lw=1, zorder=0)
        ax.text(1.0, yline - 1.4, tag, fontsize=7.2, color=MUTED, weight="bold",
                va="top")
        ax.text(1.0, yline - 3.4, label, fontsize=8.5, color=INK, weight="bold",
                va="top", linespacing=1.35)

    # --- stage 0 ------------------------------------------------------------
    box(ax, L, 84.6, 29, 6.2,
        "RNA-seq  +  replicates\nmeasures: transcript abundance", OBS, fs=8.4)
    box(ax, L + 32, 84.6, 31, 6.2,
        "Is the change cell-intrinsic?\ncell sorting / scRNA-seq controls composition",
        "#EDE9F2", tc=INK, fs=7.9)
    arrow(ax, L + 29, 87.7, L + 32, 87.7)
    ax.text(RAIL - 3.5, 87.7, "if bulk only:\ncomposition\nconfounds",
            fontsize=7.1, color=WARN, va="center", ha="right")

    # --- stage 1 ------------------------------------------------------------
    box(ax, L, 69.4, 23, 6.4,
        "WGS / WES\nmeasures: DNA sequence variation", SEQ, fs=8.1)
    box(ax, L + 25, 69.4, 25, 6.4,
        "CAGE / RAMPAGE\nmeasures: capped 5' ends → active TSS", OBS, fs=8.1)
    box(ax, L + 52, 69.4, 25, 6.4,
        "decision rule\ncoding variant → protein assay\nnoncoding → stage 2",
        "#EDE9F2", tc=INK, fs=7.6)
    arrow(ax, L + 11, 84.5, L + 11, 75.9)
    arrow(ax, L + 37, 84.5, L + 37, 75.9)
    arrow(ax, L + 50, 72.6, L + 52, 72.6)

    # --- stage 2 ------------------------------------------------------------
    lay = [
        (L,      "ATAC-seq / DNase", "transposase access\n→ open chromatin"),
        (L + 20, "ChIP-seq / CUT&Tag\nCUT&RUN", "enrichment with a\ntargeted protein/mark"),
        (L + 40, "WGBS / EM-seq", "protected cytosines\n(5mC + 5hmC together)"),
        (L + 60, "RNA Pol II / nascent", "initiation and\npausing state"),
    ]
    for x, title, sub in lay:
        box(ax, x, 52.4, 17, 6.4, title, OBS, fs=7.7, weight="bold")
        ax.text(x + 8.5, 50.6, sub, ha="center", va="top", fontsize=6.9,
                color=MUTED, linespacing=1.3)
        arrow(ax, x + 8.5, 69.3, x + 8.5, 59.1)

    ax.text(50, 45.6,
            "None of these identifies a target gene. Accessibility ≠ occupancy ≠ transcription.",
            ha="center", fontsize=8.2, color=WARN, style="italic")

    # --- stage 3 ------------------------------------------------------------
    box(ax, L, 35.4, 34, 6.4,
        "Hi-C / Micro-C / Capture-C\nmeasures: contact frequency (proximity)",
        OBS, fs=8.1)
    box(ax, L + 38, 35.4, 39, 6.4,
        "eQTL / caQTL / allele-specific reads\nmeasures: genotype ↔ molecular trait",
        OBS, fs=8.1)
    arrow(ax, L + 17, 50.0, L + 17, 42.1)
    arrow(ax, L + 57, 50.0, L + 57, 42.1)
    ax.text(50, 32.2,
            "Proximity is not regulation. A contact prioritizes a link; it does not direct it.",
            ha="center", fontsize=8.2, color=WARN, style="italic")

    # --- stage 4 ------------------------------------------------------------
    box(ax, L, 17.6, 33, 7.2,
        "MPRA / STARR-seq / luciferase\nmeasures: SEQUENCE CAPACITY\nin a construct",
        CAP, fs=8.0)
    box(ax, L + 42, 17.6, 35, 7.2,
        "CRISPRi / deletion / base editing\nmeasures: NECESSITY\nat the native locus",
        NEC, fs=8.0)
    arrow(ax, L + 33, 21.2, L + 42, 21.2, color=INK, lw=1.8)
    ax.text(L + 37.5, 25.6, "capacity\n⇒ necessity", ha="center", fontsize=7.2,
            color=INK, weight="bold", linespacing=1.25)
    arrow(ax, L + 17, 35.3, L + 16, 24.9)
    arrow(ax, L + 57, 35.3, L + 58, 24.9)

    # --- readout ------------------------------------------------------------
    box(ax, L + 6, 6.6, 62, 6.4,
        "Readout: RNA-seq + protein/activity + rescue\n"
        "If RNA rises but protein does not, the mechanism is not transcriptional.",
        "#EDE9F2", tc=INK, fs=7.9)
    arrow(ax, L + 16, 17.5, L + 28, 13.2)
    arrow(ax, L + 58, 17.5, L + 48, 13.2)

    # --- legend -------------------------------------------------------------
    handles = [
        mpatches.Patch(color=SEQ, label="sequence-level"),
        mpatches.Patch(color=OBS, label="observation (association)"),
        mpatches.Patch(color=CAP, label="capacity (construct)"),
        mpatches.Patch(color=NEC, label="necessity (endogenous)"),
    ]
    ax.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
              fontsize=8.2, bbox_to_anchor=(0.5, -0.055))

    ax.text(50, 3.2,
            "The biological question chooses the assay because each assay measures a "
            "different physical quantity,\nand only the one matching the proposed "
            "mechanism can discriminate it from the alternatives.",
            ha="center", fontsize=8.4, color=INK, style="italic", linespacing=1.4)

    fig.subplots_adjust(left=.02, right=.98, top=.98, bottom=.07)
    fig.savefig(os.path.join(FIG, "Q1_workflow.png"), dpi=200, facecolor="white")
    print("wrote figures/Q1_workflow.png")
    plt.close(fig)


# =============================================================================
# Q3 — integrated regulatory locus (schematic tracks)
# =============================================================================
def figure_q3():
    rng = np.random.default_rng(4)
    n = 1200
    x = np.linspace(0, 240, n)          # kb across the locus
    ENH, GENE_Y = 58.0, 175.0           # candidate element, Gene Y TSS

    def peak(centre, height, width):
        return height * np.exp(-0.5 * ((x - centre) / width) ** 2)

    atac = (peak(ENH, 1.0, 2.2) + peak(GENE_Y, 0.85, 2.0) +
            peak(118, 0.22, 3.0) + 0.045 * rng.random(n))
    k27 = (peak(ENH, 0.92, 5.0) + peak(GENE_Y, 0.55, 4.0) +
           0.05 * rng.random(n))
    meth = 0.82 + 0.06 * rng.standard_normal(n)
    meth -= peak(ENH, 0.62, 3.0) + peak(GENE_Y, 0.55, 2.4)
    meth = np.clip(meth, 0.02, 1.0)
    rna = peak(GENE_Y + 14, 0.9, 9.0) + 0.03 * rng.random(n)

    fig = plt.figure(figsize=(12.2, 9.6))
    gs = fig.add_gridspec(7, 1, height_ratios=[1, 1, 1, 1.25, 1, .55, 1.5],
                          hspace=.42, left=.115, right=.975, top=.915, bottom=.155)

    fig.suptitle("Q3 — Is the upstream element an enhancer of Gene Y?",
                 fontsize=14, weight="bold", color=INK, y=.975)
    fig.text(.5, .938,
             "Schematic tracks in one matched condition. Shapes illustrate analytical "
             "logic; they are not experimental data.",
             ha="center", fontsize=8.6, color=MUTED)

    def shade(ax):
        ax.axvspan(ENH - 5, ENH + 5, color=CAP, alpha=.13, zorder=0)
        ax.axvspan(GENE_Y - 5, GENE_Y + 5, color=OBS, alpha=.10, zorder=0)

    def style(ax, label, sub):
        ax.set_xlim(0, 240); ax.set_xticks([])
        ax.set_yticks([])
        for s in ("top", "right", "bottom"):
            ax.spines[s].set_visible(False)
        ax.spines["left"].set_color("#D8D3E0")
        ax.text(-.012, .5, label, transform=ax.transAxes, ha="right", va="center",
                fontsize=9, weight="bold", color=INK)
        ax.text(-.012, .12, sub, transform=ax.transAxes, ha="right", va="center",
                fontsize=6.9, color=MUTED)
        shade(ax)

    ax1 = fig.add_subplot(gs[0]); ax1.fill_between(x, atac, color=OBS, lw=0)
    style(ax1, "ATAC-seq", "transposase access")
    ax2 = fig.add_subplot(gs[1]); ax2.fill_between(x, k27, color="#8C5BA8", lw=0)
    style(ax2, "H3K27ac", "activity-associated mark")
    ax3 = fig.add_subplot(gs[2])
    ax3.vlines(x[::9], 0, meth[::9], color="#2F6F4E", lw=.8)
    ax3.set_ylim(0, 1.05); style(ax3, "DNA methylation", "WGBS/EM-seq: 5mC+5hmC")

    # Hi-C arc
    ax4 = fig.add_subplot(gs[3])
    for h, a, lw_ in ((1.0, .95, 2.6), (.72, .30, 1.4), (.5, .18, 1.0)):
        verts = [(ENH, 0), (ENH + (GENE_Y - ENH) / 2, h), (GENE_Y, 0)]
        ax4.add_patch(mpatches.PathPatch(
            Path(verts, [Path.MOVETO, Path.CURVE3, Path.CURVE3]),
            fc="none", ec=WARN, lw=lw_, alpha=a, zorder=2))
    ax4.plot([0, 240], [0, 0], color="#D8D3E0", lw=1)
    ax4.set_ylim(-.06, 1.15); style(ax4, "Hi-C / Micro-C", "contact frequency")

    ax5 = fig.add_subplot(gs[4]); ax5.fill_between(x, rna, color="#B4341F", lw=0)
    style(ax5, "RNA-seq", "transcript abundance")

    # gene model
    ax6 = fig.add_subplot(gs[5]); ax6.set_xlim(0, 240); ax6.axis("off")
    ax6.add_patch(Rectangle((ENH - 4.5, .42), 9, .3, color=CAP, zorder=3))
    ax6.text(ENH, .95, "candidate element", ha="center", fontsize=8,
             color=CAP, weight="bold")
    ax6.plot([GENE_Y, 232], [.57, .57], color=OBS, lw=1.6, zorder=2)
    for ex0, ex1 in ((GENE_Y, GENE_Y + 9), (GENE_Y + 26, GENE_Y + 38),
                     (GENE_Y + 47, 232)):
        ax6.add_patch(Rectangle((ex0, .40), ex1 - ex0, .34, color=OBS, zorder=3))
    ax6.annotate("", xy=(GENE_Y + 13, .95), xytext=(GENE_Y, .95),
                 arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.4))
    ax6.plot([GENE_Y, GENE_Y], [.57, .95], color=INK, lw=1.4)
    ax6.text(GENE_Y + 42, .95, "Gene Y", fontsize=8.6, color=OBS, weight="bold")
    ax6.set_ylim(0, 1.25)

    # perturbation panel
    ax7 = fig.add_subplot(gs[6])
    conds = ["non-targeting\ncontrol", "CRISPRi\n@ element", "CRISPRi\n@ promoter",
             "element deletion", "CTCF anchor\nperturbation"]
    vals = [1.00, 0.38, 0.12, 0.34, 0.71]
    errs = [0.07, 0.06, 0.04, 0.08, 0.09]
    cols = ["#9A9AA6", NEC, OBS, NEC, "#C58A2E"]
    b = ax7.bar(conds, vals, yerr=errs, capsize=3.5, color=cols, width=.58,
                error_kw=dict(lw=1.1, ecolor=INK))
    ax7.bar_label(b, fmt="%.2f", padding=7, fontsize=8, weight="bold", color=INK)
    ax7.axhline(1.0, color=MUTED, ls=":", lw=1)
    ax7.set_ylim(0, 1.32)
    ax7.set_ylabel("Gene Y RNA\n(relative)", fontsize=8.4, color=INK)
    ax7.tick_params(axis="x", labelsize=7.6, colors=INK)
    ax7.tick_params(axis="y", labelsize=7.6)
    for s in ("top", "right"):
        ax7.spines[s].set_visible(False)
    ax7.set_title("The perturbation panel is what converts correlation into a causal claim "
                  "— everything above is association",
                  fontsize=9, weight="bold", color=INK, pad=8)

    fig.text(.5, .052,
             "Chain:  accessibility → chromatin state → methylation → 3D contact → "
             "expression → perturbation.",
             ha="center", fontsize=8.6, color=INK, weight="bold")
    fig.text(.5, .022,
             "Alternative explanation: the element and Gene Y may both respond to a shared "
             "upstream signal, with the contact structural rather than regulatory.",
             ha="center", fontsize=8.3, color=MUTED)

    fig.savefig(os.path.join(FIG, "Q3_locus_chain.png"), dpi=200,
                facecolor="white")
    print("wrote figures/Q3_locus_chain.png")
    plt.close(fig)


if __name__ == "__main__":
    figure_q1()
    figure_q3()
