# -*- coding: utf-8 -*-
"""make_paper_figs.py —— 生成投稿用 Figure 1（研究设计）与 Figure 4（chr9q34.2 区域图）"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from scipy.stats import norm

ROOT = r"D:/GWAS/Migraine"
RES = os.path.join(ROOT, "results")
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})

# ---------------------------------------------------------------- Figure 1
def fig1():
    fig, ax = plt.subplots(figsize=(7.4, 8.4))
    ax.set_xlim(0, 10); ax.set_ylim(0, 12.6); ax.axis("off")

    def box(x, y, w, h, title, body, fc="#eef2f6", ec="#4a5a6a", fs=8.2, tfs=8.8):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08,rounding_size=0.12",
                                    fc=fc, ec=ec, lw=0.9))
        ax.text(x + w / 2, y + h - 0.22, title, ha="center", va="top",
                fontsize=tfs, fontweight="bold", color="#1f2d3a")
        ax.text(x + w / 2, y + h - 0.62, body, ha="center", va="top", fontsize=fs, color="#2b3a47")

    def arrow(x1, y1, x2, y2):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                     mutation_scale=11, lw=0.9, color="#6b7c8c"))

    box(0.35, 10.9, 4.3, 1.45, "Thyroid disease (exposure)",
        "GBMI, European ancestry\nGraves 1,866,840 · Hypothyroid 1,771,600\n"
        "Lymphocytic thyroiditis 1,405,560\nBenign nodule 1,792,740 · Cancer 2,091,280", fc="#e8eef5")
    box(5.35, 10.9, 4.3, 1.45, "Migraine (outcome)",
        "FinnGen data freeze 13\nOverall 28,504 cases / 366,556 controls\n"
        "With aura 12,743 · Without aura 10,206", fc="#eaf2ea")

    box(2.55, 9.15, 4.9, 1.1, "Instrument selection",
        "P < 5×10⁻⁸ · r² < 0.001 within 10 Mb · F > 10\n39–336 variants per exposure")
    arrow(2.5, 10.85, 4.3, 10.3)
    arrow(7.5, 10.85, 5.7, 10.3)

    box(2.55, 7.55, 4.9, 1.15, "Two-sample Mendelian randomization",
        "IVW (primary) + MR-Egger, weighted median,\nsimple and weighted mode; 15 combinations")
    arrow(5.0, 9.1, 5.0, 8.75)

    box(0.35, 5.6, 4.4, 1.5, "Sensitivity",
        "Cochran Q · MR-Egger intercept\nBonferroni α = 0.0033\n"
        "Minimum detectable OR (80% power)", fc="#f6f2e8")
    box(5.25, 5.6, 4.4, 1.5, "Sample overlap",
        "Analytic correction (Burgess)\nCross-trait LDSC intercept\n"
        "MRlap re-estimation", fc="#f6f2e8")
    arrow(4.2, 7.5, 2.55, 7.15)
    arrow(5.8, 7.5, 7.45, 7.15)

    box(2.55, 3.85, 4.9, 1.2, "Shared genetic architecture",
        "Cross-trait LD score regression\n18 trait pairs · HapMap 3, MHC excluded")
    arrow(5.0, 5.55, 5.0, 5.1)

    box(2.55, 2.35, 4.9, 1.05, "Colocalization",
        "coloc.abf · 320 loci · PP.H4 > 0.8 = high confidence")
    arrow(5.0, 3.8, 5.0, 3.45)

    box(2.55, 0.95, 4.9, 0.95, "eQTL annotation", "GTEx v8 thyroid · lead and regional genes")
    arrow(5.0, 2.3, 5.0, 1.95)

    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    plt.savefig(os.path.join(RES, "fig1_design.png"), dpi=300, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------- Figure 4
def fig4(locus=178, exposure="BNG", outcome="MigraineAura"):
    d = pd.read_csv(os.path.join(ROOT, "coloc_input", "coloc_long.tsv"), sep="\t")
    b = d[(d.locus == locus) & (d.exposure == exposure) & (d.outcome == outcome)].copy()
    pos = pd.read_csv(os.path.join(ROOT, "ldsc_cache", f"{exposure}.tsv"), sep="\t", header=None,
                      usecols=[0, 2], names=["snp", "bp"])
    b = b.merge(pos.drop_duplicates("snp"), on="snp", how="left").dropna(subset=["bp"])
    b["bp"] = b["bp"] / 1e6
    b["z_x"] = b.beta_x / b.se_x
    b["z_y"] = b.beta_y / b.se_y
    b["p_x"] = -np.log10(2 * norm.sf(np.abs(b.z_x)))
    b["p_y"] = -np.log10(2 * norm.sf(np.abs(b.z_y)))
    lead = b.loc[b.z_x.abs().idxmax()]

    fig, axes = plt.subplots(2, 1, figsize=(6.6, 4.6), sharex=True,
                             gridspec_kw={"height_ratios": [1, 1], "hspace": 0.12})
    axes[0].scatter(b.bp, b.p_x, s=7, color="#8c2d2d", alpha=0.65, lw=0)
    axes[0].axvline(lead.bp, color="#333333", ls="--", lw=0.8)
    axes[0].set_ylabel("Benign nodule\n$-\\log_{10}P$", fontsize=8.5)
    axes[0].set_title("chr9q34.2 (ABO region) · lead variant rs651007 · PP.H4 = 0.861", fontsize=9)

    axes[1].scatter(b.bp, b.p_y, s=7, color="#1f4e79", alpha=0.65, lw=0)
    axes[1].axvline(lead.bp, color="#333333", ls="--", lw=0.8)
    axes[1].set_ylabel("Migraine with aura\n$-\\log_{10}P$", fontsize=8.5)
    axes[1].set_xlabel("Position on chromosome 9 (Mb, GRCh38)", fontsize=8.5)

    for a in axes:
        a.grid(axis="y", ls=":", alpha=0.35)
    axes[1].annotate("rs651007", xy=(lead.bp, axes[1].get_ylim()[1] * 0.92),
                     xytext=(6, 0), textcoords="offset points", fontsize=7.5, color="#333333")
    plt.subplots_adjust(left=0.17, right=0.98, top=0.90, bottom=0.14)
    plt.savefig(os.path.join(RES, "fig4_regional.png"), dpi=300)
    plt.close()
    print("locus", locus, "SNPs:", len(b), "range %.2f-%.2f Mb" % (b.bp.min(), b.bp.max()),
          "lead bp %.3f" % lead.bp)


if __name__ == "__main__":
    fig1(); fig4()
    print("done")
