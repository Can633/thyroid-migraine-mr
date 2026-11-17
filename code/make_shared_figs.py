# -*- coding: utf-8 -*-
"""Figures for the shared-genetic-architecture report."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "D:/GWAS/Migraine"
RES = os.path.join(ROOT, "results")

rg = pd.read_csv(os.path.join(RES, "ldsc_rg_precise.csv"))
att = pd.read_csv(os.path.join(RES, "ldsc_attn_ratio.csv"))
blk = pd.read_csv(os.path.join(RES, "coloc_blocks.tsv"), sep="\t")

OUT_LAB = {"Migraine": "Migraine (any)", "MigraineAura": "Migraine with aura",
           "MigraineNoAura": "Migraine without aura"}
COL = {"Migraine": "#c0392b", "MigraineAura": "#e67e22", "MigraineNoAura": "#8e44ad"}

# ------------------------------------------------------------------ fig 1: rg forest
r = rg[rg.exposure != "BNG_LOO"].copy()
r["label"] = r.apply(lambda x: f"{x.exposure} → {OUT_LAB[x.outcome]}", axis=1)
r = r.sort_values("rg")
y = np.arange(len(r))
fig, ax = plt.subplots(figsize=(7.6, 6.4))
for i, (_, row) in enumerate(r.iterrows()):
    ax.plot([row.rg_lo, row.rg_hi], [i, i], color=COL[row.outcome], lw=1.6, alpha=.9)
    ax.scatter(row.rg, i, color=COL[row.outcome], s=34, zorder=3,
               edgecolor="black" if row.rg_signif else "none", linewidth=.7)
ax.axvline(0, color="#333", lw=1)
ax.axvline(0.05, color="#999", ls=":", lw=1)
ax.set_yticks(y)
ax.set_yticklabels(r.label, fontsize=8.5)
ax.set_xlabel("Genetic correlation $r_g$ (95% CI, joint block jackknife)")
ax.set_title("Positive genetic correlation between thyroid disease and migraine\n"
             "(15 non-redundant pairs; filled edge = nominal significance)", fontsize=10)
for o, c in COL.items():
    ax.scatter([], [], color=c, label=OUT_LAB[o], s=30)
ax.legend(fontsize=8, loc="lower right", frameon=False)
ax.grid(axis="x", alpha=.25)
ax.set_xlim(-0.15, 0.42)
fig.tight_layout()
fig.savefig(os.path.join(RES, "fig_rg_forest.png"), dpi=170)
plt.close(fig)

# ------------------------------------------------------------------ fig 2: attenuation
fig, ax = plt.subplots(figsize=(6.6, 4.0))
a = att.sort_values("attenuation_ratio")
cols = ["#c0392b" if t in ("HYPO",) else "#2471a3" for t in a.trait]
ax.barh(range(len(a)), a.attenuation_ratio, color=cols, alpha=.85)
ax.axvline(0.5, color="#c0392b", ls="--", lw=1.2)
ax.text(0.51, 0.3, "conventional threshold 0.5\n(above = stratification likely)",
        fontsize=7.5, color="#c0392b")
ax.set_yticks(range(len(a)))
ax.set_yticklabels(a.trait, fontsize=8.5)
ax.set_xlabel("Attenuation ratio = (intercept − 1) / (mean $\\chi^2$ − 1)")
ax.set_title("All traits are below the stratification threshold\n(no stratified LDSC correction required)", fontsize=9.5)
ax.grid(axis="x", alpha=.25)
fig.tight_layout()
fig.savefig(os.path.join(RES, "fig_attenuation.png"), dpi=170)
plt.close(fig)

# ------------------------------------------------------------------ fig 3: coloc
inf = blk[(blk.informative) & (blk.status == "ok")].copy()
fig, ax = plt.subplots(figsize=(6.8, 5.0))
for lab, c, s in [("no_signal", "#bdc3c7", 8), ("distinct", "#7f8c8d", 14),
                  ("suggestive", "#e67e22", 22), ("shared", "#c0392b", 55)]:
    d = inf[inf.label == lab]
    if len(d):
        ax.scatter(d["PP.H3"], d["PP.H4"], s=s, color=c, alpha=.75,
                   label=f"{lab} (n={len(d)})", edgecolor="none")
ax.axhline(0.8, color="#c0392b", ls="--", lw=1)
ax.axhline(0.5, color="#e67e22", ls=":", lw=1)
ax.text(0.62, 0.81, "PP.H4 = 0.8", fontsize=7.5, color="#c0392b")
ax.set_xlabel("PP.H3  (two distinct causal variants)")
ax.set_ylabel("PP.H4  (shared causal variant)")
ax.set_title("Colocalisation of thyroid disease with migraine\n"
             f"({len(inf)} informative loci; PP.H0 > 0.9 regions excluded)", fontsize=9.5)
ax.legend(fontsize=8, frameon=False, loc="upper right")
ax.grid(alpha=.2)
fig.tight_layout()
fig.savefig(os.path.join(RES, "fig_coloc.png"), dpi=170)
plt.close(fig)

print("written:")
for f in ["fig_rg_forest.png", "fig_attenuation.png", "fig_coloc.png"]:
    p = os.path.join(RES, f)
    print(" ", f, os.path.getsize(p), "bytes")

print("\n=== numbers for the report ===")
print("rg positive: %d/%d" % ((rg.rg > 0).sum(), len(rg)))
print("nominal: %d  Bonferroni(|z|>2.96): %d" % (rg.rg_signif.sum(), (rg.z_rg.abs() > 2.96).sum()))
print("rg range %+.3f to %+.3f ; median %+.3f" % (rg.rg.min(), rg.rg.max(), rg.rg.median()))
print("\nBonferroni hits:")
print(rg[rg.z_rg.abs() > 2.96][["exposure", "outcome", "rg", "se_rg", "z_rg"]]
      .to_string(index=False, float_format=lambda x: f"{x:.3f}"))
print("\nattenuation ratio range %.3f - %.3f" % (att.attenuation_ratio.min(), att.attenuation_ratio.max()))
print("\ncoloc: informative %d/%d ; shared %d ; suggestive %d ; distinct %d" %
      (len(inf), len(blk), (blk.label == "shared").sum(),
       (blk.label == "suggestive").sum(), (blk.label == "distinct").sum()))
