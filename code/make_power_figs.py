# -*- coding: utf-8 -*-
"""make_power_figs.py -- 功效曲线与“能否复现 Li 2024”对照图"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = r"D:/GWAS/Migraine/results"
D = pd.read_csv(os.path.join(RES, "power_design.csv"))
C = pd.read_csv(os.path.join(RES, "power_curves.csv"))

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#444444", "axes.linewidth": 0.8,
    "figure.dpi": 300,
})

COL = {"Graves' disease": "#8c2d2d", "Hypothyroidism": "#1f4e79",
       "Lymphocytic thyroiditis": "#c07a1e", "Benign thyroid nodule": "#4a7a4a",
       "Thyroid cancer": "#6a4a7a",
       "Migraine (overall)": "#2c3e50", "Migraine with aura": "#c0392b",
       "Migraine without aura": "#2980b9"}


def panel(ax, sub, curves, title, xmax, show_target=True):
    for key, c in sub.iterrows():
        g = curves[(curves["alpha"] == "nominal") & (curves["exposure"] == c["exposure"])
                   & (curves["outcome"] == c["outcome"])]
        ax.plot(g["OR"], g["power"], color=COL.get(c["exp_lab"], "grey"), lw=1.6,
                label=c["exp_lab"])
    ax.axhline(0.8, color="#888888", ls=":", lw=1)
    ax.text(xmax * 0.995, 0.815, "80%", ha="right", va="bottom", fontsize=7.5, color="#666666")
    if show_target:
        ax.axvline(1.066, color="#b03a2e", ls="--", lw=1)
        ax.text(1.066, 0.03, " Li 2024\n OR 1.066", fontsize=7.2, color="#b03a2e", va="bottom")
    ax.set_xlim(0.88, xmax)
    ax.set_ylim(0, 1.02)
    ax.set_xscale("log")
    xt = [0.9, 1.0, 1.1, 1.2, 1.3, 1.4] if xmax <= 1.45 else [0.9, 1.0, 1.1, 1.25, 1.5, 1.8, 2.2]
    xt = [t for t in xt if 0.88 < t <= xmax]
    ax.set_xticks(xt)
    ax.set_xticklabels(["%.2f" % t for t in xt])
    ax.minorticks_off()
    ax.set_xlabel("True causal OR (per 1-unit increase in exposure log-odds)", fontsize=8.5)
    ax.set_ylabel("Power (two-sided $\\alpha$ = 0.05)", fontsize=8.5)
    ax.set_title(title, fontsize=9.5, pad=6)
    ax.grid(axis="y", ls=":", alpha=0.35)


# ---------------- 图 1：功效曲线 ----------------
fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.6))
fwd = D[(D.direction == "forward") & (D.outcome == "Migraine")].reset_index()
rev = D[(D.direction == "reverse") & (D.outcome == "HYPO")].reset_index()

panel(axes[0], fwd, C, "A  Thyroid disease $\\rightarrow$ migraine", 1.30)
panel(axes[1], rev, C, "B  Migraine subtype $\\rightarrow$ hypothyroidism", 1.30, show_target=False)
for ax in axes:
    ax.legend(fontsize=7.0, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.22), ncol=2, handlelength=1.6,
              columnspacing=1.2, labelspacing=0.35)
plt.subplots_adjust(left=0.075, right=0.985, top=0.90, bottom=0.28, wspace=0.24)
plt.savefig(os.path.join(RES, "fig_power_curves.png"), dpi=300)
plt.close()

# ---------------- 图 2：能否复现 OR = 1.066 ----------------
d = D.copy()
d["lab"] = np.where(d.direction == "forward",
                    d.exp_lab.str.replace("Graves' disease", "Graves").str.replace("Lymphocytic thyroiditis", "Lymph. thyroiditis")
                     .str.replace("Benign thyroid nodule", "Benign nodule") + "  $\\rightarrow$ " + d.out_lab.str.replace("Migraine", "Mig."),
                    d.exp_lab.str.replace("Migraine (overall)", "Mig.").str.replace("Migraine with aura", "Mig.+A")
                     .str.replace("Migraine without aura", "Mig.-A") + "  $\\rightarrow$ " + d.out_lab.str.replace("Graves' disease", "Graves")
                     .str.replace("Lymphocytic thyroiditis", "Lymph. thyroiditis").str.replace("Benign thyroid nodule", "Benign nodule"))
d = d.sort_values(["direction", "pw_1.066"], ascending=[True, False]).reset_index(drop=True)

fig, ax = plt.subplots(figsize=(7.0, 7.2))
yy = 0
ticks, labels, seps = [], [], []
for direction, sub in d.groupby("direction", sort=True):
    for _, r in sub.iterrows():
        c = "#8c2d2d" if direction == "forward" else "#3d6a94"
        ax.barh(yy, r["pw_1.066"], height=0.62, color=c, alpha=0.9)
        ax.text(r["pw_1.066"] + 0.012, yy, "%.2f" % r["pw_1.066"], va="center", fontsize=7)
        ticks.append(yy); labels.append(r["lab"]); yy += 1
    seps.append(yy - 0.5); yy += 1.0

ax.axvline(0.8, color="#888888", ls=":", lw=1)
ax.text(0.79, -0.9, "80%", fontsize=7.5, color="#666666", ha="right")
for s in seps[:-1]:
    ax.axhline(s, color="#dddddd", lw=0.8)
ax.set_yticks(ticks); ax.set_yticklabels(labels, fontsize=7.4)
ax.invert_yaxis(); ax.set_xlim(0, 1.12)
ax.set_xlabel("Power to detect the effect size reported by Li et al. 2024 (OR 1.066)", fontsize=8.5)
ax.set_title("Power to detect OR 1.066 (Li et al. 2024)", fontsize=9.5)
plt.subplots_adjust(left=0.35, right=0.97, top=0.94, bottom=0.09)
plt.savefig(os.path.join(RES, "fig_power_li066.png"), dpi=300)
plt.close()
print("wrote fig_power_curves.png, fig_power_li066.png")
