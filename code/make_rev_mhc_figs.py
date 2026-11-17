# -*- coding: utf-8 -*-
"""
make_rev_mhc_figs.py —— 反向 MR 森林图 + MHC 敏感性图
=====================================================
fig_rev_mr.png   15 个反向组合（IVW）的 OR 森林图，按偏头痛亚型分组
fig_mhc.png      (a) 各性状全基因组显著信号落在 MHC 的比例
                 (b) r_g 在"不含 MHC / 含 MHC"两种设定下的点估计与 95% CI
标签用英文，避免 matplotlib 缺中文字体导致乱码。
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = "D:/GWAS/Migraine"
RES = os.path.join(ROOT, "results")

LAB = {"GD": "Graves' disease", "HYPO": "Hypothyroidism", "LT": "Lymphocytic thyroiditis",
       "BNG": "Benign nodular goitre", "TC": "Thyroid cancer"}
ORD = ["GD", "HYPO", "LT", "BNG", "TC"]


def fig_rev():
    d = pd.read_csv(os.path.join(RES, "rev_mr_main.csv"))
    d = d[d["method"] == "Inverse variance weighted"]
    groups = ["Migraine", "MigraineAura", "MigraineNoAura"]
    gname = {"Migraine": "Migraine (overall)", "MigraineAura": "Migraine with aura",
             "MigraineNoAura": "Migraine without aura"}
    colors = {"Migraine": "#2c3e50", "MigraineAura": "#c0392b", "MigraineNoAura": "#2980b9"}

    y, ticks, labels, seps = 0, [], [], []
    fig, ax = plt.subplots(figsize=(7.4, 5.4))
    for g in groups:
        sub = d[d["id.exposure"] == g].set_index("id.outcome").reindex(ORD).dropna(subset=["or"])
        for th in sub.index:
            r = sub.loc[th]
            ax.plot([r["or_lci95"], r["or_uci95"]], [y, y], color=colors[g], lw=1.5, alpha=.8,
                    solid_capstyle="round")
            ax.plot(r["or"], y, "o", color=colors[g], ms=6, label=gname[g] if th == ORD[0] else None)
            ticks.append(y)
            labels.append(f"{LAB[th]}  (n={int(r['nsnp'])})")
            y += 1
        seps.append(y - 0.5)
        y += 0.8

    ax.axvline(1, color="grey", ls="--", lw=1)
    for sy in seps[:-1]:
        ax.axhline(sy, color="#dddddd", lw=0.8)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(0.72, 3.4)
    ax.set_xticks([0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0])
    ax.set_xticklabels(["0.8", "1.0", "1.25", "1.5", "2.0", "2.5", "3.0"])
    ax.set_xlabel("OR for thyroid disease per genetically predicted migraine", fontsize=9)
    ax.set_title("Reverse MR: migraine \u2192 thyroid disease (IVW)", fontsize=10.5)
    ax.grid(axis="x", ls=":", alpha=.4)
    ax.legend(fontsize=8, loc="lower right", frameon=False, title="exposure")
    plt.subplots_adjust(left=0.30, right=0.97, top=0.92, bottom=0.10)
    plt.savefig(os.path.join(RES, "fig_rev_mr.png"), dpi=300)
    plt.close()


def fig_mhc():
    c = pd.read_csv(os.path.join(RES, "mhc_gws_counts.csv"))
    c = c.set_index("trait")
    order = ["GD", "HYPO", "LT", "BNG", "TC", "Migraine", "MigraineAura", "MigraineNoAura"]
    c = c.reindex(order)

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4),
                             gridspec_kw={"width_ratios": [1.05, 1]})
    ax = axes[0]
    vals = (c["frac_mhc"] * 100).values
    cols = ["#c0392b" if t in ("GD", "HYPO", "LT") else "#7f8c8d" for t in order]
    ax.bar(range(len(order)), vals, color=cols)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(["GD", "HYPO", "LT", "BNG", "TC", "Mig.", "Mig.A", "Mig.nA"],
                       fontsize=9)
    ax.set_ylabel("% of GWS SNPs in MHC", fontsize=9)
    ax.set_title("(a) MHC share of genome-wide signal", fontsize=10)
    for i, v in enumerate(vals):
        ax.text(i, v + 1.5, f"{v:.0f}%", ha="center", fontsize=8)
    ax.set_ylim(0, 105)
    ax.grid(axis="y", ls=":", alpha=.4)

    ax = axes[1]
    s = pd.read_csv(os.path.join(RES, "mhc_rg_sensitivity.csv"))
    a = s[s.variant == "A_noMHC_std"].set_index("exposure")
    b = s[s.variant == "C_MHC_selfW"].set_index("exposure")
    yy = np.arange(len(ORD))
    for i, e in enumerate(ORD):
        if e not in a.index or e not in b.index:
            continue
        ax.plot([a.loc[e, "rg"] - 1.96 * a.loc[e, "se"], a.loc[e, "rg"] + 1.96 * a.loc[e, "se"]],
                [yy[i] - .12, yy[i] - .12], color="#2c3e50", lw=1.5)
        ax.plot(a.loc[e, "rg"], yy[i] - .12, "o", color="#2c3e50", ms=6, label="without MHC" if i == 0 else "")
        ax.plot([b.loc[e, "rg"] - 1.96 * b.loc[e, "se"], b.loc[e, "rg"] + 1.96 * b.loc[e, "se"]],
                [yy[i] + .12, yy[i] + .12], color="#c0392b", lw=1.5)
        ax.plot(b.loc[e, "rg"], yy[i] + .12, "s", color="#c0392b", ms=5.5, label="with MHC" if i == 0 else "")
    ax.axvline(0, color="grey", ls="--", lw=1)
    ax.set_yticks(yy); ax.set_yticklabels([LAB[e] for e in ORD], fontsize=9)
    ax.set_xlabel(r"$r_g$ with migraine (95% CI)", fontsize=9)
    ax.set_title("(b) $r_g$ is stable when MHC is added back", fontsize=10)
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, frameon=False)
    ax.grid(axis="x", ls=":", alpha=.4)
    plt.subplots_adjust(left=0.19, right=0.98, top=0.91, bottom=0.24, wspace=0.55)
    plt.savefig(os.path.join(RES, "fig_mhc.png"), dpi=300)
    plt.close()


if __name__ == "__main__":
    fig_rev()
    fig_mhc()
    print("saved fig_rev_mr.png, fig_mhc.png")
