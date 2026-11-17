"""
mrlap_summary.py
----------------
汇总 run_mrlap.R 的 15 组合结果，并与
  * 主 MR (TwoSampleMR IVW, log OR / OR)
  * E 方案 (Burgess 解析偏倚网格)
  * F 方案 (纯 Python cross-trait LDSC 实测校正)
做对比，输出 CSV + 图。

运行: python mrlap_summary.py [default|userIV]
"""
import os
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm

ROOT = r"D:/GWAS/Migraine"
RES = os.path.join(ROOT, "results")
RUN = os.path.join(ROOT, "mrlap_run")
MODE = sys.argv[1] if len(sys.argv) > 1 else "default"

# ---------- 1. 汇总 MRlap ----------
files = sorted(glob.glob(os.path.join(RUN, MODE, "*", "row.csv")))
m = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
m = m.sort_values(["exposure", "outcome"]).reset_index(drop=True)
m.to_csv(os.path.join(RES, f"mrlap_results_{MODE}.csv"), index=False)
n_fail = int(((~m.error.isna()) & (m.error.astype(str).str.strip() != "")).sum())
print(f"[{MODE}] MRlap 组合数: {len(m)}  失败: {n_fail}")
print("  IVs:", m.m_IVs.min(), "-", m.m_IVs.max(),
      " | 校正后 p<0.05:", (m.corrected_effect_p < 0.05).sum(),
      " | p_difference<0.05:", (m.p_difference < 0.05).sum())

# ---------- 2. 主 MR / E / F ----------
main = pd.read_csv(os.path.join(RES, "mr_results.csv"))
main = main[main.method == "Inverse variance weighted"][["key", "nsnp", "b", "se", "pval"]]
main = main.rename(columns={"key": "combo", "nsnp": "main_nSNP", "b": "main_b_logOR",
                           "se": "main_se", "pval": "main_p"})
main["main_OR"] = np.exp(main.main_b_logOR)

f = pd.read_csv(os.path.join(RES, "ldsc_overlap_corrected.csv"))[
    ["exposure", "outcome", "rho_s", "OR", "OR_corrected", "p_ivw", "p_corrected"]]
f = f.rename(columns={"OR": "F_OR", "OR_corrected": "F_OR_corrected",
                      "p_ivw": "F_p", "p_corrected": "F_p_corrected"})

e = pd.read_csv(os.path.join(RES, "overlap_corrected_summary.csv"))[
    ["exposure", "outcome", "OR_corrected_max", "p_corrected_max"]]
e = e.rename(columns={"OR_corrected_max": "E_OR_extreme",
                      "p_corrected_max": "E_p_extreme"})

c = (m.merge(main, on="combo", how="left")
      .merge(f, on=["exposure", "outcome"], how="left")
      .merge(e, on=["exposure", "outcome"], how="left"))
c["main_OR_lo"] = np.exp(c.main_b_logOR - 1.96 * c.main_se)
c["main_OR_hi"] = np.exp(c.main_b_logOR + 1.96 * c.main_se)
c["mrlap_OR_approx"] = c.approx_OR
c["mrlap_OR_corr_approx"] = c.approx_OR_corrected
c.to_csv(os.path.join(RES, f"mrlap_vs_main_vs_EF_{MODE}.csv"), index=False)

# ---------- 3. 图 ----------
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3})
lbl = (c.exposure + "->" + c.outcome.str.replace("Migraine", "Mig", regex=False)).tolist()
y = np.arange(len(c))[::-1]

# (a) 观测 vs 校正 (标准化尺度)
fig, ax = plt.subplots(figsize=(9.2, 7.0))
for off, col, pre in [(-0.20, "#2c7fb8", "observed"), (0.20, "#d95f0e", "corrected")]:
    b = c[f"{pre}_effect"].to_numpy(float)
    se = c[f"{pre}_effect_se"].to_numpy(float)
    ax.errorbar(b, y + off, xerr=1.96 * se, fmt="o", ms=4.5, lw=1.1,
                color=col, label=f"MRlap {pre} (SD units)", elinewidth=1.1, capsize=2)
ax.axvline(0, color="#444", lw=0.9)
ax.set_yticks(y); ax.set_yticklabels(lbl, fontsize=8)
ax.set_xlabel("Standardised causal effect (SD of outcome per SD of exposure)")
ax.set_title(f"MRlap: observed vs corrected effect  [{MODE}]")
ax.legend(fontsize=8, loc="lower right")
fig.tight_layout()
fig.savefig(os.path.join(RES, f"mrlap_forest_{MODE}.png"), dpi=160)
plt.close(fig)

# (b) 主 MR (OR) vs MRlap 近似 OR —— 量级/方向对照
fig, ax = plt.subplots(figsize=(6.2, 5.2))
x = c.main_b_logOR.to_numpy(float)
yo = c.approx_b_logOR.to_numpy(float)
yc = c.approx_b_corrected_logOR.to_numpy(float)
ax.scatter(x, yo, s=28, color="#2c7fb8", label="MRlap observed (approx log OR)")
ax.scatter(x, yc, s=28, color="#d95f0e", marker="^", label="MRlap corrected (approx log OR)")
lim = [min(x.min(), yo.min(), yc.min()) - 0.01, max(x.max(), yo.max(), yc.max()) + 0.01]
ax.plot(lim, lim, ls="--", color="#888", lw=1, label="y = x")
ax.axhline(0, color="#bbb", lw=0.8); ax.axvline(0, color="#bbb", lw=0.8)
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("Main MR: IVW log OR (TwoSampleMR)")
ax.set_ylabel("MRlap: approximated log OR (SD -> log OR)")
r = np.corrcoef(x, yo)[0, 1]
ax.text(0.03, 0.95, f"Pearson r = {r:.3f}\n(approximate scale conversion,\nfor magnitude comparison only)",
        transform=ax.transAxes, va="top", fontsize=8, color="#444")
ax.set_title(f"Direction / magnitude cross-check [{MODE}]")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(RES, f"mrlap_vs_main_scatter_{MODE}.png"), dpi=160)
plt.close(fig)

# ---------- 4. 文字摘要 ----------
sig = c[(c.corrected_effect_p < 0.05) | (c.main_p < 0.05)]
cols = ["combo", "m_IVs", "main_OR", "main_p",
        "observed_effect", "observed_effect_p", "corrected_effect",
        "corrected_effect_p", "p_difference", "int_crosstrait",
        "int_crosstrait_se", "rho_s", "F_OR_corrected"]
print("\n=== 关键组合 (任一 p<0.05) ===")
print(sig[cols].to_string(index=False, float_format=lambda v: f"{v:.5g}"))
print("\n=== 全部 15 组合: 主MR / MRlap / F 校正 p 值 ===")
print(c[["combo", "m_IVs", "main_OR", "main_p", "observed_effect_p",
         "corrected_effect_p", "p_difference", "F_p_corrected"]]
      .to_string(index=False, float_format=lambda v: f"{v:.5g}"))
print("\n输出:", f"results/mrlap_results_{MODE}.csv",
      f"results/mrlap_vs_main_vs_EF_{MODE}.csv",
      f"results/mrlap_forest_{MODE}.png",
      f"results/mrlap_vs_main_scatter_{MODE}.png")
