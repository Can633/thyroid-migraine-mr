"""
mrlap_compare_modes.py
----------------------
把两路 MRlap 结果与主 MR / F 方案合并成最终对照表 + 图:
  default = MRlap 原生 IV 选择 (p<5e-8, 500kb 距离 pruning, 无 LD clumping)
  userIV  = 主 MR (TwoSampleMR) 完全相同的 IV 集 (已 LD clumping + F>10 + Steiger)

运行: python mrlap_compare_modes.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = r"D:/GWAS/Migraine"
RES = os.path.join(ROOT, "results")

d = pd.read_csv(os.path.join(RES, "mrlap_results_default.csv"))
u = pd.read_csv(os.path.join(RES, "mrlap_results_userIV.csv"))
main = pd.read_csv(os.path.join(RES, "mr_results.csv"))
main = main[main.method == "Inverse variance weighted"][
    ["key", "nsnp", "b", "se", "pval"]].rename(
    columns={"key": "combo", "nsnp": "main_nSNP", "b": "main_b", "se": "main_se",
             "pval": "main_p"})
f = pd.read_csv(os.path.join(RES, "ldsc_overlap_corrected.csv"))[
    ["exposure", "outcome", "rho_s", "OR_corrected", "p_corrected"]].rename(
    columns={"OR_corrected": "F_OR_corrected", "p_corrected": "F_p_corrected"})

key = ["combo", "exposure", "outcome"]
dd = d[key + ["m_IVs", "observed_effect", "observed_effect_p",
              "corrected_effect", "corrected_effect_p", "p_difference",
              "int_crosstrait", "int_crosstrait_se", "h2_exp", "rg",
              "approx_OR", "approx_OR_corrected"]].add_prefix("def_")
dd = dd.rename(columns={"def_combo": "combo", "def_exposure": "exposure",
                        "def_outcome": "outcome"})
uu = u[key + ["m_IVs", "observed_effect", "observed_effect_p",
              "corrected_effect", "corrected_effect_p", "p_difference",
              "approx_OR", "approx_OR_corrected"]].add_prefix("usr_")
uu = uu.rename(columns={"usr_combo": "combo", "usr_exposure": "exposure",
                        "usr_outcome": "outcome"})
# 避免与 dd 的 exposure/outcome 在 merge 时生成 _x/_y 后缀
uu = uu.drop(columns=["exposure", "outcome"])

c = (main.merge(dd, on="combo").merge(uu, on="combo")
     .merge(f, on=["exposure", "outcome"]))
c["main_OR"] = np.exp(c.main_b)
c["Bonferroni_p"] = 0.05 / 15
c.to_csv(os.path.join(RES, "mrlap_all_modes_comparison.csv"), index=False)

print("=== 15 组合 p 值对照 ===")
print(c[["combo", "main_nSNP", "main_OR", "main_p",
         "def_m_IVs", "def_observed_effect_p", "def_corrected_effect_p", "def_p_difference",
         "usr_m_IVs", "usr_observed_effect_p", "usr_corrected_effect_p", "usr_p_difference",
         "F_p_corrected"]].to_string(index=False, float_format=lambda v: f"{v:.5g}"))

print("\n=== 名义 p<0.05 计数 ===")
for col in ["main_p", "def_observed_effect_p", "def_corrected_effect_p",
            "usr_observed_effect_p", "usr_corrected_effect_p", "F_p_corrected"]:
    print(f"  {col:<28s}: {(c[col] < 0.05).sum():>2d} / 15   (Bonferroni: {(c[col] < 0.05/15).sum()})")
print("\n=== p_difference < 0.05 (观测 vs 校正有显著差异) ===")
print("  default:", (c.def_p_difference < 0.05).sum(),
      " userIV:", (c.usr_p_difference < 0.05).sum())

# ---------------- 图: p 值对照 ----------------
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3})
lbl = (c.exposure + "->" + c.outcome.str.replace("Migraine", "Mig", regex=False)).tolist()
y = np.arange(15)[::-1]
fig, ax = plt.subplots(figsize=(8.6, 6.8))
for off, col, mark, pcol, name in [
    (-0.24, "main_p", "o", "#34495e", "Main MR (TwoSampleMR IVW)"),
    (0.00, "def_corrected_effect_p", "s", "#c0392b", "MRlap default (native IVs)"),
    (0.24, "usr_corrected_effect_p", "^", "#2471a3", "MRlap userIV = main MR IVs"),
]:
    ax.scatter(c[col], y + off, s=30, marker=mark, color=pcol,
               label=name, zorder=3)
ax.axvline(0.05, color="#e67e22", ls="--", lw=1.2, label="nominal p = 0.05")
ax.axvline(0.05/15, color="#7d3c98", ls=":", lw=1.2, label="Bonferroni p = 0.0033")
ax.set_xscale("log")
ax.set_xlim(1e-9, 1.2)
ax.set_yticks(y); ax.set_yticklabels(lbl, fontsize=8)
ax.set_xlabel("p-value (log scale)  [MRlap values = corrected effect]")
ax.set_title("Three MR pipelines: significance comparison (15 combinations)")
ax.legend(fontsize=8, loc="lower left")
fig.tight_layout()
fig.savefig(os.path.join(RES, "mrlap_modes_pvalues.png"), dpi=160)
plt.close(fig)

# ---------------- 图: 观测 -> 校正效应变化 ----------------
fig, ax = plt.subplots(figsize=(8.6, 6.8))
for off, pfx, col, name in [(-0.18, "def_", "#c0392b", "default"),
                            (0.18, "usr_", "#2471a3", "userIV")]:
    b = c[pfx + "corrected_effect"].to_numpy(float)
    # se 从各模式原始结果按 combo 对齐补入
    src = d if pfx == "def_" else u
    se = src.set_index("combo").loc[c.combo, "corrected_effect_se"].to_numpy(float)
    ax.errorbar(b, y + off, xerr=1.96 * se, fmt="o", ms=4, lw=1,
                color=col, label=f"MRlap corrected [{name}]", capsize=2)
ax.axvline(0, color="#444", lw=0.9)
ax.set_yticks(y); ax.set_yticklabels(lbl, fontsize=8)
ax.set_xlabel("Corrected standardised causal effect (SD/SD)")
ax.set_title("MRlap corrected effects: native IVs vs main-MR IVs")
ax.legend(fontsize=8, loc="lower right")
fig.tight_layout()
fig.savefig(os.path.join(RES, "mrlap_modes_forest.png"), dpi=160)
plt.close(fig)

print("\n输出: results/mrlap_all_modes_comparison.csv, "
      "results/mrlap_modes_pvalues.png, results/mrlap_modes_forest.png")
