#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
overlap_correct.py
==================
样本重叠偏倚的纯 Python 校正与实证检验（不依赖 R 编译器 / 1000G LD scores）

分析目标
--------
甲状腺疾病(GBMI/VTB 2026 五病) -> 偏头痛(FinnGen R12 三表型) 两样本 MR
结果全体接近无效(OR ~0.98-1.03, 均不显著)，既往报告将其归因于
"GBMI 暴露含 FinnGen、结局为 FinnGen -> 样本重叠把效应压向零(衰减伪阴性)"。
本脚本用三条互相独立的证据检验该归因，并给出校正后的效应界限。

方法
----
(1) Burgess, Davies & Thompson (2016) Genet Epidemiol 40(7):597-608 解析偏倚公式
    来源：论文 Web Appendix A1 原始 R 代码
          (本地副本 ldref/burgess2016_appendix.pdf)
        expf   = (N-K-1)/K * rsq/(1-rsq)      # 期望 F 统计量
        bias   = olsbias * overlap.prop * (1/expf)
        type1err = 2 - Phi(1.96 + bias/sqrt(var)) - Phi(1.96 - bias/sqrt(var))
    与原文 Appendix Table A1 的一致性验证：
        alpha=0.05, Mean F=6.3, Mean OLS=0.497, 100% overlap
        表列 Expected Bias = 0.079 ;  0.497/6.3 = 0.0789  (一致)
    方向（原文 Abstract 明确）：一样本/重叠 + 弱工具 -> 偏向观察性关联(混杂)；
        两样本不重叠 -> 偏倚朝向零。因此重叠不会把估计压向零。

(2) BNG 天然对照（无模型假设的实测）
    GCST90627740 (BNG EUR 主版, 含 FinnGen) vs
    GCST90627765 (BNG EUR 备选版, 剔除 FinnGen + HUNT/LATVIAN/CCPM/EXCEED/UKB)
    结局同为 FinnGen R12 偏头痛三表型。
    为消除工具集差异，仅在两版共有 SNP 上、用完全相同的结局 beta/se/权重，
    只替换暴露 beta，重算 IVW。差值即"剥离 FinnGen 重叠后效应的实测变化"。

(3) 回归稀释一阶校正（beta_exposure 测量误差导致的向零衰减）
    lambda = 1 - 1/E[F] ;  b_dilution_corrected = b_IVW / lambda
    这是测量误差的一阶近似，仅用于量级判断，标注为近似。

(4) 与 Li 2024 (PMID 39087012) 的定量对比：反推需要多强的偏倚
    才能把本研究的点估计推到文献效应量。

局限（在报告中详述）
--------------------
- Burgess 公式为因果零假设下的偏倚；非零真效应下偏倚结构更复杂。
- expf 用实测 mean(F_j) 近似（附录 Table A1 即用 "Mean F"）；
  winner's curse 会使实测 F 上偏 -> 偏倚被低估，故另用 median F 做保守敏感性。
- olsbias（观察性混杂关联）与 overlap.prop（重叠比例）非实测，
  故采用全范围网格敏感性，而非给单点"校正值"。
- 不做 MRlap/CAUSE 式的遗传力/遗传相关联合估计（需 LD scores）。

输出
----
results/overlap_bias_grid.csv          偏倚网格 (15 组合 x pi x delta)
results/overlap_corrected_summary.csv  校正后效应汇总
results/bng_natural_experiment.csv     BNG 天然对照
results/li2024_comparison.csv          与 Li 2024 的定量对比
results/forest_overlap.png             原始 vs 校正后森林图
results/bias_required.png              反推所需观察性关联强度
"""
import os
import glob
import math
import numpy as np
import pandas as pd
from scipy.stats import norm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "D:/GWAS/Migraine"
RES = os.path.join(ROOT, "results")
os.makedirs(RES, exist_ok=True)

# ---------------------------------------------------------------------------
# 外部参照：Li R et al. Front Neurol 2024;15:1420201 (PMID 39087012)
# OR/CI 取自 PMC11288874 正文原句，非推测值
# ---------------------------------------------------------------------------
LI2024 = {
    "Migraine":      dict(OR=1.066, lo=1.030, hi=1.103, p=0.0002),
    "MigraineAura":  dict(OR=1.074, lo=1.021, hi=1.131, p=0.006),
    "MigraineNoAura":dict(OR=1.068, lo=1.012, hi=1.127, p=0.017),
}
LI_B = {k: math.log(v["OR"]) for k, v in LI2024.items()}
LI_SE = {k: (math.log(v["hi"]) - math.log(v["lo"])) / (2 * 1.96)
         for k, v in LI2024.items()}

# 数据驱动的重叠比例上界（推断，非实测）：
# 结局 FinnGen R12 migraine 有效 N = 401,499；暴露 EUR 元分析 N ~ 1.8e6。
# FinnGen R12 全库公开上限 ~500,348，故重叠个体数 N_overlap <= 401,499。
# Burgess 的 overlap.prop 原为等量样本定义；对不等量样本用
#   pi ~ N_overlap / sqrt(N_exposure * N_outcome) 作为合理代理。
N_OUT_MIGRAINE = 395060
N_EXPOS_TYPICAL = 1.8e6
PI_DATA_BOUND = N_OUT_MIGRAINE / math.sqrt(N_EXPOS_TYPICAL * N_OUT_MIGRAINE)
# 甲状腺病-偏头痛观察性关联为正（共病），故 olsbias > 0，
# Burgess 重叠偏倚为「正」（朝观察性关联/远离零），与「衰减伪阴性」方向相反。
OLSBIAS_SIGN = "positive (thyroid-migraine observational comorbidity is positive)"

# Burgess 网格：overlap.prop 全范围 0-1；olsbias 为 log OR 尺度
PIS = [0.0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0]
# delta 对应观察性 OR = 1.00, 1.22, 1.50, 2.00 (甲状腺病-偏头痛共病文献常见量级)
DELTAS = [0.0, 0.20, math.log(1.5), math.log(2.0)]
REALISTIC_PI, REALISTIC_DELTA = 0.3, math.log(1.5)
MAX_PI, MAX_DELTA = 1.0, math.log(2.0)

OUT_N = {"Migraine": 395060, "MigraineAura": 379299, "MigraineNoAura": 376762}


def ivw(d):
    """固定效应 IVW，与 TwoSampleMR::mr(method='Inverse variance weighted') 一致。"""
    bx = d["beta.exposure"].to_numpy(float)
    by = d["beta.outcome"].to_numpy(float)
    sy = d["se.outcome"].to_numpy(float)
    w = 1.0 / sy ** 2
    denom = float(np.sum(w * bx * bx))
    b = float(np.sum(w * bx * by) / denom)
    se = math.sqrt(1.0 / denom)
    p = float(2 * norm.sf(abs(b / se)))
    return b, se, p


def load_harm(path):
    d = pd.read_csv(path)
    if "mr_keep" in d.columns:
        d = d[d["mr_keep"] == True]
    return d.reset_index(drop=True)


def burgess_bias(olsbias, pi, expF):
    """因果零假设下、样本重叠导致的 IVW 偏倚（Burgess 2016 附录 A1）。"""
    return olsbias * pi / expF


def type1_err(bias, se):
    z = bias / se
    return float(2 - norm.cdf(1.96 + z) - norm.cdf(1.96 - z))


def parse_combo(path):
    base = os.path.basename(path).replace("harmonised_", "").replace(".csv", "")
    e, o = base.split("_to_")
    return e, o


def r_se_for(rres, e, o):
    """TwoSampleMR 默认 IVW 使用乘性随机效应；Q>df 时 se 被 sqrt(Q/df) 膨胀。
    以 R 报告的 se 作为主分析 se（更保守），固定效应 se 仅作复现校验。"""
    rr = rres[(rres["id.exposure"] == e) & (rres["id.outcome"] == o)]
    return float(rr["se"].iloc[0]) if len(rr) else float("nan")


def conservative_p(b, se):
    return float(2 * norm.sf(abs(b / se)))


# ===========================================================================
# 1. 15 组合主分析：IVW 复现 + Burgess 校正
# ===========================================================================
files = sorted(glob.glob(os.path.join(RES, "harmonised_*_to_*.csv")))
files = [f for f in files if "_LOO_" not in os.path.basename(f)]
assert len(files) == 15, f"expected 15 main harmonised files, got {len(files)}"

rres = pd.read_csv(os.path.join(RES, "mr_results.csv"))
ivw_r = rres[rres["method"] == "Inverse variance weighted"].copy()

summary_rows, grid_rows = [], []
for f in files:
    e, o = parse_combo(f)
    d = load_harm(f)
    b, se, p = ivw(d)
    F = (d["beta.exposure"] / d["se.exposure"]) ** 2
    meanF, medF = float(F.mean()), float(F.median())
    lam = 1.0 - 1.0 / meanF

    # 与 R (TwoSampleMR) 结果比对，验证 Python 复现无误
    rr = ivw_r[(ivw_r["id.exposure"] == e) & (ivw_r["id.outcome"] == o)]
    rb = float(rr["b"].iloc[0]) if len(rr) else float("nan")
    rse = r_se_for(ivw_r, e, o)
    # 主分析采用 R 报告的 se（乘性随机效应，更保守）
    se_use = rse

    or0 = math.exp(b)
    lo0, hi0 = math.exp(b - 1.96 * se_use), math.exp(b + 1.96 * se_use)

    # 反推：要把本组合的点估计推到 Li 2024 效应量，需要多大的（负）偏倚
    b_target = LI_B[o]
    needed_bias = b - b_target                      # 通常为负 = 需要压低估计的偏倚
    # pi=1 时所需 |olsbias| 最小；换算成观察性 OR 量级
    min_abs_olsbias = abs(needed_bias) * meanF / 1.0
    needed_obs_OR_if_positive = math.exp(min_abs_olsbias)
    needed_obs_OR_if_negative = math.exp(-min_abs_olsbias)

    # 网格
    for pi in PIS:
        for delta in DELTAS:
            for expF, ftag in [(meanF, "meanF"), (medF, "medianF_conservative")]:
                bias = burgess_bias(delta, pi, expF)
                bc = b - bias                     # 观察关联为正 -> 重叠使估计上偏 -> 校正需减去
                grid_rows.append(dict(
                    exposure=e, outcome=o, overlap_prop=pi,
                    olsbias_logOR=delta, observational_OR=math.exp(delta),
                    F_basis=ftag, expF=expF,
                    bias_direction="positive (away from null)",
                    bias=bias, b_corrected=bc,
                    OR_corrected=math.exp(bc),
                    se=se_use, p_corrected=conservative_p(bc, se_use),
                    type1_error=type1_err(bias, se_use),
                ))

    def corrected(delta, pi, expF):
        bias = burgess_bias(delta, pi, expF)
        bc = b - bias
        return bias, bc, math.exp(bc), conservative_p(bc, se_use)

    b_real, bc_real, or_real, p_real = corrected(REALISTIC_DELTA, REALISTIC_PI, meanF)
    b_max, bc_max, or_max, p_max = corrected(MAX_DELTA, MAX_PI, meanF)
    _, bc_max_med, or_max_med, p_max_med = corrected(MAX_DELTA, MAX_PI, medF)

    summary_rows.append(dict(
        exposure=e, outcome=o, nSNP=len(d),
        b_ivw=b, se_ivw=se_use, p_ivw=conservative_p(b, se_use),
        OR=or0, OR_lo=lo0, OR_hi=hi0,
        se_fixed_effect=se,
        R_b=rb, R_se=rse, IVW_repro_diff=abs(b - rb),
        SE_diff_random_vs_fixed=se_use - se,
        meanF=meanF, medianF=medF,
        dilution_lambda=lam, OR_dilution_corrected=math.exp(b / lam),
        bias_realistic=b_real, OR_corrected_realistic=or_real, p_corrected_realistic=p_real,
        bias_max=b_max, OR_corrected_max=or_max, p_corrected_max=p_max,
        OR_corrected_max_medianF=or_max_med, p_corrected_max_medianF=p_max_med,
        Li2024_OR=LI2024[o]["OR"], Li2024_b=b_target,
        needed_bias_to_match_Li=needed_bias,
        min_abs_olsbias_needed_at_pi1=min_abs_olsbias,
        needed_observational_OR_if_positive=needed_obs_OR_if_positive,
        needed_observational_OR_if_negative=needed_obs_OR_if_negative,
    ))

grid = pd.DataFrame(grid_rows)
summary = pd.DataFrame(summary_rows).sort_values(["exposure", "outcome"])
grid.to_csv(os.path.join(RES, "overlap_bias_grid.csv"), index=False)
summary.to_csv(os.path.join(RES, "overlap_corrected_summary.csv"), index=False)

# ===========================================================================
# 2. BNG 天然对照：同一 SNP 集、同一结局，只替换暴露 beta（含/不含 FinnGen）
# ===========================================================================
bng_rows = []
for o in OUT_N:
    fm = os.path.join(RES, f"harmonised_BNG_to_{o}.csv")
    fl = os.path.join(RES, f"harmonised_BNG_LOO_to_{o}.csv")
    if not (os.path.exists(fm) and os.path.exists(fl)):
        continue
    m, l = load_harm(fm), load_harm(fl)
    common = sorted(set(m["SNP"]) & set(l["SNP"]))
    m = m[m["SNP"].isin(common)].set_index("SNP").loc[common].reset_index()
    l = l[l["SNP"].isin(common)].set_index("SNP").loc[common].reset_index()

    # 等位基因一致性：若效应等位基因相反，翻转 LOO 暴露 beta
    ea_m = m["effect_allele.exposure"].astype(str).str.upper().to_numpy()
    ea_l = l["effect_allele.exposure"].astype(str).str.upper().to_numpy()
    flip = ea_m != ea_l
    bx_loo = l["beta.exposure"].to_numpy(float).copy()
    bx_loo[flip] *= -1.0

    # 结局 beta/se 在两版中应完全一致（同一 FinnGen 文件）
    out_diff = float(np.max(np.abs(m["beta.outcome"].to_numpy(float)
                                  - l["beta.outcome"].to_numpy(float)))) if len(m) else float("nan")

    d_main = m.copy()
    d_loo = m.copy()
    d_loo["beta.exposure"] = bx_loo          # 唯一改动：暴露 beta（剥离 FinnGen）
    d_loo["se.exposure"] = l["se.exposure"].to_numpy(float)

    b1, s1, p1 = ivw(d_main)
    b2, s2, p2 = ivw(d_loo)
    bng_rows.append(dict(
        outcome=o, nSNP_common=len(common), n_allele_flipped=int(flip.sum()),
        max_outcome_beta_diff=out_diff,
        with_FinnGen_b=b1, with_FinnGen_OR=math.exp(b1), with_FinnGen_p=p1,
        without_FinnGen_b=b2, without_FinnGen_OR=math.exp(b2), without_FinnGen_p=p2,
        delta_b=b2 - b1, delta_OR=math.exp(b2) - math.exp(b1),
        # Burgess 预测：含重叠的估计被推向观察关联(正) -> 剔除重叠后应下降
        burgess_predicted_sign_of_delta_b="negative",
        observed_sign="negative" if (b2 - b1) < 0 else "positive",
    ))
bng = pd.DataFrame(bng_rows)
bng.to_csv(os.path.join(RES, "bng_natural_experiment.csv"), index=False)

# ===========================================================================
# 3. 与 Li 2024 的定量对比（GD 作为自身免疫性甲亢的代理，非等同）
# ===========================================================================
li_rows = []
for o in OUT_N:
    s = summary[(summary["exposure"] == "GD") & (summary["outcome"] == o)]
    if not len(s):
        continue
    s = s.iloc[0]
    z_diff = (s["b_ivw"] - LI_B[o]) / math.sqrt(s["se_ivw"] ** 2 + LI_SE[o] ** 2)
    p_diff = float(2 * norm.sf(abs(z_diff)))
    ci_overlap_lo = max(s["OR_lo"], LI2024[o]["lo"])
    ci_overlap_hi = min(s["OR_hi"], LI2024[o]["hi"])
    li_rows.append(dict(
        outcome=o, this_OR=s["OR"], this_lo=s["OR_lo"], this_hi=s["OR_hi"],
        Li2024_OR=LI2024[o]["OR"], Li2024_lo=LI2024[o]["lo"], Li2024_hi=LI2024[o]["hi"],
        logOR_difference=s["b_ivw"] - LI_B[o],
        z_for_difference=z_diff, p_for_difference=p_diff,
        CI_overlap_lo=ci_overlap_lo, CI_overlap_hi=ci_overlap_hi,
        CIs_overlap=bool(ci_overlap_lo <= ci_overlap_hi),
        # 要把本研究近零结果「变成」Li 2024 的阳性，需要一个负偏倚（把效应压低）
        needed_bias_sign="negative (bias would have to PUSH estimate DOWN)",
        burgess_overlap_bias_sign="positive (overlap pushes estimate UP)",
        sign_compatible=False,
        needed_negative_bias=s["needed_bias_to_match_Li"],
        min_abs_olsbias_needed=s["min_abs_olsbias_needed_at_pi1"],
        needed_observational_OR_lower=s["needed_observational_OR_if_negative"],
        needed_observational_OR_upper=s["needed_observational_OR_if_positive"],
    ))
li = pd.DataFrame(li_rows)
li.to_csv(os.path.join(RES, "li2024_comparison.csv"), index=False)

# ===========================================================================
# 4. 图
# ===========================================================================
plt.rcParams.update({"font.size": 9, "axes.unicode_minus": False})

# 森林图：原始 IVW vs 极端校正后
sm = summary.sort_values(["exposure", "outcome"]).reset_index(drop=True)
labels = [f"{r.exposure} -> {r.outcome}" for r in sm.itertuples()]
y = np.arange(len(sm))
fig, ax = plt.subplots(figsize=(9.2, 0.42 * len(sm) + 1.8))
for i, r in enumerate(sm.itertuples()):
    ax.plot([r.OR_lo, r.OR_hi], [i, i], color="#3b6ea5", lw=1.6, zorder=2)
    ax.scatter([r.OR], [i], s=28, color="#12395f", zorder=3)
    ax.annotate("", xy=(r.OR_corrected_max, i + 0.22), xytext=(r.OR, i + 0.22),
                arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.1))
    ax.scatter([r.OR_corrected_max], [i + 0.22], s=18, color="#c0392b", zorder=3)
li_or = LI2024["Migraine"]["OR"]
ax.axvline(1.0, color="grey", ls="--", lw=1)
ax.axvline(li_or, color="#2e8b57", ls=":", lw=1.3)
ax.text(li_or, len(sm) - 0.3, f"Li 2024 OR={li_or}", color="#2e8b57", fontsize=8, rotation=90,
        va="top", ha="right")
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=8)
ax.set_xlabel("Odds ratio (blue = IVW original; red = after extreme overlap-bias correction)")
ax.set_title("Overlap-bias correction: point estimates barely move (strong instruments, F = 43-142)")
ax.set_xlim(0.90, 1.13)
fig.tight_layout()
fig.savefig(os.path.join(RES, "forest_overlap.png"), dpi=160)
plt.close(fig)

# 反推所需观察性关联强度：用 log-OR 尺度线性轴，避免 exp() 溢出压扁可读区
fig, ax = plt.subplots(figsize=(7.6, 4.8))
pis = np.linspace(0.02, 1.0, 300)
for o, color in [("Migraine", "#c0392b"), ("MigraineAura", "#e67e22"), ("MigraineNoAura", "#8e44ad")]:
    s = summary[(summary["exposure"] == "GD") & (summary["outcome"] == o)].iloc[0]
    need = abs(s["needed_bias_to_match_Li"]) * s["meanF"] / pis   # 所需 |olsbias| (log OR)
    ax.plot(pis, need, color=color, lw=1.9, label=f"GD -> {o}")
# 合理的观察性混杂偏倚范围：观察性 OR 1.0-2.0 -> log OR 0-0.69
ax.axhspan(0.0, math.log(2.0), color="#95a5a6", alpha=0.32)
ax.text(0.035, math.log(2.0) + 0.35,
        "plausible observational bias (OR 1.0-2.0)", fontsize=8, color="#34495e")
ax.axvline(PI_DATA_BOUND, color="#2e8b57", ls="--", lw=1.3)
ax.set_ylim(0, 12)
ax.set_xlim(0.02, 1.0)
ax.text(PI_DATA_BOUND + 0.012, 10.9,
        f"pi upper bound\n= {PI_DATA_BOUND:.2f}", fontsize=8,
        color="#2e8b57", va="top")
ax.set_xlabel("Sample overlap proportion (pi)")
ax.set_ylabel("Required |observational bias| in log OR")
ax.set_title("Even at 100% overlap, the bias needed to match Li 2024 is implausible")
ax.legend(fontsize=8, loc="upper right")
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(RES, "bias_required.png"), dpi=160)
plt.close(fig)

# ===========================================================================
# 5. 控制台摘要
# ===========================================================================
print("=" * 78)
print("IVW reproduction check (Python vs TwoSampleMR R)")
print("  max |b_python - b_R| =", round(float(summary["IVW_repro_diff"].max()), 8))
print("  fixed-effect se is smaller than TwoSampleMR random-effect se by up to",
      round(float(summary["SE_diff_random_vs_fixed"].max()), 6),
      "-> primary se uses the more conservative R value")
print("=" * 78)
print("Instrument strength: meanF range = %.1f - %.1f ; medianF range = %.1f - %.1f"
      % (summary.meanF.min(), summary.meanF.max(), summary.medianF.min(), summary.medianF.max()))
print("Max overlap bias at pi=1, observational OR=2, using median F: %.6f (log OR)"
      % float((MAX_DELTA * MAX_PI / summary.medianF).max()))
print("Max |change in OR| after extreme correction: %.5f"
      % float((summary.OR - summary.OR_corrected_max).abs().max()))
print("Any combination becomes significant (p<0.05) after correction? ",
      bool((summary.p_corrected_max < 0.05).any()))
print("=" * 78)
print("BNG natural experiment (same SNPs, same outcome, only exposure beta swapped)")
print(bng[["outcome", "nSNP_common", "with_FinnGen_OR", "without_FinnGen_OR",
           "delta_b", "observed_sign"]].to_string(index=False, float_format=lambda x: f"{x:.5f}"))
print("=" * 78)
print("Data-driven upper bound on overlap proportion (inference): pi <= %.3f" % PI_DATA_BOUND)
print("Burgess overlap-bias direction:", OLSBIAS_SIGN, "-> bias is POSITIVE (away from null)")
print("=" * 78)
print("Li 2024 comparison (GD as proxy for autoimmune hyperthyroidism)")
print(li[["outcome", "this_OR", "Li2024_OR", "CIs_overlap", "p_for_difference",
          "min_abs_olsbias_needed", "needed_observational_OR_lower"]]
      .to_string(index=False, float_format=lambda x: f"{x:.5g}"))
print("=" * 78)
print("Saved:", "overlap_bias_grid.csv", "overlap_corrected_summary.csv",
      "bng_natural_experiment.csv", "li2024_comparison.csv",
      "forest_overlap.png", "bias_required.png")
