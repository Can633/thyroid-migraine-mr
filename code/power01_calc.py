# -*- coding: utf-8 -*-
"""
power01_calc.py -- 双样本 MR 功效 / 最小可检出效应（MDE）计算

方向 1（前向）：甲状腺病（GBMI）-> 偏头痛（FinnGen R12），15 个组合
方向 2（反向）：偏头痛（FinnGen R12）-> 甲状腺病（GBMI），15 个组合

做法：
  1. 从 harmonised 数据集（前向）与 IV/结局原始文件（反向，重新 harmonise 并校验）
     取实际入模 SNP 集，计算每个组合的
       - nSNP、总解释方差 R2 = sum(2f(1-f)b^2)（liability 近似）、F 均值/最小
  2. 以实测 IVW 标准误 SE 为精度基准（不依赖任何分布假设）：
       MDE(80%, a) = exp((z_{1-a/2} + z_0.8) * SE)
       power(OR_true) = Phi(|log OR_true| / SE - z_{1-a/2})
  3. 用解析式交叉验证：SE_hat = 1 / sqrt(sum(b_j^2 / se_Yj^2))
  4. 输出设计表、功效曲线表、两幅图

以二分类暴露的 MR，效应量单位为“暴露 log-odds 每增加 1 单位”，
与 Li 2024 等同类研究（同为二分类暴露）可直接比较。
"""
import os
import math
import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = r"D:/GWAS/Migraine"
RES = os.path.join(ROOT, "results")
REV = os.path.join(ROOT, "rev")

N_OUT_MIG = {"Migraine": 395060, "MigraineAura": 379299, "MigraineNoAura": 376762}  # FinnGen R13

EXPOS = ["GD", "HYPO", "LT", "BNG", "TC"]
EXP_LAB = {"GD": "Graves' disease", "HYPO": "Hypothyroidism",
           "LT": "Lymphocytic thyroiditis", "BNG": "Benign thyroid nodule",
           "TC": "Thyroid cancer"}
OUTS = ["Migraine", "MigraineAura", "MigraineNoAura"]
OUT_LAB = {"Migraine": "Migraine (overall)", "MigraineAura": "Migraine with aura",
           "MigraineNoAura": "Migraine without aura"}

ALPHA = 0.05
N_TEST = 15                      # 每方向 15 个组合 -> Bonferroni
ALPHA_B = ALPHA / N_TEST
TARGET_ORS = [1.02, 1.03, 1.05, 1.066, 1.10, 1.15, 1.20, 1.25, 1.30]


def z_crit(a):
    return norm.ppf(1 - a / 2)


def power_at(or_true, se, alpha):
    return float(norm.cdf(abs(math.log(or_true)) / se - z_crit(alpha)))


def mde(se, alpha, power=0.80):
    return float(math.exp((z_crit(alpha) + norm.ppf(power)) * se))


def is_palindromic(a1, a2):
    return (a1, a2) in {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}


# ---------------------------------------------------------------- 前向
def forward_design():
    """返回前向 15 组合的设计指标 + 实测 IVW 估计。"""
    # 注意：mr_summary_or.csv 为 11-12 旧表（HYPO 修复前），以 mr_results.csv 为准
    mr = pd.read_csv(os.path.join(RES, "mr_results.csv"))
    mr = mr[mr["method"] == "Inverse variance weighted"].copy()
    mr["OR"] = np.exp(mr["b"])
    mr["OR_lo"] = np.exp(mr["b"] - 1.96 * mr["se"])
    mr["OR_hi"] = np.exp(mr["b"] + 1.96 * mr["se"])
    rows = []
    for _, r in mr.iterrows():
        e, o = r["id.exposure"], r["id.outcome"]
        f = os.path.join(RES, "harmonised_%s_to_%s.csv" % (e, o))
        if not os.path.exists(f):
            print("  [skip] missing", f)
            continue
        d = pd.read_csv(f)
        if "mr_keep" in d.columns:
            d = d[d["mr_keep"].astype(str).str.upper().isin(["TRUE", "1", "T"])]
        d = d.dropna(subset=["beta.exposure", "se.exposure", "eaf.exposure"])
        b = d["beta.exposure"].values.astype(float)
        se = d["se.exposure"].values.astype(float)
        f_ = d["eaf.exposure"].values.astype(float)
        f_ = np.clip(f_, 1e-4, 1 - 1e-4)
        rsq_j = 2 * f_ * (1 - f_) * b ** 2          # liability 近似
        F_j = (b / se) ** 2
        # 解析式复核：IVW 标准误
        seY = d["se.outcome"].values.astype(float)
        bY = d["beta.outcome"].values.astype(float)
        w = b ** 2 / seY ** 2
        se_hat = 1 / math.sqrt(w.sum())          # 固定效应（同质因果效应）
        b_ivw = (w * (bY / b)).sum() / w.sum() if np.all(b != 0) else np.nan
        rows.append(dict(
            direction="forward", exposure=e, outcome=o,
            exp_lab=EXP_LAB[e], out_lab=OUT_LAB[o],
            n_snp=len(d), n_snp_mr=int(r["nsnp"]),
            rsq_logodds=rsq_j.sum(), F_mean=F_j.mean(), F_min=F_j.min(),
            N_out=N_OUT_MIG[o],
            OR=float(r["OR"]), OR_lo=float(r["OR_lo"]), OR_hi=float(r["OR_hi"]),
            b=float(r["b"]), se=float(r["se"]), p=float(r["pval"]),
            se_analytic=se_hat, b_analytic=b_ivw,
        ))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- 反向
def harmonise_reverse(iv, ou):
    """把偏头痛 IV 与甲状腺结局对齐（action=2 规则：回文且 MAF>0.42 者丢弃）。"""
    m = iv.merge(ou, on="rsid", how="inner", suffixes=("_x", "_y"))
    m = m.dropna(subset=["beta_x", "se_x", "eaf_x", "beta_y", "se_y"])
    m = m[~m["rsid"].duplicated()]

    ea_x = m["alt"].str.upper()
    oa_x = m["ref"].str.upper()
    ea_y = m["ea"].str.upper()
    oa_y = m["oa"].str.upper()

    same = (ea_x == ea_y) & (oa_x == oa_y)
    swap = (ea_x == oa_y) & (oa_x == ea_y)
    keep = same | swap

    maf_x = np.minimum(m["eaf_x"], 1 - m["eaf_x"])
    pal = np.array([is_palindromic(a, b) for a, b in zip(ea_x, oa_x)], dtype=bool)
    amb = pal & (maf_x > 0.42)             # 无法判定链方向
    keep = keep & ~amb

    m = m[keep].copy()
    sw = swap[keep].values
    m.loc[sw, "beta_y"] = -m.loc[sw, "beta_y"]
    m.loc[sw, "eaf_y"] = 1 - m.loc[sw, "eaf_y"]
    return m


def reverse_design():
    mr = pd.read_csv(os.path.join(RES, "rev_mr_main.csv"))
    mr = mr[mr["method"] == "Inverse variance weighted"].copy()
    rows = []
    for _, r in mr.iterrows():
        e, o = r["id.exposure"], r["id.outcome"]      # e = migraine subtype
        iv = pd.read_csv(os.path.join(REV, "iv_mig_%s.tsv" % e), sep="\t")
        ou = pd.read_csv(os.path.join(REV, "out_%s.tsv" % o), sep="\t")
        m = harmonise_reverse(iv, ou)
        b = m["beta_x"].values.astype(float)
        se = m["se_x"].values.astype(float)
        f_ = np.clip(m["eaf_x"].values.astype(float), 1e-4, 1 - 1e-4)
        rsq_j = 2 * f_ * (1 - f_) * b ** 2
        F_j = (b / se) ** 2
        seY = m["se_y"].values.astype(float)
        bY = m["beta_y"].values.astype(float)
        w = b ** 2 / seY ** 2
        se_hat = 1 / math.sqrt(w.sum())
        b_ivw = (w * (bY / b)).sum() / w.sum()
        rows.append(dict(
            direction="reverse", exposure=e, outcome=o,
            exp_lab=OUT_LAB[e], out_lab=EXP_LAB[o],
            n_snp=len(m), n_snp_mr=int(r["nsnp"]),
            rsq_logodds=rsq_j.sum(), F_mean=F_j.mean(), F_min=F_j.min(),
            N_out=float(m["n"].median()) if "n" in m else np.nan,
            OR=float(r["or"]), OR_lo=float(r["or_lci95"]), OR_hi=float(r["or_uci95"]),
            b=float(r["b"]), se=float(r["se"]), p=float(r["pval"]),
            se_analytic=se_hat, b_analytic=b_ivw,
        ))
    return pd.DataFrame(rows)


def main():
    print("== 前向设计 ==")
    fw = forward_design()
    print("== 反向设计 ==")
    rv = reverse_design()
    d = pd.concat([fw, rv], ignore_index=True)

    # se 为论文实际报告的标准误（TwoSampleMR 的乘性随机效应 MRE）；
    # se_analytic 为固定效应（同质因果效应）标准误，作为设计功效的参照
    d["se_fe"] = d["se_analytic"]
    d["se_inflation"] = d["se"] / d["se_fe"]
    d["se_rel_err"] = (d["se_analytic"] - d["se"]) / d["se"]
    d["mde80"] = [mde(s, ALPHA) for s in d["se"]]                 # 基于报告 CI
    d["mde80_bonf"] = [mde(s, ALPHA_B) for s in d["se"]]
    d["mde80_fe"] = [mde(s, ALPHA) for s in d["se_fe"]]           # 基于设计（同质效应）
    for t in TARGET_ORS:
        d["pw_%.3f" % t] = [power_at(t, s, ALPHA) for s in d["se"]]
    d["pw_bonf_1.10"] = [power_at(1.10, s, ALPHA_B) for s in d["se"]]
    # 观测 CI 已经排除到哪
    d["ci_excl_above"] = d["OR_hi"]
    d["ci_excl_below"] = d["OR_lo"]

    cols = ["direction", "exposure", "outcome", "exp_lab", "out_lab", "n_snp",
            "rsq_logodds", "F_mean", "F_min", "N_out", "OR", "OR_lo", "OR_hi",
            "b", "se", "se_fe", "se_inflation", "se_rel_err", "p",
            "mde80", "mde80_bonf", "mde80_fe", "ci_excl_above", "ci_excl_below"]
    out = d[cols + ["pw_%.3f" % t for t in TARGET_ORS] + ["pw_bonf_1.10"]]
    out.to_csv(os.path.join(RES, "power_design.csv"), index=False)

    # 功效曲线
    grid = np.exp(np.linspace(math.log(0.85), math.log(1.6), 300))
    curves = []
    for _, r in d.iterrows():
        for a_lab, a in (("nominal", ALPHA), ("bonferroni", ALPHA_B)):
            for orv in grid:
                curves.append(dict(direction=r["direction"], exposure=r["exposure"],
                                   outcome=r["outcome"], exp_lab=r["exp_lab"],
                                   out_lab=r["out_lab"], alpha=a_lab,
                                   OR=float(orv), power=power_at(orv, r["se"], a)))
    pd.DataFrame(curves).to_csv(os.path.join(RES, "power_curves.csv"), index=False)

    # ---- 控制台汇总 ----
    pd.set_option("display.width", 250)
    print("\n===== 前向（甲状腺 -> 偏头痛）=====")
    a = d[d.direction == "forward"].copy()
    print(a[["exposure", "outcome", "n_snp", "rsq_logodds", "F_mean",
             "OR", "OR_lo", "OR_hi", "se", "se_inflation", "mde80", "mde80_fe",
             "mde80_bonf", "pw_1.066", "pw_1.100"]].to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    print("\n===== 反向（偏头痛 -> 甲状腺）=====")
    b_ = d[d.direction == "reverse"].copy()
    print(b_[["exposure", "outcome", "n_snp", "rsq_logodds", "F_mean",
              "OR", "OR_lo", "OR_hi", "se", "se_inflation", "mde80", "mde80_fe",
              "mde80_bonf", "pw_1.066", "pw_1.100"]].to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    print("\n===== 管线校验：解析式 SE vs 实测 SE 的相对误差 =====")
    print(d[["direction", "exposure", "outcome", "se", "se_analytic", "se_rel_err"]]
          .to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    return d


if __name__ == "__main__":
    main()
