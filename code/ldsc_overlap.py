#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ldsc_overlap.py  —— F 方案：用本地 1000G EUR LD scores 实测样本重叠并重新校正 MR
================================================================================
背景
----
E 方案（overlap_correct.py）因当时无 LD scores，只能用 Burgess 2016 解析公式
配合假设的 overlap.prop 网格做偏倚界限分析。用户随后提供
D:/GWAS/LD core 下的 1000G Phase3 EUR 参考文件，使真正的
cross-trait LD Score regression 成为可能。

本脚本做三件事
--------------
(1) 单性状 LDSC：估计 9 个 GWAS（6 暴露 + 3 结局）的截距 I11/I22 与观测尺度 h2。
(2) cross-trait LDSC：18 组（含 BNG 主版/LOO 版对照）估计
      截距 I12、归一化采样误差相关 rho_s = I12/sqrt(I11*I22)、遗传协方差与 rg。
    理论（Bulik-Sullivan 2015 Nat Genet; Genomic SEM 2019）：
      E[z1*z2] = sqrt(N1*N2)*rho_g/M * l_j + rho*Ns/sqrt(N1*N2)
    其中截距项即样本重叠；斜率不受重叠影响。
(3) 用实测 rho_s 直接校正 15 组合 IVW 的样本重叠偏倚：
      Cov(e_Xj, e_Yj) = rho_s * se_Xj * se_Yj
      bias(IVW | null) = sum_j w_j*Cov_j / sum_j w_j*E[beta_Xj^2],  w_j = 1/se_Yj^2
    并与 E 方案的假设网格校正结果对比。

与 E 方案的关键差别
--------------------
E: overlap.prop 与 olsbias 均为假设网格（0-1、OR 1-2），只能给偏倚"上界"。
F: rho_s 由全基因组实测截距给出，无需假设重叠比例与观察性关联强度。

参考文件（用户本地提供，1000G Phase3 EUR, GRCh37）
-------------------------------------------------
  LDscore/LDscore.{1..22}.l2.ldscore.gz     ref-ld
  LDscore/LDscore.{1..22}.l2.M              M = 9,997,231
  1000G_Phase3_weights_hm3_no_MHC/...       w-ld (regression weights)
  1000G_Phase3_frq/1000G.EUR.QC.*.frq      MAF
  hm3_no_MHC.list.txt                      1,217,311 HapMap3 SNPs

输出（results/）
---------------
  ldsc_univariate.csv              单性状 LDSC（截距、h2、SE）
  ldsc_crosstrait.csv              cross-trait LDSC（I12, rho_s, rg, Ns 反推）
  ldsc_overlap_corrected.csv      实测 rho_s 校正后的 15 组合 IVW
  ldsc_vs_E_comparison.csv        F vs E 校正幅度对比
  ldsc_intercepts.png             实测截距 / rho_s 图（含 BNG 天然对照）
  forest_ldsc_corrected.png       原始 vs F 校正 vs E 极端校正 森林图
  ldsc_rg.png                     遗传相关 rg 图
"""
import os
import math
import warnings
import numpy as np
import pandas as pd
from scipy.stats import norm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = "D:/GWAS/Migraine"
LDDIR = "D:/GWAS/LD core/extracted"
HM3 = "D:/GWAS/LD core/hm3_no_MHC.list.txt"
CACHE = os.path.join(ROOT, "ldsc_cache")
RES = os.path.join(ROOT, "results")
os.makedirs(RES, exist_ok=True)

M_TOTAL = 9997231          # sum of LDscore.*.l2.M (1000G Phase3 EUR)
N_BLOCKS = 200             # block jackknife
CHI2_CUT = 100.0           # |z| > 10 剔除（主分析），敏感性保留全部

EXPS = ["GD", "HYPO", "LT", "BNG", "TC", "BNG_LOO"]
OUTS = ["Migraine", "MigraineAura", "MigraineNoAura"]
OUT_N = {"Migraine": 395060, "MigraineAura": 379299, "MigraineNoAura": 376762}  # FinnGen R13
ACC = {"GD": "GCST90627744", "HYPO": "GCST90627749", "LT": "GCST90627754",
       "BNG": "GCST90627740", "TC": "GCST90627763", "BNG_LOO": "GCST90627765"}
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

# 方法学正对照：同一 FinnGen 队列内的偏头痛亚型之间（已知 100% 样本重叠、
# 表型高度相关）应给出明显 >0 的截距；GWAS 与自身则理论 rho_s = 1。
# 若这些方法能测到强信号而甲状腺×偏头痛测不到，才说明「未检出重叠」不是方法失灵。
POSITIVE_PAIRS = [
    ("Migraine", "MigraineAura"), ("Migraine", "MigraineNoAura"),
    ("MigraineAura", "MigraineNoAura"),
    ("Migraine", "Migraine"), ("GD", "GD"),
]

COLS = ["SNP", "CHR", "BP", "EA", "OA", "BETA", "SE", "EAF", "N", "TAG"]


# =============================================================================
# 1. 参考面板
# =============================================================================
def load_reference():
    print("[1/6] 载入 1000G EUR LD score 参考面板 ...", flush=True)
    hm3 = set(pd.read_csv(HM3, header=None, names=["SNP"])["SNP"].astype(str))
    print(f"      HapMap3 (no MHC) SNPs: {len(hm3):,}", flush=True)

    lds, wts, frqs = [], [], []
    for c in range(1, 23):
        p = os.path.join(LDDIR, "LDscore", f"LDscore.{c}.l2.ldscore.gz")
        d = pd.read_csv(p, sep="\t", usecols=["CHR", "SNP", "BP", "L2"])
        lds.append(d[d["SNP"].isin(hm3)])

        p = os.path.join(LDDIR, "1000G_Phase3_weights_hm3_no_MHC",
                         f"weights.hm3_noMHC.{c}.l2.ldscore.gz")
        d = pd.read_csv(p, sep="\t", usecols=["SNP", "BP", "L2"])
        wts.append(d.rename(columns={"L2": "W2", "BP": "WBP"}))

        p = os.path.join(LDDIR, "1000G_Phase3_frq", f"1000G.EUR.QC.{c}.frq")
        d = pd.read_csv(p, sep=r"\s+")
        d.columns = [x.strip() for x in d.columns]
        frqs.append(d[["SNP", "MAF"]])

    ref = pd.concat(lds, ignore_index=True)
    w = pd.concat(wts, ignore_index=True)
    f = pd.concat(frqs, ignore_index=True).drop_duplicates("SNP")
    ref = ref.merge(w[["SNP", "W2"]], on="SNP", how="inner")
    ref = ref.merge(f, on="SNP", how="left")
    ref = ref.drop_duplicates("SNP").sort_values(["CHR", "BP"]).reset_index(drop=True)
    ref = ref[(ref["MAF"].isna()) | ((ref["MAF"] >= 0.01) & (ref["MAF"] <= 0.99))]
    ref = ref.reset_index(drop=True)
    # 连续 SNP 分块（标准 ldsc 200 blocks）
    edges = np.linspace(0, len(ref), N_BLOCKS + 1).astype(int)
    blk = np.zeros(len(ref), dtype=int)
    for k in range(N_BLOCKS):
        blk[edges[k]:edges[k + 1]] = k
    ref["block"] = blk
    print(f"      回归 SNP: {len(ref):,}  (ref-ld ∩ w-ld ∩ hm3, MAF 0.01-0.99)", flush=True)
    print(f"      L2 mean={ref.L2.mean():.2f}  W2 mean={ref.W2.mean():.2f}  M={M_TOTAL:,}", flush=True)
    return ref


# =============================================================================
# 2. 载入 hm3 过滤后的 GWAS
# =============================================================================
def load_gwas(name):
    p = os.path.join(CACHE, f"{name}.tsv")
    d = pd.read_csv(p, sep="\t", header=None, names=COLS,
                    dtype={"SNP": str, "EA": str, "OA": str, "TAG": str})
    for c in ["BETA", "SE", "EAF", "N", "CHR", "BP"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.dropna(subset=["BETA", "SE"])
    d = d[(d["SE"] > 0) & np.isfinite(d["BETA"]) & np.isfinite(d["SE"])]
    d = d.drop_duplicates("SNP", keep="first")
    d["z"] = d["BETA"] / d["SE"]
    d = d[np.isfinite(d["z"])]
    return d.set_index("SNP")


# =============================================================================
# 3. LD score 回归（WLS + 200 block jackknife）
# =============================================================================
def _point(X, y, w):
    wx = w[:, None] * X
    A = X.T @ wx
    b = X.T @ (w * y)
    try:
        return np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        return np.linalg.lstsq(A, b, rcond=None)[0]


def _jk(X, y, w):
    """delete-one-block jackknife; O(n + B) via cumulative sums."""
    n = len(y)
    comp = np.column_stack([
        w * X[:, 0] * X[:, 0], w * X[:, 0] * X[:, 1], w * X[:, 1] * X[:, 1],
        w * X[:, 0] * y,       w * X[:, 1] * y,
    ])
    cs = np.vstack([np.zeros(5), np.cumsum(comp, axis=0)])
    A = np.array([[cs[-1, 0], cs[-1, 1]], [cs[-1, 1], cs[-1, 2]]])
    b = cs[-1, 3:5]
    full = np.linalg.solve(A, b)
    edges = np.linspace(0, n, N_BLOCKS + 1).astype(int)
    est = np.empty((N_BLOCKS, 2))
    for k in range(N_BLOCKS):
        s = cs[edges[k + 1]] - cs[edges[k]]
        Ak = A - np.array([[s[0], s[1]], [s[1], s[2]]])
        bk = b - s[3:5]
        est[k] = np.linalg.solve(Ak, bk)
    var = (N_BLOCKS - 1) / N_BLOCKS * np.sum((est - est.mean(0)) ** 2, axis=0)
    return full, np.sqrt(var)


def ldsc_fit(L2, W2, y, n_iter=3):
    """ldsc 风格迭代异方差 WLS：w = 1/(1+c*L2)^2 / max(W2,1)。"""
    X = np.column_stack([np.ones(len(L2)), L2])
    c = 0.0
    beta = None
    for _ in range(n_iter):
        w = 1.0 / (1.0 + max(c, 0.0) * L2) ** 2 / np.maximum(W2, 1.0)
        beta = _point(X, y, w)
        c = max(float(beta[1]), 0.0)
    w = 1.0 / (1.0 + max(c, 0.0) * L2) ** 2 / np.maximum(W2, 1.0)
    full, se = _jk(X, y, w)
    # OLS 对照（检验权重形式对点估计的影响）
    ols = _point(X, y, np.ones(len(y)))
    return dict(intercept=float(full[0]), se_intercept=float(se[0]),
                slope=float(full[1]), se_slope=float(se[1]),
                intercept_ols=float(ols[0]), slope_ols=float(ols[1]),
                n_snp=len(y))


# =============================================================================
# 4. 等位基因对齐
# =============================================================================
def align_z(x, y, ref_l2):
    """把暴露 x 的 z 对齐到结局 y 的效应等位基因；返回 merged DataFrame。"""
    a = x[["z", "EA", "OA", "CHR", "BP"]].reset_index()
    b = y[["z", "EA", "OA", "CHR", "BP"]].reset_index()
    m = a.merge(b, on="SNP", suffixes=("_x", "_y"))
    ok = (m.EA_x.str.len() == 1) & (m.OA_x.str.len() == 1) & \
         (m.EA_y.str.len() == 1) & (m.OA_y.str.len() == 1)
    m = m[ok]
    same = (m.EA_x == m.EA_y) & (m.OA_x == m.OA_y)
    flip = (m.EA_x == m.OA_y) & (m.OA_x == m.EA_y)
    cx, cox = m.EA_x.map(COMP), m.OA_x.map(COMP)
    same_c = (cx == m.EA_y) & (cox == m.OA_y)
    flip_c = (cx == m.OA_y) & (cox == m.EA_y)
    sign_all = np.where(same | same_c, 1, np.where(flip | flip_c, -1, 0))
    keep = sign_all != 0
    m = m[keep].copy()
    sign = sign_all[keep]
    m["zx"] = m["z_x"].to_numpy(float) * sign
    m["zy"] = m["z_y"].to_numpy(float)
    m = m.merge(ref_l2, on="SNP", how="inner")
    m["chi_x"] = m["zx"] ** 2
    m["chi_y"] = m["zy"] ** 2
    return m


# =============================================================================
# main
# =============================================================================
def main():
    ref = load_reference()
    ref_l2 = ref[["SNP", "L2", "W2"]].copy()

    print("[2/6] 载入 9 个 GWAS 的 hm3 z-score ...", flush=True)
    g = {}
    for nm in EXPS + OUTS:
        g[nm] = load_gwas(nm)
        nval = g[nm]["N"].median() if nm in EXPS else OUT_N[nm]
        print(f"      {nm:<15s} nSNP={len(g[nm]):>9,}  median N={nval:>10,.0f}", flush=True)
    N1 = {nm: float(g[nm]["N"].median()) for nm in EXPS}
    N2 = OUT_N

    # ---------- 单性状 ----------
    print("[3/6] 单性状 LDSC ...", flush=True)
    uni_rows = []
    single = {}
    for nm in EXPS + OUTS:
        d = g[nm].reset_index().merge(ref_l2, on="SNP", how="inner")
        d = d[np.isfinite(d["z"]) & (d["z"] ** 2 <= CHI2_CUT)]
        r = ldsc_fit(d["L2"].to_numpy(float), d["W2"].to_numpy(float),
                     (d["z"].to_numpy(float)) ** 2)
        N = N1[nm] if nm in EXPS else N2[nm]
        r.update(trait=nm, N=N, acc=ACC.get(nm, ""),
                 h2_observed=r["slope"] * M_TOTAL / N,
                 se_h2=r["se_slope"] * M_TOTAL / N)
        uni_rows.append(r)
        single[nm] = r
        print(f"      {nm:<15s} intercept={r['intercept']:.4f} (se {r['se_intercept']:.4f})  "
              f"h2_obs={r['h2_observed']:.4f}  OLS_int={r['intercept_ols']:.4f}", flush=True)
    uni = pd.DataFrame(uni_rows)
    uni.to_csv(os.path.join(RES, "ldsc_univariate.csv"), index=False)

    # ---------- cross-trait ----------
    print("[4/6] cross-trait LDSC (18 组) ...", flush=True)
    ct_rows = []
    cross = {}
    for e in EXPS:
        for o in OUTS:
            m = align_z(g[e], g[o], ref_l2)
            n_raw = len(m)
            m = m[(m.chi_x <= CHI2_CUT) & (m.chi_y <= CHI2_CUT)]
            r = ldsc_fit(m["L2"].to_numpy(float), m["W2"].to_numpy(float),
                         (m["zx"] * m["zy"]).to_numpy(float))
            I11, I22 = single[e]["intercept"], single[o]["intercept"]
            I12 = r["intercept"]
            rho_s = I12 / math.sqrt(max(I11 * I22, 1e-12))
            # rg = slope12 / sqrt(slope1*slope2)（M 与 N 相消，避免尺度假设）
            s1, s2, s12 = single[e]["slope"], single[o]["slope"], r["slope"]
            rg = s12 / math.sqrt(s1 * s2) if (s1 > 0 and s2 > 0) else float("nan")
            rho_g = r["slope"] * M_TOTAL / math.sqrt(N1[e] * N2[o])
            # 截距反推重叠个体数（需假设重叠样本中表型相关 rho）
            sqrtN = math.sqrt(N1[e] * N2[o])
            ns = {f"Ns_rho{int(q*100):02d}": I12 * sqrtN / q for q in (0.02, 0.05, 0.10, 0.20)}
            row = dict(exposure=e, outcome=o, acc_exposure=ACC[e],
                       N_exposure=N1[e], N_outcome=N2[o],
                       n_snp_raw=n_raw, n_snp_used=r["n_snp"],
                       I11=I11, I22=I22, I12=I12, se_I12=r["se_intercept"],
                       I12_ols=r["intercept_ols"], rho_s=rho_s,
                       rho_g=rho_g, slope=r["slope"], se_slope=r["se_slope"],
                       rg=rg, **ns)
            ct_rows.append(row)
            cross[(e, o)] = row
            print(f"      {e:<9s}->{o:<15s} n={r['n_snp']:>8,}  I12={I12:+.4f} "
                  f"(se {r['se_intercept']:.4f})  rho_s={rho_s:+.4f}  rg={rg:+.3f}", flush=True)
    ct = pd.DataFrame(ct_rows)
    ct.to_csv(os.path.join(RES, "ldsc_crosstrait.csv"), index=False)

    # ---------- 方法学正对照 ----------
    print("[4b] 方法学正对照（已知重叠应被检出）...", flush=True)
    pos_rows = []
    for a, b in POSITIVE_PAIRS:
        m = align_z(g[a], g[b], ref_l2)
        m = m[(m.chi_x <= CHI2_CUT) & (m.chi_y <= CHI2_CUT)]
        r = ldsc_fit(m["L2"].to_numpy(float), m["W2"].to_numpy(float),
                     (m["zx"] * m["zy"]).to_numpy(float))
        Ia, Ib = single[a]["intercept"], single[b]["intercept"]
        rho_s = r["intercept"] / math.sqrt(max(Ia * Ib, 1e-12))
        pos_rows.append(dict(trait_a=a, trait_b=b, n_snp=r["n_snp"],
                             I12=r["intercept"], se_I12=r["se_intercept"],
                             I12_ols=r["intercept_ols"], rho_s=rho_s,
                             z_vs_zero=r["intercept"] / r["se_intercept"] if r["se_intercept"] > 0 else float("nan"),
                             expected="strong overlap" if a != b else "self (rho_s=1)"))
        print(f"      {a:<15s} x {b:<15s} I12={r['intercept']:+.4f} "
              f"(se {r['se_intercept']:.4f})  rho_s={rho_s:+.3f}", flush=True)
    pos = pd.DataFrame(pos_rows)
    pos.to_csv(os.path.join(RES, "ldsc_positive_controls.csv"), index=False)

    # ---------- 用实测 rho_s 校正 IVW ----------
    print("[5/6] 用实测 rho_s 校正 15 组合 IVW ...", flush=True)
    mr = pd.read_csv(os.path.join(RES, "mr_results.csv"))
    ivw_r = mr[mr["method"].astype(str).str.contains("Inverse variance weighted", case=False, na=False)]
    corr_rows = []
    for e in ["GD", "HYPO", "LT", "BNG", "TC"]:
        for o in OUTS:
            hp = os.path.join(RES, f"harmonised_{e}_to_{o}.csv")
            h = pd.read_csv(hp)
            h = h[h["mr_keep"] == True] if "mr_keep" in h else h
            bx = h["beta.exposure"].to_numpy(float)
            by = h["beta.outcome"].to_numpy(float)
            sx = h["se.exposure"].to_numpy(float)
            sy = h["se.outcome"].to_numpy(float)
            w = 1.0 / sy ** 2
            num = float(np.sum(w * bx * by))
            den = float(np.sum(w * bx ** 2))
            b = num / den
            se_fixed = math.sqrt(1.0 / den)
            rr = ivw_r[(ivw_r["id.exposure"] == e) & (ivw_r["id.outcome"] == o)]
            r_se = float(rr["se"].iloc[0]) if len(rr) else se_fixed
            r_b = float(rr["b"].iloc[0]) if len(rr) else b

            rs = cross[(e, o)]["rho_s"]
            # 保守敏感性：OLS 截距（不加权）通常更大，用它做偏倚上界
            rs_ols = cross[(e, o)]["I12_ols"] / math.sqrt(
                max(cross[(e, o)]["I11"] * cross[(e, o)]["I22"], 1e-12))
            cov = rs * sx * sy
            cov_ols = rs_ols * sx * sy
            bias_num = float(np.sum(w * cov))
            bias_num_ols = float(np.sum(w * cov_ols))
            bias = bias_num / den
            bias_ols = bias_num_ols / den
            den_unb = float(np.sum(w * np.maximum(bx ** 2 - sx ** 2, 1e-12)))
            bias_unb = bias_num / den_unb
            bc = r_b - bias
            bc_unb = r_b - bias_unb
            F = (bx / sx) ** 2
            row = dict(exposure=e, outcome=o, nSNP=len(h),
                       rho_s=rs, rho_s_ols=rs_ols, I12=cross[(e, o)]["I12"],
                       I12_ols=cross[(e, o)]["I12_ols"],
                       b_ivw=r_b, se_ivw=r_se,
                       OR=math.exp(r_b), OR_lo=math.exp(r_b - 1.96 * r_se),
                       OR_hi=math.exp(r_b + 1.96 * r_se),
                       p_ivw=float(2 * norm.sf(abs(r_b / r_se))),
                       bias_ldsc=bias, b_corrected=bc, OR_corrected=math.exp(bc),
                       p_corrected=float(2 * norm.sf(abs(bc / r_se))),
                       bias_ols=bias_ols, b_corrected_ols=r_b - bias_ols,
                       OR_corrected_ols=math.exp(r_b - bias_ols),
                       p_corrected_ols=float(2 * norm.sf(abs((r_b - bias_ols) / r_se))),
                       bias_unbiased_den=bias_unb, b_corrected_unb=bc_unb,
                       OR_corrected_unb=math.exp(bc_unb),
                       p_corrected_unb=float(2 * norm.sf(abs(bc_unb / r_se))),
                       meanF=float(F.mean()), medianF=float(np.median(F)),
                       N_exposure=cross[(e, o)]["N_exposure"],
                       N_outcome=N2[o])
            corr_rows.append(row)
    corr = pd.DataFrame(corr_rows)
    corr.to_csv(os.path.join(RES, "ldsc_overlap_corrected.csv"), index=False)

    # F vs E
    ep = os.path.join(RES, "overlap_corrected_summary.csv")
    cmp_rows = []
    if os.path.exists(ep):
        E = pd.read_csv(ep)
        for _, r in corr.iterrows():
            er = E[(E["exposure"] == r["exposure"]) & (E["outcome"] == r["outcome"])]
            if len(er) == 0:
                continue
            er = er.iloc[0]
            cmp_rows.append(dict(
                exposure=r["exposure"], outcome=r["outcome"],
                OR_raw=r["OR"], OR_F_ldsc=r["OR_corrected"],
                OR_E_extreme=er["OR_corrected_max_medianF"],
                shift_F=abs(r["OR"] - r["OR_corrected"]),
                shift_E=abs(r["OR"] - er["OR_corrected_max_medianF"]),
                p_raw=r["p_ivw"], p_F=r["p_corrected"],
                rho_s=r["rho_s"],
                pi_E_grid_max=1.0,
                note="F uses LDSC-measured rho_s; E uses assumed pi<=1 & OR_obs<=2"))
    cmp = pd.DataFrame(cmp_rows)
    cmp.to_csv(os.path.join(RES, "ldsc_vs_E_comparison.csv"), index=False)

    # ---------- 图 ----------
    print("[6/6] 绘图 ...", flush=True)
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})

    # 4.1 截距与 rho_s
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2))
    lab = [f"{r.exposure}→{r.outcome.replace('Migraine','Mig')}" for _, r in ct.iterrows()]
    xpos = np.arange(len(ct))
    ax = axes[0]
    ax.bar(xpos, ct["I12"], yerr=1.96 * ct["se_I12"], color="#4a7ebb", alpha=.85,
           error_kw=dict(ecolor="#2c3e50", capsize=2, lw=.8))
    ax.axhline(0, color="#c0392b", lw=1, ls="--")
    ax.set_xticks(xpos); ax.set_xticklabels(lab, rotation=75, ha="right", fontsize=6)
    ax.set_ylabel("cross-trait LDSC intercept  $I_{12}$")
    ax.set_title("Measured sample-overlap signal (Bulik-Sullivan 2015)")
    ax.grid(axis="y", alpha=.3)
    ax = axes[1]
    ax.bar(xpos, ct["rho_s"], color="#8e6bbf", alpha=.85)
    ax.axhline(0, color="#c0392b", lw=1, ls="--")
    ax.set_xticks(xpos); ax.set_xticklabels(lab, rotation=75, ha="right", fontsize=6)
    ax.set_ylabel(r"$\rho_s = I_{12}/\sqrt{I_{11}I_{22}}$")
    ax.set_title("Normalised sampling-error correlation")
    ax.grid(axis="y", alpha=.3)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "ldsc_intercepts.png"), dpi=160)
    plt.close(fig)

    # 4.2 森林图
    fig, ax = plt.subplots(figsize=(8.6, 8.2))
    order = corr.sort_values(["outcome", "exposure"]).reset_index(drop=True)
    ypos = np.arange(len(order))[::-1]
    for i, (_, r) in zip(ypos, order.iterrows()):
        ax.plot([r.OR_lo, r.OR_hi], [i, i], color="#b0b8c0", lw=2.4, zorder=1)
        ax.scatter([r.OR], [i + .20], s=34, color="#2c3e50", zorder=3, label="Raw IVW" if i == ypos[0] else None)
        ax.scatter([r.OR_corrected], [i], s=38, color="#c0392b", zorder=3, label="F: LDSC-measured correction" if i == ypos[0] else None)
        ax.scatter([r.OR_corrected_unb], [i - .20], s=24, marker="^", color="#e67e22", zorder=3,
                   label="F: unbiased-denominator variant" if i == ypos[0] else None)
    ax.axvline(1.0, color="#7f8c8d", ls="--", lw=1)
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{r.exposure} → {r.outcome}" for _, r in order.iterrows()], fontsize=7.5)
    ax.set_xlabel("Odds ratio (95% CI for raw IVW)")
    ax.set_title("15 combinations: raw IVW vs LDSC-measured overlap correction")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(axis="x", alpha=.3)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "forest_ldsc_corrected.png"), dpi=160)
    plt.close(fig)

    # 4.3 rg
    fig, ax = plt.subplots(figsize=(9.5, 4.4))
    ax.bar(xpos, ct["rg"], color="#3a9d8f", alpha=.85)
    ax.axhline(0, color="#2c3e50", lw=1)
    ax.set_xticks(xpos); ax.set_xticklabels(lab, rotation=75, ha="right", fontsize=6)
    ax.set_ylabel("genetic correlation $r_g$ (observed scale)")
    ax.set_title("Thyroid disease × migraine genetic correlation from LDSC slope")
    ax.grid(axis="y", alpha=.3)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "ldsc_rg.png"), dpi=160)
    plt.close(fig)

    # 4.4 方法学正对照：已知重叠应能被测出
    fig, ax = plt.subplots(figsize=(9.0, 4.4))
    xp = np.arange(len(pos))
    ax.bar(xp, pos["I12"], yerr=1.96 * pos["se_I12"], color="#c0392b", alpha=.85,
           error_kw=dict(ecolor="#2c3e50", capsize=3, lw=.9))
    lo, hi = float(ct["I12"].min()), float(ct["I12"].max())
    ax.axhspan(lo, hi, color="#7f8c8d", alpha=.30)
    ax.text(0.01, hi, f"  thyroid × migraine I12 range [{lo:+.4f}, {hi:+.4f}]",
            fontsize=8, color="#34495e", va="bottom")
    ax.axhline(0, color="#2c3e50", lw=1)
    ax.set_xticks(xp)
    ax.set_xticklabels([f"{r.trait_a} × {r.trait_b}" for _, r in pos.iterrows()],
                       rotation=20, ha="right", fontsize=7.5)
    ax.set_ylabel("cross-trait LDSC intercept  $I_{12}$")
    ax.set_title("Positive controls: the method does detect known sample overlap")
    ax.grid(axis="y", alpha=.3)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "ldsc_positive_controls.png"), dpi=160)
    plt.close(fig)

    # ---------- 汇总输出 ----------
    print("\n" + "=" * 80)
    print("F 方案汇总")
    print("=" * 80)
    print(f"单性状截距范围: {uni.intercept.min():.4f} – {uni.intercept.max():.4f}")
    print(f"cross-trait I12 范围: {ct.I12.min():+.4f} – {ct.I12.max():+.4f}")
    print(f"rho_s 范围: {ct.rho_s.min():+.4f} – {ct.rho_s.max():+.4f}")
    print(f"rg 范围: {ct.rg.min():+.3f} – {ct.rg.max():+.3f}")
    print(f"\n实测校正后 OR 变化: 最大 {np.abs(corr.OR - corr.OR_corrected).max():.6f} "
          f"(无偏分母变体 {np.abs(corr.OR - corr.OR_corrected_unb).max():.6f})")
    print(f"保守 OLS 截距校正 OR 变化: 最大 {np.abs(corr.OR - corr.OR_corrected_ols).max():.6f}")
    print(f"校正后 p<0.05 个数: {(corr.p_corrected < .05).sum()} / {len(corr)} "
          f"(OLS 保守 {(corr.p_corrected_ols < .05).sum()} / {len(corr)})")
    print(f"校正后 p<0.05 个数: {(corr.p_corrected < .05).sum()} / {len(corr)}")
    print(f"原始 p<0.05 个数:   {(corr.p_ivw < .05).sum()} / {len(corr)}")
    if len(cmp):
        print(f"E 方案极端网格 OR 最大位移: {cmp.shift_E.max():.6f} ; F 实测: {cmp.shift_F.max():.6f}")
    print(f"\n正对照 I12: {pos.I12.min():+.4f} – {pos.I12.max():+.4f} "
          f"(自对照 rho_s={pos[pos.trait_a == pos.trait_b].rho_s.iloc[0]:.3f})")
    print(f"甲状腺×偏头痛 I12: {ct.I12.min():+.4f} – {ct.I12.max():+.4f}")
    bng = ct[(ct.exposure == "BNG") & (ct.outcome == "Migraine")]["I12"].iloc[0]
    bngl = ct[(ct.exposure == "BNG_LOO") & (ct.outcome == "Migraine")]["I12"].iloc[0]
    print(f"\nBNG 天然对照 (Migraine): I12 含 FinnGen {bng:+.4f} vs 剔 FinnGen {bngl:+.4f} "
          f"(Δ={bng - bngl:+.4f})")
    print("输出文件: results/ldsc_*.csv, results/ldsc_*.png, results/forest_ldsc_corrected.png")


if __name__ == "__main__":
    main()
