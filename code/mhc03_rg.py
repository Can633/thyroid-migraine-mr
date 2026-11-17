# -*- coding: utf-8 -*-
"""
mhc03_rg.py —— MHC 敏感性分析第三步：把 MHC 加回 LD score 回归，重估 r_g
=========================================================================
主分析（ldsc_rg_precise.py）用的 HapMap3 参考面板 **去除了 MHC**
（hm3_no_MHC.list.txt），因此 GD 这类遗传信号高度集中于 HLA 的性状，
其 r_g 只反映非 MHC 的多基因成分。本脚本做 2×2 敏感性：

  变体 A  无 MHC + 标准权重（= 主分析，用于复现基准）
  变体 B  无 MHC + 均匀权重（W2≡1，隔离"权重形式"本身的影响）
  变体 C  含 MHC + MHC SNP 用自身 L2 作权重（近似 down-weight 高 LD 区）
  变体 D  含 MHC + MHC SNP 用 W2=1（不惩罚，等价于让 MHC 全权重参与）

读 A/B 与 C/D 的差，可把"加入 MHC"的效应与"换权重"的效应分开。

权重说明（必须写进方法学限制）
------------------------------
官方 1000G Phase3 weights 文件本身就是 hm3_noMHC，MHC SNP 没有对应的
w-ld，因此 C/D 只能用近似权重。MHC 内长程 LD 也违反 LDSC 的独立性假设，
所以含 MHC 的结果只能作为**敏感性参考**，不能替代主分析。

输出：results/mhc_rg_sensitivity.csv
"""
import math
import os

import numpy as np
import pandas as pd

import ldsc_overlap as L
from ldsc_rg_precise import fit_with_blocks

ROOT = "D:/GWAS/Migraine"
MHC = os.path.join(ROOT, "mhc")
RES = os.path.join(ROOT, "results")
LDDIR = "D:/GWAS/LD core/extracted/LDscore"
LO, HI = 25_000_000, 34_000_000
EXPS = ["GD", "HYPO", "LT", "BNG", "TC"]
OUT = "Migraine"


def load_mhc_ref():
    d = pd.read_csv(os.path.join(LDDIR, "LDscore.6.l2.ldscore.gz"),
                    sep="\t", usecols=["CHR", "SNP", "BP", "L2"])
    return d[(d.BP >= LO) & (d.BP <= HI)].copy()


def mhc_gwas(name):
    p = os.path.join(MHC, f"mhc_{name}.tsv")
    d = pd.read_csv(p, sep="\t", dtype={"a1": str, "a2": str})
    d = d.dropna(subset=["beta", "se"])
    d = d[(d.se > 0) & np.isfinite(d.beta)].drop_duplicates("rsid")
    d["z"] = d.beta / d.se
    d = d.rename(columns={"rsid": "SNP", "a1": "EA", "a2": "OA",
                          "pos": "BP", "chr": "CHR"})
    return d.set_index("SNP")


def rg_fit(gx, gy, ref_l2, N1, N2, chi2_cut=100.0):
    m = L.align_z(gx, gy, ref_l2)
    m = m[(m.chi_x <= chi2_cut) & (m.chi_y <= chi2_cut)]
    if len(m) < 100:
        return None
    L2 = m["L2"].to_numpy(float)
    W2 = m["W2"].to_numpy(float)
    zx = m["zx"].to_numpy(float)
    zy = m["zy"].to_numpy(float)
    chr_col = "CHR_x" if "CHR_x" in m.columns else "CHR"
    bp_col = "BP_x" if "BP_x" in m.columns else "BP"
    n_mhc = int(((m[chr_col] == 6) & (m[bp_col] >= LO) & (m[bp_col] <= HI)).sum())

    fx, ex, _ = fit_with_blocks(L2, W2, zx ** 2)
    fy, ey, _ = fit_with_blocks(L2, W2, zy ** 2)
    fc, ec, _ = fit_with_blocks(L2, W2, zx * zy)
    s1, s2, s12 = float(fx[1]), float(fy[1]), float(fc[1])
    if s1 <= 0 or s2 <= 0:
        return None
    rg = s12 / math.sqrt(s1 * s2)
    rgk = np.array([ec[k, 1] / math.sqrt(ex[k, 1] * ey[k, 1])
                    for k in range(len(ec)) if ex[k, 1] > 0 and ey[k, 1] > 0])
    se = (math.sqrt((len(rgk) - 1) / len(rgk) *
                    float(np.sum((rgk - rgk.mean()) ** 2))) if len(rgk) > 10 else float("nan"))
    return dict(rg=rg, se=se, z=(rg / se if se > 0 else float("nan")),
                n_snp=len(m), n_mhc_snp=n_mhc,
                h2_x=s1 * L.M_TOTAL / N1, h2_y=s2 * L.M_TOTAL / N2)


def main():
    print("[1/3] 载入参考与 GWAS ...", flush=True)
    ref = L.load_reference()
    ref_base = ref[["SNP", "L2", "W2"]].copy()

    mhc_ref = load_mhc_ref()
    mhc_ref = mhc_ref[~mhc_ref.SNP.isin(ref_base.SNP)].copy()
    print(f"      MHC 窗口内 LD score SNP: {len(mhc_ref)}", flush=True)

    gy = mhc_gwas(OUT)
    base_y = L.load_gwas(OUT)
    N2 = L.OUT_N[OUT]
    print(f"      {OUT}: MHC 窗口 GWAS SNP {len(gy)}", flush=True)

    rows = []
    for e in EXPS:
        gx = mhc_gwas(e)
        base_x = L.load_gwas(e)
        N1 = float(base_x["N"].median())

        gx_in = pd.concat([base_x, gx])
        gy_in = pd.concat([base_y, gy])
        variants = [
            # label, 加MHC, 均匀权重, chi2 截断
            ("A_noMHC_std",       False, False, 100.0),
            ("B_noMHC_unif",      False, True,  100.0),
            ("C_MHC_selfW",       True,  False, 100.0),
            ("D_MHC_unitW",       True,  True,  100.0),
            ("E_noMHC_nocut",     False, False, 1e9),
            ("F_MHC_selfW_nocut", True,  False, 1e9),
        ]
        for label, add_mhc, uniform, chi2 in variants:
            gx_ = gx_in if add_mhc else base_x
            gy_ = gy_in if add_mhc else base_y
            r_in = ref_base
            if add_mhc:
                add = mhc_ref[["SNP", "L2"]].copy()
                add["W2"] = add["L2"] if label.endswith("selfW") or "selfW" in label else 1.0
                r_in = pd.concat([ref_base, add], ignore_index=True)
            if uniform:
                r_in = r_in.copy(); r_in["W2"] = 1.0
            r = rg_fit(gx_, gy_, r_in, N1, N2, chi2_cut=chi2)
            if r is None:
                print(f"      {e} {label}: fit failed", flush=True); continue
            rows.append(dict(exposure=e, outcome=OUT, variant=label, **r))
            print(f"      {e:<5s} {label:<13s} rg={r['rg']:+.4f} (se {r['se']:.4f}) "
                  f"z={r['z']:+.2f}  nSNP={r['n_snp']:,} (MHC {r['n_mhc_snp']}) "
                  f"h2x={r['h2_x']:.4f}", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "mhc_rg_sensitivity.csv"), index=False)
    print("\nsaved -> results/mhc_rg_sensitivity.csv")


if __name__ == "__main__":
    main()
