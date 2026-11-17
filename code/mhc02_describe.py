# -*- coding: utf-8 -*-
"""
mhc02_describe.py —— MHC 敏感性分析第二步：描述性证据
=====================================================
回答三个问题，都是"加不加 MHC"能不能改变结论的前提：

 Q1 各性状的全基因组显著信号有多少落在 MHC？MHC 是不是它的主要信号源？
 Q2 甲状腺性状在 MHC 内的显著信号，在偏头痛里有信号吗？方向是否一致？
 Q3 整个 MHC 区间内，两侧 z 的相关性有多大（粗略的"区域遗传重叠"）？

输出：results/mhc_signal_summary.csv     每性状的 MHC 信号占比
      results/mhc_shared_snps.csv        MHC 内两侧都有强信号的 SNP
      results/mhc_z_correlation.csv      MHC 区间 z 相关（含分层：严格 MHC / 扩展）
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = "D:/GWAS/Migraine"
MHC = os.path.join(ROOT, "mhc")
RES = os.path.join(ROOT, "results")
CACHE = os.path.join(ROOT, "ldsc_cache")
LO, HI = 25_000_000, 34_000_000
MHC_STRICT = (28_000_000, 34_000_000)   # GRCh38 MHC 本体约 28.5–33.5 Mb，取宽一点的 28–34
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
THY = ["GD", "HYPO", "LT", "BNG", "TC"]
MIG = ["Migraine", "MigraineAura", "MigraineNoAura"]


def load_mhc(name):
    d = pd.read_csv(os.path.join(MHC, f"mhc_{name}.tsv"), sep="\t",
                    dtype={"a1": str, "a2": str})
    d = d.dropna(subset=["beta", "se"])
    d = d[(d.se > 0) & np.isfinite(d.beta)]
    d["z"] = d.beta / d.se
    d["p_calc"] = 2 * norm.sf(np.abs(d.z))
    d = d.drop_duplicates("rsid")
    return d


def count_gws_nomhc(name):
    """非 MHC 区（hm3 缓存）的 GWS 个数，用于算 MHC 占比。"""
    p = os.path.join(CACHE, f"{name}.tsv")
    d = pd.read_csv(p, sep="\t", header=None,
                    names=["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "EAF", "N", "TAG"],
                    usecols=[5, 6], dtype={5: float, 6: float})
    z = (d.iloc[:, 0] / d.iloc[:, 1]).to_numpy()
    z = z[np.isfinite(z)]
    return int(np.sum(2 * norm.sf(np.abs(z)) < 5e-8))


def align(a, b):
    """a = 甲状腺（保持方向），b = 偏头痛（对齐到 a 的效应等位基因）。"""
    m = a.merge(b, on="rsid", suffixes=("_t", "_m"))
    same = (m.a1_t == m.a1_m) & (m.a2_t == m.a2_m)
    flip = (m.a1_t == m.a2_m) & (m.a2_t == m.a1_m)
    ct, cot = m.a1_t.map(COMP), m.a2_t.map(COMP)
    same_c = (ct == m.a1_m) & (cot == m.a2_m)
    flip_c = (ct == m.a2_m) & (cot == m.a1_m)
    sign = np.where(same | same_c, 1.0, np.where(flip | flip_c, -1.0, 0.0))
    m = m[sign != 0].copy()
    m["z_mig_aligned"] = m.z_m.to_numpy() * sign[sign != 0]
    return m


def main():
    rows, shared_all, corr_rows = [], [], []
    for t in THY:
        d = load_mhc(t)
        gws = d[d.p_calc < 5e-8]
        n_nomhc = count_gws_nomhc(t)
        n_mhc = len(gws)
        top = gws.nsmallest(1, "p")
        rows.append(dict(
            trait=t, n_gws_nomhc=n_nomhc, n_gws_mhc=n_mhc,
            frac_gws_mhc=(n_mhc / (n_mhc + n_nomhc) if (n_mhc + n_nomhc) else np.nan),
            top_mhc_rsid=(top.rsid.iloc[0] if len(top) else ""),
            top_mhc_pos=(int(top.pos.iloc[0]) if len(top) else np.nan),
            top_mhc_p=(float(top.p.iloc[0]) if len(top) else np.nan),
            mhc_h2_like=float(np.sum(gws.z ** 2)) if n_mhc else 0.0))

        for o in MIG:
            g = load_mhc(o)
            m = align(d, g)
            # Q2：甲状腺 GWS（MHC 内）对应的偏头痛 z
            if len(gws):
                sub = m[m.rsid.isin(gws.rsid)]
                if len(sub):
                    shared_all.append(dict(
                        thyroid=t, outcome=o, n_thy_mhc_gws=len(sub),
                        n_mig_z_gt2=int((sub.z_mig_aligned.abs() > 2).sum()),
                        n_same_sign=int((np.sign(sub.z_mig_aligned) == np.sign(sub.z_t)).sum()),
                        mean_z_mig=float(sub.z_mig_aligned.mean()),
                        min_p_mig=float(2 * norm.sf(sub.z_mig_aligned.abs().max())),
                        top_rsid=sub.loc[sub.z_mig_aligned.abs().idxmax(), "rsid"]))
            # Q3：整个窗口的 z 相关
            for lbl, (a, b) in {"extended_25_34": (LO, HI),
                                "strict_28_34": MHC_STRICT}.items():
                s = m[(m.pos_t >= a) & (m.pos_t <= b)]
                s = s[np.isfinite(s.z_t) & np.isfinite(s.z_mig_aligned)]
                if len(s) > 500:
                    corr_rows.append(dict(
                        thyroid=t, outcome=o, window=lbl, n_snp=len(s),
                        r=np.corrcoef(s.z_t, s.z_mig_aligned)[0, 1],
                        r_strong=(np.corrcoef(
                            s[(s.z_t.abs() > 2) & (s.z_mig_aligned.abs() > 1)].z_t,
                            s[(s.z_t.abs() > 2) & (s.z_mig_aligned.abs() > 1)].z_mig_aligned)[0, 1]
                            if ((s.z_t.abs() > 2) & (s.z_mig_aligned.abs() > 1)).sum() > 30
                            else np.nan)))

    # 偏头痛自身的 MHC 信号（用于对照）
    for o in MIG:
        g = load_mhc(o)
        gws = g[g.p_calc < 5e-8]
        top = gws.nsmallest(1, "p")
        rows.append(dict(
            trait=o, n_gws_nomhc=count_gws_nomhc(o), n_gws_mhc=len(gws),
            frac_gws_mhc=(len(gws) / (len(gws) + count_gws_nomhc(o))
                          if (len(gws) + count_gws_nomhc(o)) else np.nan),
            top_mhc_rsid=(top.rsid.iloc[0] if len(top) else ""),
            top_mhc_pos=(int(top.pos.iloc[0]) if len(top) else np.nan),
            top_mhc_p=(float(top.p.iloc[0]) if len(top) else np.nan),
            mhc_h2_like=float(np.sum(gws.z ** 2)) if len(gws) else 0.0))

    s1 = pd.DataFrame(rows)
    s1.to_csv(os.path.join(RES, "mhc_signal_summary.csv"), index=False)
    s2 = pd.DataFrame(shared_all)
    s2.to_csv(os.path.join(RES, "mhc_shared_snps.csv"), index=False)
    s3 = pd.DataFrame(corr_rows)
    s3.to_csv(os.path.join(RES, "mhc_z_correlation.csv"), index=False)

    pd.set_option("display.width", 200)
    print("=== Q1 各性状 MHC 内 GWS 信号 ===\n")
    print(s1.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    print("\n=== Q2 甲状腺 MHC 信号在偏头痛中的表现 ===\n")
    print(s2.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print("\n=== Q3 MHC 区间 z 相关 ===\n")
    print(s3.to_string(index=False, float_format=lambda x: f"{x:.4f}"))


if __name__ == "__main__":
    main()
