"""
mrlap_prepare.py
----------------
把 ldsc_cache/*.tsv（按 1000G hm3 过滤后的全基因组 summary stats）
转换成 MRlap::MRlap() 需要的输入格式：

    rsid  chr  pos  alt  ref  beta  se  N

要点
  * MRlap 的 tidy_inputGWAS() 用 Z = beta/se、std_beta = Z/sqrt(N)，
    所以 N 必须是每个 SNP 的【总样本量】(README 明确: total, not effective)。
    暴露 (GBMI 五病) 的 per-SNP n 直接来自 GWAS Catalog harmonised 文件；
    结局 (FinnGen R12) 原文件无 N 列，用 FinnGen API 核实的病例+对照总数。
  * alt = 效应等位基因 (GWAS Catalog: effect_allele; FinnGen: alt),
    ref = 非效应等位基因。MRlap/GenomicSEM::munge 会按 hm3 三列参考再次对齐。
  * 同时输出 scale_factors.csv: median(se*sqrt(N))，用于把 MRlap 的
    标准化效应近似换算回 log OR 量级 (仅供量级对照，不作为正式结论)。

运行: python mrlap_prepare.py
"""
import os
import numpy as np
import pandas as pd

ROOT = r"D:/GWAS/Migraine"
CACHE = os.path.join(ROOT, "ldsc_cache")
OUTDIR = os.path.join(ROOT, "mrlap_input")
HM3 = r"D:/GWAS/MRlap_ref/hm3_eur_3col.txt"

EXPS = ["GD", "HYPO", "LT", "BNG", "TC", "BNG_LOO"]
OUTS = ["Migraine", "MigraineAura", "MigraineNoAura"]
OUT_N = {"Migraine": 395060, "MigraineAura": 379299, "MigraineNoAura": 376762}

COLS = ["SNP", "CHR", "BP", "EA", "OA", "BETA", "SE", "EAF", "N", "TAG"]
VALID_ALLELES = {"A", "C", "G", "T"}

os.makedirs(OUTDIR, exist_ok=True)
hm3 = pd.read_csv(HM3, sep="\t", dtype=str)
print("hm3 三列参考:", len(hm3), "行, 列:", list(hm3.columns))
print("  重复 SNP:", hm3.SNP.duplicated().sum())

rows = []
for name in EXPS + OUTS:
    p = os.path.join(CACHE, name + ".tsv")
    d = pd.read_csv(p, sep="\t", header=None, names=COLS,
                    dtype={"SNP": str, "EA": str, "OA": str, "TAG": str})
    n0 = len(d)
    for c in ["CHR", "BP", "BETA", "SE", "N"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")

    d = d[d.SNP.astype(str).str.startswith("rs")]
    d = d[d.CHR.between(1, 22) & (d.BP > 0)]
    d = d[d.EA.isin(VALID_ALLELES) & d.OA.isin(VALID_ALLELES)]
    d = d.dropna(subset=["BETA", "SE", "N"])
    d = d[(d.SE > 0) & np.isfinite(d.BETA) & np.isfinite(d.SE)]
    d = d.drop_duplicates("SNP", keep="first")

    if name in OUT_N:
        d["N"] = OUT_N[name]
    d = d[d.N > 0]

    out = pd.DataFrame({
        "rsid": d.SNP.values,
        "chr": d.CHR.astype(int).values,
        "pos": d.BP.astype(np.int64).values,
        "alt": d.EA.values,
        "ref": d.OA.values,
        "beta": d.BETA.values,
        "se": d.SE.values,
        "N": d.N.values,
    })
    fp = os.path.join(OUTDIR, name + ".tsv")
    out.to_csv(fp, sep="\t", index=False)

    # ---- 参考面板预检：与真实三列 hm3 的匹配率 / 等位基因方向 ----
    mg = out.merge(hm3, left_on="rsid", right_on="SNP", how="inner")
    same = (mg.alt == mg.A1) & (mg.ref == mg.A2)
    flip = (mg.alt == mg.A2) & (mg.ref == mg.A1)
    h = float(np.median(out.se * np.sqrt(out.N)))

    print(f"{name:<15s} raw={n0:>9,}  clean={len(out):>9,}  "
          f"N={out.N.median():>9,.0f}  hm3_intersect={len(mg):>9,} "
          f"(same={same.sum():,}, flip={flip.sum():,}, other={len(mg)-same.sum()-flip.sum():,})  "
          f"median(se*sqrt(N))={h:.3f}")

    rows.append(dict(trait=name, n_snp=len(out), N_median=float(out.N.median()),
                     hm3_intersect=len(mg), allele_same=int(same.sum()),
                     allele_flip=int(flip.sum()),
                     allele_other=int(len(mg) - same.sum() - flip.sum()),
                     se_sqrtN_median=h,
                     file=fp))

pd.DataFrame(rows).to_csv(os.path.join(OUTDIR, "scale_factors.csv"), index=False)
print("\n输入文件写入:", OUTDIR)
