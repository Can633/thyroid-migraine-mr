# -*- coding: utf-8 -*-
"""
rev01_extract_ivs.py —— 反向 MR 第一步：从 FinnGen R12 偏头痛 GWAS 提取候选 IV
=============================================================================
输入：finngen_R13_G6_MIGRAINE{,_WITH_AURA,_NO_AURA}.gz（GRCh38，#chrom pos ref alt rsids ...）
输出：rev/iv_mig_<trait>.raw.tsv  列 = rsid chr pos ref alt pval beta se eaf
过滤：p < 5e-8、F=(beta/se)^2 > 10、MAF 可用（af 在 0.01–0.99）
说明：FinnGen 的 beta 以 alt 为效应等位基因，下游统一按 A1=alt 处理。
"""
import gzip
import os
import sys

ROOT = "D:/GWAS/Migraine"
OUT = os.path.join(ROOT, "rev")
os.makedirs(OUT, exist_ok=True)

TRAITS = {
    "Migraine": "finngen_R13_G6_MIGRAINE.gz",
    "MigraineAura": "finngen_R13_G6_MIGRAINE_WITH_AURA.gz",
    "MigraineNoAura": "finngen_R13_G6_MIGRAINE_NO_AURA.gz",
}
P_CUT = 5e-8
F_CUT = 10.0


def run(trait, fn):
    path = os.path.join(ROOT, fn)
    rows = []
    n = 0
    with gzip.open(path, "rt") as fh:
        header = fh.readline()
        for line in fh:
            n += 1
            p = line.split("\t", 8)
            # 0 #chrom 1 pos 2 ref 3 alt 4 rsids 5 gene 6 pval 7 mlogp 8 其余
            try:
                pval = float(p[6])
            except ValueError:
                continue
            if pval >= P_CUT:
                continue
            rest = p[8].split("\t")
            try:
                beta, se, af = float(rest[0]), float(rest[1]), float(rest[2])
            except (ValueError, IndexError):
                continue
            if se <= 0 or (beta / se) ** 2 <= F_CUT:
                continue
            if not (0.01 <= af <= 0.99):
                continue
            rsid = p[4].split(",")[0].strip()
            if not rsid.startswith("rs"):
                continue
            rows.append((rsid, p[0], p[1], p[2].upper(), p[3].upper(),
                         pval, beta, se, af))
    outp = os.path.join(OUT, f"iv_mig_{trait}.raw.tsv")
    with open(outp, "w") as fh:
        fh.write("rsid\tchr\tpos\tref\talt\tpval\tbeta\tse\teaf\n")
        for r in rows:
            fh.write("%s\t%s\t%s\t%s\t%s\t%.6g\t%.6g\t%.6g\t%.6g\n" % r)
    print(f"{trait}: scanned {n:,} rows -> {len(rows)} candidate SNPs -> {outp}", flush=True)
    return len(rows)


if __name__ == "__main__":
    total = 0
    for t, f in TRAITS.items():
        total += run(t, f)
    print("total candidates:", total)
