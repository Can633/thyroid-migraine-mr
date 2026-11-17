# -*- coding: utf-8 -*-
"""
r13_build.py -- 用 FinnGen R13 替换结局数据（单次流式解压，同时产出两套下游文件）

对每个 R13 文件做一遍扫描，输出：
  1. outcome_<code>.tsv         -- 仅含 iv_rsids.txt 的 693 个 rsID，供 TwoSampleMR 使用
  2. ldsc_cache/<Trait>.tsv     -- HapMap3(no MHC) 过滤后的全基因组缓存，供 LDSC / coloc 使用
     （格式与 R12 缓存一致：SNP CHR BP A1 A2 BETA SE EAF N TAG，无表头）

R13 列（与 R12 相同）：
  0 #chrom 1 pos 2 ref 3 alt 4 rsids 5 nearest_genes 6 pval 7 mlogp
  8 beta 9 sebeta 10 af_alt 11 af_alt_cases 12 af_alt_controls
"""
import os
import gzip
import time

ROOT = r"D:/GWAS/Migraine"
HM3 = r"D:/GWAS/LD core/hm3_no_MHC.list.txt"

FILES = {
    "G6_MIGRAINE": "Migraine",
    "G6_MIGRAINE_WITH_AURA": "MigraineAura",
    "G6_MIGRAINE_NO_AURA": "MigraineNoAura",
}

IVS = set(open(os.path.join(ROOT, "iv_rsids.txt"), encoding="utf-8").read().split())
print(f"IV rsIDs: {len(IVS)}", flush=True)

hm3 = set(open(HM3, encoding="utf-8").read().split())
print(f"HapMap3 (no MHC): {len(hm3)}", flush=True)

CACHE = os.path.join(ROOT, "ldsc_cache")


def build(code, trait):
    src = os.path.join(ROOT, f"finngen_R13_{code}.gz")
    out_small = os.path.join(ROOT, f"outcome_{code}.tsv")
    out_cache = os.path.join(CACHE, f"{trait}.tsv")

    n_line = n_small = n_cache = n_bad = 0
    t0 = time.time()
    with gzip.open(src, "rt", encoding="utf-8", errors="replace") as f, \
            open(out_small, "w", encoding="utf-8") as fs, \
            open(out_cache, "w", encoding="utf-8") as fc:
        fs.write("rsid\tref\talt\tpval_out\tbeta_out\tse_out\teaf_out\n")
        header = f.readline().rstrip("\n").split("\t")
        if header[0].lstrip("#") != "chrom":
            raise SystemExit(f"unexpected header in {src}: {header[:4]}")
        for line in f:
            n_line += 1
            p = line.rstrip("\n").split("\t")
            if len(p) < 11:
                n_bad += 1
                continue
            rs_field = p[4]
            if not rs_field:
                continue
            ref, alt = p[2], p[3]
            beta, se, af = p[8], p[9], p[10]
            pval = p[6]
            if beta in ("NA", "") or se in ("NA", "") or af in ("NA", ""):
                n_bad += 1
                continue

            hit_iv = ""
            hit_hm3 = ""
            ids = rs_field.split(",") if "," in rs_field else (rs_field,)
            for rs in ids:
                if not hit_iv and rs in IVS:
                    hit_iv = rs
                if not hit_hm3 and rs in hm3:
                    hit_hm3 = rs
                if hit_iv and hit_hm3:
                    break

            if hit_iv:
                fs.write("\t".join([hit_iv, ref, alt, pval, beta, se, af]) + "\n")
                n_small += 1
            if hit_hm3:
                # SNP CHR BP A1 A2 BETA SE EAF N TAG  (A1=alt, A2=ref，与 R12 一致)
                fc.write("\t".join([hit_hm3, p[0], p[1], alt, ref, beta, se, af, "0", "OUT"]) + "\n")
                n_cache += 1

    print(f"{code:26s} -> {trait:15s} lines={n_line:,} IV={n_small} hm3={n_cache:,} "
          f"skipped={n_bad} ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    for code, trait in FILES.items():
        build(code, trait)
    print("DONE", flush=True)
