# -*- coding: utf-8 -*-
"""
ukb_build.py —— 从 UK Biobank WGS 偏头痛 GWAS (GCST90473326) 构建下游文件

该文件无 rsID 列，ID 为 chr:pos:ea:oa。因此：
  1) 用 FinnGen R13（同为 GRCh38）建立 HapMap3 的 (chr,pos)->rsid 映射
  2) 工具变量本身带 chr/base_pair_location，直接用坐标匹配
输出：
  outcome_ukb.tsv                 rsid ref alt pval_out beta_out se_out eaf_out
  ldsc_cache/MigraineUKB.tsv      SNP CHR BP A1 A2 BETA SE EAF N TAG（无表头）
"""
import os
import gzip

ROOT = r"D:/GWAS/Migraine"
HM3 = r"D:/GWAS/LD core/hm3_no_MHC.list.txt"
FINNGEN = os.path.join(ROOT, "finngen_R13_G6_MIGRAINE.gz")
UKB = os.path.join(ROOT, "GCST90473326.tsv.gz")
UKB_N = 458440

print("[1] 读取 HapMap3 rsID", flush=True)
hm3 = set(open(HM3, encoding="utf-8").read().split())
print("    hm3:", len(hm3), flush=True)

print("[2] 从 FinnGen R13 建立 hm3 位点映射", flush=True)
pos2rs = {}
with gzip.open(FINNGEN, "rt", encoding="utf-8", errors="replace") as f:
    f.readline()
    for line in f:
        p = line.split("\t", 6)
        rs = p[4]
        if not rs:
            continue
        ids = rs.split(",")
        for r in ids:
            if r in hm3:
                pos2rs[(p[0], p[1])] = r
                break
print("    hm3 位点映射:", len(pos2rs), flush=True)

print("[3] 读取工具变量坐标", flush=True)
iv = {}
for k in ["GD", "HYPO", "LT", "BNG", "TC"]:
    with open(os.path.join(ROOT, f"iv_{k}.tsv"), encoding="utf-8") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ci, pi = hdr.index("chr"), hdr.index("base_pair_location")
        ri, ei, oi = hdr.index("rsid"), hdr.index("effect_allele"), hdr.index("other_allele")
        for line in f:
            p = line.rstrip("\n").split("\t")
            iv[(p[ci], p[pi])] = (p[ri], p[ei].upper(), p[oi].upper())
print("    工具变量位点数:", len(iv), flush=True)

print("[4] 扫描 UKB 文件", flush=True)
out_small = open(os.path.join(ROOT, "outcome_ukb.tsv"), "w", encoding="utf-8")
out_small.write("rsid\tref\talt\tpval_out\tbeta_out\tse_out\teaf_out\n")
out_cache = open(os.path.join(ROOT, "ldsc_cache", "MigraineUKB.tsv"), "w", encoding="utf-8")

n_line = n_iv = n_cache = n_amb = 0
with gzip.open(UKB, "rt", encoding="utf-8", errors="replace") as f:
    f.readline()
    for line in f:
        n_line += 1
        p = line.rstrip("\n").split("\t")
        if len(p) < 8:
            continue
        ch, bp, ea, oa, beta, se, eaf, pval = p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]
        if beta in ("NA", "") or se in ("NA", ""):
            continue
        key = (ch, bp)
        rec = iv.get(key)
        if rec is not None:
            rsid, iea, ioa = rec
            if {ea.upper(), oa.upper()} == {iea, ioa}:
                out_small.write("\t".join([rsid, oa, ea, pval, beta, se, eaf]) + "\n")
                n_iv += 1
        rs_hm3 = pos2rs.get(key)
        if rs_hm3:
            out_cache.write("\t".join([rs_hm3, ch, bp, ea, oa, beta, se, eaf, str(UKB_N), "OUT"]) + "\n")
            n_cache += 1

out_small.close()
out_cache.close()
print(f"    {n_line:,} 行 | 命中工具变量 {n_iv} | hm3 缓存 {n_cache:,}", flush=True)
print("DONE", flush=True)
