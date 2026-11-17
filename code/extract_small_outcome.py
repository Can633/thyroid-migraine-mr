#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从南亚 Genes&Health 偏头痛 GWAS (GCST90727115, 外显子测序) 中提取 5 个甲状腺暴露的
工具变量，并按位置 (chr,pos) 跨种族对齐等位基因，输出每个暴露的 harmonized 文件。

关键约束:
- 该结局 variant_id 为 chr_pos_ref_alt 形式, 无 rsID -> 用位置匹配
- 暴露 IV 来自欧洲 GBMI GWAS, 含 chr/base_pair_location/alleles
- 跨种族: 仅做位置匹配 + 等位基因对齐, 不假设 LD 一致
"""
import csv, os, sys

OUTCOME = "D:/GWAS/Migraine/GCST90727115.tsv"
EXPOS   = {
    "GD":   "D:/GWAS/Migraine/iv_GD.tsv",
    "LT":   "D:/GWAS/Migraine/iv_LT.tsv",
    "BNG":  "D:/GWAS/Migraine/iv_BNG.tsv",
    "HYPO": "D:/GWAS/Migraine/iv_HYPO.tsv",
    "TC":   "D:/GWAS/Migraine/iv_TC.tsv",
}
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

def comp(a):
    return COMP.get(a.upper(), a.upper())

# ---- 读入各暴露 IV, 以 (chr,pos) 为键 ----
expos_map = {}          # (chr,pos) -> list of (exp_tag, exp_info)
for tag, path in EXPOS.items():
    with open(path, encoding="utf-8") as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            try:
                key = (int(r["chr"]), int(r["base_pair_location"]))
            except Exception:
                continue
            expos_map.setdefault(key, []).append({
                "tag": tag,
                "rsid": r["rsid"],
                "A1": r["effect_allele"].upper(),   # exposure effect
                "A2": r["other_allele"].upper(),     # exposure other
                "b_exp": float(r["beta"]),
                "se_exp": float(r["se"]),
                "eaf_exp": float(r["eaf"]),
            })

# ---- 流式扫描结局, 匹配并对齐 ----
out_files = {t: open(f"D:/GWAS/Migraine/harm_small_{t}.csv", "w", newline="", encoding="utf-8") for t in EXPOS}
writers = {t: csv.writer(out_files[t]) for t in EXPOS}
for t in EXPOS:
    writers[t].writerow(["rsid","chr","pos","A1_exp","A2_exp","b_exp","se_exp","eaf_exp",
                         "A1_out","A2_out","b_out","se_out","eaf_out","p_out","orient"])

total = 0
matched = 0
palindromic = 0
with open(OUTCOME, encoding="utf-8") as f:
    header = f.readline().rstrip("\n").split("\t")
    # chromosome base_pair_location effect_allele other_allele beta standard_error effect_allele_frequency p_value variant_id n TEST CHISQ EXTRA
    idx = {h: i for i, h in enumerate(header)}
    for line in f:
        total += 1
        if total % 5000000 == 0:
            sys.stderr.write(f"  scanned {total} lines, matched {matched}\n")
        p = line.rstrip("\n").split("\t")
        if len(p) <= idx["p_value"]:
            continue
        try:
            chr_i = int(p[idx["chromosome"]])
            pos_i = int(p[idx["base_pair_location"]])
        except Exception:
            continue
        key = (chr_i, pos_i)
        if key not in expos_map:
            continue
        O1 = p[idx["effect_allele"]].upper()      # outcome effect
        O2 = p[idx["other_allele"]].upper()       # outcome other
        try:
            b_out = float(p[idx["beta"]]); se_out = float(p[idx["standard_error"]])
            eaf_out = float(p[idx["effect_allele_frequency"]]); pv = float(p[idx["p_value"]])
        except Exception:
            continue
        for e in expos_map[key]:
            A1, A2 = e["A1"], e["A2"]
            # 排除回文(模糊) SNP
            if A1 == comp(A2) or O1 == comp(O2):
                palindromic += 1
                continue
            use_b = None; orient = None
            if O1 == A1 and O2 == A2:
                use_b, orient = b_out, "+"
            elif O1 == A2 and O2 == A1:
                use_b, orient = -b_out, "-"
            elif O1 == comp(A1) and O2 == comp(A2):
                use_b, orient = b_out, "comp+"      # 链翻转不改变 beta
            elif O1 == comp(A2) and O2 == comp(A1):
                use_b, orient = -b_out, "comp-"
            else:
                continue   # 无法对齐 -> 丢弃
            writers[e["tag"]].writerow([
                e["rsid"], chr_i, pos_i, A1, A2, f"{e['b_exp']:.6g}", f"{e['se_exp']:.6g}", f"{e['eaf_exp']:.5f}",
                O1, O2, f"{use_b:.6g}", f"{se_out:.6g}", f"{eaf_out:.5f}", f"{pv:.6g}", orient])
            matched += 1

for t in EXPOS:
    out_files[t].close()

# ---- 报告每个暴露匹配数 ----
sys.stderr.write(f"TOTAL outcome lines scanned: {total}\n")
sys.stderr.write(f"TOTAL matched SNP-instrument rows: {matched} (palindromic skipped: {palindromic})\n")
for t in EXPOS:
    n = sum(1 for _ in open(f"D:/GWAS/Migraine/harm_small_{t}.csv")) - 1
    sys.stderr.write(f"  {t}: {n} harmonized instruments -> harm_small_{t}.csv\n")
