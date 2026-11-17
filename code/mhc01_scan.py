# -*- coding: utf-8 -*-
"""
mhc01_scan.py —— MHC 敏感性分析第一步：抽取 6 号染色体 MHC 区间的全部 SNP
==========================================================================
为什么需要单独扫
----------------
主分析用的 ldsc_cache/*.tsv 只保留 HapMap3 **去 MHC** 的 SNP，MHC 区被整体
丢弃，因此无法直接回答"MHC 是否影响 r_g"。

本脚本直接从 8 个原始 GWAS 中抽取 chr6:25–34 Mb（LDSC 惯例窗口，GRCh37 坐标；
对 GRCh38 文件为近似，MHC 本体约 28.5–33.5 Mb，故同时记录 GRCh38 严格区间）
的全部记录，供后续两用：
  (1) 描述性：各性状在 MHC 内的全基因组显著信号有多少、是否与偏头痛同向；
  (2) 回归：把 MHC SNP 加回 LD score 回归，重估 r_g。

输出：mhc/mhc_<trait>.tsv  列 rsid chr pos a1 a2 beta se eaf p n build
"""
import glob
import gzip
import os
import sys
from concurrent.futures import ProcessPoolExecutor

ROOT = "D:/GWAS/Migraine"
OUT = os.path.join(ROOT, "mhc")
os.makedirs(OUT, exist_ok=True)
THYR = "D:/GWAS/TH+SVD/thyroid/harmonised"

LO, HI = 25_000_000, 34_000_000

THY = {"GD": "GCST90627744", "HYPO": "GCST90627749", "LT": "GCST90627754",
       "BNG": "GCST90627740", "TC": "GCST90627763"}
MIG = {"Migraine": "finngen_R13_G6_MIGRAINE.gz",
       "MigraineAura": "finngen_R13_G6_MIGRAINE_WITH_AURA.gz",
       "MigraineNoAura": "finngen_R13_G6_MIGRAINE_NO_AURA.gz"}


def scan_thy(item):
    name, acc = item
    path = os.path.join(THYR, acc, f"{acc}.h.tsv.gz")
    rows = []
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {k: i for i, k in enumerate(header)}
        ci = idx.get("rsid", 9)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p[0] != "6":
                continue
            try:
                pos = int(p[1])
            except ValueError:
                continue
            if not (LO <= pos <= HI):
                continue
            if len(p) <= ci or not p[ci].startswith("rs"):
                continue
            try:
                rows.append((p[ci], p[0], pos, p[2].upper(), p[3].upper(),
                             float(p[4]), float(p[5]), float(p[6]), float(p[7]),
                             float(p[10]) if len(p) > 10 and p[10] not in ("", "NA") else float("nan")))
            except (ValueError, IndexError):
                continue
    op = os.path.join(OUT, f"mhc_{name}.tsv")
    with open(op, "w") as fh:
        fh.write("rsid\tchr\tpos\ta1\ta2\tbeta\tse\teaf\tp\tn\tbuild\n")
        for r in rows:
            fh.write("%s\t%s\t%d\t%s\t%s\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\tGRCh38\n" % r)
    return name, len(rows), op


def scan_mig(item):
    name, fn = item
    path = os.path.join(ROOT, fn)
    rows = []
    with gzip.open(path, "rt") as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p[0] != "6":
                continue
            try:
                pos = int(p[1])
            except ValueError:
                continue
            if not (LO <= pos <= HI):
                continue
            rsid = p[4].split(",")[0].strip()
            if not rsid.startswith("rs"):
                continue
            try:
                rows.append((rsid, p[0], pos, p[3].upper(), p[2].upper(),
                             float(p[8]), float(p[9]), float(p[10]), float(p[6])))
            except (ValueError, IndexError):
                continue
    op = os.path.join(OUT, f"mhc_{name}.tsv")
    with open(op, "w") as fh:
        fh.write("rsid\tchr\tpos\ta1\ta2\tbeta\tse\teaf\tp\tn\tbuild\n")
        for r in rows:
            fh.write("%s\t%s\t%d\t%s\t%s\t%.6g\t%.6g\t%.6g\t%.6g\t\tGRCh38\n" % r)
    return name, len(rows), op


if __name__ == "__main__":
    jobs = [("thy", n, a) for n, a in THY.items()] + [("mig", n, f) for n, f in MIG.items()]
    with ProcessPoolExecutor(max_workers=4) as ex:
        futs = []
        for kind, n, v in jobs:
            if kind == "thy":
                futs.append(ex.submit(scan_thy, (n, v)))
            else:
                futs.append(ex.submit(scan_mig, (n, v)))
        for f in futs:
            name, k, op = f.result()
            print(f"{name}: {k} SNPs in chr6:{LO/1e6:.0f}-{HI/1e6:.0f}Mb -> {op}", flush=True)
