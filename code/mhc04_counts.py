# -*- coding: utf-8 -*-
"""
mhc04_counts.py —— 全基因组显著信号的分布（MHC vs 非 MHC）
=========================================================
为什么不能用 ldsc_cache 计数
----------------------------
缓存只含 HapMap3 SNP，且 FinnGen R12 的 rsID 与 HapMap3 存在大量版本不匹配
（实测：偏头痛 381 个显著 SNP 中仅 65 个能对上 hm3），用它作分母会**系统性
低估**非 MHC 区的信号数，从而夸大"MHC 占比"。

本脚本直接扫原始文件，按位置切分 chr6:25–34 Mb（LDSC 惯例窗口）与非 MHC 区，
统计 p<5e-8 的 SNP 数，并在窗口内再报 chr6:28–34 Mb（GRCh38 MHC 本体附近）。

输出：results/mhc_gws_counts.csv
"""
import gzip
import os
from concurrent.futures import ProcessPoolExecutor

ROOT = "D:/GWAS/Migraine"
RES = os.path.join(ROOT, "results")
THYR = "D:/GWAS/TH+SVD/thyroid/harmonised"
P_CUT = 5e-8
LO, HI = 25_000_000, 34_000_000
SLO, SHI = 28_000_000, 34_000_000

THY = {"GD": "GCST90627744", "HYPO": "GCST90627749", "LT": "GCST90627754",
       "BNG": "GCST90627740", "TC": "GCST90627763"}
MIG = {"Migraine": "finngen_R13_G6_MIGRAINE.gz",
       "MigraineAura": "finngen_R13_G6_MIGRAINE_WITH_AURA.gz",
       "MigraineNoAura": "finngen_R13_G6_MIGRAINE_NO_AURA.gz"}


def scan_thy(item):
    name, acc = item
    pcol, ccol, poscol = 7, 0, 1
    path = os.path.join(THYR, acc, f"{acc}.h.tsv.gz")
    tot = mhc = strict = nomhc = 0
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {k: i for i, k in enumerate(header)}
        ccol, poscol, pcol = idx["chromosome"], idx["base_pair_location"], idx["p_value"]
        for line in fh:
            p = line.split("\t")
            if len(p) <= max(ccol, poscol, pcol):
                continue
            try:
                if float(p[pcol]) >= P_CUT:
                    continue
                pos = int(p[poscol])
            except ValueError:
                continue
            tot += 1
            if p[ccol] == "6" and LO <= pos <= HI:
                mhc += 1
                if SLO <= pos <= SHI:
                    strict += 1
            else:
                nomhc += 1
    return name, dict(trait=name, n_gws_total=tot, n_gws_mhc=mhc,
                      n_gws_mhc_strict=strict, n_gws_nomhc=nomhc,
                      frac_mhc=(mhc / tot if tot else 0.0))


def scan_mig(item):
    name, fn = item
    path = os.path.join(ROOT, fn)
    tot = mhc = strict = nomhc = 0
    with gzip.open(path, "rt") as fh:
        fh.readline()
        for line in fh:
            p = line.split("\t")
            if len(p) < 8:
                continue
            try:
                if float(p[6]) >= P_CUT:
                    continue
                pos = int(p[1])
            except ValueError:
                continue
            tot += 1
            if p[0] == "6" and LO <= pos <= HI:
                mhc += 1
                if SLO <= pos <= SHI:
                    strict += 1
            else:
                nomhc += 1
    return name, dict(trait=name, n_gws_total=tot, n_gws_mhc=mhc,
                      n_gws_mhc_strict=strict, n_gws_nomhc=nomhc,
                      frac_mhc=(mhc / tot if tot else 0.0))


if __name__ == "__main__":
    import pandas as pd
    rows = []
    with ProcessPoolExecutor(max_workers=4) as ex:
        futs = [ex.submit(scan_thy, (n, a)) for n, a in THY.items()] + \
               [ex.submit(scan_mig, (n, f)) for n, f in MIG.items()]
        for f in futs:
            name, d = f.result()
            rows.append(d)
            print(f"{name:15s} total={d['n_gws_total']:8,d}  MHC={d['n_gws_mhc']:7,d} "
                  f"(strict {d['n_gws_mhc_strict']:7,d})  non-MHC={d['n_gws_nomhc']:8,d}  "
                  f"frac={d['frac_mhc']:.3f}", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "mhc_gws_counts.csv"), index=False)
    print("\nsaved -> results/mhc_gws_counts.csv")
