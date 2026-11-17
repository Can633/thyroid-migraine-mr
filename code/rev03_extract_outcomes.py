# -*- coding: utf-8 -*-
"""
rev03_extract_outcomes.py —— 反向 MR 第二步：从 5 个甲状腺病 GWAS 取结局 SNP
============================================================================
输入：rev/iv_mig_*.raw.tsv 的 rsID 并集（含 clumping 后可能被保留的全部候选）
      甲状腺原始文件 D:/GWAS/TH+SVD/thyroid/harmonised/<ACC>/<ACC>.h.tsv.gz
      列约定（GWAS-SSF harmonised）：
        1 chromosome, 2 base_pair_location, 3 effect_allele, 4 other_allele,
        5 beta, 6 standard_error, 7 effect_allele_frequency, 8 p_value,
        9 variant_id, 10 rsid, 11 n
输出：rev/out_<NAME>.tsv
      列 = rsid chr pos ea oa beta se eaf p n
说明：5 个文件并行流式扫描（不整表读入），每文件约 1–3 分钟。
"""
import glob
import gzip
import os
from concurrent.futures import ProcessPoolExecutor

ROOT = "D:/GWAS/Migraine"
REV = os.path.join(ROOT, "rev")
THYR = "D:/GWAS/TH+SVD/thyroid/harmonised"

ACCS = {
    "GD": "GCST90627744",
    "HYPO": "GCST90627749",
    "LT": "GCST90627754",
    "BNG": "GCST90627740",
    "TC": "GCST90627763",
}


def load_rsids():
    s = set()
    for f in glob.glob(os.path.join(REV, "iv_mig_*.raw.tsv")):
        with open(f) as fh:
            next(fh)
            for line in fh:
                s.add(line.split("\t", 1)[0])
    return s


def scan(item):
    name, acc, keep = item
    path = os.path.join(THYR, acc, f"{acc}.h.tsv.gz")
    out = []
    n = 0
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {k: i for i, k in enumerate(header)}
        ci = idx.get("rsid", 9)
        for line in fh:
            n += 1
            p = line.rstrip("\n").split("\t")
            if len(p) <= ci:
                continue
            rsid = p[ci]
            if rsid not in keep:
                continue
            try:
                out.append((rsid, p[0], p[1], p[2].upper(), p[3].upper(),
                            float(p[4]), float(p[5]),
                            float(p[6]) if p[6] not in ("", "NA") else float("nan"),
                            float(p[7]),
                            float(p[10]) if len(p) > 10 and p[10] not in ("", "NA") else float("nan")))
            except (ValueError, IndexError):
                continue
    op = os.path.join(REV, f"out_{name}.tsv")
    with open(op, "w") as fh:
        fh.write("rsid\tchr\tpos\tea\toa\tbeta\tse\teaf\tp\tn\n")
        for r in out:
            fh.write("%s\t%s\t%s\t%s\t%s\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\n" % r)
    return name, n, len(out), op


if __name__ == "__main__":
    keep = load_rsids()
    print(f"unique candidate rsIDs: {len(keep)}", flush=True)
    jobs = [(n, a, keep) for n, a in ACCS.items()]
    with ProcessPoolExecutor(max_workers=5) as ex:
        for name, n, k, op in ex.map(scan, jobs):
            print(f"{name}: scanned {n:,} rows -> matched {k} -> {op}", flush=True)
