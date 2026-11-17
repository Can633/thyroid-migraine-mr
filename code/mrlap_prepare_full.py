"""
mrlap_prepare_full.py
---------------------
为 userIV 模式补齐输入：主 MR (TwoSampleMR) 的工具变量多数不在 HapMap3
面板内，而 mrlap_input/*.tsv 只保留了 hm3 SNP，导致 userIV 只能用到
17/49 个 IV。这里把主 MR 的 IV 并集 (iv_rsids.txt, 693 个) 从原始
GWAS 补进一份新的输入目录 mrlap_input_full/。

注意
  * 补齐的非 hm3 SNP 会被 GenomicSEM::munge 剔除，因此**不影响 LDSC**
    (h2 / 截距 / 遗传相关)；但 MRlap 的 run_MR 用未经 munge 的完整输入，
    所以 userIV 模式能用上全部主 MR 工具变量。
  * 单独目录，避免干扰已在运行的 default 批量。

运行: python mrlap_prepare_full.py
"""
import os
import shutil
import numpy as np
import pandas as pd

ROOT = r"D:/GWAS/Migraine"
SRC = os.path.join(ROOT, "mrlap_input")
DST = os.path.join(ROOT, "mrlap_input_full")
EXPS = ["GD", "HYPO", "LT", "BNG", "TC"]
OUTS = ["Migraine", "MigraineAura", "MigraineNoAura"]
OUT_N = {"Migraine": 395060, "MigraineAura": 379299, "MigraineNoAura": 376762}
CODES = dict(zip(OUTS, ["G6_MIGRAINE", "G6_MIGRAINE_WITH_AURA", "G6_MIGRAINE_NO_AURA"]))
COLS = ["rsid", "chr", "pos", "alt", "ref", "beta", "se", "N"]

os.makedirs(DST, exist_ok=True)
for f in os.listdir(SRC):
    if f.endswith(".tsv"):
        shutil.copy2(os.path.join(SRC, f), os.path.join(DST, f))

# ---- 暴露侧: iv_<e>.tsv 已是主 MR 的 IV (含 chr/pos/等位基因/beta/se/n) ----
for e in EXPS:
    iv = pd.read_csv(os.path.join(ROOT, f"iv_{e}.tsv"), sep="\t")
    add = pd.DataFrame({
        "rsid": iv.rsid.astype(str),
        "chr": pd.to_numeric(iv.chr, errors="coerce"),
        "pos": pd.to_numeric(iv.base_pair_location, errors="coerce"),
        "alt": iv.effect_allele.astype(str).str.upper(),
        "ref": iv.other_allele.astype(str).str.upper(),
        "beta": pd.to_numeric(iv.beta, errors="coerce"),
        "se": pd.to_numeric(iv.se, errors="coerce"),
        "N": pd.to_numeric(iv.n, errors="coerce"),
    })
    fp = os.path.join(DST, e + ".tsv")
    cur = pd.read_csv(fp, sep="\t")
    add = add[~add.rsid.isin(cur.rsid)]
    add = add.dropna().query("se > 0 and N > 0")
    add = add[add.chr.between(1, 22)]
    add = add.drop_duplicates("rsid")
    pd.concat([cur, add[COLS]], ignore_index=True).to_csv(fp, sep="\t", index=False)
    print(f"{e:<8s} +{len(add):>4d} IV 行 -> 共 {len(cur)+len(add):,}")

# ---- 结局侧: outcome_full_<code>.tsv 含 chr/pos/等位基因/beta/se ----
for o in OUTS:
    od = pd.read_csv(os.path.join(ROOT, f"outcome_full_{CODES[o]}.tsv"), sep="\t")
    add = pd.DataFrame({
        "rsid": od.rsid.astype(str),
        "chr": pd.to_numeric(od.chr, errors="coerce"),
        "pos": pd.to_numeric(od.pos, errors="coerce"),
        "alt": od.alt.astype(str).str.upper(),
        "ref": od.ref.astype(str).str.upper(),
        "beta": pd.to_numeric(od.beta_out, errors="coerce"),
        "se": pd.to_numeric(od.se_out, errors="coerce"),
        "N": OUT_N[o],
    })
    fp = os.path.join(DST, o + ".tsv")
    cur = pd.read_csv(fp, sep="\t")
    add = add[~add.rsid.isin(cur.rsid)]
    add = add.dropna().query("se > 0 and N > 0")
    add = add[add.chr.between(1, 22)]
    add = add.drop_duplicates("rsid")
    pd.concat([cur, add[COLS]], ignore_index=True).to_csv(fp, sep="\t", index=False)
    print(f"{o:<16s} +{len(add):>4d} IV 行 -> 共 {len(cur)+len(add):,}")

# ---- 覆盖率校验: 15 组合 mr_keep IV 是否全部在 full 输入中 ----
print("\n=== userIV 覆盖率校验 ===")
for e in EXPS:
    ep = set(pd.read_csv(os.path.join(DST, e + ".tsv"), sep="\t", usecols=["rsid"]).rsid)
    for o in OUTS:
        op = set(pd.read_csv(os.path.join(DST, o + ".tsv"), sep="\t", usecols=["rsid"]).rsid)
        h = pd.read_csv(os.path.join(ROOT, "results", f"harmonised_{e}_to_{o}.csv"))
        h = h[h.mr_keep == True]
        need = set(h.SNP.astype(str))
        print(f"  {e:<6s}->{o:<16s} mr_keep={len(need):>3d}  两侧齐全={len(need & ep & op):>3d}")
