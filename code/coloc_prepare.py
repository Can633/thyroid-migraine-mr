# -*- coding: utf-8 -*-
"""
Coloc data preparation for thyroid -> migraine shared causal variant analysis.

Strategy
--------
1. Pool all IVs (lead SNPs) from the 5 thyroid exposures, cluster them into
   independent loci using a 1 Mb greedy window (coloc standard practice).
2. For every locus, take lead_pos +/- WIN bp and pull beta/se for every
   exposure and every migraine outcome from the hm3 caches built earlier.
3. Align alleles by rsID (same-strand / swapped; everything else dropped).
4. Attach 1000G EUR MAF (same reference used for LDSC) for coloc.abf.
5. Emit a long table: locus, exposure, outcome, snp, beta_x, se_x,
   beta_y, se_y, maf.

Outputs (coloc_input/)
---------------------
loci.tsv              one row per independent locus
coloc_long.tsv        long table consumed by run_coloc.R
prepare_report.txt    counts / diagnostics
"""
import os
import sys
import numpy as np
import pandas as pd

ROOT = "D:/GWAS/Migraine"
CACHE = os.path.join(ROOT, "ldsc_cache")
OUTDIR = os.path.join(ROOT, "coloc_input")
FRQ = "D:/GWAS/LD core/extracted/1000G_Phase3_frq"

EXPOSURES = ["GD", "HYPO", "LT", "BNG", "TC"]
OUTCOMES = ["Migraine", "MigraineAura", "MigraineNoAura"]

WIN = 500_000          # half-window around lead SNP
CLUSTER_GAP = 1_000_000  # max distance merging two IVs into one locus
MIN_SNP = 30           # minimum aligned SNPs to attempt coloc

COLS = ["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N", "TAG"]
DTYPES = {0: "str", 1: "int16", 2: "int32", 3: "category",
          4: "category", 5: "float64", 6: "float64"}

os.makedirs(OUTDIR, exist_ok=True)
log_lines = []


def log(msg):
    print(msg)
    log_lines.append(str(msg))


# ---------------------------------------------------------------- 1. loci
ivs = []
for e in EXPOSURES:
    f = os.path.join(ROOT, f"iv_{e}.tsv")
    if not os.path.exists(f):
        log(f"[WARN] missing {f}")
        continue
    d = pd.read_csv(f, sep="\t")
    d = d.rename(columns={"rsid": "SNP", "chr": "CHR", "base_pair_location": "BP"})
    d["exposure"] = e
    ivs.append(d[["SNP", "CHR", "BP", "pval", "exposure"]])
ivs = pd.concat(ivs, ignore_index=True)
log(f"[1] pooled IVs: {len(ivs)} from {len(EXPOSURES)} exposures")

ivs = ivs.sort_values(["CHR", "BP"]).reset_index(drop=True)
locus_id = 0
prev_chr, prev_pos = None, None
assign = []
for chr_, pos in zip(ivs.CHR.values, ivs.BP.values):
    if chr_ != prev_chr or (pos - prev_pos) > CLUSTER_GAP:
        locus_id += 1
    assign.append(locus_id)
    prev_chr, prev_pos = chr_, pos
ivs["locus"] = assign

loci = (ivs.sort_values("pval")
           .groupby("locus", as_index=False)
           .agg(chr=("CHR", "first"), lead_snp=("SNP", "first"),
                lead_pos=("BP", "first"), n_iv=("SNP", "size"),
                exposures=("exposure", lambda s: ",".join(sorted(set(s))))))
loci["start"] = np.maximum(loci.lead_pos - WIN, 1)
loci["end"] = loci.lead_pos + WIN
log(f"[2] independent loci after {CLUSTER_GAP/1e6:.1f} Mb clustering: {len(loci)}")
log("    loci per chromosome: " + str(loci.chr.value_counts().sort_index().to_dict()))
log("    IVs per locus: median %.1f  max %d" % (loci.n_iv.median(), loci.n_iv.max()))
log("    multi-exposure loci: %d" % (loci.exposures.str.contains(",")).sum())

loci.to_csv(os.path.join(OUTDIR, "loci.tsv"), sep="\t", index=False)

# ------------------------------------------------- 2. load GWAS and slice
chroms_needed = sorted(loci.chr.unique())
meta = {}
for t in EXPOSURES + OUTCOMES:
    f = os.path.join(CACHE, f"{t}.tsv")
    if not os.path.exists(f):
        log(f"[WARN] missing cache {f}")
        continue
    d = pd.read_csv(f, sep="\t", header=None, names=COLS, usecols=[0, 1, 2, 3, 4, 5, 6],
                    dtype=DTYPES)
    d = d[d.CHR.isin(chroms_needed)]
    d = d.dropna(subset=["BETA", "SE"])
    d = d[(d.SE > 0) & np.isfinite(d.BETA) & np.isfinite(d.SE)]
    d["A1"] = d.A1.astype(str).str.upper()
    d["A2"] = d.A2.astype(str).str.upper()
    meta[t] = d
    log(f"[3] loaded {t}: {len(d)} SNPs on target chromosomes")

# ------------------------------------------------- 3. slice per locus
records = []
diag = []
for _, L in loci.iterrows():
    sub = {t: d[(d.CHR == L.chr) & (d.BP >= L.start) & (d.BP <= L.end)]
           for t, d in meta.items()}
    n_ref = len(sub[OUTCOMES[0]]) if OUTCOMES[0] in sub else 0
    row = dict(locus=int(L.locus), chr=int(L.chr), lead_snp=L.lead_snp,
               n_snp_window=n_ref, n_exposures_in_locus=L.exposures)
    for e in EXPOSURES:
        for o in OUTCOMES:
            if e not in sub or o not in sub:
                continue
            x = sub[e].drop_duplicates("SNP").set_index("SNP")
            y = sub[o].drop_duplicates("SNP").set_index("SNP")
            common = x.index.intersection(y.index)
            common = common[common.notna()]
            if len(common) < MIN_SNP:
                diag.append(dict(locus=int(L.locus), exposure=e, outcome=o,
                                 n_aligned=len(common), status="too_few_snps"))
                continue
            a = x.loc[common]
            b = y.loc[common]
            # allele alignment
            same = (a.A1.values == b.A1.values) & (a.A2.values == b.A2.values)
            flip = (a.A1.values == b.A2.values) & (a.A2.values == b.A1.values)
            keep = same | flip
            if keep.sum() < MIN_SNP:
                diag.append(dict(locus=int(L.locus), exposure=e, outcome=o,
                                 n_aligned=int(keep.sum()), status="too_few_after_align"))
                continue
            snp = common[keep]
            bx = a.BETA.values[keep].astype(float)
            sx = a.SE.values[keep].astype(float)
            by = b.BETA.values[keep].astype(float)
            sy = b.SE.values[keep].astype(float)
            by = np.where(flip[keep], -by, by)
            records.append(pd.DataFrame(dict(
                locus=int(L.locus), chr=int(L.chr), lead_snp=L.lead_snp,
                exposure=e, outcome=o, snp=snp,
                beta_x=bx, se_x=sx, beta_y=by, se_y=sy)))
            diag.append(dict(locus=int(L.locus), exposure=e, outcome=o,
                             n_aligned=int(keep.sum()), status="ok"))
    records.append(None)  # placeholder keeps loop cheap
    if len(records) % 20 == 0:
        log(f"    ... processed {len(records)} locus entries")

recs = [r for r in records if r is not None]
long = pd.concat(recs, ignore_index=True) if recs else pd.DataFrame()
log(f"[4] candidate (locus, exposure, outcome) blocks: {long.groupby(['locus','exposure','outcome']).ngroups if len(long) else 0}")

# ------------------------------------------------- 4. attach 1000G MAF
needed = set(long.snp.unique()) if len(long) else set()
maf_map = {}
for c in chroms_needed:
    f = os.path.join(FRQ, f"1000G.EUR.QC.{c}.frq")
    if not os.path.exists(f):
        continue
    d = pd.read_csv(f, sep=r"\s+")
    d = d[d.SNP.isin(needed)]
    if len(d):
        maf_map.update(dict(zip(d.SNP, d.MAF.astype(float))))
log(f"[5] MAF resolved for {len(maf_map)}/{len(needed)} SNPs "
    f"({100*len(maf_map)/max(len(needed),1):.1f}%)")
long["maf"] = long.snp.map(maf_map)
before = len(long)
long = long.dropna(subset=["maf"])
log(f"    dropped {before-len(long)} SNP rows without MAF")

# loci must retain enough SNPs after MAF merge
cnt = long.groupby(["locus", "exposure", "outcome"]).size().rename("n").reset_index()
good = cnt[cnt.n >= MIN_SNP][["locus", "exposure", "outcome"]]
long = long.merge(good, on=["locus", "exposure", "outcome"], how="inner")
log(f"[6] final blocks: {len(good)}   final SNP rows: {len(long)}")
log("    block size: median %.0f  min %d  max %d" %
    (good.merge(cnt, on=["locus","exposure","outcome"]).n.median(),
     good.merge(cnt, on=["locus","exposure","outcome"]).n.min(),
     good.merge(cnt, on=["locus","exposure","outcome"]).n.max()))

long.to_csv(os.path.join(OUTDIR, "coloc_long.tsv"), sep="\t", index=False)
pd.DataFrame(diag).to_csv(os.path.join(OUTDIR, "prepare_diagnostics.tsv"), sep="\t", index=False)

with open(os.path.join(OUTDIR, "prepare_report.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(log_lines))
log("[DONE] written to " + OUTDIR)
