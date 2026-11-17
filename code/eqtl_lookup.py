# -*- coding: utf-8 -*-
"""
Annotate colocalising loci with GTEx v8 Thyroid eQTLs.

For every locus with coloc PP.H4 >= 0.5, pull the eQTL records in the region
via remote tabix, then report:
  * whether the coloc top SNP / locus lead SNP is itself a significant eQTL
  * which genes are the strongest eQTL targets in the region
  * the eQTL direction at the shared SNP

Outputs
-------
eqtl/thyroid_locus_eqtl.tsv    all eQTL rows in candidate regions
results/eqtl_snp_annotation.csv  per-SNP annotation for reported SNPs
results/eqtl_gene_summary.csv    per-locus top eQTL genes
"""
import os
import sys
import pandas as pd

import tabix_query as T

ROOT = "D:/GWAS/Migraine"
URL = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/Thyroid.tsv.gz"
IDX = os.path.join(ROOT, "eqtl", "Thyroid.tsv.gz.tbi")

COLS = ["variant", "r2", "pvalue", "molecular_trait_object_id", "molecular_trait_id",
        "maf", "gene_id", "median_tpm", "beta", "se", "an", "ac",
        "chromosome", "position", "ref", "alt", "type", "rsid"]

P_SIG = 1e-5   # nominal eQTL threshold reported alongside (GTEx cis-eQTL)

loci = pd.read_csv(os.path.join(ROOT, "coloc_input", "loci.tsv"), sep="\t")
hits = pd.read_csv(os.path.join(ROOT, "results", "coloc_hits.csv"))  # comma-separated

# loci that reached PP.H4 >= 0.5 in any (exposure, outcome) pair
cand_ids = sorted(hits.locus.unique())
print(f"candidate loci from coloc (PP.H4>=0.5): {len(cand_ids)} -> {cand_ids}")

idx = T.open_index(IDX)
rows = []
for lid in cand_ids:
    r = loci[loci.locus == lid].iloc[0]
    chrom, start, end = str(int(r.chr)), int(r.start), int(r.end)
    try:
        got = T.fetch(URL, idx, chrom, start, end)
    except Exception as e:
        print(f"  locus {lid} chr{chrom}:{start}-{end} FAILED: {e}")
        continue
    print(f"  locus {lid} chr{chrom}:{start}-{end} -> {len(got)} eQTL rows")
    for _, line in got:
        f = line.split("\t")
        if len(f) >= len(COLS):
            rows.append(f[:len(COLS)])

if not rows:
    print("ERROR: no eQTL rows retrieved")
    sys.exit(1)

e = pd.DataFrame(rows, columns=COLS)
for c in ["pvalue", "beta", "se", "maf", "position", "an", "median_tpm"]:
    e[c] = pd.to_numeric(e[c], errors="coerce")
e["locus"] = None
for lid in cand_ids:
    r = loci[loci.locus == lid].iloc[0]
    m = (e.chromosome.astype(str) == str(int(r.chr))) & \
        (e.position >= r.start) & (e.position <= r.end)
    e.loc[m, "locus"] = lid
e = e.dropna(subset=["locus"])
e["locus"] = e.locus.astype(int)

out_eqtl = os.path.join(ROOT, "eqtl", "thyroid_locus_eqtl.tsv")
e.to_csv(out_eqtl, sep="\t", index=False)
print(f"\n[eQTL] {len(e)} rows, {e.gene_id.nunique()} genes, {e.locus.nunique()} loci")
print(f"       nominal p<{P_SIG}: {(e.pvalue < P_SIG).sum()} rows, "
      f"{(e[e.pvalue < P_SIG]).gene_id.nunique()} genes")

# ---------------------------------------------------------------- SNP level
snp_rows = []
for lid in cand_ids:
    sub = hits[hits.locus == lid]
    locus = loci[loci.locus == lid].iloc[0]
    snps = set(sub.top_snp.dropna().astype(str)) | {str(locus.lead_snp)}
    reg = e[e.locus == lid]
    for s in snps:
        ss = reg[(reg.rsid == s) | (reg.variant == s)]
        if len(ss) == 0:
            snp_rows.append(dict(locus=lid, snp=s, is_lead=(s == locus.lead_snp),
                                 n_eqtl_genes=0, best_gene=None, best_p=None,
                                 best_beta=None, sig=False,
                                 note="not in thyroid eQTL panel"))
            continue
        ss = ss.sort_values("pvalue")
        best = ss.iloc[0]
        snp_rows.append(dict(
            locus=lid, snp=s, is_lead=(s == locus.lead_snp),
            n_eqtl_genes=len(ss),
            n_sig_genes=int((ss.pvalue < P_SIG).sum()),
            best_gene=best.gene_id, best_p=best.pvalue, best_beta=best.beta,
            sig=bool(best.pvalue < P_SIG),
            note=""))
snpdf = pd.DataFrame(snp_rows)
snpdf.to_csv(os.path.join(ROOT, "results", "eqtl_snp_annotation.csv"), index=False)
print("\n=== per-SNP thyroid eQTL annotation ===")
print(snpdf.to_string(index=False, float_format=lambda x: f"{x:.3g}"))

# ---------------------------------------------------------------- gene level
gs = (e.sort_values("pvalue")
        .groupby(["locus", "gene_id"], as_index=False)
        .first()[["locus", "gene_id", "rsid", "variant", "beta", "se", "pvalue", "maf"]])
gs = gs.sort_values(["locus", "pvalue"])
short = (gs[gs.pvalue < P_SIG]
           .groupby("locus")
           .head(5)
           .sort_values(["locus", "pvalue"]))
short.to_csv(os.path.join(ROOT, "results", "eqtl_gene_summary.csv"), index=False)
print("\n=== top thyroid eQTL genes per locus (p<%.0e) ===" % P_SIG)
print(short.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
