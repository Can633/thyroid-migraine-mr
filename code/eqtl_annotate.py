# -*- coding: utf-8 -*-
"""Map Ensembl gene IDs in the eQTL outputs to HGNC symbols and label the
shared loci with their candidate effector genes."""
import os
import pandas as pd

ROOT = "D:/GWAS/Migraine"
RES = os.path.join(ROOT, "results")
HGNC = os.path.join(ROOT, "eqtl", "hgnc_complete_set.txt")

h = pd.read_csv(HGNC, sep="\t", low_memory=False)
h = h[["ensembl_gene_id", "symbol", "name", "location"]].dropna(subset=["ensembl_gene_id"])
h["ensembl_gene_id"] = h.ensembl_gene_id.astype(str).str.split(".").str[0]
sym = dict(zip(h.ensembl_gene_id, h.symbol))
full = dict(zip(h.ensembl_gene_id, h.name))
loc = dict(zip(h.ensembl_gene_id, h.location))
print(f"HGNC map: {len(sym)} Ensembl IDs")


def add(df, col="gene_id"):
    df = df.copy()
    gid = df[col].astype(str).str.split(".").str[0]
    df["gene_symbol"] = gid.map(sym)
    df["gene_name"] = gid.map(full)
    df["gene_location"] = gid.map(loc)
    df[col] = gid
    return df


# --- per-SNP annotation
snp = pd.read_csv(os.path.join(RES, "eqtl_snp_annotation.csv"))
# 幂等：首次运行时列为 best_gene，重跑后已改名为 best_ensg
_col = "best_gene" if "best_gene" in snp.columns else "best_ensg"
snp = add(snp, _col)
if _col == "best_gene":
    snp = snp.rename(columns={"best_gene": "best_ensg"})
snp.to_csv(os.path.join(RES, "eqtl_snp_annotation.csv"), index=False)
print("\n=== per-SNP annotation (with symbols) ===")
print(snp[["locus", "snp", "is_lead", "n_eqtl_genes", "n_sig_genes",
           "best_ensg", "gene_symbol", "best_p", "best_beta", "sig"]]
      .to_string(index=False, float_format=lambda x: f"{x:.3g}"))

# --- per-gene summary
gs = pd.read_csv(os.path.join(RES, "eqtl_gene_summary.csv"))
gs = add(gs, "gene_id")
gs.to_csv(os.path.join(RES, "eqtl_gene_summary.csv"), index=False)
print("\n=== top thyroid eQTL genes per locus ===")
print(gs[["locus", "gene_symbol", "gene_name", "rsid", "beta", "pvalue"]]
      .to_string(index=False, float_format=lambda x: f"{x:.3g}"))

# --- locus-level one-line summary for the report
print("\n=== locus summary ===")
for lid in sorted(gs.locus.unique()):
    s = gs[gs.locus == lid].head(3)
    genes = ", ".join(f"{r.gene_symbol or r.gene_id}({r.pvalue:.1e})" for r in s.itertuples())
    print(f"  locus {lid}: {genes}")
