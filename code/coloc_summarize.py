# -*- coding: utf-8 -*-
"""
Summarise coloc.abf output.

Key distinction
---------------
A high PP.H0 is only meaningful if BOTH traits actually have a signal in the
region. Otherwise it is "no data", not "no sharing". So every block is labelled:

  informative   min_z_exp >= 3 and min_z_out >= 3
  shared        informative and PP.H4 >= 0.8      (strong shared causal variant)
  suggestive    informative and 0.5 <= PP.H4 < 0.8
  distinct      informative and PP.H3 >= 0.8      (two independent signals)
  no_signal     not informative -> excluded from the denominator

Outputs
-------
results/coloc_blocks.tsv       per-block table with labels
results/coloc_by_exposure.csv  per-exposure summary
results/coloc_hits.csv         loci with PP.H4 >= 0.5 (for reporting)
"""
import os
import numpy as np
import pandas as pd

ROOT = "D:/GWAS/Migraine"
RES = os.path.join(ROOT, "results")
IN = os.path.join(ROOT, "coloc_results", "coloc_main.tsv")

Z_SIG = 3.0
H4_STRONG, H4_SUG, H3_STRONG = 0.8, 0.5, 0.8

d = pd.read_csv(IN, sep="\t")
d = d[d.status == "ok"].copy()
print(f"blocks with status=ok: {len(d)}")

d["informative"] = (d.min_z_exp >= Z_SIG) & (d.min_z_out >= Z_SIG)
d["label"] = "no_signal"
inf = d.informative
d.loc[inf & (d["PP.H4"] >= H4_STRONG), "label"] = "shared"
d.loc[inf & (d["PP.H4"] >= H4_SUG) & (d["PP.H4"] < H4_STRONG), "label"] = "suggestive"
d.loc[inf & (d["PP.H3"] >= H3_STRONG) & (d["PP.H4"] < H4_SUG), "label"] = "distinct"

d.to_csv(os.path.join(RES, "coloc_blocks.tsv"), sep="\t", index=False)

print("\n=== label counts (all blocks) ===")
print(d.label.value_counts().to_string())
print(f"\ninformative blocks: {int(inf.sum())} / {len(d)} "
      f"({100*inf.mean():.1f}%)")

sub = d[inf]
print(f"\n=== among informative blocks (n={len(sub)}) ===")
print(sub.label.value_counts().to_string())
print("\nPP.H4 distribution among informative:")
q = sub["PP.H4"].quantile([0.5, 0.75, 0.9, 0.95, 1.0]).round(4)
print(q.to_string())

by = (sub.groupby("exposure")
         .agg(n_informative=("locus", "size"),
              n_shared=("label", lambda s: (s == "shared").sum()),
              n_suggestive=("label", lambda s: (s == "suggestive").sum()),
              n_distinct=("label", lambda s: (s == "distinct").sum()),
              max_H4=("PP.H4", "max"), median_H4=("PP.H4", "median"))
         .reset_index())
by["pct_shared"] = (100 * by.n_shared / by.n_informative).round(1)
by.to_csv(os.path.join(RES, "coloc_by_exposure.csv"), index=False)
print("\n=== by exposure (informative blocks only) ===")
print(by.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

hits = d[d["PP.H4"] >= H4_SUG].sort_values("PP.H4", ascending=False)
hits.to_csv(os.path.join(RES, "coloc_hits.csv"), index=False)
print(f"\n=== loci with PP.H4 >= {H4_SUG}: {len(hits)} ===")
if len(hits):
    cols = ["locus", "exposure", "outcome", "nsnps", "min_z_exp", "min_z_out",
            "PP.H3", "PP.H4", "top_snp", "top_snp_pp", "label"]
    print(hits[cols].head(20).to_string(index=False, float_format=lambda x: f"{x:.4f}"))

# 每个 (exposure, outcome) 组合的最高 H4 locus -- useful for reporting
best = (d[inf].sort_values("PP.H4", ascending=False)
        .groupby(["exposure", "outcome"], as_index=False).first())
best.to_csv(os.path.join(RES, "coloc_best_per_combo.csv"), index=False)
print("\n=== best locus per (exposure, outcome) ===")
print(best[["exposure", "outcome", "locus", "PP.H4", "PP.H3", "min_z_exp",
            "min_z_out", "label"]].to_string(index=False, float_format=lambda x: f"{x:.3f}"))
