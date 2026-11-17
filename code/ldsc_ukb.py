# -*- coding: utf-8 -*-
"""ldsc_ukb.py —— 复用 ldsc_overlap 的回归实现，把结局换成 UK Biobank 偏头痛"""
import os, sys
sys.path.insert(0, r"D:/GWAS/Migraine")
import ldsc_overlap as L

L.OUTS = ["MigraineUKB"]
L.OUT_N = {"MigraineUKB": 458440}
L.ACC = {**L.ACC, "MigraineUKB": "GCST90473326"}
L.POSITIVE_PAIRS = [("MigraineUKB", "MigraineUKB")]
L.RES = os.path.join(L.ROOT, "results_ukb")
os.makedirs(L.RES, exist_ok=True)

if __name__ == "__main__":
    L.main()
