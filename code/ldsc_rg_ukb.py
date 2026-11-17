import os, sys
sys.path.insert(0, r"D:/GWAS/Migraine")
import ldsc_overlap as L
L.OUTS = ["MigraineUKB"]; L.OUT_N = {"MigraineUKB": 458440}
L.ACC = {**L.ACC, "MigraineUKB": "GCST90473326"}
L.RES = os.path.join(L.ROOT, "results_ukb")
import ldsc_rg_precise as P
P.RES = L.RES
if __name__ == "__main__":
    P.main()
