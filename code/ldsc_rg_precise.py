# -*- coding: utf-8 -*-
"""
Precise rg estimation with joint block jackknife.

Why a separate script
---------------------
ldsc_overlap.py estimates rg as slope12 / sqrt(slope1*slope2) using h2 slopes
taken from the *univariate* runs, i.e. on different SNP sets, and reports no
SE for rg. For a shared-genetic-architecture paper we need:
  (a) all three regressions on the SAME SNP set (aligned, chi2-filtered),
  (b) a joint delete-one-block jackknife SE for rg,
  (c) attenuation ratio (intercept-1)/(mean chi2-1) to justify skipping
      stratified LDSC.

Outputs
-------
results/ldsc_rg_precise.csv    rg with SE/CI, h2_x, h2_y, rho_g, attenuation
results/ldsc_attn_ratio.csv    per-trait intercept diagnostics
"""
import os
import math
import numpy as np
import pandas as pd

import ldsc_overlap as L

RES = L.RES
N_B = L.N_BLOCKS
CHI2 = L.CHI2_CUT


def fit_with_blocks(L2, W2, y, n_iter=3):
    """Return (full beta, per-block beta matrix) with the same WLS weighting
    scheme used in ldsc_overlap.ldsc_fit."""
    X = np.column_stack([np.ones(len(L2)), L2])
    c = 0.0
    for _ in range(n_iter):
        w = 1.0 / (1.0 + max(c, 0.0) * L2) ** 2 / np.maximum(W2, 1.0)
        wx = w[:, None] * X
        A = X.T @ wx
        b = X.T @ (w * y)
        try:
            beta = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            beta = np.linalg.lstsq(A, b, rcond=None)[0]
        c = max(float(beta[1]), 0.0)
    w = 1.0 / (1.0 + max(c, 0.0) * L2) ** 2 / np.maximum(W2, 1.0)

    comp = np.column_stack([
        w * X[:, 0] * X[:, 0], w * X[:, 0] * X[:, 1], w * X[:, 1] * X[:, 1],
        w * X[:, 0] * y, w * X[:, 1] * y,
    ])
    cs = np.vstack([np.zeros(5), np.cumsum(comp, axis=0)])
    A = np.array([[cs[-1, 0], cs[-1, 1]], [cs[-1, 1], cs[-1, 2]]])
    b = cs[-1, 3:5]
    full = np.linalg.solve(A, b)

    edges = np.linspace(0, len(y), N_B + 1).astype(int)
    est = np.empty((N_B, 2))
    for k in range(N_B):
        s = cs[edges[k + 1]] - cs[edges[k]]
        Ak = A - np.array([[s[0], s[1]], [s[1], s[2]]])
        bk = b - s[3:5]
        try:
            est[k] = np.linalg.solve(Ak, bk)
        except np.linalg.LinAlgError:
            est[k] = np.linalg.lstsq(Ak, bk, rcond=None)[0]
    return full, est, w


def jk_se(est):
    return math.sqrt((N_B - 1) / N_B * float(np.sum((est - est.mean()) ** 2)))


def main():
    ref = L.load_reference()
    ref_l2 = ref[["SNP", "L2", "W2"]].copy()

    print("[1/3] loading GWAS ...", flush=True)
    g = {}
    for nm in L.EXPS + L.OUTS:
        g[nm] = L.load_gwas(nm)
    N1 = {nm: float(g[nm]["N"].median()) for nm in L.EXPS}
    N2 = L.OUT_N

    # ---------- attenuation ratio per trait ----------
    print("[2/3] attenuation ratio per trait (stratification check) ...", flush=True)
    att = []
    for nm in L.EXPS + L.OUTS:
        d = g[nm].reset_index().merge(ref_l2, on="SNP", how="inner")
        d = d[np.isfinite(d["z"]) & (d["z"] ** 2 <= CHI2)]
        mean_chi2 = float(np.mean(d["z"].to_numpy(float) ** 2))
        full, est, _ = fit_with_blocks(d["L2"].to_numpy(float),
                                       d["W2"].to_numpy(float),
                                       d["z"].to_numpy(float) ** 2)
        intercept, se_int = float(full[0]), jk_se(est[:, 0])
        ratio = (intercept - 1.0) / (mean_chi2 - 1.0) if mean_chi2 > 1 else float("nan")
        att.append(dict(trait=nm, mean_chi2=mean_chi2, intercept=intercept,
                        se_intercept=se_int, attenuation_ratio=ratio,
                        n_snp=len(d)))
        print(f"      {nm:<15s} mean_chi2={mean_chi2:7.4f}  intercept={intercept:.4f} "
              f"(se {se_int:.4f})  attn_ratio={ratio:+.3f}", flush=True)
    pd.DataFrame(att).to_csv(os.path.join(RES, "ldsc_attn_ratio.csv"), index=False)

    # ---------- rg with joint jackknife ----------
    print("[3/3] rg with joint block jackknife (18 combos) ...", flush=True)
    rows = []
    for e in L.EXPS:
        for o in L.OUTS:
            m = L.align_z(g[e], g[o], ref_l2)
            m = m[(m.chi_x <= CHI2) & (m.chi_y <= CHI2)]
            L2 = m["L2"].to_numpy(float)
            W2 = m["W2"].to_numpy(float)
            zx = m["zx"].to_numpy(float)
            zy = m["zy"].to_numpy(float)

            fx, ex, _ = fit_with_blocks(L2, W2, zx ** 2)
            fy, ey, _ = fit_with_blocks(L2, W2, zy ** 2)
            fc, ec, _ = fit_with_blocks(L2, W2, zx * zy)

            s1, s2, s12 = float(fx[1]), float(fy[1]), float(fc[1])
            rg = s12 / math.sqrt(s1 * s2) if (s1 > 0 and s2 > 0) else float("nan")

            # joint delete-one-block: recompute rg without block k
            rg_k = []
            for k in range(N_B):
                a, b, c_ = ex[k, 1], ey[k, 1], ec[k, 1]
                if a > 0 and b > 0:
                    rg_k.append(c_ / math.sqrt(a * b))
            rg_k = np.array(rg_k)
            se_rg = math.sqrt((len(rg_k) - 1) / len(rg_k) *
                              float(np.sum((rg_k - rg_k.mean()) ** 2))) if len(rg_k) > 10 else float("nan")
            lo, hi = rg - 1.96 * se_rg, rg + 1.96 * se_rg
            z = rg / se_rg if se_rg and se_rg > 0 else float("nan")

            h2x = s1 * L.M_TOTAL / N1[e]
            h2y = s2 * L.M_TOTAL / N2[o]
            rho_g = s12 * L.M_TOTAL / math.sqrt(N1[e] * N2[o])

            rows.append(dict(exposure=e, outcome=o, n_snp=len(m),
                             rg=rg, se_rg=se_rg, rg_lo=lo, rg_hi=hi, z_rg=z,
                             p_rg=float(2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2)))))
                             if np.isfinite(z) else float("nan"),
                             rg_signif=bool(np.isfinite(z) and abs(z) > 1.96),
                             h2_exposure=h2x, h2_outcome=h2y,
                             rho_g=rho_g, slope12=s12,
                             intercept_cross=float(fc[0]),
                             se_intercept_cross=jk_se(ec[:, 0])))
            print(f"      {e:<9s}->{o:<15s} rg={rg:+.4f} (se {se_rg:.4f}) "
                  f"z={z:+.2f}  h2x={h2x:.4f} h2y={h2y:.4f}", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "ldsc_rg_precise.csv"), index=False)

    print("\n=== summary ===")
    print("rg range: %+.3f to %+.3f" % (out.rg.min(), out.rg.max()))
    print("n positive: %d / %d" % ((out.rg > 0).sum(), len(out)))
    print("n nominal (|z|>1.96): %d" % out.rg_signif.sum())
    print("n Bonferroni (|z|>2.96, 18 tests): %d" % (out.z_rg.abs() > 2.96).sum())
    print("\ntop by |rg|:")
    print(out.reindex(out.rg.abs().sort_values(ascending=False).index)
          [["exposure", "outcome", "rg", "se_rg", "z_rg", "h2_exposure", "h2_outcome"]]
          .head(8).to_string(index=False, float_format=lambda x: f"{x:.4f}"))


if __name__ == "__main__":
    main()
