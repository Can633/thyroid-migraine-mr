# -*- coding: utf-8 -*-
"""秩一（单共享因子）加权拟合 + r_g 随机效应合并。"""
import numpy as np, pandas as pd
from scipy.optimize import least_squares
from scipy.stats import norm

EXPS = ["GD", "HYPO", "LT", "BNG", "TC", "BNG_LOO"]
OUTS = ["Migraine", "MigraineAura", "MigraineNoAura"]


def pool(y, se):
    w = 1.0 / se ** 2
    b = (w * y).sum() / w.sum()
    Q = (w * (y - b) ** 2).sum()
    df = len(y) - 1
    C = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max(0.0, (Q - df) / C)
    wr = 1.0 / (se ** 2 + tau2)
    br = (wr * y).sum() / wr.sum()
    ser = np.sqrt(1.0 / wr.sum())
    i2 = max(0.0, (Q - df) / Q * 100)
    return br, br - 1.96 * ser, br + 1.96 * ser, 2 * (1 - norm.cdf(abs(br / ser))), Q, i2, tau2


def rank1(R, SE):
    W = 1.0 / SE ** 2
    W = W / W.max()          # 归一化，避免溢出
    ne, no = R.shape
    x0 = np.concatenate([np.full(ne, 0.4), np.full(no, 0.4)])

    def res(p):
        lam, mu = p[:ne], p[ne:]
        return ((R - np.outer(lam, mu)) * np.sqrt(W)).ravel()

    sol = least_squares(res, x0, method="lm", max_nfev=20000)
    lam, mu = sol.x[:ne], sol.x[ne:]
    fit = np.outer(lam, mu)
    ssr = np.sum(W * (R - fit) ** 2)
    sst = np.sum(W * (R - np.average(R, weights=W)) ** 2)
    return lam, mu, 1 - ssr / sst, np.max(np.abs(R - fit))


d = pd.read_csv("results/ldsc_rg_precise.csv")
R = d.pivot(index="exposure", columns="outcome", values="rg").loc[EXPS, OUTS].values
SE = d.pivot(index="exposure", columns="outcome", values="se_rg").loc[EXPS, OUTS].values
lam, mu, r2, maxres = rank1(R, SE)
print("秩一（单共享因子）模型  加权 R² = %.3f  最大残差 = %.3f" % (r2, maxres))
print("  暴露载荷:", {k: round(float(v), 3) for k, v in zip(EXPS, lam)})
print("  结局载荷:", {k: round(float(v), 3) for k, v in zip(OUTS, mu)})
print("  18 对中 r_g<0 的个数:", int((R < 0).sum()))

b, lo, hi, p, Q, i2, tau2 = pool(d.rg.values, d.se_rg.values)
print("\n18 对合并 r_g = %.3f (95%% CI %.3f–%.3f), P = %.2g; Q = %.1f, I² = %.0f%%, tau² = %.4f"
      % (b, lo, hi, p, Q, i2, tau2))
d15 = d[d.exposure != "BNG_LOO"]
b2, lo2, hi2, p2, Q2, i22, _ = pool(d15.rg.values, d15.se_rg.values)
print("15 对（去 LOO）合并 r_g = %.3f (%.3f–%.3f), P = %.2g; I² = %.0f%%" % (b2, lo2, hi2, p2, i22))

u = pd.read_csv("results_ukb/ldsc_rg_precise.csv")
print("\nUKB r_g:")
print(u[["exposure", "rg", "se_rg", "p_rg"]].to_string(index=False, float_format=lambda x: f"{x:.3f}"))
b3, lo3, hi3, p3, Q3, i23, _ = pool(u.rg.values, u.se_rg.values)
print("UKB 6 对合并 r_g = %.3f (%.3f–%.3f), P = %.2g; I² = %.0f%%" % (b3, lo3, hi3, p3, i23))
print("FinnGen − UKB 差值 = %.3f" % (b - b3))
