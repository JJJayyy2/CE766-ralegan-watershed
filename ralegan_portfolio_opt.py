"""
CE 766 Project - Topic 1: Multi-objective optimisation of watershed intervention
portfolios, Ralegan Siddhi watershed (Parner, Ahmednagar, Maharashtra), 2000-2023.

Run:  python ralegan_portfolio_opt.py            (uses PLACEHOLDER inputs - demo only)
      python ralegan_portfolio_opt.py rain.csv   (uses real daily rainfall)

rain.csv columns: date (YYYY-MM-DD), rain_mm   <- extract from IMD 0.25 deg gridded
rainfall for the grid cell(s) covering the watershed.

!! Every value in the PARAMS / ZONES blocks marked [ASSUMED] is a placeholder.
   Replace with values derived from DEM, Bhuvan LULC, ISRIC soils, literature.
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pymoo.core.problem import Problem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize

rng = np.random.default_rng(42)
YEARS = np.arange(2000, 2024)

# ---------------------------------------------------------------- RAINFALL
def load_rain(path=None):
    """Return array (n_years, 365) of daily rainfall (mm)."""
    if path:
        df = pd.read_csv(path, parse_dates=["date"])
        df = df[(df.date.dt.year >= 2000) & (df.date.dt.year <= 2023)]
        df = df[~((df.date.dt.month == 2) & (df.date.dt.day == 29))]
        return df.groupby(df.date.dt.year).rain_mm.apply(lambda s: s.values[:365]).tolist() and \
            np.vstack([g.rain_mm.values[:365] for _, g in df.groupby(df.date.dt.year)])
    # PLACEHOLDER synthetic rainfall (semi-arid Deccan, mean ~500 mm, CV ~0.3)
    ann = rng.gamma(shape=1 / 0.3**2, scale=500 * 0.3**2, size=len(YEARS))
    R = np.zeros((len(YEARS), 365))
    for i, a in enumerate(ann):
        wet = rng.random(365) < np.r_[np.full(151, 0.05), np.full(122, 0.35), np.full(92, 0.04)]
        d = rng.exponential(1.0, 365) * wet
        R[i] = d / d.sum() * a
    return R

# ----------------------------------------------------------------- ZONES
# 5 sub-units. area_ha, CN (AMC II), K, LS, C0 : all [ASSUMED] placeholders
ZONES = pd.DataFrame({
    "zone": ["Z1-upper", "Z2-upper-mid", "Z3-mid", "Z4-lower-mid", "Z5-lower"],
    "area_ha": [220, 200, 210, 180, 172],           # total ~ 982 ha [ASSUMED]
    "CN":     [86, 84, 82, 80, 78],
    "K":      [0.30, 0.28, 0.26, 0.24, 0.22],
    "LS":     [3.5, 2.6, 1.8, 1.2, 0.8],
    "C0":     [0.30, 0.28, 0.25, 0.22, 0.20],
})
NZ = len(ZONES)

PARAMS = dict(                                      # all [ASSUMED]
    dCN_bund=4.0, dCN_aff=10.0,                     # CN reduction at 100 % treatment
    cap_cd=2500.0, fills_cd=3.0, rech_cd=0.60,      # check dam: m3/dam, fills/yr, recharge frac
    cap_fp=400.0, fills_fp=2.0, rech_fp=0.30,       # farm pond
    bund_capture=0.15,                              # extra runoff infiltrated by bunds (x frac)
    P_bund=0.5, C_forest=0.02, SDR=0.35, trap_eff=0.8,
    base_rech_frac=0.06,                            # natural recharge, frac of rainfall
    irr_demand=4000.0,                              # m3/ha supplemental irrigation
    cost_cd=450_000, cost_fp=60_000, cost_bund_ha=8_000, cost_aff_ha=25_000,   # INR
    budget=25_000_000,
)
MAXN = dict(cd=4, fp=20, bund=0.6, aff=0.3)

# ----------------------------------------------------------------- MODEL
def simulate(x, R):
    """x = (NZ x 4): n_checkdam, n_farmpond, frac_bund, frac_aff. Returns per-year metrics."""
    p = PARAMS
    ann_rain = R.sum(1)
    rech = np.zeros(len(R)); sed = np.zeros(len(R)); irr_ha = np.zeros(len(R)); cost = 0.0
    for z in range(NZ):
        A = ZONES.area_ha[z]; a_m2 = A * 1e4
        ncd, nfp, fb, fa = x[z]
        cn = np.clip(ZONES.CN[z] - p["dCN_bund"] * fb / MAXN["bund"] - p["dCN_aff"] * fa / MAXN["aff"], 40, 98)
        S = 25400 / cn - 254
        Q = np.where(R > 0.2 * S, (R - 0.2 * S) ** 2 / (R + 0.8 * S + 1e-9), 0).sum(1)   # mm/yr
        V = Q / 1000 * a_m2                                                              # m3/yr
        cap_cd = ncd * p["cap_cd"] * p["fills_cd"]; cap_fp = nfp * p["cap_fp"] * p["fills_fp"]
        cd_in = np.minimum(V, cap_cd); fp_in = np.minimum(V - cd_in, cap_fp)
        bund_in = fb * p["bund_capture"] * (V - cd_in - fp_in)
        rech_z = (p["base_rech_frac"] * ann_rain / 1000 * a_m2 * (1 - fa * 0.2)
                  + cd_in * p["rech_cd"] + fp_in * p["rech_fp"] + bund_in)
        rech += rech_z
        # RUSLE
        Rf = 79 + 0.363 * ann_rain
        C = ZONES.C0[z] * (1 - fa / MAXN["aff"] * (fa / MAXN["aff"]) * 0 ) if False else \
            ZONES.C0[z] * (1 - fa) + p["C_forest"] * fa
        Pf = 1 - (1 - p["P_bund"]) * fb / MAXN["bund"] * (fb / MAXN["bund"] > 0)
        loss = Rf * ZONES.K[z] * ZONES.LS[z] * C * Pf * A                               # t/yr
        trap = p["trap_eff"] * np.minimum(1, (cd_in + 0.5 * fp_in) / np.maximum(V, 1))
        sed += loss * p["SDR"] * (1 - trap)
        # livelihood: irrigable area supported (ha), minus cropland lost to afforestation
        irr_ha += (cd_in * p["rech_cd"] + fp_in) / p["irr_demand"] - fa * A * 0.3
        cost += ncd * p["cost_cd"] + nfp * p["cost_fp"] + fb * A * p["cost_bund_ha"] + fa * A * p["cost_aff_ha"]
    return cost, rech, sed, irr_ha

class PortfolioProblem(Problem):
    def __init__(self, R):
        xu = np.tile([MAXN["cd"], MAXN["fp"], MAXN["bund"], MAXN["aff"]], NZ)
        super().__init__(n_var=4 * NZ, n_obj=4, n_ieq_constr=1, xl=np.zeros(4 * NZ), xu=xu)
        self.R = R
    def _evaluate(self, X, out, *a, **k):
        F = np.zeros((len(X), 4)); G = np.zeros((len(X), 1))
        for i, row in enumerate(X):
            x = row.reshape(NZ, 4).copy(); x[:, :2] = np.round(x[:, :2])
            c, r, s, l = simulate(x, self.R)
            F[i] = [c / 1e6, -r.mean() / 1e3, s.mean(), -l.mean()]   # min cost, max recharge, min sed, max irr ha
            G[i] = (c - PARAMS["budget"]) / 1e6
        out["F"], out["G"] = F, G

# ------------------------------------------------------------------ RUN
if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else None
    R = load_rain(path)
    print("Rainfall mean annual (mm): %.0f | min %.0f | max %.0f" % (R.sum(1).mean(), R.sum(1).min(), R.sum(1).max()))
    # baseline: no interventions
    base = simulate(np.zeros((NZ, 4)), R)
    res = minimize(PortfolioProblem(R), NSGA2(pop_size=120), ("n_gen", 150), seed=1, verbose=False)
    F, X = res.F, res.X
    df = pd.DataFrame({"cost_MINR": F[:, 0], "recharge_TCM_per_yr": -F[:, 1],
                       "sediment_t_per_yr": F[:, 2], "irrigated_ha": -F[:, 3]})
    for z in range(NZ):
        for j, n in enumerate(["checkdams", "farmponds", "bund_frac", "aff_frac"]):
            df[f"Z{z+1}_{n}"] = X[:, 4 * z + j] if j >= 2 else np.round(X[:, 4 * z + j])
    df.to_csv("pareto_portfolios.csv", index=False)

    # robustness: recharge in the driest 25 % of years for each Pareto solution
    dry = np.argsort(R.sum(1))[: len(R) // 4]
    dry_rech = []
    for row in X:
        x = row.reshape(NZ, 4).copy(); x[:, :2] = np.round(x[:, :2])
        dry_rech.append(simulate(x, R)[1][dry].mean() / 1e3)
    df["dry_year_recharge_TCM"] = dry_rech

    # ---------- figures
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    sc = ax[0].scatter(df.cost_MINR, df.recharge_TCM_per_yr, c=df.sediment_t_per_yr, cmap="viridis_r", s=22)
    ax[0].scatter(0, base[1].mean() / 1e3, marker="*", s=200, c="red", label="No intervention")
    ax[0].set(xlabel="Cost (million INR)", ylabel="Mean recharge (thousand m3/yr)", title="Cost vs recharge")
    plt.colorbar(sc, ax=ax[0], label="Sediment yield (t/yr)"); ax[0].legend()
    sc2 = ax[1].scatter(df.recharge_TCM_per_yr, df.sediment_t_per_yr, c=df.cost_MINR, cmap="plasma", s=22)
    ax[1].scatter(base[1].mean() / 1e3, base[2].mean(), marker="*", s=200, c="red")
    ax[1].set(xlabel="Mean recharge (thousand m3/yr)", ylabel="Sediment yield (t/yr)", title="Recharge vs erosion")
    plt.colorbar(sc2, ax=ax[1], label="Cost (M INR)")
    freq = pd.DataFrame({
        "Check dams": [(df[f"Z{z+1}_checkdams"] > 0).mean() for z in range(NZ)],
        "Farm ponds": [(df[f"Z{z+1}_farmponds"] > 0).mean() for z in range(NZ)],
        "Contour bunds": [(df[f"Z{z+1}_bund_frac"] > 0.05).mean() for z in range(NZ)],
        "Afforestation": [(df[f"Z{z+1}_aff_frac"] > 0.02).mean() for z in range(NZ)]},
        index=[f"Z{i+1}" for i in range(NZ)])
    freq.plot.bar(ax=ax[2], rot=0); ax[2].set(ylabel="Share of Pareto portfolios using it", title="Selection frequency by zone")
    plt.tight_layout(); plt.savefig("preliminary_results.png", dpi=170)

    print("\nBaseline (no interventions): recharge %.1f TCM/yr, sediment %.0f t/yr" % (base[1].mean() / 1e3, base[2].mean()))
    print("Pareto solutions:", len(df))
    print(df[["cost_MINR", "recharge_TCM_per_yr", "sediment_t_per_yr", "irrigated_ha", "dry_year_recharge_TCM"]].describe().round(1))
    print("\nSelection frequency:\n", freq.round(2))
    df.to_csv("pareto_portfolios.csv", index=False)
