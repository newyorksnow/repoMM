import os
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
import numpy as np, pandas as pd, itertools, json
rng = np.random.default_rng(7)

P = pd.read_csv(os.path.join(HERE, "poisson_permeal.csv"))
se_fisher = P.se_beta.mean(); sig_pair = np.sqrt(2) * se_fisher
h = lambda t: np.where(t > 1e-12, (1 - np.exp(-t)) / np.maximum(t, 1e-12), 1.0)

def inv_h(v):
    if v >= 1: return 0.0
    if v <= 0: return np.inf
    lo, hi = 1e-9, 1e4
    for _ in range(200):
        m = .5 * (lo + hi)
        if h(m) > v: lo = m
        else: hi = m
    return .5 * (lo + hi)

# ---------- 1. sigma_D floor from counting alone ----------
Nbar = P.N.mean()
print("=== sigma_D ===")
print(f"mean bites per meal N = {Nbar:.1f}")
print(f"Poisson floor sigma_D = sqrt(N_i+N_j) = {np.sqrt(2*Nbar):.2f} bites")
print(f"  (irreducible: holds under perfect measurement and a perfect model)")
print(f"holdout-RMSE route            = {np.sqrt(2)*6.79:.2f} bites")
print(f"  -> below the Poisson floor, because pooled RMSE over all time points")
print(f"     understates end-of-meal error; treat as optimistic")
print(f"split-half forward prediction = 19.52 bites")

# ---------- 2. pseudo-dyads under Poisson fits ----------
PD = pd.DataFrame([dict(G=abs(P.beta[b] - P.beta[a]),
                        S=float(P.beta[a] + P.beta[b]),
                        Tm=float(.5 * (P["T"][a] + P["T"][b])))
                   for a, b in itertools.combinations(range(21), 2)])
PD["AGT"] = PD.G * PD.Tm; PD["nu"] = PD.G / sig_pair
print(f"\n=== pseudo-dyads (Poisson fits), n={len(PD)} ===")
print(f"nu   median={PD.nu.median():.3f} max={PD.nu.max():.3f} "
      f"frac>=3: {(PD.nu>=3).mean():.3f}  frac>=2: {(PD.nu>=2).mean():.3f}")
print(f"AGT  median={PD.AGT.median():.2f} max={PD.AGT.max():.2f}")

# ---------- 3. sensitivity surface: theta_max over (sigma_D, nu_crit) ----------
print("\n=== theta_max as a function of sigma_D and the criterion ===")
print(f"{'sigma_D':>8} " + "".join(f"{('nu*=%g'%k):>10}" for k in (1.5, 2, 3, 4)))
AGT = PD.AGT.median()
for sD in (2.6, 5, 9.6, 11.2, 15, 19.5):
    row = f"{sD:8.1f} "
    for k in (1.5, 2, 3, 4):
        hm = k * sD / AGT
        row += f"{'empty' if hm >= 1 else f'{inv_h(hm):.2f}':>10}"
    print(row)
print(f"(at the median pseudo-dyad gap A*Gamma*T = {AGT:.1f} bites)")

# break-even sigma_D at each criterion
print("\nlargest sigma_D admitting any theta, at the median gap:")
for k in (1.5, 2, 3, 4):
    print(f"  nu*={k}: sigma_D < {AGT/k:.2f} bites")

# ---------- 4. simulation study ----------
def simulate(nu_target, theta, m_solo, T=740.0, S=0.18, R=4000):
    """Beta_i,beta_j set to hit nu_target given m_solo meals per diner."""
    lam_solo = S / 2 * T                       # expected bites per solo meal
    se_b = np.sqrt(lam_solo / m_solo) / T      # SE of beta_hat from m meals
    Gam = nu_target * np.sqrt(2) * se_b
    bi, bj = (S - Gam) / 2, (S + Gam) / 2
    if bi <= 0: return None
    hth = h(theta)
    ki = rng.poisson(S * T / 2 - Gam * T * hth / 2, R)
    kj = rng.poisson(S * T / 2 + Gam * T * hth / 2, R)
    Ni = rng.poisson(bi * T * m_solo, R); Nj = rng.poisson(bj * T * m_solo, R)
    bih, bjh = Ni / (T * m_solo), Nj / (T * m_solo)
    Gh, Sh = bjh - bih, bjh + bih
    with np.errstate(divide="ignore", invalid="ignore"):
        hh = (kj - ki) / (kj + ki) * Sh / Gh
    ok = np.isfinite(hh)
    th = np.array([inv_h(v) if 0 < v < 1 else np.nan for v in hh[ok]])
    good = np.isfinite(th)
    # Fieller-style interval via the delta variance, for coverage
    varh = hh[ok] ** 2 * (1 / nu_target ** 2 * (1 + (Gam / S) ** 2)
                          + (S * T) / (Gam * T * hth) ** 2)
    lo, hi = hh[ok] - 1.96 * np.sqrt(varh), hh[ok] + 1.96 * np.sqrt(varh)
    cov = np.mean((lo <= hth) & (hth <= hi))
    return dict(nu=nu_target, theta=theta, m=m_solo,
                frac_usable=good.mean(),
                median_theta=np.nanmedian(th),
                bias=np.nanmedian(th) - theta,
                iqr=np.nanpercentile(th[good], 75) - np.nanpercentile(th[good], 25)
                if good.sum() > 10 else np.nan,
                coverage=cov)

print("\n=== simulation: recovery of theta ===")
print(f"{'nu':>5}{'theta':>7}{'m':>4}{'usable':>8}{'median th':>11}"
      f"{'bias':>8}{'IQR':>8}{'cover':>8}")
res = []
for nu in (1, 2, 3, 5, 8):
    for th in (0.5, 1.0, 2.0):
        r = simulate(nu, th, m_solo=4)
        if r:
            res.append(r)
            print(f"{r['nu']:5.0f}{r['theta']:7.1f}{r['m']:4d}"
                  f"{r['frac_usable']:8.2f}{r['median_theta']:11.2f}"
                  f"{r['bias']:+8.2f}{r['iqr']:8.2f}{r['coverage']:8.2f}")
pd.DataFrame(res).to_csv(os.path.join(HERE, "sim.csv"), index=False)

print("\n=== effect of solo replication at theta=1 ===")
print(f"{'nu':>5}{'m':>4}{'usable':>8}{'bias':>8}{'IQR':>8}{'cover':>8}")
res2 = []
for nu in (2, 3, 5):
    for m in (1, 4, 16):
        r = simulate(nu, 1.0, m_solo=m)
        if r:
            res2.append(r)
            print(f"{r['nu']:5.0f}{m:4d}{r['frac_usable']:8.2f}"
                  f"{r['bias']:+8.2f}{r['iqr']:8.2f}{r['coverage']:8.2f}")
pd.DataFrame(res2).to_csv(os.path.join(HERE, "sim_m.csv"), index=False)

json.dump(dict(sig_pair=float(sig_pair), AGT_med=float(AGT),
               nu_max=float(PD.nu.max()), nu_med=float(PD.nu.median()),
               sigmaD_floor=float(np.sqrt(2 * Nbar))),
          open(os.path.join(HERE, "nums2.json"), "w"))
