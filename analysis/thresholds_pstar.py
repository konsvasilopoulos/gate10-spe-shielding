"""
thresholds_pstar.py - Analytic explanation of the material comparison with the NIST PSTAR tables (no Monte Carlo).
For each material and areal density t:
  E_th : minimum energy for a proton to cross the shield (CSDA range = t)
  Phi(>E_th) for each spectrum (transmitted protons)
and comparison of the ratios to Al with the skin- and mean-dose ratios of the simulations.
Usage (from the repository root):  python analysis/thresholds_pstar.py
"""
import os, sys, csv
import numpy as np
sys.path.insert(0, os.getcwd())
from spectra import SPECTRA

MATS = ["LITHIUM_HYDRIDE", "POLYETHYLENE", "WATER", "PLEXIGLASS", "Al", "Ti", "Fe", "Cu", "Pb"]
AREALS = [5, 10, 20]
RUNS = {"webber1956": "spe_runs/materials/materials_summary.csv",
        "king1972": "spe_runs/materials_king1972/materials_summary.csv",
        "kold1956": "spe_runs/materials_kold1956/materials_summary.csv",
        "kold1972": "spe_runs/materials_kold1972/materials_summary.csv"}
SKIN, MEAN = "D_cGy_skin_0-0.3", "D_cGy_mean_0-30"


def read_pstar(m):
    rows = [l.split() for l in open(f"pstar/G4_{m}.txt") if len(l.split()) == 7]
    a = np.array([[float(x) for x in r] for r in rows])
    return a[:, 0], a[:, 3], a[:, 4]        # E (MeV), total stopping power (MeV cm2/g), CSDA range (g/cm2)


def e_of_r(R, E, Rc):
    return float(np.exp(np.interp(np.log(R), np.log(Rc), np.log(E))))


sim = {}
for sp, f in RUNS.items():
    for r in csv.DictReader(open(f)):
        sim[(sp, r["material"], float(r["areal_g_cm2"]))] = (float(r[SKIN]), float(r[MEAN]))

tab = {}
for m in MATS:
    E, S, Rc = read_pstar(m)
    S100 = float(np.interp(100, E, S))
    for t in AREALS:
        tab[(m, t)] = dict(Eth=e_of_r(t, E, Rc), S100=S100)

out = []
print("Minimum transmission energy E_th (MeV) from PSTAR (CSDA) and stopping power at 100 MeV (MeV cm2/g)\n")
print(f"{'Material':<14}{'S(100 MeV)':>11}" + "".join(f"{f'E_th {t}':>10}" for t in AREALS))
for m in MATS:
    print(f"{m:<14}{tab[(m,5)]['S100']:>11.3f}" + "".join(f"{tab[(m,t)]['Eth']:>10.1f}" for t in AREALS))

for sp in RUNS:
    J = SPECTRA[sp]
    print(f"\n=== {sp}: ratio to Al - Phi(>E_th) (analytic) vs skin dose and mean dose (GATE) ===")
    print(f"{'Material':<14}" + "".join(f"{f'{t} g/cm2: Phi / skin / mean':>32}" for t in AREALS))
    for m in MATS:
        line = f"{m:<14}"
        for t in AREALS:
            phi = J(tab[(m, t)]["Eth"]) / J(tab[("Al", t)]["Eth"])
            sk = sim[(sp, m, t)][0] / sim[(sp, "Al", t)][0]
            mn = sim[(sp, m, t)][1] / sim[(sp, "Al", t)][1]
            line += f"{phi:>14.2f} /{sk:>6.2f} /{mn:>6.2f}"
            out.append([sp, m, t, tab[(m, t)]["Eth"], J(tab[(m, t)]["Eth"]), phi, sk, mn])
        print(line)

with open("analysis/thresholds_pstar.csv", "w", newline="") as fo:
    w = csv.writer(fo)
    w.writerow(["spectrum", "material", "areal_g_cm2", "E_th_MeV", "fluence_above_Eth_cm2",
                "ratio_fluence_to_Al", "ratio_skin_dose_to_Al_GATE", "ratio_mean_dose_to_Al_GATE"])
    w.writerows(out)
r = np.array([[x[5], x[6], x[7]] for x in out])
print(f"\nCorrelation of log(Phi ratio) with log(skin-dose ratio): r = {np.corrcoef(np.log(r[:,0]), np.log(r[:,1]))[0,1]:.3f}"
      f" | with log(mean-dose ratio): r = {np.corrcoef(np.log(r[:,0]), np.log(r[:,2]))[0,1]:.3f}")
print("Written: analysis/thresholds_pstar.csv")


# --- Analytic dose (CSDA, straight-ahead, no nuclear reactions) ----------------------------------
# For each proton with E > E_th: residual range R_mat(E) - t -> E_out from the material, then into water.
# Dose in a water layer [z1, z2] (g/cm2) = int phi(E) [E(z1) - E(z2)] dE / (z2 - z1)   (MeV/g per source cm2)
Ew, Sw, Rw = read_pstar("WATER")
def dose_layer(m, t, sp, z1, z2, n=4000):
    E, S, Rc = read_pstar(m)
    J = SPECTRA[sp]
    Eg = np.geomspace(max(tab[(m, t)]["Eth"] * 1.0005, 10), 2000, n)
    phi = -np.gradient(J(Eg), Eg)                                   # dPhi/dE
    Rres = np.exp(np.interp(np.log(Eg), np.log(E), np.log(Rc))) - t   # residual g/cm2 of material
    ok = Rres > Rc.min()
    Eout = np.zeros_like(Eg)
    Eout[ok] = np.exp(np.interp(np.log(Rres[ok]), np.log(Rc), np.log(E)))
    Rw_out = np.where(ok, np.exp(np.interp(np.log(np.maximum(Eout, 1e-3)), np.log(Ew), np.log(Rw))), 0)
    def E_at(z):
        rr = Rw_out - z
        return np.where(rr > Rw.min(), np.exp(np.interp(np.log(np.maximum(rr, Rw.min())), np.log(Rw), np.log(Ew))), 0)
    dep = E_at(z1) - E_at(z2)
    return (np.trapezoid if hasattr(np,'trapezoid') else np.trapz)(phi * dep, Eg) / (z2 - z1)

print("\n=== Analytic dose (CSDA, no nuclear reactions) as ratio to Al, vs GATE ===")
avg = []
for sp in RUNS:
    print(f"\n{sp}: skin 0-0.3 g/cm2 (analytic / GATE) | mean 0-30 (analytic / GATE)")
    for m in MATS:
        line = f"{m:<14}"
        for t in AREALS:
            a_s = dose_layer(m, t, sp, 0, 0.3) / dose_layer("Al", t, sp, 0, 0.3)
            a_m = dose_layer(m, t, sp, 0, 30) / dose_layer("Al", t, sp, 0, 30)
            g_s = sim[(sp, m, t)][0] / sim[(sp, "Al", t)][0]
            g_m = sim[(sp, m, t)][1] / sim[(sp, "Al", t)][1]
            line += f"  {t:>2}: {a_s:5.2f}/{g_s:5.2f} · {a_m:5.2f}/{g_m:5.2f}"
            avg.append([sp, m, t, round(tab[(m, t)]["Eth"], 2), round(a_s, 4), round(g_s, 4), round(a_m, 4), round(g_m, 4),
                        dose_layer(m, t, sp, 0, 30) * 1.602176634e-8, sim[(sp, m, t)][1]])
        print(line)

with open("analysis/analytic_vs_gate.csv", "w", newline="") as fo:
    w = csv.writer(fo)
    w.writerow(["spectrum", "material", "areal_g_cm2", "E_th_MeV", "skin_ratio_analytic", "skin_ratio_GATE",
                "mean_ratio_analytic", "mean_ratio_GATE", "mean_D_analytic_cGy", "mean_D_GATE_cGy"])
    w.writerows(avg)
a = np.array([[r[6], r[7]] for r in avg]); k = np.array([[r[4], r[5]] for r in avg])
print(f"\nMean dose: r(log) = {np.corrcoef(np.log(a[:,0]), np.log(a[:,1]))[0,1]:.4f}, max |GATE - analytic| = {np.max(np.abs(a[:,1]-a[:,0])):.3f}")
print(f"Skin:      r(log) = {np.corrcoef(np.log(k[:,0]), np.log(k[:,1]))[0,1]:.4f}, max |GATE - analytic| = {np.max(np.abs(k[:,1]-k[:,0])):.3f}")
print("Written: analysis/analytic_vs_gate.csv")
