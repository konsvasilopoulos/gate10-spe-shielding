"""spectral_uncertainty.py - Propagates the uncertainty of the Koldobskiy et al. (2021) reconstructions
through the response matrix (no new simulations).

Koldobskiy et al. give, for each GLE, the average width Delta of the 68 % envelope of the fit (their Table 2,
column Delta, averaged from 230 MV to the highest neutron-monitor rigidity). The envelope is represented by
two perturbations of the integral spectrum Phi(>E) of each GLE:
  norm : Phi -> Phi * (1 +- Delta)                         (changes absolute doses, not the ratios)
  tilt : Phi -> Phi * (1 +- Delta * x(E)),  x = ln(E/30 MeV) / ln(1000/30), clipped to [0, 1]
         (no change at 30 MeV, +-Delta above 1 GeV: a harder or softer spectrum)
For the 1972 and 1989 sums all GLEs are perturbed in the same direction (fully correlated, conservative).
Reports the PE/Al and Pb/Al ratios of the mean D and H, the BFO dose behind Al and the minimum areal density
of PE and Al for the NASA 30-day BFO limit (250 mGy-Eq, RBE 1.5).
Usage (from the repository root):  python analysis/spectral_uncertainty.py   -> analysis/spectral_uncertainty.csv
"""
import os, sys, csv
import numpy as np
sys.path.insert(0, os.getcwd()); sys.path.insert(0, "analysis")
from spectra import koldobskiy_gle_integral
import response_matrix as rm

DELTA = {5: 0.210, 24: 0.267, 25: 0.215, 42: 0.096, 43: 0.135, 44: 0.156, 45: 0.130}   # Koldobskiy et al. 2021, Table 2
EVENTS = {"kold1956": [5], "kold1972": [24, 25], "kold1989": [42, 43, 44, 45]}
AREALS = [5, 10, 20]
RBE, LIMIT = 1.5, 250.0


def x_tilt(E):
    return np.clip(np.log(np.asarray(E, float) / 30.0) / np.log(1000.0 / 30.0), 0.0, 1.0)


def make_spectrum(gles, mode, sign):
    def F(E):
        tot = 0.0
        for g in gles:
            f = koldobskiy_gle_integral(E, g)
            if mode == "norm":
                f = f * (1 + sign * DELTA[g])
            elif mode == "tilt":
                f = f * (1 + sign * DELTA[g] * x_tilt(E))
            tot = tot + f
        return tot
    return F


def min_areal(bfo):
    """Areal density where the BFO dose reaches the limit (log-linear interpolation, as in nasa_limits.py)."""
    t = np.array(AREALS, float); y = np.log(np.array(bfo))
    if bfo[0] <= LIMIT: return "<=5"
    if bfo[-1] > LIMIT: return ">20"
    i = np.where(np.array(bfo) <= LIMIT)[0][0]
    return f"{np.interp(np.log(LIMIT), [y[i], y[i - 1]], [t[i], t[i - 1]]):.1f}"


rows = []
for ev, gles in EVENTS.items():
    print(f"\n=== {ev} (GLE {', '.join(map(str, gles))}; Delta = {', '.join(f'{100*DELTA[g]:.0f} %' for g in gles)}) ===")
    print(f"{'case':<10}" + "".join(f"{f'{t} g/cm2: PE/Al D,H  Pb/Al D,H':>34}" for t in AREALS) + f"{'BFO Al20':>10}{'PE min':>8}{'Al min':>8}")
    for mode, sign in (("central", 0), ("norm", +1), ("norm", -1), ("tilt", +1), ("tilt", -1)):
        F = make_spectrum(gles, mode, sign)
        lab = "central" if mode == "central" else f"{mode}{'+' if sign > 0 else '-'}"
        line = f"{lab:<10}"; bfo = {}
        for t in AREALS:
            al = {q: rm.fold(F, "Al", t, q) for q in "DH"}
            pe = [rm.fold(F, "POLYETHYLENE", t, q) / al[q] for q in "DH"]
            pb = [rm.fold(F, "Pb", t, q) / al[q] for q in "DH"]
            line += f"{pe[0]:>14.3f}{pe[1]:>6.3f}{pb[0]:>8.3f}{pb[1]:>6.3f}"
            for m in ("Al", "POLYETHYLENE"):
                bfo.setdefault(m, []).append(RBE * 10 * rm.fold(F, m, t, "D", "BFO"))
            rows.append([ev, lab, t, pe[0], pe[1], pb[0], pb[1], al["D"], al["H"]])
        line += f"{bfo['Al'][-1]:>10.0f}{min_areal(bfo['POLYETHYLENE']):>8}{min_areal(bfo['Al']):>8}"
        rows[-1] += [bfo["Al"][-1], min_areal(bfo["POLYETHYLENE"]), min_areal(bfo["Al"])]
        print(line)

with open("analysis/spectral_uncertainty.csv", "w", newline="") as fo:
    w = csv.writer(fo)
    w.writerow(["event", "case", "areal_g_cm2", "PE_Al_D", "PE_Al_H", "Pb_Al_D", "Pb_Al_H", "Al_mean_D_cGy", "Al_mean_H_cSv",
                "BFO_Al20_mGyEq", "min_areal_PE", "min_areal_Al"])
    for r in rows: w.writerow(r + [""] * (12 - len(r)))
print("\nWritten: analysis/spectral_uncertainty.csv")
