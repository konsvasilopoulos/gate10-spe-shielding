"""
spectrum_sensitivity.py - How much the results depend on the energy range / extrapolation of the spectra.
(a) Analytic (CSDA, PSTAR): mean and skin dose behind Al and ratios to Al, with the spectrum
    truncated at different limits (e.g. King 1972 only within its fit range 30-200 MeV).
(b) From GATE: share of the mean dose from each energy group (spe_dose_depth.csv per configuration).
Usage:  python analysis/spectrum_sensitivity.py
"""
import os, sys, csv, glob
import numpy as np
sys.path.insert(0, os.getcwd())
from spectra import SPECTRA

MATS = ["POLYETHYLENE", "PLEXIGLASS", "Al", "Ti", "Fe", "Cu", "Pb"]
AREALS = [5, 10, 20]

def read_pstar(m):
    a = np.array([[float(x) for x in l.split()] for l in open(f"pstar/G4_{m}.txt") if len(l.split()) == 7])
    return a[:, 0], a[:, 4]
Ew, Rw = read_pstar("WATER")
TAB = {m: read_pstar(m) for m in MATS}

def dose(m, t, sp, z1, z2, emin, emax, n=4000):
    E, Rc = TAB[m]; J = SPECTRA[sp]
    Eth = float(np.exp(np.interp(np.log(t), np.log(Rc), np.log(E))))
    lo = max(Eth * 1.0005, emin)
    if lo >= emax: return 0.0
    Eg = np.geomspace(lo, emax, n)
    phi = -np.gradient(J(Eg), Eg)
    Rres = np.exp(np.interp(np.log(Eg), np.log(E), np.log(Rc))) - t
    ok = Rres > Rc.min()
    Eout = np.where(ok, np.exp(np.interp(np.log(np.maximum(Rres, Rc.min())), np.log(Rc), np.log(E))), 0)
    Rwo = np.where(ok, np.exp(np.interp(np.log(np.maximum(Eout, Ew[0])), np.log(Ew), np.log(Rw))), 0)
    Eat = lambda z: np.where(Rwo - z > Rw.min(), np.exp(np.interp(np.log(np.maximum(Rwo - z, Rw.min())), np.log(Rw), np.log(Ew))), 0)
    return float(np.trapezoid(phi * (Eat(z1) - Eat(z2)), Eg) / (z2 - z1))

CASES = {
    "webber1956": [("base 10-2000", 10, 2000), ("without <30 MeV", 30, 2000), ("without >500 MeV", 10, 500), ("without >1000 MeV", 10, 1000)],
    "king1972":   [("base 10-1000", 10, 1000), ("fit range only 30-200", 30, 200), ("without >200 MeV", 10, 200),
                   ("without >400 MeV", 10, 400), ("without <30 MeV", 30, 1000)],
}
out = []
for sp, cases in CASES.items():
    print(f"\n=== {sp} (analytic, CSDA) ===")
    base = {}
    for lab, lo, hi in cases:
        line = f"{lab:<32}"
        for t in AREALS:
            al_m = dose("Al", t, sp, 0, 30, lo, hi); al_s = dose("Al", t, sp, 0, 0.3, lo, hi)
            pe = dose("POLYETHYLENE", t, sp, 0, 30, lo, hi) / al_m
            pb = dose("Pb", t, sp, 0, 30, lo, hi) / al_m
            if lab.startswith("base"): base[t] = (al_m, al_s, pe, pb)
            b = base[t]
            line += f" | {t:>2} g/cm2: Al mean {100*al_m/b[0]:5.1f}% skin {100*al_s/b[1]:5.1f}% | PE/Al {pe:.2f} Pb/Al {pb:.2f}"
            out.append([sp, lab, lo, hi, t, al_m / b[0], al_s / b[1], pe, pb])
        print(line)

print("\n=== GATE: share of the mean dose (0-30 g/cm2) per energy group ===")
for sp, base in (("webber1956", "spe_runs/materials"), ("king1972", "spe_runs/materials_king1972")):
    for m in ("POLYETHYLENE", "Al", "Pb"):
        for t in AREALS:
            f = f"{base}/{m}_{t}/spe_dose_depth.csv"
            if not os.path.exists(f): continue
            a = np.genfromtxt(f, delimiter=",", names=True)
            gcols = [c for c in a.dtype.names if c.endswith("MeV_cGy")]
            tot = sum(a[c].mean() for c in gcols)
            print(f"{sp:<11}{m:<14}{t:>3} g/cm2: " + " | ".join(f"{c.replace('MeV_cGy','')} {100*a[c].mean()/tot:4.1f}%" for c in gcols))

with open("analysis/spectrum_sensitivity.csv", "w", newline="") as fo:
    w = csv.writer(fo)
    w.writerow(["spectrum", "case", "Emin_MeV", "Emax_MeV", "areal_g_cm2", "Al_mean_dose_rel_to_base", "Al_skin_dose_rel_to_base", "PE_over_Al", "Pb_over_Al"])
    w.writerows(out)
