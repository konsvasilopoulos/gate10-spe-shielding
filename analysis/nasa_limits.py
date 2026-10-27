"""nasa_limits.py - Comparison with the NASA 30-day limits (NASA-STD-3001 Vol. 1 Rev. C, Tables 4.8-1, 4.8-2).
Slab approximation, normal incidence, total event fluence (as in NASA/TP-2009-215560):
  skin ~ mean dose 0-0.3 g/cm2, lens ~ same layer (0.3 cm), BFO ~ dose at 5 g/cm2 of water.
  mGy-Eq = 1.5 x D (RBE of protons > 2 MeV; protons give > 90 % of D) - conservative for e-/gamma.
Writes analysis/nasa_limits.csv. Usage: python analysis/nasa_limits.py"""
import csv, os
import numpy as np
LIM = {"skin": 1500, "lens": 1000, "BFO": 250}            # mGy-Eq, 30 days
RBE = 1.5
BASES = {"kold1956": "spe_runs/materials_kold1956", "kold1972": "spe_runs/materials_kold1972",
         "kold1989": "spe_runs/materials_kold1989", "webber1956": "spe_runs/materials", "king1972": "spe_runs/materials_king1972"}
MATS = ["POLYETHYLENE", "PLEXIGLASS", "Al", "Ti", "Fe", "Cu", "Pb"]
rows = []; need = []
for sp, b in BASES.items():
    f = f"{b}/materials_summary.csv"
    if not os.path.exists(f): continue
    d = {(r["material"], int(float(r["areal_g_cm2"]))): r for r in csv.DictReader(open(f))}
    print(f"\n=== {sp}: mGy-Eq (30-day limits: skin 1500, lens 1000, BFO 250) ===")
    print(f"{'Material':<14}" + "".join(f"{f'{t} g/cm2: skin / BFO':>26}" for t in (5, 10, 20)))
    for m in MATS:
        line = f"{m:<14}"; bfo = []
        for t in (5, 10, 20):
            r = d[(m, t)]
            sk = RBE * float(r["D_cGy_skin_0-0.3"]) * 10; bf = RBE * float(r["D_cGy_5_gcm2"]) * 10
            bfo.append(bf)
            flag = ("!" if sk > LIM["lens"] else "") + ("*" if bf > LIM["BFO"] else "")
            line += f"{sk:>12.0f} / {bf:>7.0f}{flag:<3}"
            rows.append([sp, m, t, sk, sk, bf, sk > LIM["skin"], sk > LIM["lens"], bf > LIM["BFO"]])
        # minimum areal density for BFO < 250 mGy-Eq (logarithmic interpolation)
        t = np.array([5, 10, 20.]); y = np.log(np.array(bfo))
        if bfo[0] <= LIM["BFO"]: tmin = "≤ 5"
        elif bfo[-1] > LIM["BFO"]: tmin = "> 20"
        else:
            i = np.where(np.array(bfo) <= LIM["BFO"])[0][0]
            tmin = f"{np.interp(np.log(LIM['BFO']), [y[i], y[i-1]], [t[i], t[i-1]]):.1f}"
        need.append([sp, m, tmin]); print(line + f"  -> BFO limit met from {tmin} g/cm2")
print("\n* = above the BFO limit (250 mGy-Eq) | ! = above the lens limit (1000 mGy-Eq)")
with open("analysis/nasa_limits.csv", "w", newline="") as fo:
    w = csv.writer(fo); w.writerow(["spectrum", "material", "areal_g_cm2", "skin_mGyEq", "lens_mGyEq", "BFO_mGyEq",
                                    "over_skin_1500", "over_lens_1000", "over_BFO_250"]); w.writerows(rows)
with open("analysis/nasa_limits_min_areal.csv", "w", newline="") as fo:
    w = csv.writer(fo); w.writerow(["spectrum", "material", "min_areal_for_BFO_250_g_cm2"]); w.writerows(need)
print("Written: analysis/nasa_limits.csv, analysis/nasa_limits_min_areal.csv")
