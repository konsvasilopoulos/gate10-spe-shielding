"""iso_test.py - Results of run_iso_test.sh: PE/Al and Pb/Al ratios (mean D 0-30 and BFO 4-6 g/cm2)
for each angle, for isotropic incidence (mean of the 4 cosine-law angles) and for 0 deg.
All with the energy groups >= 120 MeV. Writes analysis/iso_test.csv."""
import os, glob, json, csv
import numpy as np
ANG = ["0", "20.7", "37.8", "52.2", "69.3"]; ISO = ANG[1:]
MATS = ["Al", "POLYETHYLENE", "Pb"]

def dose(d):
    """Depth-dose (Gy) summed over the energy groups, each scaled by Phi_group x A_source / N."""
    tot = None; n = 0
    for g in sorted(glob.glob(f"{d}/g*")):
        sj = f"{g}/setup.json"
        if not os.path.exists(sj): continue
        s = json.load(open(sj)); f = s["group_fluence_cm2"] * s["source_area_cm2"] / s["n"]
        mhd = [x for x in glob.glob(f"{g}/dose_*_dose.mhd") if "uncertainty" not in x]
        a = np.fromfile(mhd[0][:-4] + ".raw", dtype=np.float64) * f
        tot = a if tot is None else tot + a; n += 1
    return tot, n

rows = []
for sp in ("kold1972", "kold1956"):
    base = f"spe_runs/iso_{sp}"
    if not os.path.isdir(base): continue
    D = {}
    for m in MATS:
        for a in ANG:
            v, n = dose(f"{base}/{m}_20_t{a}")
            if v is not None and n == 5: D[(m, a)] = v
    z = (np.arange(300) + 0.5) * 0.1
    W = {"mean": (z >= 0) & (z <= 30), "BFO": (z >= 4) & (z <= 6)}
    print(f"\n=== {sp}, 20 g/cm2, groups >= 120 MeV ===")
    print(f"{'angle':<8}{'PE/Al mean':>12}{'Pb/Al mean':>12}{'PE/Al BFO':>11}{'Pb/Al BFO':>11}{'Al mean/0deg':>14}")
    cases = [(a, [a]) for a in ANG] + [("iso", ISO)]
    ref = None
    for lab, angs in cases:
        if not all((m, a) in D for m in MATS for a in angs): print(f"{lab:<8} runs missing"); continue
        v = {m: np.mean([D[(m, a)] for a in angs], axis=0) for m in MATS}
        r = {w: {m: v[m][sel].mean() for m in MATS} for w, sel in W.items()}
        if lab == "0": ref = r["mean"]["Al"]
        rel = r["mean"]["Al"] / ref if ref else float("nan")
        print(f"{lab:<8}{r['mean']['POLYETHYLENE']/r['mean']['Al']:>12.3f}{r['mean']['Pb']/r['mean']['Al']:>12.3f}"
              f"{r['BFO']['POLYETHYLENE']/r['BFO']['Al']:>11.3f}{r['BFO']['Pb']/r['BFO']['Al']:>11.3f}{rel:>14.3f}")
        rows.append([sp, lab, r["mean"]["POLYETHYLENE"] / r["mean"]["Al"], r["mean"]["Pb"] / r["mean"]["Al"],
                     r["BFO"]["POLYETHYLENE"] / r["BFO"]["Al"], r["BFO"]["Pb"] / r["BFO"]["Al"], rel])
with open("analysis/iso_test.csv", "w", newline="") as fo:
    w = csv.writer(fo); w.writerow(["spectrum", "angle", "PE_Al_mean", "Pb_Al_mean", "PE_Al_BFO", "Pb_Al_BFO", "Al_mean_rel_0deg"]); w.writerows(rows)
print("\nWritten: analysis/iso_test.csv")
