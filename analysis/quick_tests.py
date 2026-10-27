"""quick_tests.py - Results of run_quick_tests.sh.
A. Lateral geometry 120/200/200 vs 80/120/120 (Koldobskiy 1956, 20 g/cm2): ratios of D and H.
B. Layered shields (Koldobskiy 1972, 20 g/cm2): Al10+PE10 and PE10+Al10 vs Al 20 and PE 20."""
import csv, os
import numpy as np
def prof(d):
    f = f"{d}/dose_eq_depth.csv"
    return np.genfromtxt(f, delimiter=",", names=True) if os.path.exists(f) else None
def v(a, k, lo, hi):
    m = (a["depth_g_cm2"] >= lo) & (a["depth_g_cm2"] <= hi); return a[k][m].mean()
W = {"skin 0-0.3": (0, 0.3), "5 g/cm2": (4, 6), "25-30": (25, 30), "mean 0-30": (0, 30)}
out = []
print("A. Wider geometry / standard (Koldobskiy 1956, 20 g/cm2)")
print(f"{'Material':<14}" + "".join(f"{'D ' + w:>16}{'H':>7}" for w in W))
for m in ("Al", "POLYETHYLENE", "Pb"):
    a, b = prof(f"spe_runs/lateral_kold1956/{m}_20"), prof(f"spe_runs/materials_kold1956/{m}_20")
    if a is None or b is None: print(f"{m:<14} data missing"); continue
    line = f"{m:<14}"
    for w, (lo, hi) in W.items():
        rd, rh = v(a, "D_cGy", lo, hi) / v(b, "D_cGy", lo, hi), v(a, "H_cSv", lo, hi) / v(b, "H_cSv", lo, hi)
        line += f"{rd:>16.3f}{rh:>7.3f}"; out.append(["lateral", m, w, rd, rh])
    print(line)
print("\nB. Layered shields 20 g/cm2 (Koldobskiy 1972): ratio to Al 20")
print(f"{'Shield':<26}" + "".join(f"{'D ' + w:>16}{'H':>7}" for w in W))
ref = prof("spe_runs/materials_kold1972/Al_20")
for lab, d in (("PE 20", "spe_runs/materials_kold1972/POLYETHYLENE_20"), ("Al 10 + PE 10 (PE inside)", "spe_runs/layered_kold1972/Al10_PE10_20"),
               ("PE 10 + Al 10 (Al inside)", "spe_runs/layered_kold1972/PE10_Al10_20")):
    a = prof(d)
    if a is None or ref is None: print(f"{lab:<26} data missing"); continue
    line = f"{lab:<26}"
    for w, (lo, hi) in W.items():
        rd, rh = v(a, "D_cGy", lo, hi) / v(ref, "D_cGy", lo, hi), v(a, "H_cSv", lo, hi) / v(ref, "H_cSv", lo, hi)
        line += f"{rd:>16.3f}{rh:>7.3f}"; out.append(["layered", lab, w, rd, rh])
    print(line)
with open("analysis/quick_tests.csv", "w", newline="") as fo:
    wr = csv.writer(fo); wr.writerow(["test", "case", "window", "ratio_D", "ratio_H"]); wr.writerows(out)
print("\nWritten: analysis/quick_tests.csv")
