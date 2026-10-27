"""
decompose_em_nuclear.py - Ratio of the mean dose / H (0-30 g/cm2 water) to Al at 20 g/cm2, at three levels:
  analytic CSDA (PSTAR, no nuclear) | GATE EM-only | GATE full physics.
EM-only - analytic ~ scattering/straggling | full - EM-only = nuclear part.
"""
import csv, os
import os
_ALL = {"kold1956": "spe_runs/materials_kold1956/materials_summary.csv",
        "kold1972": "spe_runs/materials_kold1972/materials_summary.csv",
        "webber1956": "spe_runs/materials/materials_summary.csv",
        "king1972": "spe_runs/materials_king1972/materials_summary.csv"}
EMO = {sp: f"spe_runs/materials_EMonly_{sp}/materials_summary.csv" for sp in _ALL}
FULL = {sp: f for sp, f in _ALL.items() if os.path.exists(EMO[sp])}   # only spectra with EM-only runs
AN = "analysis/thresholds_pstar.csv"
kD, kH = "D_cGy_mean_0-30", "H_cSv_mean_0-30"
MATS = ["POLYETHYLENE", "PLEXIGLASS", "Al", "Ti", "Fe", "Cu", "Pb"]

def load(f):
    return {(r["material"], float(r["areal_g_cm2"])): r for r in csv.DictReader(open(f))} if os.path.exists(f) else {}

an = {(r["spectrum"], r["material"], float(r["areal_g_cm2"])): float(r["ratio_mean_dose_to_Al_GATE"]) for r in []}
# analytic mean-dose ratios: recomputed with thresholds_pstar (function dose_layer)
import importlib.util, sys
sys.argv = ["x"]
spec = importlib.util.spec_from_file_location("tp", "analysis/thresholds_pstar.py")
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    tp = importlib.util.module_from_spec(spec); spec.loader.exec_module(tp)

rows_out = []
for sp in FULL:
    F, E = load(FULL[sp]), load(EMO[sp])
    print(f"\n=== {sp}, 20 g/cm2: ratio to Al (mean 0-30) ===")
    print(f"{'Material':<14}{'anal. D':>10}{'EM-only D':>11}{'full D':>10}{'nucl. dD':>11}{'EM-only H':>11}{'full H':>10}{'nucl. dH':>11}")
    for m in MATS:
        a = tp.dose_layer(m, 20, sp, 0, 30) / tp.dose_layer("Al", 20, sp, 0, 30)
        f = lambda T, k: float(T[(m, 20.0)][k]) / float(T[("Al", 20.0)][k]) if (m, 20.0) in T else float("nan")
        eD, fD, eH, fH = f(E, kD), f(F, kD), f(E, kH), f(F, kH)
        print(f"{m:<14}{a:>10.2f}{eD:>11.2f}{fD:>10.2f}{fD - eD:>+11.2f}{eH:>11.2f}{fH:>10.2f}{fH - eH:>+11.2f}")
        rows_out.append([sp, m, round(a, 4), round(eD, 4), round(fD, 4), round(eH, 4), round(fH, 4),
                         float(E[(m, 20.0)][kD]), float(F[(m, 20.0)][kD]), float(E[(m, 20.0)][kH]), float(F[(m, 20.0)][kH])])
    if ("Al", 20.0) in E and ("Al", 20.0) in F:
        print(f"Al absolute: D EM-only {float(E[('Al',20.0)][kD]):.3f} vs full {float(F[('Al',20.0)][kD]):.3f} cGy | "
              f"H {float(E[('Al',20.0)][kH]):.3f} vs {float(F[('Al',20.0)][kH]):.3f} cSv")

with open("analysis/decompose_em_nuclear.csv", "w", newline="") as fo:
    w = csv.writer(fo)
    w.writerow(["spectrum", "material", "ratio_D_analytic", "ratio_D_EMonly", "ratio_D_full", "ratio_H_EMonly", "ratio_H_full",
                "D_EMonly_cGy", "D_full_cGy", "H_EMonly_cSv", "H_full_cSv"])
    w.writerows(rows_out)
print("Written: analysis/decompose_em_nuclear.csv")
