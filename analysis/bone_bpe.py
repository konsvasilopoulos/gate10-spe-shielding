"""bone_bpe.py - Two sensitivity tests of the ratios to Al (November 2026):
  (a) bone layer: 2 g/cm2 of cortical bone (G4_BONE_CORTICAL_ICRP, 1.92 g/cm3) between a 20 g/cm2 shield (Al, PE, Pb)
      and the water phantom (run_set8_lih_water_rm.sh -> spe_runs/bone_<spectrum>), compared with the runs without bone
      (spe_runs/materials_<spectrum>, same spectrum and shield);
  (b) borated polyethylene: PE with 5 % boron by mass (1.04 g/cm3, I = 58.35 eV from Bragg additivity), run together
      with PE and Al at 20 g/cm2 with opengate 10.1.1 / Geant4 11.4.2 (spe_runs/paper2_custom_<spectrum>);
      only ratios within that release are used.
Writes analysis/bone_bpe.csv (ratios for mean 0-30 g/cm2 and BFO (5 g/cm2) dose and dose equivalent).
Usage (from the repository root):  python analysis/bone_bpe.py"""
import csv, os

Q = [("D mean", "D_cGy_mean_0-30"), ("H mean", "H_cSv_mean_0-30"), ("D BFO", "D_cGy_5_gcm2"), ("H BFO", "H_cSv_5_gcm2")]


def read(path):
    return {r["material"]: r for r in csv.DictReader(open(path))}


def main():
    rows = [["test", "spectrum", "ratio"] + [f"{q} with" for q, _ in Q] + [f"{q} without" for q, _ in Q] + [f"{q} change %" for q, _ in Q]]
    for sp in ("kold1956", "kold1972"):
        fb, fn = f"spe_runs/bone_{sp}/materials_summary.csv", f"spe_runs/materials_{sp}/materials_summary.csv"
        if os.path.exists(fb) and os.path.exists(fn):
            b = read(fb); n = {r["material"]: r for r in csv.DictReader(open(fn)) if float(r["areal_g_cm2"]) == 20}
            al_b = next(v for k, v in b.items() if k.startswith("Al"))
            for key, mat, lab in (("PE20", "POLYETHYLENE", "PE/Al"), ("Pb20", "Pb", "Pb/Al")):
                rb = next(v for k, v in b.items() if k.startswith(key))
                w = [float(rb[c]) / float(al_b[c]) for _, c in Q]
                wo = [float(n[mat][c]) / float(n["Al"][c]) for _, c in Q]
                rows.append(["bone 2 g/cm2", sp, lab] + [round(x, 4) for x in w + wo] + [round(100 * (a / c - 1), 2) for a, c in zip(w, wo)])
            rows.append(["bone 2 g/cm2", sp, "Al with/without bone"] + [round(float(al_b[c]) / float(n["Al"][c]), 4) for _, c in Q] + [""] * 8)
        fp = f"spe_runs/paper2_custom_{sp}/materials_summary.csv"
        if os.path.exists(fp):
            p = read(fp)
            rows.append(["borated PE (GATE 10.1.1)", sp, "BPE5/PE"] + [round(float(p["My_BPE5"][c]) / float(p["POLYETHYLENE"][c]), 4) for _, c in Q] + [""] * 8)
    if len(rows) == 1:
        print("No run summaries found (spe_runs/ is part of the Zenodo dataset); analysis/bone_bpe.csv kept."); return
    for r in rows[1:]:
        print("  " + "  ".join(str(x) for x in r[:7]))
    with open("analysis/bone_bpe.csv", "w", newline="") as fo:
        csv.writer(fo).writerows(rows)
    print("Written: analysis/bone_bpe.csv")


if __name__ == "__main__":
    main()
