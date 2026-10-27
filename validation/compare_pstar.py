import os
import re
import sys
import csv
import json
import glob
import itk
import numpy as np
import uproot

print(">>> compare_pstar.py")

# Usage:  python validation/compare_pstar.py                      (all materials in the list)
#         python validation/compare_pstar.py G4_Al G4_Ti          (only these)
#
# PSTAR tables (NIST): one file per material in the folder pstar/, named after the Geant4 material:
#   pstar/G4_WATER.txt, pstar/G4_Al.txt, pstar/G4_Ti.txt ...
PSTAR_DIR = "pstar"
REF = "G4_AIR"                 # run without shield (reference for the R80 shift)
DEFAULT = ["G4_POLYETHYLENE", "G4_PLEXIGLASS", "G4_Al", "G4_Ti", "G4_Fe", "G4_Cu", "G4_Pb"]
# Geant4 densities (g/cm3). CHECK them against the GATE/Geant4 version you use.
RHO = {"G4_POLYETHYLENE": 0.94, "G4_PLEXIGLASS": 1.19, "G4_Al": 2.699, "G4_Ti": 4.54,
       "G4_Fe": 7.874, "G4_Cu": 8.96, "G4_Pb": 11.35}
NZ = 200                       # voxels along z (dose.size[2])
E_CUT_FRAC = 0.3               # exit spectrum: protons with E > 30% of the initial energy

materials = sys.argv[1:] if len(sys.argv) > 1 else DEFAULT

NUM = re.compile(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?")


def read_pstar(path):
    """Reads a PSTAR table (7 columns: E, electronic, nuclear, total stopping power,
    CSDA range, projected range, detour factor). Returns (E MeV, CSDA g/cm2)."""
    rows = []
    for line in open(path, encoding="utf-8", errors="ignore"):
        nums = NUM.findall(line)
        if len(nums) == 7 and re.match(r"^\s*[-+]?\d", line):
            rows.append([float(x) for x in nums])
    if len(rows) < 10:
        raise ValueError(f"{path}: only {len(rows)} lines with 7 numbers found. "
                         "Was the PSTAR table copied correctly?")
    a = np.array(sorted(rows, key=lambda r: r[0]))
    a = a[np.concatenate([[True], np.diff(a[:, 0]) > 0])]        # remove duplicate energies
    return a[:, 0], a[:, 4]


def r_of_e(E, tabE, tabR):
    return np.exp(np.interp(np.log(E), np.log(tabE), np.log(tabR)))      # log-log interpolation


def e_of_r(R, tabE, tabR):
    return np.exp(np.interp(np.log(R), np.log(tabR), np.log(tabE)))


def read_mhd(material, kind):
    exact = f"output_{material}/dose_{material}_{kind}.mhd"
    files = [exact] if os.path.exists(exact) else glob.glob(f"output_{material}/dose_{material}*{kind}.mhd")
    if not files:
        return None, None
    img = itk.imread(files[0])
    raw = np.array(itk.array_view_from_image(img))
    if raw.size != NZ:
        print(f"[!] {files[0]}: shape {raw.shape}, expected {NZ} voxels. Skipped.")
        return None, None
    return raw.reshape(-1), img.GetSpacing()[2]


def r80(depth, prof):
    i = int(np.argmax(prof))
    if i == 0:
        return np.nan
    target = 0.8 * prof[i]
    for j in range(i, len(prof) - 1):
        if prof[j] >= target > prof[j + 1]:
            return depth[j] + (prof[j] - target) / (prof[j] - prof[j + 1]) * (depth[j + 1] - depth[j])
    return np.nan


def load_run(material):
    """Returns a dict with the simulation results, or None."""
    sj = f"output_{material}/setup.json"
    if not os.path.exists(sj):
        print(f"[!] {sj} does not exist. Run run_material.py again for {material}.")
        return None
    setup = json.load(open(sj))
    prof, dz = read_mhd(material, "dose")
    if prof is None:
        print(f"[!] No dose file found for {material}")
        return None
    depth = (np.arange(NZ) + 0.5) * dz / 10.0
    out = dict(setup=setup, r80=r80(depth, prof), e_out=np.nan)
    pf = glob.glob(f"output_{material}/phsp_{material}*.root")
    if pf:
        with uproot.open(pf[0]) as f:
            t = f[f.keys()[0]]
            ph = t.arrays(["ParticleName", "KineticEnergy", "Direction_Z"], library="np")
        sel = (ph["ParticleName"] == "proton") & (ph["Direction_Z"] > 0) \
              & (ph["KineticEnergy"] > E_CUT_FRAC * setup["energy_MeV"])
        if sel.sum() > 0:
            out["e_out"] = float(np.median(ph["KineticEnergy"][sel]))
    return out


# --- Reference (no shield) and water ---
ref = load_run(REF)
if ref is None:
    sys.exit(f"The {REF} run (no shield) is needed for the R80 shift.")
E0 = ref["setup"]["energy_MeV"]

water_file = os.path.join(PSTAR_DIR, "G4_WATER.txt")
if not os.path.exists(water_file):
    sys.exit(f"PSTAR table for water missing: {water_file}")
wE, wR = read_pstar(water_file)
Rw0 = float(r_of_e(E0, wE, wR))
print(f"E0 = {E0:.0f} MeV | CSDA range in water (PSTAR) = {Rw0:.3f} g/cm2 | "
      f"reference R80 (simulation) = {ref['r80']:.2f} cm")

rows = []
for m in materials:
    pf = os.path.join(PSTAR_DIR, f"{m}.txt")
    if not os.path.exists(pf):
        print(f"[!] PSTAR table missing: {pf}. Skipping {m}.")
        continue
    run = load_run(m)
    if run is None:
        continue
    t, n = run["setup"]["thickness_cm"], run["setup"]["n"]
    rho = RHO[m]
    mE, mR = read_pstar(pf)
    R_left = float(r_of_e(E0, mE, mR)) - rho * t              # g/cm2 left after the slab
    if R_left <= mR.min():
        print(f"[!] {m}: according to PSTAR the protons stop inside the slab.")
        continue
    e_pstar = float(e_of_r(R_left, mE, mR))
    dR_pstar = Rw0 - float(r_of_e(e_pstar, wE, wR))           # cm of water (density 1)
    dR_sim = run["r80"] - ref["r80"]
    rows.append(dict(material=m, t_cm=t, rho=rho, n=n,
                     e_sim=run["e_out"], e_pstar=e_pstar,
                     e_diff_pct=100 * (run["e_out"] - e_pstar) / e_pstar,
                     dR_sim_cm=-dR_sim, dR_pstar_cm=dR_pstar,
                     dR_diff_mm=10 * (-dR_sim - dR_pstar)))

if not rows:
    sys.exit("No results to compare.")

# --- Final table ---
print(f"\n{'Material':<16}{'t cm':>6}{'rho g/cm3':>9}{'E_out sim.':>15}{'E_out PSTAR':>13}{'D %':>7}"
      f"{'dR80 sim.':>14}{'dR PSTAR':>10}{'D mm':>7}")
for r in rows:
    print(f"{r['material']:<16}{r['t_cm']:>6.2f}{r['rho']:>9.3f}{r['e_sim']:>15.1f}{r['e_pstar']:>13.1f}"
          f"{r['e_diff_pct']:>7.2f}{r['dR_sim_cm']:>14.2f}{r['dR_pstar_cm']:>10.2f}{r['dR_diff_mm']:>7.1f}")
print("\nE_out sim. = median kinetic energy of protons crossing the plane forward "
      f"(E > {E_CUT_FRAC:.0%} of the initial energy)")
print("E_out PSTAR  = exit energy from the CSDA range: R(E_out) = R(E0) - rho*t")
print("dR80 sim. = reduction of R80 relative to the run without shield | dR PSTAR = R_water(E0) - R_water(E_out PSTAR)")

with open("comparison_pstar.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print("--> Saved as 'comparison_pstar.csv'")
