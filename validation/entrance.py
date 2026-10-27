import re
import os
import sys
import csv
import json
import glob

import itk
import numpy as np

# Usage:
#   python validation/entrance.py                          -> all output_G4_* folders
#   python validation/entrance.py output_G4_Al output_G4_Ti -> only these
#
# Entrance energy per targeted proton (MeV) = edep(first voxel) / (N x 0.5625).
# Reference: electronic stopping power of water at the energy reaching the water (PSTAR) x 0.1 cm.
#   - air: 150 MeV (no slab)
#   - slab materials: simulated E_out from comparison_pstar.csv

PSTAR_WATER = "pstar/G4_WATER.txt"
NUM = re.compile(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?")
E_AIR = 150.0
THICK_WATER_CM = 0.1


def load_water():
    rows = []
    for line in open(PSTAR_WATER, encoding="utf-8", errors="ignore"):
        nums = NUM.findall(line)
        if len(nums) == 7 and re.match(r"^\s*[-+]?\d", line):
            rows.append([float(x) for x in nums])
    if len(rows) < 10:
        sys.exit(f"Table {PSTAR_WATER} was not read correctly ({len(rows)} lines).")
    a = np.array(sorted(rows, key=lambda r: r[0]))
    return a[:, 0], a[:, 1]          # E (MeV), electronic stopping power (MeV cm2/g)


def stopping(E, tabE, tabS):
    return float(np.exp(np.interp(np.log(E), np.log(tabE), np.log(tabS))))


def sim_eout():
    out = {}
    if os.path.exists("comparison_pstar.csv"):
        for row in csv.DictReader(open("comparison_pstar.csv", encoding="utf-8")):
            out[row["material"]] = float(row["e_sim"])
    else:
        print("[i] comparison_pstar.csv not found: slab materials will have no reference.")
    return out


tabE, tabS = load_water()
eout = sim_eout()
folders = sys.argv[1:] or sorted(glob.glob("output_G4_*"))

print(f"{'Folder':<22}{'E at water MeV':>16}{'Entrance MeV/p':>16}{'Reference MeV/p':>16}{'D %':>8}")
for folder in folders:
    sj = os.path.join(folder, "setup.json")
    if not os.path.exists(sj):
        print(f"[!] {folder}: no setup.json, skipped")
        continue
    setup = json.load(open(sj))
    mat, n = setup["material"], setup["n"]
    files = sorted(glob.glob(f"{folder}/*edep*.mhd"))
    if not files:
        print(f"[!] {folder}: no edep file found, skipped")
        continue
    arr = np.array(itk.array_view_from_image(itk.imread(files[0]))).reshape(-1)
    entr = arr[0] / (n * 0.5625)

    E = E_AIR if mat == "G4_AIR" else eout.get(mat, np.nan)
    if np.isnan(E):
        ref, dev = np.nan, np.nan
    else:
        ref = stopping(E, tabE, tabS) * THICK_WATER_CM
        dev = 100 * (entr - ref) / ref
    print(f"{folder:<22}{E:>16.2f}{entr:>16.4f}{ref:>16.4f}{dev:>8.2f}")
