"""
seed_stats.py - Mean and spread of dR80 and E_out over the independent-seed runs
(run_validation_150MeV.sh), and comparison with PSTAR.

Usage:   python validation/seed_stats.py            (reads runs/s*/output_*)
         python validation/seed_stats.py runs_wide  (another base folder)

For each seed dR80 = R80(air) - R80(material), using the air run
OF THE SAME seed. Writes seed_stats_per_run.csv and seed_stats_summary.csv.

The .mhd files are read directly (no itk). E_out needs uproot.
"""
import os
import sys
import csv
import glob
import json
import numpy as np

BASE = sys.argv[1] if len(sys.argv) > 1 else "runs"
# CSVs of the base folder "runs" keep the plain names; other folders get a prefix.
PREFIX = "" if BASE.rstrip("/") == "runs" else os.path.basename(BASE.rstrip("/")) + "_"
REF = "G4_AIR"
E_CUT_FRAC = 0.3          # E_out: forward protons with E > 30% of E0

# PSTAR values (from compare_pstar.py) for the slab thicknesses of run_material.py.
# Valid ONLY for these thicknesses; if a thickness changes, run compare_pstar.py again.
PSTAR = {  # material: (thickness cm, E_out MeV, dR cm)
    "G4_POLYETHYLENE": (2.0, 138.85, 1.992),
    "G4_Al":           (2.0, 125.61, 4.212),
    "G4_Ti":           (2.0, 111.34, 6.433),
    "G4_Cu":           (1.0, 114.21, 6.002),
    "G4_PLEXIGLASS":   (4.0, 123.01, 4.629),
    "G4_Pb":           (0.4, 137.34, 2.253),
    "G4_Fe":           (1.1, 114.04, 6.027),
}


def read_mhd(path):
    """Reads a MetaImage (.mhd + .raw). Returns (1D depth profile, z step in mm)."""
    hdr = {}
    with open(path) as f:
        for line in f:
            if "=" in line:
                k, v = line.split("=", 1)
                hdr[k.strip()] = v.strip()
    types = {"MET_DOUBLE": "f8", "MET_FLOAT": "f4", "MET_INT": "i4", "MET_UINT": "u4"}
    dt = np.dtype(types[hdr["ElementType"]])
    msb = hdr.get("BinaryDataByteOrderMSB", hdr.get("ElementByteOrderMSB", "False"))
    dt = dt.newbyteorder(">" if msb.lower() == "true" else "<")
    dims = [int(x) for x in hdr["DimSize"].split()]
    spacing = [float(x) for x in hdr["ElementSpacing"].split()]
    raw = os.path.join(os.path.dirname(path), hdr["ElementDataFile"])
    arr = np.fromfile(raw, dtype=dt).astype(float)
    if arr.size != int(np.prod(dims)):
        raise ValueError(f"{path}: {arr.size} values, expected {np.prod(dims)}")
    # One voxel in X-Y: the whole array is the depth profile.
    if dims[0] != 1 or dims[1] != 1:
        arr = arr.reshape(dims[2], dims[1], dims[0]).sum(axis=(1, 2))
    return arr, spacing[2]


def r80(depth, prof):
    """Depth (cm) at the distal 80% of the maximum dose, with linear interpolation."""
    i = int(np.argmax(prof))
    lvl = 0.8 * prof[i]
    for j in range(i + 1, len(prof)):
        if prof[j] < lvl:
            return depth[j - 1] + (lvl - prof[j - 1]) * (depth[j] - depth[j - 1]) / (prof[j] - prof[j - 1])
    return np.nan


def e_out(folder, mat, e0):
    files = glob.glob(os.path.join(folder, f"phsp_{mat}*.root"))
    if not files:
        return np.nan
    import uproot
    with uproot.open(files[0]) as f:
        t = f[f.keys()[0]]
        p = t.arrays(["ParticleName", "KineticEnergy", "Direction_Z"], library="np")
    sel = (p["ParticleName"] == "proton") & (p["Direction_Z"] > 0) & (p["KineticEnergy"] > E_CUT_FRAC * e0)
    return float(np.median(p["KineticEnergy"][sel])) if sel.any() else np.nan


def load(folder):
    sj = os.path.join(folder, "setup.json")
    if not os.path.exists(sj):
        return None
    setup = json.load(open(sj))
    mat = setup["material"]
    mhd = os.path.join(folder, f"dose_{mat}_dose.mhd")
    if not os.path.exists(mhd):
        print(f"[!] Missing {mhd}")
        return None
    prof, dz = read_mhd(mhd)
    depth = (np.arange(len(prof)) + 0.5) * dz / 10.0
    return dict(setup=setup, mat=mat, r80=r80(depth, prof), e_out=e_out(folder, mat, setup["energy_MeV"]))


# --- Load all runs ---
runs = {}   # (seed, material) -> results
for sdir in sorted(glob.glob(os.path.join(BASE, "s*"))):
    for folder in sorted(glob.glob(os.path.join(sdir, "output_*"))):
        r = load(folder)
        if r is None:
            continue
        seed = r["setup"].get("seed", os.path.basename(sdir)[1:])
        cut = r["setup"].get("cut_mm")
        if cut is not None and abs(float(cut) - 0.1) > 1e-9:
            print(f"[!] {folder}: cut_mm = {cut}, not 0.1")
        runs[(str(seed), r["mat"])] = r

seeds = sorted({s for s, _ in runs})
if not seeds:
    sys.exit(f"No runs found in {BASE}/s*/output_*")

per_run = []
for s in seeds:
    ref = runs.get((s, REF))
    if ref is None:
        print(f"[!] seed {s}: air run missing, skipped")
        continue
    for (s2, m), r in runs.items():
        if s2 != s or m == REF:
            continue
        per_run.append(dict(seed=s, material=m, thickness_cm=r["setup"]["thickness_cm"],
                            r80_cm=r["r80"], dR80_cm=ref["r80"] - r["r80"], e_out_MeV=r["e_out"]))

with open(PREFIX + "seed_stats_per_run.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(per_run[0].keys()))
    w.writeheader()
    w.writerows(per_run)

# --- Summary per material ---
print(f"\nSeeds: {', '.join(seeds)}   (+- = standard deviation between seeds)\n")
hdr = (f"{'Material':<16}{'n':>3}{'dR80 cm':>16}{'dR PSTAR':>10}{'Diff. mm':>15}"
       f"{'E_out MeV':>18}{'E PSTAR':>9}{'Diff. %':>9}")
print(hdr)
print("-" * len(hdr))
summary = []
for m in sorted({r["material"] for r in per_run}):
    rows = [r for r in per_run if r["material"] == m]
    dr = np.array([r["dR80_cm"] for r in rows])
    eo = np.array([r["e_out_MeV"] for r in rows])
    n = len(rows)
    sd = (lambda x: float(np.std(x, ddof=1)) if len(x) > 1 else np.nan)
    dr_m, dr_s = float(np.mean(dr)), sd(dr)
    eo_m, eo_s = float(np.nanmean(eo)), sd(eo[~np.isnan(eo)])
    t, e_ps, dr_ps = PSTAR.get(m, (np.nan, np.nan, np.nan))
    if m in PSTAR and abs(rows[0]["thickness_cm"] - t) > 1e-6:
        print(f"[!] {m}: thickness {rows[0]['thickness_cm']} cm != {t} cm of the PSTAR table; comparison not valid")
        e_ps = dr_ps = np.nan
    d_mm = 10 * (dr_m - dr_ps)
    d_pct = 100 * (eo_m - e_ps) / e_ps
    print(f"{m:<16}{n:>3}{dr_m:>9.3f} ±{dr_s:5.3f}{dr_ps:>10.3f}{d_mm:>+8.2f} ±{10*dr_s:4.2f}"
          f"{eo_m:>10.2f} ±{eo_s:5.2f}{e_ps:>9.2f}{d_pct:>+9.2f}")
    summary.append(dict(material=m, n_seeds=n, dR80_mean_cm=dr_m, dR80_sd_cm=dr_s, dR_pstar_cm=dr_ps,
                        dR_diff_mm=d_mm, dR_diff_sd_mm=10 * dr_s, e_out_mean_MeV=eo_m,
                        e_out_sd_MeV=eo_s, e_pstar_MeV=e_ps, e_diff_pct=d_pct))

with open(PREFIX + "seed_stats_summary.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
    w.writeheader()
    w.writerows(summary)

print("\nThe +- mm is the spread of dR80 between seeds (statistical uncertainty of one run).")
print(f"Written: {PREFIX}seed_stats_per_run.csv, {PREFIX}seed_stats_summary.csv")
