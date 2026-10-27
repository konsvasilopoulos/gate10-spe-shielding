"""
analyze_spectrum.py - Combines the energy groups of an SPE run into the dose for the whole event.

Usage:  python analyze_spectrum.py spe_runs/webber1956_Al20

Reads every subfolder <base>/g*/ (one run per energy group, from run_materials_spe.sh).
Each run is multiplied by  Phi_group x A_source / N  (from setup.json), so that the dose
corresponds to the whole event (Gy; printed here in cGy as in NASA TP-2009-215560).

Writes to the base folder:
  spe_dose_depth.csv     dose vs depth: total, per group, per particle type (if --split was used)
  spe_fluence_plane.csv  fluence of all species at the water entrance (forward/backward), per MeV
  spe_summary.png        plots
"""
import os
import sys
import csv
import glob
import json
import numpy as np

BASE = sys.argv[1] if len(sys.argv) > 1 else "spe_runs/webber1956_Al20"
PLANE_AREA = 30.0 * 30.0          # cm2, scoring plane (DETECTOR_XY^2)
REPORT_DEPTHS = [0, 10, 20, 30]   # g/cm2 of water, as in the benchmark
TAGS = ["proton", "neutron", "photon_electron"]


def read_mhd(path):
    hdr = {}
    for line in open(path):
        if "=" in line:
            k, v = line.split("=", 1)
            hdr[k.strip()] = v.strip()
    dt = np.dtype({"MET_DOUBLE": "f8", "MET_FLOAT": "f4"}[hdr["ElementType"]])
    dt = dt.newbyteorder(">" if hdr.get("BinaryDataByteOrderMSB", "False").lower() == "true" else "<")
    dims = [int(x) for x in hdr["DimSize"].split()]
    dz = float(hdr["ElementSpacing"].split()[2])
    arr = np.fromfile(os.path.join(os.path.dirname(path), hdr["ElementDataFile"]), dtype=dt).astype(float)
    # One voxel in X-Y: the data are the depth profile. GATE writes the uncertainty file
    # with DimSize "N 1 1" instead of "1 1 N", so we do not reshape.
    if sorted(dims)[:2] != [1, 1]:
        arr = arr.reshape(dims[2], dims[1], dims[0]).sum(axis=(1, 2))
    return arr, dz


groups = []
for d in sorted(glob.glob(os.path.join(BASE, "g*")), key=lambda p: float(os.path.basename(p)[1:].split("-")[0])):
    sj = os.path.join(d, "setup.json")
    if not os.path.exists(sj):
        print(f"[!] {d}: setup.json missing (run not finished?) - skipped")
        continue
    s = json.load(open(sj))
    mat = s["material"]
    f = s["group_fluence_cm2"] * s["source_area_cm2"] / s["n"]      # event protons per simulated proton
    dose, dz = read_mhd(os.path.join(d, f"dose_{mat}_dose.mhd"))
    unc_path = os.path.join(d, f"dose_{mat}_dose_uncertainty.mhd")
    unc = read_mhd(unc_path)[0] if os.path.exists(unc_path) else np.zeros_like(dose)
    parts = {}
    for t in TAGS:
        p = os.path.join(d, f"dose_{mat}_{t}_dose.mhd")
        if os.path.exists(p):
            parts[t] = read_mhd(p)[0] * f
    groups.append(dict(dir=d, s=s, mat=mat, f=f, dose=dose * f, sig=np.nan_to_num(unc) * dose * f,
                       parts=parts, dz=dz))
    print(f"  {os.path.basename(d):<14} {s['emin_MeV']:>6.0f}–{s['emax_MeV']:<6.0f} MeV  "
          f"Phi = {s['group_fluence_cm2']:.3e} cm^-2  N = {s['n']}")

if not groups:
    sys.exit(f"No runs found in {BASE}/g*/")

# consistency checks
ref = groups[0]["s"]
for g in groups[1:]:
    for k in ("material", "areal_g_cm2", "spectrum", "phantom_cm", "dz_mm", "physics", "cut_mm"):
        if g["s"][k] != ref[k]:
            sys.exit(f"[!] Inconsistent {k}: {g['dir']} has {g['s'][k]}, {groups[0]['dir']} has {ref[k]}")
edges = sorted((g["s"]["emin_MeV"], g["s"]["emax_MeV"]) for g in groups)
for (a, b), (c, d) in zip(edges, edges[1:]):
    if abs(b - c) > 1e-9:
        print(f"[!] Groups are not contiguous: {a}-{b} and {c}-{d}")

nz = len(groups[0]["dose"])
dz = groups[0]["dz"]
depth = (np.arange(nz) + 0.5) * dz / 10.0                    # cm = g/cm2 in water
total = sum(g["dose"] for g in groups) * 100                  # cGy
sigma = np.sqrt(sum(g["sig"] ** 2 for g in groups)) * 100
has_split = all(len(g["parts"]) == len(TAGS) for g in groups)
parts = {t: sum(g["parts"][t] for g in groups) * 100 for t in TAGS} if has_split else {}

print(f"\n{ref['spectrum']} | {ref['material']} {ref['areal_g_cm2']} g/cm2 | water {ref['phantom_cm']} cm | {ref['physics']}")
hdr = f"{'Depth (g/cm2)':>14}{'Dose (cGy)':>13}{'+- (cGy)':>10}"
if has_split:
    hdr += "".join(f"{t[:8] + ' %':>11}" for t in TAGS) + f"{'ions etc %':>13}"
print(hdr)
for z in REPORT_DEPTHS:
    i = int(np.clip(np.searchsorted(depth, z), 0, nz - 1))
    line = f"{depth[i]:>14.2f}{total[i]:>13.4f}{sigma[i]:>10.4f}"
    if has_split:
        line += "".join(f"{100 * parts[t][i] / total[i]:>11.1f}" for t in TAGS)
        line += f"{100 - sum(100 * parts[t][i] / total[i] for t in TAGS):>13.1f}"
    print(line)

if has_split:
    print("(Percentages refer to the particle that DEPOSITS the energy. Neutrons deposit through recoil")
    print(" protons and ions, so their contribution appears under 'proton' and 'ions etc'.)")
print("\nContribution of each energy group to the entrance dose (first voxel):")
for g in sorted(groups, key=lambda g: g["s"]["emin_MeV"]):
    print(f"  {g['s']['emin_MeV']:>6.0f}–{g['s']['emax_MeV']:<6.0f} MeV: {100 * g['dose'][0] * 100 / total[0]:6.2f} %")

with open(os.path.join(BASE, "spe_dose_depth.csv"), "w", newline="") as fo:
    w = csv.writer(fo)
    gl = [f"{g['s']['emin_MeV']:.0f}-{g['s']['emax_MeV']:.0f}MeV_cGy" for g in groups]
    w.writerow(["depth_g_cm2", "dose_cGy", "sigma_cGy"] + [f"{t}_cGy" for t in parts] + gl)
    for i in range(nz):
        w.writerow([f"{depth[i]:.3f}", total[i], sigma[i]] + [parts[t][i] for t in parts]
                   + [g["dose"][i] * 100 for g in groups])

# --- Fluence at the scoring plane (water entrance) ---
ebins = np.geomspace(0.01, 2000, 101)   # MeV (for ions: total kinetic energy)
# Fluence per species (forward = towards the water, backward = from the water). All species are kept,
# because the secondary fluence is the input for later folding (e.g. a Geant4-DNA look-up table).
SPECIES = ["proton", "neutron", "alpha", "deuteron", "triton", "He3", "e-", "gamma"]
flu = {f"{sp}_{dr}": np.zeros(len(ebins) - 1) for sp in SPECIES + ["ions_Z3plus"] for dr in ("fwd", "bwd")}
try:
    import uproot
    for g in groups:
        files = glob.glob(os.path.join(g["dir"], "phsp_*.root"))
        if not files:
            continue
        with uproot.open(files[0]) as fr:
            t = fr[fr.keys()[0]]
            p = t.arrays(["ParticleName", "KineticEnergy", "Direction_Z"], library="np")
        cosz = np.clip(np.abs(p["Direction_Z"]), 0.05, None)    # fluence = sum 1/|cos(theta)| / A
        names = p["ParticleName"].astype(str)
        known = np.isin(names, SPECIES + ["e+"])
        for dr, sel_dir in (("fwd", p["Direction_Z"] > 0), ("bwd", p["Direction_Z"] < 0)):
            for sp in SPECIES + ["ions_Z3plus"]:
                sel = (~known if sp == "ions_Z3plus" else names == sp) & sel_dir
                h, _ = np.histogram(p["KineticEnergy"][sel], bins=ebins, weights=1.0 / cosz[sel])
                flu[f"{sp}_{dr}"] += h * g["f"] / PLANE_AREA
    de = np.diff(ebins)
    with open(os.path.join(BASE, "spe_fluence_plane.csv"), "w", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(["E_low_MeV", "E_high_MeV"] + [f"{k}_per_cm2_MeV" for k in flu])
        for i in range(len(de)):
            w.writerow([ebins[i], ebins[i + 1]] + [flu[k][i] / de[i] for k in flu])
    have_flu = True
except ImportError:
    print("[!] uproot not found; plane fluence skipped.")
    have_flu = False

try:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3 if have_flu else 2, figsize=(17 if have_flu else 12, 4.8))
    ax[0].plot(depth, total, "k-", label="Total")
    ax[0].fill_between(depth, total - sigma, total + sigma, color="k", alpha=0.2)
    for t in parts:
        ax[0].plot(depth, parts[t], label=t)
    ax[0].set_xlabel("Depth in water (g/cm$^2$)"); ax[0].set_ylabel("Dose (cGy)"); ax[0].set_yscale("log")
    ax[0].set_title("Dose vs depth"); ax[0].legend(); ax[0].grid(alpha=0.3)
    for g in groups:
        ax[1].plot(depth, g["dose"] * 100, label=f"{g['s']['emin_MeV']:.0f}–{g['s']['emax_MeV']:.0f} MeV")
    ax[1].set_xlabel("Depth in water (g/cm$^2$)"); ax[1].set_ylabel("Dose (cGy)"); ax[1].set_yscale("log")
    ax[1].set_title("Contribution per energy group"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    if have_flu:
        ec = np.sqrt(ebins[:-1] * ebins[1:]); de = np.diff(ebins)
        for k, lab in (("proton_fwd", "p forward"), ("neutron_fwd", "n forward"), ("neutron_bwd", "n backward")):
            ax[2].loglog(ec, flu[k] / de, drawstyle="steps-mid", label=lab)
        ax[2].set_xlabel("E (MeV)"); ax[2].set_ylabel("Fluence (cm$^{-2}$ MeV$^{-1}$)")
        ax[2].set_title("Fluence at the water entrance"); ax[2].legend(); ax[2].grid(alpha=0.3, which="both")
    fig.suptitle(f"{ref['spectrum']} · {ref['material']} {ref['areal_g_cm2']} g/cm²")
    fig.tight_layout(); fig.savefig(os.path.join(BASE, "spe_summary.png"), dpi=150)
    print(f"\nWritten: {BASE}/spe_dose_depth.csv, spe_fluence_plane.csv, spe_summary.png")
except ImportError:
    pass
