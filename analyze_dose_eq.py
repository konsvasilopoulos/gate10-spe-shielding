"""
analyze_dose_eq.py - Dose equivalent H (ICRP-60) from run_spectrum.py runs made with --dose-eq,
and comparison with NASA/TP-2009-215560 (their Fig. 4).

Usage:  python analyze_dose_eq.py spe_runs/webber1956_Al20_wide_H

Method: each LETActor (particle x energy bin) gives the numerator N = sum Edep*L and the
denominator E = sum Edep per voxel. The dose-averaged LET of the bin is Lbar = N/E (keV/um, since
1 MeV/mm = 1 keV/um). With E_tot the total Edep of the voxel (DoseActor):
    H = D * [ 1 + sum_bins E_bin * (Q(Lbar_bin) - 1) / E_tot ]
Where no actor exists (e.g. protons > 8 MeV, electrons) L < 10 keV/um, so Q = 1.
"""
import os, sys, glob, json, csv
import numpy as np

BASE = sys.argv[1] if len(sys.argv) > 1 else "spe_runs/webber1956_Al20_wide_H"
REF = "benchmark/nasa_tp2009_215560_digitized.csv"
DEPTHS = [0, 5, 10, 15, 20, 25, 30]


def Q_icrp60(L):
    L = np.asarray(L, dtype=float)
    q = np.ones_like(L)
    m = (L >= 10) & (L <= 100)
    q[m] = 0.32 * L[m] - 2.2
    m = L > 100
    q[m] = 300.0 / np.sqrt(L[m])
    return q


def read_mhd(path):
    hdr = {}
    for line in open(path):
        if "=" in line:
            k, v = line.split("=", 1)
            hdr[k.strip()] = v.strip()
    dt = np.dtype({"MET_DOUBLE": "f8", "MET_FLOAT": "f4"}[hdr["ElementType"]])
    dt = dt.newbyteorder(">" if hdr.get("BinaryDataByteOrderMSB", "False").lower() == "true" else "<")
    dims = [int(x) for x in hdr["DimSize"].split()]
    arr = np.fromfile(os.path.join(os.path.dirname(path), hdr["ElementDataFile"]), dtype=dt).astype(float)
    if sorted(dims)[:2] != [1, 1]:
        arr = arr.reshape(dims[2], dims[1], dims[0]).sum(axis=(1, 2))
    return arr, float(hdr["ElementSpacing"].split()[2])


def species(name):
    b = os.path.basename(name)            # let_<particle>_<bin>... or let_rest_all...
    return b.split("_")[1]


D = H = None
dH_by = {}
for d in sorted(glob.glob(os.path.join(BASE, "g*"))):
    sj = os.path.join(d, "setup.json")
    if not os.path.exists(sj):
        print(f"[!] {d}: setup.json missing"); continue
    s = json.load(open(sj))
    if not s.get("dose_eq"):
        print(f"[!] {d}: run was not made with --dose-eq"); continue
    mat = s["material"]
    f = s["group_fluence_cm2"] * s["source_area_cm2"] / s["n"] * 100          # -> cGy (cSv)
    dose, dz = read_mhd(os.path.join(d, f"dose_{mat}_dose.mhd"))
    edep, _ = read_mhd(os.path.join(d, f"dose_{mat}_edep.mhd"))
    nums = sorted(glob.glob(os.path.join(d, "let_*numerator*.mhd")))
    if not nums:
        print(f"[!] {d}: no let_*numerator*.mhd files found"); continue
    extra = np.zeros_like(dose)                  # sum E_bin (Q - 1)
    e_scored = np.zeros_like(dose)
    for npath in nums:
        dpath = npath.replace("numerator", "denominator")
        num, _ = read_mhd(npath); den, _ = read_mhd(dpath)
        with np.errstate(divide="ignore", invalid="ignore"):
            Lbar = np.where(den > 0, num / den, 0.0)
        term = den * (Q_icrp60(Lbar) - 1.0)
        extra += term
        e_scored += den
        sp = species(npath)
        with np.errstate(divide="ignore", invalid="ignore"):
            dH_by[sp] = dH_by.get(sp, 0) + np.where(edep > 0, dose * term / edep, 0) * f
    if np.any(e_scored > edep * 1.001 + 1e-30):
        print(f"[!] {d}: sum Edep of the LET actors > total Edep in some voxels (overlapping bins?)")
    with np.errstate(divide="ignore", invalid="ignore"):
        h = dose * (1 + np.where(edep > 0, extra / edep, 0))
    D = dose * f if D is None else D + dose * f
    H = h * f if H is None else H + h * f
    print(f"  {os.path.basename(d):<12} {s['emin_MeV']:>5.0f}–{s['emax_MeV']:<5.0f} MeV  "
          f"{len(nums)} LET actors")

if D is None:
    sys.exit("No runs made with --dose-eq found")
z = (np.arange(len(D)) + 0.5) * dz / 10

ref = None
if os.path.exists(REF):
    ref = np.genfromtxt(REF, delimiter=",", names=True, skip_header=3)

# H in a single 1 mm voxel is noisy (rare high-Q events), so we report averages over a depth
# window: 0-0.3 g/cm2 for '0', 29-30 for '30' and +-1 g/cm2 elsewhere. The +- is the standard
# error of the mean from the spread of the voxels.
def window(zc):
    if zc == 0:
        return z < 0.3
    if zc >= 30:
        return z > 29.0
    return np.abs(z - zc) <= 1.0


print(f"\n{'Depth':>6}{'D (cGy)':>10}{'H (cSv)':>16}{'Qbar':>7}", end="")
if ref is not None:
    print(f"{'H FLUKA':>9}{'H HETC*':>9}{'H HZETRN':>10}{'GATE/FLUKA':>12}{'GATE/HETC*':>12}{'GATE/HZETRN':>13}")
else:
    print()
rows = []
for k, zr in enumerate(DEPTHS):
    m = window(zr)
    d_, h_ = D[m].mean(), H[m].mean()
    se = H[m].std(ddof=1) / np.sqrt(m.sum()) if m.sum() > 1 else np.nan
    line = f"{zr:>6}{d_:>10.3f}{h_:>9.3f}±{se:<6.3f}{h_/d_:>7.2f}"
    row = [zr, d_, h_, se, h_ / d_]
    if ref is not None:
        fl, he, hz = ref["H_FLUKA_cSv"][k], ref["H_HETC_corr_cSv"][k], ref["H_HZETRN_cSv"][k]
        line += f"{fl:>9.2f}{he:>9.2f}{hz:>10.2f}{h_/fl:>12.3f}{h_/he:>12.3f}{h_/hz:>13.3f}"
        row += [fl, he, hz]
    print(line); rows.append(row)
with open(os.path.join(BASE, "dose_eq_compare_nasa.csv"), "w", newline="") as fo:
    w = csv.writer(fo)
    w.writerow(["depth_g_cm2", "D_cGy", "H_cSv", "H_se_cSv", "Qbar"] + (["H_FLUKA", "H_HETC_corr", "H_HZETRN"] if ref is not None else []))
    w.writerows(rows)
print("(HETC* = HETC-HEDS with the corrected neutron coefficients, 'Corrected' in their Fig. 4)")

print("\nContribution to (H - D) per particle, %:")
tot = H - D
print(f"{'Depth':>6}" + "".join(f"{sp:>10}" for sp in dH_by))
for zr in DEPTHS:
    m = window(zr)
    print(f"{zr:>6}" + "".join(f"{100 * dH_by[sp][m].sum() / tot[m].sum():>10.1f}" for sp in dH_by))

with open(os.path.join(BASE, "dose_eq_depth.csv"), "w", newline="") as fo:
    w = csv.writer(fo)
    w.writerow(["depth_g_cm2", "D_cGy", "H_cSv", "Qbar"] + [f"dH_{sp}_cSv" for sp in dH_by])
    for i in range(len(z)):
        w.writerow([f"{z[i]:.3f}", D[i], H[i], H[i] / D[i] if D[i] > 0 else 0] + [dH_by[sp][i] for sp in dH_by])
print(f"\nWritten: {BASE}/dose_eq_depth.csv")
