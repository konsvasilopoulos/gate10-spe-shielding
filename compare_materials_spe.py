"""
compare_materials_spe.py - Comparison of shielding materials for one SPE spectrum (wide geometry).
For each (material, g/cm2) it computes D and H (ICRP-60) with analyze_dose_eq.py and reports:
  - D, H and mean Q at the skin (0-0.3 g/cm2), at 5 (BFO, 4-6) and 10 (9-11) g/cm2 of water,
    and averaged over the whole phantom (0-30 g/cm2)
  - the ratio to Al at the same areal density
Usage:  python compare_materials_spe.py spe_runs/materials_kold1956  -> materials_summary.csv, materials_DH.png
"""
import os, glob, csv, subprocess, sys
import numpy as np

BASE = sys.argv[1] if len(sys.argv) > 1 else "spe_runs/materials"
ORDER = ["POLYETHYLENE", "PLEXIGLASS", "Al", "Ti", "Fe", "Cu", "Pb"]
POINTS = {"skin_0-0.3": (0, 0.3), "5_gcm2": (4, 6), "10_gcm2": (9, 11), "mean_0-30": (0, 30)}

res = {}
for d in sorted(glob.glob(os.path.join(BASE, "*_*"))):
    name = os.path.basename(d)
    if name == "logs" or not glob.glob(os.path.join(d, "g*/setup.json")):
        continue
    mat, ar = name.rsplit("_", 1)
    n = len(glob.glob(os.path.join(d, "g*/setup.json")))
    if n < 6:
        print(f"[!] {name}: {n}/6 groups - skipped"); continue
    subprocess.run([sys.executable, "analyze_dose_eq.py", d], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    f = os.path.join(d, "dose_eq_depth.csv")
    if not os.path.exists(f):
        print(f"[!] {name}: analysis failed"); continue
    a = np.genfromtxt(f, delimiter=",", names=True)
    res[(mat, float(ar))] = (a["depth_g_cm2"], a["D_cGy"], a["H_cSv"])

if not res:
    sys.exit("No complete runs found in " + BASE)

def val(z, y, lo, hi):
    m = (z >= lo) & (z <= hi)
    return y[m].mean()

areals = sorted({k[1] for k in res})
mats = [m for m in ORDER if any(k[0] == m for k in res)] + sorted({k[0] for k in res} - set(ORDER))
rows = []
for ar in areals:
    print(f"\n=== {ar:g} g/cm2 ===")
    print(f"{'Material':<14}" + "".join(f"{'D ' + p:>22}{'H':>9}{'Qbar':>6}" for p in POINTS))
    for m in mats:
        if (m, ar) not in res: continue
        z, D, H = res[(m, ar)]
        line = f"{m:<14}"; row = [m, ar]
        for p, (lo, hi) in POINTS.items():
            dv, hv = val(z, D, lo, hi), val(z, H, lo, hi)
            line += f"{dv:>22.3f}{hv:>9.3f}{hv/dv:>6.2f}"; row += [dv, hv]
        print(line); rows.append(row)
    if ("Al", ar) in res:
        za, Da, Ha = res[("Al", ar)]
        lo, hi = POINTS["mean_0-30"]
        parts = []
        for m in mats:
            if (m, ar) not in res: continue
            z, D, H = res[(m, ar)]
            parts.append(f"{m}: D {val(z, D, lo, hi) / val(za, Da, lo, hi):.2f}, H {val(z, H, lo, hi) / val(za, Ha, lo, hi):.2f}")
        print("  Ratio to Al (mean 0-30): " + " | ".join(parts))

hdr = ["material", "areal_g_cm2"] + [f"{q}_{p}" for p in POINTS for q in ("D_cGy", "H_cSv")]
with open(os.path.join(BASE, "materials_summary.csv"), "w", newline="") as fo:
    w = csv.writer(fo); w.writerow(hdr); w.writerows(rows)

try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    for k, (p, (lo, hi)) in enumerate([("skin_0-0.3", POINTS["skin_0-0.3"]), ("mean_0-30", POINTS["mean_0-30"])]):
        for m in mats:
            xs = [ar for ar in areals if (m, ar) in res]
            ax[0, k].plot(xs, [val(res[(m, a)][0], res[(m, a)][1], lo, hi) for a in xs], "o-", label=m)
            ax[1, k].plot(xs, [val(res[(m, a)][0], res[(m, a)][2], lo, hi) for a in xs], "o-", label=m)
        ax[0, k].set_title(f"Absorbed dose - {p}"); ax[1, k].set_title(f"Dose equivalent - {p}")
        ax[0, k].set_ylabel("D (cGy)"); ax[1, k].set_ylabel("H (cSv)")
        for a in (ax[0, k], ax[1, k]):
            a.set_xlabel("Shield areal density (g/cm$^2$)"); a.set_yscale("log"); a.grid(alpha=.3, which="both"); a.legend(fontsize=8)
    fig.suptitle(f"{os.path.basename(BASE)}: shielding materials (GATE 10, QGSP_BIC_EMZ)")
    fig.tight_layout(); fig.savefig(os.path.join(BASE, "materials_DH.png"), dpi=150)
    print(f"\nWritten: {BASE}/materials_summary.csv, materials_DH.png")
except ImportError:
    pass
