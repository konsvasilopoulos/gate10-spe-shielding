"""
compare_nasa.py - Compares the GATE dose (Webber 1956, 20 g/cm2 Al) with NASA/TP-2009-215560.
Usage:  python compare_nasa.py spe_runs/webber1956_Al20
Writes: <base>/compare_nasa.csv and compare_nasa.png
"""
import os, sys, csv
import numpy as np
BASE = sys.argv[1] if len(sys.argv) > 1 else "spe_runs/webber1956_Al20"
REF = "benchmark/nasa_tp2009_215560_digitized.csv"

ref = np.genfromtxt(REF, delimiter=",", names=True, skip_header=3)  # 3 comment lines
sim = np.genfromtxt(os.path.join(BASE, "spe_dose_depth.csv"), delimiter=",", names=True)
z, d, s = sim["depth_g_cm2"], sim["dose_cGy"], sim["sigma_cGy"]

rows = []
print(f"{'Depth':>6}{'GATE (cGy)':>15}{'FLUKA':>8}{'HETC':>8}{'HZETRN':>8}{'GATE/FLUKA':>12}{'GATE/HETC':>11}{'GATE/HZETRN':>13}")
for r in ref:
    zr = r["depth_g_cm2"]
    i = int(np.argmin(np.abs(z - max(zr, z[0]))))        # for depth 0 use the first voxel (0-1 mm)
    g, gs = d[i], s[i]
    fl, he, hz = r["dose_FLUKA_cGy"], r["dose_HETC_cGy"], r["dose_HZETRN_cGy"]
    print(f"{zr:>6.0f}{g:>9.3f}±{gs:.3f}{fl:>8.2f}{he:>8.2f}{hz:>8.2f}{g/fl:>12.3f}{g/he:>11.3f}{g/hz:>13.3f}")
    rows.append([zr, g, gs, fl, he, hz, g / fl, g / he, g / hz])

with open(os.path.join(BASE, "compare_nasa.csv"), "w", newline="") as fo:
    w = csv.writer(fo)
    w.writerow(["depth_g_cm2", "GATE_cGy", "GATE_sigma_cGy", "FLUKA_cGy", "HETC_cGy", "HZETRN_cGy",
                "GATE_over_FLUKA", "GATE_over_HETC", "GATE_over_HZETRN"])
    w.writerows(rows)
r = np.array(rows)
print(f"\nGATE/FLUKA: {r[:,6].min():.3f}–{r[:,6].max():.3f} · GATE/HETC: {r[:,7].min():.3f}–{r[:,7].max():.3f} · "
      f"GATE/HZETRN: {r[:,8].min():.3f}–{r[:,8].max():.3f}")

try:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"width_ratios": [1.4, 1]})
    ax[0].plot(z, d, "k-", lw=1.5, label="GATE 10 (QGSP_BIC_EMZ)")
    ax[0].fill_between(z, d - s, d + s, color="k", alpha=0.2)
    for col, lab, mk in (("dose_FLUKA_cGy", "FLUKA", "o"), ("dose_HETC_cGy", "HETC-HEDS", "^"), ("dose_HZETRN_cGy", "HZETRN", "x")):
        ax[0].plot(ref["depth_g_cm2"], ref[col], mk, ms=7, label=f"{lab} (NASA TP-2009-215560)")
    ax[0].set_xlabel("Depth in water (g/cm$^2$)"); ax[0].set_ylabel("Dose (cGy)")
    ax[0].set_title("Feb 1956 SPE (Webber) behind 20 g/cm$^2$ Al"); ax[0].legend(); ax[0].grid(alpha=0.3)
    for k, lab, mk in ((6, "GATE/FLUKA", "o"), (7, "GATE/HETC-HEDS", "^"), (8, "GATE/HZETRN", "x")):
        ax[1].plot(r[:, 0], r[:, k], mk + "-", label=lab)
    ax[1].axhline(1, color="k", lw=0.8); ax[1].axhspan(0.95, 1.05, color="g", alpha=0.1, label="±5 %")
    ax[1].set_xlabel("Depth in water (g/cm$^2$)"); ax[1].set_ylabel("Dose ratio"); ax[1].legend(); ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(BASE, "compare_nasa.png"), dpi=150)
    print(f"Written: {BASE}/compare_nasa.csv, compare_nasa.png")
except ImportError:
    pass
