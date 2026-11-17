"""fig_response.py - Response functions: mean D and H (0-30 g/cm2 water) per unit fluence vs proton energy."""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FormatStrFormatter, NullFormatter
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "savefig.bbox": "tight"})
M = np.load("analysis/response_matrix.npz"); e = M["edges"]; ec = np.sqrt(e[:-1] * e[1:]); z = M["z"]
MATS = [("POLYETHYLENE", "PE"), ("PLEXIGLASS", "PMMA"), ("Al", "Al"), ("Ti", "Ti"), ("Fe", "Fe"), ("Cu", "Cu"), ("Pb", "Pb")]
fig, ax = plt.subplots(1, 3, figsize=(7.48, 2.8), sharey=True)
for a, t in zip(ax, (5, 10, 20)):
    for i, (m, lab) in enumerate(MATS):
        D = M[f"D_{m}_{t}"].mean(axis=1) * 1e9; H = M[f"H_{m}_{t}"].mean(axis=1) * 1e9      # nGy·cm² (nSv·cm²)
        a.loglog(ec, D, "-", color=f"C{i}", lw=1.2, label=lab); a.loglog(ec, H, ":", color=f"C{i}", lw=1)
    a.set_title(f"{t} g/cm² shield"); a.set_xlabel("Proton energy (MeV)"); a.grid(alpha=0.3, which="both")
    a.xaxis.set_major_locator(FixedLocator([30, 100, 300, 1000, 3000])); a.xaxis.set_major_formatter(FormatStrFormatter("%g")); a.xaxis.set_minor_formatter(NullFormatter())
    a.set_xlim(30, 5000)
ax[0].set_ylabel("Mean dose in 0–30 g/cm² water\nper unit fluence (nGy cm²)"); ax[0].legend(frameon=False, fontsize=6.5, ncol=2, loc="lower right")
ax[2].text(0.97, 0.05, "solid: D, dotted: H (nSv cm²)", transform=ax[2].transAxes, ha="right", fontsize=6.5)
for a, t in zip(ax, ("(a)", "(b)", "(c)")): a.text(0.03, 0.92, t, transform=a.transAxes, fontweight="bold")
fig.tight_layout()
for x in ("png", "pdf"): fig.savefig(f"figures/response_functions.{x}", dpi=300)
print("Written: figures/response_functions")
