"""make_fig_geometry.py — Fig. 1: simulation geometry. (a) side view x–z to scale (Al 20 g/cm² shown),
(b) front view x–y. Values from run_spectrum.py (wide geometry).
Usage: python analysis/make_fig_geometry.py  → figures/geometry.png/.pdf"""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.lines import Line2D
plt.rcParams.update({"font.size": 8, "savefig.bbox": "tight"})
T = 7.41                      # cm, Al 20 g/cm²
ZS, ZW = -T - 20, 30.0
C = dict(src="#d62728", slab="#7f7f7f", water="#c6dbef", score="#6baed6", plane="#f2a900",
         sec="#7b3294", n="#1a9850", g="#ff7f00")

fig = plt.figure(figsize=(7.2, 4.4))
ax = fig.add_axes([0.0, 0.07, 0.69, 0.90]); bx = fig.add_axes([0.72, 0.50, 0.26, 0.40])

# ---------- (a) side view ----------
ax.add_patch(Rectangle((-60, -T), 120, T, fc=C["slab"], ec="k", lw=0.8))
ax.add_patch(Rectangle((-60, 0.0), 120, ZW, fc=C["water"], ec="k", lw=0.8))
ax.add_patch(Rectangle((-15, 0.0), 30, ZW, fc=C["score"], ec="#08519c", lw=0.9, alpha=0.55))
for z in np.arange(1, ZW, 1.0): ax.plot([-15, 15], [z, z], color="#08519c", lw=0.2, alpha=0.6)
ax.add_patch(Rectangle((-15, -0.35), 30, 0.7, fc=C["plane"], ec="k", lw=0.5, zorder=5))
ax.plot([-40, 40], [ZS, ZS], color=C["src"], lw=3.2, solid_capstyle="butt")
for x in np.linspace(-36, 36, 9):
    ax.add_patch(FancyArrowPatch((x, ZS + 0.6), (x, -T - 0.4), arrowstyle="-|>", mutation_scale=7, color=C["src"], lw=0.9))
# illustrative tracks (left half of the column, away from labels)
ax.plot([-9, -8.6, -8.0, -7.3], [-T, 0, 12, 24], color=C["src"], lw=1.3)                       # primary p
ax.plot([3, 3.2], [-T, -3.0], color=C["src"], lw=1.3); ax.plot(3.2, -3.0, "o", ms=2.5, color="k", zorder=6)
ax.plot([3.2, 6.5, 9.0], [-3.0, 6, 13], color=C["sec"], lw=1.1)                                 # secondary p
ax.plot([3.2, 26], [-3.0, 27], color=C["n"], lw=1.1, ls="--")                                 # neutron
t = np.linspace(0, 1, 120); gx = 3.2 + 36 * t; gz = -3.0 + 14 * t + 0.5 * np.sin(t * 45)
ax.plot(gx, gz, color=C["g"], lw=0.9)                                                           # gamma
# dimensions (left)
def dim(z0, z1, txt):
    ax.annotate("", xy=(-64, z0), xytext=(-64, z1), arrowprops=dict(arrowstyle="<->", lw=0.6))
    ax.text(-66, (z0 + z1) / 2, txt, ha="right", va="center", fontsize=7)
dim(ZS, -T, "20 cm air"); dim(-T, 0, "5, 10, 20\ng/cm²"); dim(0, ZW, "30 g/cm²\nwater")
# labels (right) with leader lines
def lab(xy, xytext, txt, color="k"):
    ax.annotate(txt, xy=xy, xytext=xytext, fontsize=7.3, color=color, va="center",
                arrowprops=dict(arrowstyle="-", lw=0.5, color="0.3"))
lab((40, ZS), (64, ZS), "Proton source 80 × 80 cm\nSPE spectrum, normal incidence", C["src"])
lab((60, -T / 2), (64, -5.5), "Shield 120 × 120 cm\nPE, LiH, water, PMMA,\nAl, Ti, Fe, Cu, Pb")
lab((15, 0), (64, 3.0), "Phase-space plane, 1 mm air,\n30 × 30 cm: fluence of all species")
lab((15, 16), (64, 13.5), "Scoring column 30 × 30 cm,\n1 mm voxels: D and H\n(LET-binned, ICRP-60 Q)", "#08306b")
lab((60, 25), (64, 25), "Water 120 × 120 × 30 cm")
for z, t_ in ((0.3, "skin 0–0.3"), (5, "5 (BFO)"), (10, "10 g/cm²")):
    ax.plot([-18, -15], [z, z], color="k", lw=0.6); ax.text(-18.8, z, t_, va="center", ha="right", fontsize=6.5)
ax.annotate("", xy=(-52, 24), xytext=(-52, 16), arrowprops=dict(arrowstyle="-|>", lw=0.8)); ax.text(-50.5, 20, "z", fontsize=8)
ax.set_xlim(-84, 120); ax.set_ylim(ZW + 1.5, ZS - 4.5); ax.set_aspect(2.0); ax.axis("off")
ax.text(-84, ZS - 3.3, "(a) Side view (Al 20 g/cm² shown; z drawn 2× x)", fontsize=8, fontweight="bold")
hl = [Line2D([], [], color=C["src"], lw=1.3), Line2D([], [], color=C["sec"], lw=1.1),
      Line2D([], [], color=C["n"], lw=1.1, ls="--"), Line2D([], [], color=C["g"], lw=0.9)]
fig.legend(hl, ["primary proton", "secondary proton", "neutron", "photon"], ncol=4, frameon=False,
           loc="lower left", bbox_to_anchor=(0.06, 0.04), fontsize=7.3)

# ---------- (b) front view ----------
bx.add_patch(Rectangle((-60, -60), 120, 120, fc=C["water"], ec="k", lw=0.8))
bx.add_patch(Rectangle((-40, -40), 80, 80, fc="none", ec=C["src"], lw=1.6))
bx.add_patch(Rectangle((-15, -15), 30, 30, fc=C["score"], ec="#08519c", lw=0.9, alpha=0.8))
bx.text(0, 50, "shield, water: 120 cm", ha="center", fontsize=6.8)
bx.text(0, 28, "source: 80 cm", ha="center", color=C["src"], fontsize=6.8)
bx.text(0, 0, "scoring\n30 cm", ha="center", va="center", fontsize=6.8)
bx.set_xlim(-62, 62); bx.set_ylim(-62, 62); bx.set_aspect("equal"); bx.axis("off")
bx.text(-62, 70, "(b) Front view (x–y)", fontsize=8, fontweight="bold")
fig.text(0.72, 0.44, "World: air, 3 × 3 × 3 m\nGATE 10.0.0, QGSP_BIC_EMZ\nProduction cut 0.1 mm\nSource 20 cm upstream\nof the shield",
         fontsize=7, va="top", linespacing=1.4)
for ext in ("png", "pdf"): fig.savefig(f"figures/geometry.{ext}", dpi=300)
print("Written: figures/geometry (.png + .pdf)")
