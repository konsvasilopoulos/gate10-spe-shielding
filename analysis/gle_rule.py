"""gle_rule.py - Threshold rule: the dose ratio of a material to Al depends on the slope gamma of the integral
spectrum at the transmission threshold:  D_m/D_Al ~ (E_th,m / E_th,Al)^(-gamma(E*)),  E* = sqrt(E_th,m * E_th,Al),
gamma(E) = -d ln Phi(>E) / d ln E. Checked (a) against folding for 58 GLEs, (b) against the direct GATE runs,
(c) against folding of the OLTARIS spectra Feb 1956 (LaRC) and Aug 1972 (LaRC) (oltaris_larc.py).
The direct runs are read from spe_runs/*/materials_summary.csv; without them (repository without the Zenodo
dataset) their rows are taken from the existing analysis/gle_rule.csv.
Writes analysis/gle_rule.csv and figures/gle_rule.png/.pdf."""
import os, sys, csv
import numpy as np
sys.path.insert(0, os.getcwd()); sys.path.insert(0, "analysis")
from gle_table import F, T2
from spectra import SPECTRA
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FormatStrFormatter, NullFormatter, NullLocator
ETH = {}
for r in csv.DictReader(open("analysis/thresholds_pstar.csv")):
    ETH[(r["material"], int(float(r["areal_g_cm2"])))] = float(r["E_th_MeV"])
def gam(Ff, E): return -(np.log(Ff(E * 1.05)) - np.log(Ff(E / 1.05))) / np.log(1.05 ** 2)
def rule(Ff, m, t):
    ea, em = ETH[("Al", t)], ETH[(m, t)]; return (em / ea) ** (-gam(Ff, np.sqrt(ea * em)))
scan = list(csv.DictReader(open("analysis/fold_gle_scan.csv")))
out = []; dev = []
for r in scan:
    Ff = lambda E, g=int(r["GLE"]): F(E, g)
    for m, key in (("POLYETHYLENE", "PE"), ("Pb", "Pb")):
        for t in (5, 10, 20):
            fo = float(r[f"{key}_Al_{t}"]); ru = float(rule(Ff, m, t)); dev.append(fo / ru)
            out.append(["GLE" + r["GLE"], m, t, float(gam(Ff, np.sqrt(ETH[("Al", t)] * ETH[(m, t)]))), ru, fo, "fold"])
dev = np.array(dev); print(f"Folding / rule (58 GLEs x PE, Pb x 3 thicknesses): median {np.median(dev):.3f}, range {dev.min():.3f}-{dev.max():.3f}")
# direct GATE runs
PREV = [r for r in csv.DictReader(open("analysis/gle_rule.csv")) if r["source"] == "GATE"] if os.path.exists("analysis/gle_rule.csv") else []
DIR = {"kold1956": "spe_runs/materials_kold1956", "kold1972": "spe_runs/materials_kold1972", "kold1989": "spe_runs/materials_kold1989",
       "webber1956": "spe_runs/materials", "king1972": "spe_runs/materials_king1972"}
dd = []
for k, b in DIR.items():
    f = f"{b}/materials_summary.csv"
    if not os.path.exists(f):
        for r in (r for r in PREV if r["event"] == k):
            g_, ru = float(r["ratio"]), float(r["rule_ratio"]); dd.append((k, r["material"], int(r["areal"]), g_ / ru))
            out.append([k, r["material"], int(r["areal"]), float(r["gamma_at_Estar"]), ru, g_, "GATE"])
        continue
    s = {(x["material"], int(float(x["areal_g_cm2"]))): float(x["D_cGy_mean_0-30"]) for x in csv.DictReader(open(f))}
    for m in ("POLYETHYLENE", "PLEXIGLASS", "Ti", "Fe", "Cu", "Pb"):
        for t in (5, 10, 20):
            g_ = s[(m, t)] / s[("Al", t)]; ru = float(rule(SPECTRA[k], m, t)); dd.append((k, m, t, g_ / ru))
            out.append([k, m, t, float(gam(SPECTRA[k], np.sqrt(ETH[("Al", t)] * ETH[(m, t)]))), ru, g_, "GATE"])
dd_r = np.array([x[3] for x in dd])
print(f"GATE / rule ({len(dd)} cases, all materials): median {np.median(dd_r):.3f}, range {dd_r.min():.3f}-{dd_r.max():.3f}")
for k in DIR:
    v = [x[3] for x in dd if x[0] == k]
    if v: print(f"   {k:<11} {np.median(v):.3f} ({min(v):.3f}–{max(v):.3f})")
# OLTARIS LaRC spectra folded with the response matrix
import response_matrix as rm
from oltaris_larc import flux_integral
LARC = {"larc1956": flux_integral("LaRC56")[0], "larc1972": flux_integral("LaRC72")[0]}
for k, Ff in LARC.items():
    for m in ("POLYETHYLENE", "Pb"):
        for t in (5, 10, 20):
            fo = rm.fold(Ff, m, t, "D", "mean") / rm.fold(Ff, "Al", t, "D", "mean"); ru = float(rule(Ff, m, t))
            out.append([k, m, t, float(gam(Ff, np.sqrt(ETH[("Al", t)] * ETH[(m, t)]))), ru, fo, "fold-OLTARIS"])
            print(f"   {k} {m:<13}{t:>3} g/cm2: gamma {out[-1][3]:.2f}, fold {fo:.3f}, rule {ru:.3f}")
with open("analysis/gle_rule.csv", "w", newline="") as fo:
    w = csv.writer(fo); w.writerow(["event", "material", "areal", "gamma_at_Estar", "rule_ratio", "ratio", "source"]); w.writerows(out)
# figure
plt.rcParams.update({"font.size": 9, "savefig.bbox": "tight"})
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.3), sharey=True)
gg = np.linspace(0.5, 6.5, 100)
for a, t in zip(ax, (5, 20)):
    for m, c in (("POLYETHYLENE", "C0"), ("Pb", "C3")):
        x = [o[3] for o in out if o[1] == m and o[2] == t and o[6] == "fold"]; y = [o[5] for o in out if o[1] == m and o[2] == t and o[6] == "fold"]
        a.plot(x, y, ".", color=c, alpha=0.45, ms=5)
        ratio = ETH[(m, t)] / ETH[("Al", t)]; a.plot(gg, ratio ** (-gg), "-", color=c, lw=1, label=f"{'PE' if m.startswith('P') and m != 'Pb' else 'Pb'}: (E_th,m/E_th,Al)^−γ")
        for k, mk in (("kold1956", "o"), ("kold1972", "s"), ("kold1989", "D"), ("webber1956", "^"), ("king1972", "v"), ("larc1956", "P"), ("larc1972", "X")):
            pts = [o for o in out if o[0] == k and o[1] == m and o[2] == t]
            if pts: a.plot(pts[0][3], pts[0][5], mk, mfc="white", mec="k", ms=6, zorder=5)
    a.set_yscale("log"); a.yaxis.set_major_locator(FixedLocator([0.3, 0.5, 0.7, 1, 1.5, 2, 3, 4])); a.yaxis.set_major_formatter(FormatStrFormatter("%g"))
    a.yaxis.set_minor_locator(NullLocator()); a.set_ylim(0.27, 4.2); a.set_xlim(0.4, 6.5); a.axhline(1, color="0.5", lw=0.6); a.grid(alpha=0.3)
    ea = ETH[("Al", t)]; a.set_xlabel(f"Integral spectral index γ near E_th (≈{ea:.0f} MeV)"); a.set_title(f"{t} g/cm² shield", fontsize=9)
ax[0].set_ylabel("Mean dose relative to Al")
from matplotlib.lines import Line2D
h = [Line2D([], [], ls="", marker=".", color="0.4", ms=7), Line2D([], [], color="0.3", lw=1)] + \
    [Line2D([], [], ls="", marker=mk, mfc="white", mec="k") for mk in ("o", "s", "D", "^", "v", "P", "X")]
ax[1].legend(h, ["58 GLEs (folded)", "threshold rule", "Feb 1956 K.", "Aug 1972 K.", "Oct 1989 K.", "Webber 1956", "King 1972", "Feb 1956 LaRC", "Aug 1972 LaRC"],
             frameon=False, fontsize=6.5, loc="center right")
ax[0].text(0.5, 0.31, "PE (blue), Pb (red)", fontsize=7)
for a, t in zip(ax, ("(a)", "(b)")): a.text(0.02, 0.93, t, transform=a.transAxes, fontweight="bold")
fig.tight_layout()
for e in ("png", "pdf"): fig.savefig(f"figures/gle_rule.{e}", dpi=300)
print("Written: analysis/gle_rule.csv, figures/gle_rule.png/.pdf")
