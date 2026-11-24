"""make_figures.py — Paper 1 figures (English), PNG (300 dpi) + PDF (vector), in figures/.
  benchmark_D_H          absorbed dose and dose equivalent vs FLUKA, HETC-HEDS, HZETRN (Webber 1956, 20 g/cm² Al)
  materials_ratio_4spectra   dose ratio to Al per material and areal density, four SPE spectra
  pe_pb_vs_spectrum      PE/Al and Pb/Al mean dose vs shield areal density, four spectra
Usage (from the repository root): python analysis/make_figures.py"""
import csv, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FormatStrFormatter, NullFormatter, NullLocator

# Sizes are the final print sizes (Elsevier: 190 mm double column, 140 mm 1.5 column); no text below 6.5 pt.
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "legend.fontsize": 7, "xtick.labelsize": 7,
                     "ytick.labelsize": 7, "axes.linewidth": 0.7, "lines.linewidth": 1.2, "savefig.bbox": "tight"})
W2, W15 = 7.48, 5.51          # inches: double column (190 mm), 1.5 column (140 mm)
os.makedirs("figures", exist_ok=True)

def log_ticks(ax, ticks):
    """Log y-axis with plain-number labels at the given values (no 10^n)."""
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
    ax.yaxis.set_minor_locator(NullLocator()); ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_ylim(min(ticks) * 0.9, max(ticks) * 1.1)

def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(f"figures/{name}.{ext}", dpi=300)
    plt.close(fig)

# ---------------- benchmark ----------------
B = "spe_runs/webber1956_Al20_wide_H"
nasa = np.genfromtxt("benchmark/nasa_tp2009_215560_digitized.csv", delimiter=",", names=True, skip_header=3)
dd = np.genfromtxt(f"{B}/spe_dose_depth.csv", delimiter=",", names=True)
he = np.genfromtxt(f"{B}/dose_eq_depth.csv", delimiter=",", names=True)
fig, ax = plt.subplots(1, 2, figsize=(W2, 3.0))
ax[0].plot(dd["depth_g_cm2"], dd["dose_cGy"], "k-", lw=1.4, label="GATE 10 (this work)")
ax[0].fill_between(dd["depth_g_cm2"], dd["dose_cGy"] - dd["sigma_cGy"], dd["dose_cGy"] + dd["sigma_cGy"], color="k", alpha=0.2, lw=0)
for k, lab, m in (("dose_FLUKA_cGy", "FLUKA", "o"), ("dose_HETC_cGy", "HETC-HEDS", "s"), ("dose_HZETRN_cGy", "HZETRN", "^")):
    ax[0].plot(nasa["depth_g_cm2"], nasa[k], m, mfc="none", ms=4.5, label=lab)
ax[0].set_ylabel("Absorbed dose (cGy)"); log_ticks(ax[0], [0.4, 0.5, 0.7, 1, 2, 3, 5, 7, 10]); ax[0].set_ylim(0.42, 9)
Hs = np.convolve(he["H_cSv"], np.ones(21) / 21, mode="same"); ok = (he["depth_g_cm2"] > 1) & (he["depth_g_cm2"] < 29)
ax[1].plot(he["depth_g_cm2"], he["H_cSv"], "-", color="0.75", lw=0.6, label="GATE 10, 1 mm bins")
ax[1].plot(he["depth_g_cm2"][ok], Hs[ok], "k-", lw=1.4, label="GATE 10, ±1 g/cm² mean")
for k, lab, m in (("H_FLUKA_cSv", "FLUKA", "o"), ("H_HETC_corr_cSv", "HETC-HEDS (corr.)", "s"), ("H_HZETRN_cSv", "HZETRN", "^")):
    ax[1].plot(nasa["depth_g_cm2"], nasa[k], m, mfc="none", ms=4.5, label=lab)
ax[1].set_ylabel("Dose equivalent (cSv)"); log_ticks(ax[1], [0.7, 1, 2, 3, 5, 7, 10, 15]); ax[1].set_ylim(0.65, 14)
for a, t in zip(ax, ("(a)", "(b)")):
    a.set_xlabel("Depth in water (g/cm²)"); a.set_xlim(-0.5, 30.5); a.grid(alpha=0.3); a.legend(frameon=False)
    a.text(0.02, 0.03, t, transform=a.transAxes, fontweight="bold")
save(fig, "benchmark_D_H")

# ---------------- materials ----------------
SP = [("kold1956", "Feb 1956 — Koldobskiy", "spe_runs/materials_kold1956"),
      ("kold1972", "Aug 1972 — Koldobskiy", "spe_runs/materials_kold1972"),
      ("webber1956", "Feb 1956 — Webber", "spe_runs/materials"),
      ("king1972", "Aug 1972 — King", "spe_runs/materials_king1972"),
      ("kold1989", "Sep–Oct 1989 — Koldobskiy", "spe_runs/materials_kold1989")]
SP = [s_ for s_ in SP if os.path.exists(f"{s_[2]}/materials_summary.csv")]
MATS = ["LITHIUM_HYDRIDE", "POLYETHYLENE", "WATER", "PLEXIGLASS", "Al", "Ti", "Fe", "Cu", "Pb"]
LAB = {"POLYETHYLENE": "PE", "PLEXIGLASS": "PMMA", "LITHIUM_HYDRIDE": "LiH", "WATER": "H₂O"}
AR = [5, 10, 20]
def load(base):
    return {(r["material"], float(r["areal_g_cm2"])): (float(r["D_cGy_mean_0-30"]), float(r["H_cSv_mean_0-30"]))
            for r in csv.DictReader(open(f"{base}/materials_summary.csv"))}
data = {k: load(b) for k, _, b in SP}
RT = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1, 1.25, 1.5, 2, 3, 4]
ncol = 3 if len(SP) > 4 else 2
fig, axs = plt.subplots(2, ncol, figsize=(W2, 4.9), sharex=True, sharey=True)
order = ["kold1956", "kold1972", "kold1989", "webber1956", "king1972"] if ncol == 3 else ["kold1956", "kold1972", "webber1956", "king1972"]
SP = sorted(SP, key=lambda s_: order.index(s_[0]))
cols = ["#1b9e77", "#d95f02", "#7570b3"]
for a, (k, lab, _), t in zip(axs.flat, SP, ("(a)", "(b)", "(c)", "(d)", "(e)")):
    x = np.arange(len(MATS))
    for j, ar in enumerate(AR):
        D = [data[k][(m, ar)][0] / data[k][("Al", ar)][0] for m in MATS]
        H = [data[k][(m, ar)][1] / data[k][("Al", ar)][1] for m in MATS]
        a.plot(x, H, "-o", color=cols[j], ms=3.5, lw=1.1, label=f"{ar} g/cm², H")
        a.plot(x, D, "--", color=cols[j], lw=1, label=f"{ar} g/cm², D")
    a.axhline(1, color="0.5", lw=0.7); a.set_xticks(x); a.set_xticklabels([LAB.get(m, m) for m in MATS], rotation=45, ha="right", rotation_mode="anchor")
    log_ticks(a, RT); a.set_ylim(0.27, 3.8); a.grid(alpha=0.3)
    a.set_title(lab); a.text(0.03, 0.92, t, transform=a.transAxes, fontweight="bold")
for a in axs[:, 0]: a.set_ylabel("Ratio to Al (mean, 0–30 g/cm² water)")
h_, l_ = axs[0, 0].get_legend_handles_labels()
for i_, a in enumerate(list(axs.flat)[len(SP):]):
    a.axis("off"); axs[0, ncol - 1 - i_].tick_params(labelbottom=True)
    if i_ == 0: a.legend(h_, l_, frameon=False, ncol=1, loc="center")
if len(SP) == len(axs.flat): axs[0, 0].legend(frameon=False, ncol=1, loc="lower right")
for a in axs[1, :]: a.tick_params(labelbottom=True)
fig.tight_layout(); save(fig, "materials_ratio_4spectra")

fig, ax = plt.subplots(figsize=(4.6, 3.1))
styles = {"kold1956": "C0", "kold1972": "C1", "webber1956": "C2", "king1972": "C3", "kold1989": "C4"}
for k, lab, _ in SP:
    ax.plot(AR, [data[k][("POLYETHYLENE", t)][0] / data[k][("Al", t)][0] for t in AR], "-o", color=styles[k], ms=4, label=f"PE, {lab}")
    ax.plot(AR, [data[k][("Pb", t)][0] / data[k][("Al", t)][0] for t in AR], "--s", color=styles[k], ms=4, label=f"Pb, {lab}")
ax.axhline(1, color="0.5", lw=0.7); log_ticks(ax, RT); ax.set_ylim(0.27, 3.8); ax.set_xticks(AR)
ax.set_xlabel("Shield areal density (g/cm²)"); ax.set_ylabel("Mean dose relative to Al"); ax.grid(alpha=0.3)
ax.legend(frameon=False, fontsize=7, ncol=1, loc="center left", bbox_to_anchor=(1.01, 0.5))
save(fig, "pe_pb_vs_spectrum")
print("Written: figures/benchmark_D_H, materials_ratio_4spectra, pe_pb_vs_spectrum (.png + .pdf)")

# ---------------- BFO dose vs NASA 30-day limit ----------------
SPL = [("kold1956", "Feb 1956 — Koldobskiy", "spe_runs/materials_kold1956"), ("webber1956", "Feb 1956 — Webber", "spe_runs/materials"),
       ("kold1972", "Aug 1972 — Koldobskiy", "spe_runs/materials_kold1972"), ("king1972", "Aug 1972 — King", "spe_runs/materials_king1972"),
       ("kold1989", "Sep–Oct 1989 — Koldobskiy", "spe_runs/materials_kold1989")]
SPL = [s for s in SPL if os.path.exists(f"{s[2]}/materials_summary.csv")]
fig, axs = plt.subplots(1, len(SPL), figsize=(W2, 2.6), sharey=True)
for a, (k, lab, b) in zip(np.atleast_1d(axs), SPL):
    d = {(r["material"], float(r["areal_g_cm2"])): float(r["D_cGy_5_gcm2"]) * 10 * 1.5 for r in csv.DictReader(open(f"{b}/materials_summary.csv"))}
    fo = {r["material"]: [float(r[f"BFO_{t}_mGyEq"]) for t in (5, 10, 20, 30, 40)]
          for r in csv.DictReader(open("analysis/nasa_limits_thick.csv")) if r["spectrum"] == k}
    for m, c, mk in (("POLYETHYLENE", "C0", "o"), ("Al", "C7", "s"), ("Pb", "C3", "^")):
        a.plot([5, 10, 20, 30, 40], fo[m], "-", color=c, lw=1.2)                     # response matrix 5-40
        a.plot(AR, [d[(m, t)] for t in AR], mk, color=c, ms=4, label=LAB.get(m, m))  # direct simulations
    a.axhline(250, color="k", ls=":", lw=1); a.set_title(lab.replace(" — ", "\n"), fontsize=7.5); a.set_xticks([5, 10, 20, 30, 40])
    log_ticks(a, [2, 5, 10, 20, 50, 100, 250, 500, 1000, 2000]); a.set_ylim(1.5, 2600); a.grid(alpha=0.3)
    a.set_xlabel("Areal density (g/cm²)"); a.tick_params(axis="x", labelsize=6.5)
np.atleast_1d(axs)[0].set_ylabel("BFO dose at 5 g/cm² (mGy-Eq)"); np.atleast_1d(axs)[0].legend(frameon=False, loc="lower left")
np.atleast_1d(axs)[min(1, len(SPL) - 1)].text(16, 300, "NASA 30-day limit", fontsize=6.5)
fig.tight_layout(); save(fig, "bfo_vs_limit")
print("Written: figures/bfo_vs_limit (.png + .pdf)")

# ---------------- Validation at 150 MeV: depth-dose, ΔR80 vs PSTAR, nuclear attenuation ----------------
def read_mhd(path):
    hdr = {l.split("=")[0].strip(): l.split("=", 1)[1].strip() for l in open(path) if "=" in l}
    dims = [int(x) for x in hdr["DimSize"].split()]
    a = np.fromfile(os.path.join(os.path.dirname(path), hdr["ElementDataFile"]), "<f8")
    if sorted(dims)[:2] != [1, 1]: a = a.reshape(dims[2], dims[1], dims[0]).sum(axis=(1, 2))
    return a, float(hdr["ElementSpacing"].split()[2])
ML = {"G4_AIR": "no shield", "G4_POLYETHYLENE": "PE", "G4_PLEXIGLASS": "PMMA", "G4_Al": "Al", "G4_Ti": "Ti",
      "G4_Fe": "Fe", "G4_Cu": "Cu", "G4_Pb": "Pb"}
fig, ax = plt.subplots(1, 3, figsize=(W2, 2.6))
for i, (m, lab) in enumerate(ML.items()):
    f = f"runs/s1234/output_{m}/dose_{m}_dose.mhd"
    if not os.path.exists(f): continue
    d, dz = read_mhd(f); z = (np.arange(len(d)) + 0.5) * dz / 10
    ax[0].plot(z, d / d.max(), lw=1.1, color="k" if m == "G4_AIR" else f"C{i}", label=lab)
ax[0].set_xlabel("Depth in water (cm)"); ax[0].set_ylabel("Relative dose"); ax[0].set_xlim(0, 17); ax[0].grid(alpha=0.3)
ax[0].legend(frameon=False, fontsize=6.5, ncol=2, loc="center left", bbox_to_anchor=(0.0, 0.62))
ss = list(csv.DictReader(open("validation/results/seed_stats_summary.csv")))
ss.sort(key=lambda r: float(r["dR_pstar_cm"]))
xs = np.arange(len(ss)); labs = [ML[r["material"]] for r in ss]
h1 = ax[1].errorbar(xs, [float(r["dR_diff_mm"]) for r in ss], yerr=[float(r["dR_diff_sd_mm"]) for r in ss], fmt="o", color="C0", ms=4, label="ΔR80 (left axis)")
ax[1].set_ylabel("GATE − PSTAR, ΔR80 (mm)", color="C0"); ax[1].set_ylim(0, 0.8)
a2 = ax[1].twinx(); h2, = a2.plot(xs, [float(r["e_diff_pct"]) for r in ss], "s", mfc="none", color="C3", ms=5, label="E_out (right axis)")
ax[1].legend([h1, h2], ["ΔR80 (left axis)", "E_out (right axis)"], frameon=False, fontsize=6.5, loc="lower right")
a2.set_ylabel("GATE − PSTAR, E_out (%)", color="C3"); a2.set_ylim(-0.4, 0.4); a2.axhline(0, color="C3", lw=0.6, ls=":")
ax[1].set_xticks(xs); ax[1].set_xticklabels(labs); ax[1].grid(alpha=0.3)
na = list(csv.DictReader(open("analysis/nuclear_attenuation.csv")))
pe = [100 * float(r["P_int_expected"]) for r in na]; pg = [100 * float(r["P_int_GATE"]) for r in na]
er = [100 * float(r["P_int_GATE_se"]) for r in na]
ax[2].errorbar(pe, pg, yerr=er, fmt="o", ms=4, color="C2")
OFF = {"G4_POLYETHYLENE": (6, -4), "G4_Pb": (-14, 6), "G4_Al": (6, -2), "G4_PLEXIGLASS": (-30, 2), "G4_Ti": (6, -6), "G4_Fe": (6, -8), "G4_Cu": (-18, 4)}
for r, x, y in zip(na, pe, pg): ax[2].annotate(ML[r["material"]], (x, y), textcoords="offset points", xytext=OFF.get(r["material"], (4, -9)), fontsize=7)
lim2 = [0, 8]; ax[2].plot(lim2, lim2, "k-", lw=0.7); ax[2].fill_between(lim2, [0.95 * v for v in lim2], [1.05 * v for v in lim2], color="0.85")
ax[2].set_xlim(lim2); ax[2].set_ylim(lim2)
ax[2].set_xlabel("Expected from measured σ_R (%)"); ax[2].set_ylabel("GATE: primaries interacting (%)"); ax[2].grid(alpha=0.3)
for a, t in zip(ax, ("(a)", "(b)", "(c)")): a.text(0.02, 0.92, t, transform=a.transAxes, fontweight="bold")
fig.tight_layout(); save(fig, "validation_150MeV")

# ---------------- SPE spectra ----------------
import sys; sys.path.insert(0, os.getcwd())
from spectra import SPECTRA
fig, ax = plt.subplots(figsize=(4.2, 3.0))
E = np.geomspace(10, 5000, 400)
for k, lab, ls, c in (("webber1956", "Feb 1956 — Webber", "--", "C0"), ("kold1956", "Feb 1956 — Koldobskiy", "-", "C0"),
                      ("king1972", "Aug 1972 — King", "--", "C1"), ("kold1972", "Aug 1972 — Koldobskiy", "-", "C1"),
                      ("kold1989", "Sep–Oct 1989 — Koldobskiy", "-", "C2")):
    e = E[E >= 30] if k.startswith("kold") else E
    ax.loglog(e, SPECTRA[k](e), ls, color=c, lw=1.4, label=lab)
ax.set_xlabel("Proton kinetic energy E (MeV)"); ax.set_ylabel("Integral fluence Φ(>E) (cm⁻²)")
ax.set_ylim(1e3, 3e10); ax.set_xlim(10, 5000); ax.grid(alpha=0.3, which="both"); ax.legend(frameon=False, fontsize=7)
ax.xaxis.set_major_locator(FixedLocator([10, 30, 100, 300, 1000, 3000])); ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
ax.xaxis.set_minor_formatter(NullFormatter())
save(fig, "spe_spectra")
print("Written: figures/validation_150MeV, spe_spectra (.png + .pdf)")

# ---------------- Fig. 5: mechanisms ----------------
SPC = {"kold1956": ("C0", "o", "Feb 1956 — Koldobskiy"), "kold1972": ("C1", "o", "Aug 1972 — Koldobskiy"),
       "webber1956": ("C0", "^", "Feb 1956 — Webber"), "king1972": ("C1", "^", "Aug 1972 — King")}
fig, ax = plt.subplots(1, 3, figsize=(W2, 2.8))
av = [r for r in csv.DictReader(open("analysis/analytic_vs_gate.csv")) if r["material"] != "Al"]
for k, (c, mk, lab) in SPC.items():
    sel = [r for r in av if r["spectrum"] == k]
    ax[0].plot([float(r["mean_ratio_analytic"]) for r in sel], [float(r["mean_ratio_GATE"]) for r in sel], mk,
               color=c, mfc=c if mk == "o" else "none", ms=4.5, label=lab)
x = np.array([float(r["mean_ratio_analytic"]) for r in av]); y = np.array([float(r["mean_ratio_GATE"]) for r in av])
rr = np.corrcoef(np.log(x), np.log(y))[0, 1]
ax[0].plot([0.27, 3.8], [0.27, 3.8], "k-", lw=0.7)
log_ticks(ax[0], [0.3, 0.5, 0.7, 1, 1.5, 2, 3]); ax[0].set_xscale("log"); ax[0].xaxis.set_major_locator(FixedLocator([0.3, 0.5, 0.7, 1, 1.5, 2, 3]))
ax[0].xaxis.set_major_formatter(FormatStrFormatter("%g")); ax[0].xaxis.set_minor_formatter(NullFormatter())
ax[0].set_xlim(0.27, 3.8); ax[0].set_ylim(0.27, 3.8); ax[0].set_aspect("equal")
ax[0].set_xlabel("Analytic CSDA model, ratio to Al"); ax[0].set_ylabel("Full simulation, ratio to Al (mean dose)")
ax[0].text(0.05, 0.80, f"r = {rr:.3f}\n({len(av)} cases)", transform=ax[0].transAxes, fontsize=7)
ax[0].legend(frameon=False, fontsize=6.5, loc="lower right", handletextpad=0.2, borderaxespad=0.2); ax[0].grid(alpha=0.3)
# (b) nuclear part of the ratio to Al, 20 g/cm²
de = [r for r in csv.DictReader(open("analysis/decompose_em_nuclear.csv")) if r["material"] != "Al"]
mats = ["POLYETHYLENE", "PLEXIGLASS", "Ti", "Fe", "Cu", "Pb"]; xs = np.arange(len(mats)); w = 0.2
for j, (k, q, hatch, lab) in enumerate((("kold1956", "D", "", "1956, D"), ("kold1956", "H", "///", "1956, H"),
                                        ("kold1972", "D", "", "1972, D"), ("kold1972", "H", "///", "1972, H"))):
    v = [float(next(r for r in de if r["spectrum"] == k and r["material"] == m)[f"ratio_{q}_full"]) -
         float(next(r for r in de if r["spectrum"] == k and r["material"] == m)[f"ratio_{q}_EMonly"]) for m in mats]
    ax[1].bar(xs + (j - 1.5) * w, v, w, color=SPC[k][0], alpha=0.85 if q == "D" else 0.45, hatch=hatch, edgecolor="k", lw=0.4, label=lab)
ax[1].axhline(0, color="k", lw=0.7); ax[1].set_xticks(xs); ax[1].set_xticklabels([LAB.get(m, m) for m in mats])
ax[1].set_ylabel("Nuclear part of the ratio to Al\n(full − EM-only), 20 g/cm²"); ax[1].grid(alpha=0.3, axis="y")
ax[1].set_ylim(-0.075, 0.1); ax[1].legend(frameon=False, fontsize=6.5, ncol=2, loc="upper center", columnspacing=0.8, handlelength=1.4)
# (c) contribution of each particle to H − D, Al 20 g/cm², Koldobskiy 1972
dq = np.genfromtxt("spe_runs/materials_kold1972/Al_20/dose_eq_depth.csv", delimiter=",", names=True)
z = dq["depth_g_cm2"]; sp = [("proton", "p"), ("alpha", "α"), ("rest", "Z ≥ 3"), ("deuteron", "d"), ("He3", "³He"), ("triton", "t")]
win = np.ones(31) / 31
parts = np.array([np.convolve(dq[f"dH_{k}_cSv"], win, mode="same") for k, _ in sp])
tot = parts.sum(axis=0); ok = (z > 1.5) & (z < 28.5)
ax[2].stackplot(z[ok], 100 * parts[:, ok] / tot[ok], labels=[l for _, l in sp],
                colors=["#d62728", "#2ca02c", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"], alpha=0.85)
ax[2].set_xlim(1.5, 28.5); ax[2].set_ylim(0, 100); ax[2].set_xlabel("Depth in water (g/cm²)")
ax[2].set_ylabel("Share of H − D (%)"); ax[2].legend(frameon=False, fontsize=6.5, loc="upper right", ncol=2)
for a, t in zip(ax, ("(a)", "(b)", "(c)")): a.text(-0.02, 1.03, t, transform=a.transAxes, fontweight="bold")
fig.tight_layout(); save(fig, "mechanisms")
print("Written: figures/mechanisms (.png + .pdf)")
