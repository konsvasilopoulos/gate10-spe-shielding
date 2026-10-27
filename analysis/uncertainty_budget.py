"""uncertainty_budget.py - Quantification of the uncertainties for the paper (no new runs).
Writes analysis/uncertainty_budget.csv (and the console output can be saved with tee).
  A. Statistics: relative sigma of the mean dose 0-30 and of the skin dose, and of the ratios to Al (all material runs)
  B. LET binning for H: merging of neighbouring LET bins (benchmark)
  C. Spectrum sampling: 100 logarithmic bins per group, flat inside each bin (analytic CSDA)
  D. Spectrum limits: fluence < 30 MeV (Koldobskiy) and > 5 GeV / > 2 GeV (analytic CSDA)
  E. Shield stopping power (PSTAR) +-2 %: effect on the ratios to Al (analytic CSDA)
  F. Digitisation of the NASA curves: depth shift x slope of the curve
"""
import os, sys, glob, json, csv
import numpy as np
sys.path.insert(0, os.getcwd())
from spectra import SPECTRA
trapz = getattr(np, "trapezoid", None) or np.trapz
out = []          # (category, quantity, range, comment)

def read_mhd(path):
    hdr = {l.split("=")[0].strip(): l.split("=", 1)[1].strip() for l in open(path) if "=" in l}
    dims = [int(x) for x in hdr["DimSize"].split()]
    a = np.fromfile(os.path.join(os.path.dirname(path), hdr["ElementDataFile"]), "<f8")
    if sorted(dims)[:2] != [1, 1]: a = a.reshape(dims[2], dims[1], dims[0]).sum(axis=(1, 2))
    return a

# ---------------- A. Statistics ----------------
print("A. Statistical uncertainty (relative sigma, %): upper bound = fully correlated voxels (sum sigma / sum D)")
BASES = {"kold1956": "spe_runs/materials_kold1956", "kold1972": "spe_runs/materials_kold1972",
         "webber1956": "spe_runs/materials", "king1972": "spe_runs/materials_king1972"}
stat = {}
for sp, b in BASES.items():
    for f in glob.glob(f"{b}/*_*/spe_dose_depth.csv"):
        name = os.path.basename(os.path.dirname(f)); m, t = name.rsplit("_", 1)
        d = np.genfromtxt(f, delimiter=",", names=True)
        D, s = d["dose_cGy"], d["sigma_cGy"]
        stat[(sp, m, int(t))] = dict(mean_corr=s.sum() / D.sum(), mean_unc=np.sqrt((s**2).sum()) / D.sum(),
                                     skin=s[:3].sum() / D[:3].sum(), deep=s[-10:].sum() / D[-10:].sum())
for k in ("mean_corr", "mean_unc", "skin", "deep"):
    v = np.array([x[k] for x in stat.values()]) * 100
    print(f"   {k:<10}: median {np.median(v):.2f} %, max {v.max():.2f} %  ({max(stat, key=lambda q: stat[q][k])})")
rat = [np.hypot(stat[(sp, m, t)]["mean_corr"], stat[(sp, "Al", t)]["mean_corr"]) for (sp, m, t) in stat if m != "Al"]
print(f"   ratio to Al (mean dose): median {100*np.median(rat):.2f} %, max {100*max(rat):.2f} %")
v = np.array([x["mean_corr"] for x in stat.values()]) * 100
out.append(["Statistics (mean dose 0-30)", "D", f"{np.median(v):.1f} % (max {v.max():.1f} %)", "upper bound, fully correlated voxels"])
v = np.array([x["skin"] for x in stat.values()]) * 100
out.append(["Statistics (skin 0-0.3)", "D", f"{np.median(v):.1f} % (max {v.max():.1f} %)", ""])
out.append(["Statistics (ratio to Al)", "D/D_Al", f"{100*np.median(rat):.1f} % (max {100*max(rat):.1f} %)", ""])

# ---------------- B. LET binning ----------------
print("\nB. LET binning (benchmark Webber 1956, Al 20): H with neighbouring energy bins merged per particle")
def Q(L):
    L = np.asarray(L, float); q = np.ones_like(L)
    m = (L >= 10) & (L <= 100); q[m] = 0.32 * L[m] - 2.2
    m = L > 100; q[m] = 300 / np.sqrt(L[m]); return q
B = "spe_runs/webber1956_Al20_wide_H"
def H_profile(merge):
    Dt = Ht = 0
    for d in sorted(glob.glob(f"{B}/g*")):
        s = json.load(open(f"{d}/setup.json")); mat = s["material"]
        f = s["group_fluence_cm2"] * s["source_area_cm2"] / s["n"] * 100
        dose = read_mhd(f"{d}/dose_{mat}_dose.mhd"); edep = read_mhd(f"{d}/dose_{mat}_edep.mhd")
        by = {}
        for npth in glob.glob(f"{d}/let_*numerator*.mhd"):
            b = os.path.basename(npth); spc = b.split("_")[1]
            lo = 0.0 if spc == "rest" else float(b.split("_")[2].split("-")[0])
            by.setdefault(spc, []).append((lo, read_mhd(npth), read_mhd(npth.replace("numerator", "denominator"))))
        extra = np.zeros_like(dose)
        for spc, lst in by.items():
            lst.sort(key=lambda x: x[0])
            groups = [lst] if merge == "all" else [lst[i:i + merge] for i in range(0, len(lst), merge)]
            for g in groups:
                num = sum(x[1] for x in g); den = sum(x[2] for x in g)
                with np.errstate(divide="ignore", invalid="ignore"):
                    L = np.where(den > 0, num / den, 0)
                extra += den * (Q(L) - 1)
        with np.errstate(divide="ignore", invalid="ignore"):
            h = dose * (1 + np.where(edep > 0, extra / edep, 0))
        Dt = Dt + dose * f; Ht = Ht + h * f
    return Dt, Ht
z = None; res = {}
for mg in (1, 2, 3, "all"):
    Dd, Hh = H_profile(mg); z = (np.arange(len(Dd)) + 0.5) / 10; res[mg] = Hh
ref = res[1]
for zc in (0, 10, 20, 29.5):
    w = (np.abs(z - zc) <= 1.0) if zc > 0 else (z < 0.3)
    print(f"   depth {zc:>4}: " + " | ".join(f"merge {mg}: {res[mg][w].mean()/ref[w].mean():.3f}" for mg in (2, 3, "all")))
dev2 = max(abs(res[2][(np.abs(z-zc)<=1)].mean()/ref[(np.abs(z-zc)<=1)].mean()-1) for zc in (5, 10, 20, 29))
devA = max(abs(res["all"][(np.abs(z-zc)<=1)].mean()/ref[(np.abs(z-zc)<=1)].mean()-1) for zc in (5, 10, 20, 29))
out.append(["LET binning / Q(Lbar)", "H", f"{100*dev2:.1f} % (bins 2x coarser)", f"one bin per particle: {100*devA:.0f} %"])

# ---------------- C, D, E: analytic CSDA ----------------
def pst(m):
    a = np.array([[float(x) for x in l.split()] for l in open(f"pstar/G4_{m}.txt") if len(l.split()) == 7]); return a[:, 0], a[:, 4]
Ew, Rw = pst("WATER")
Rf = lambda E: np.exp(np.interp(np.log(np.maximum(E, Ew[0])), np.log(Ew), np.log(Rw)))
Ef = lambda r: np.exp(np.interp(np.log(np.maximum(r, Rw[0])), np.log(Rw), np.log(Ew)))
def E_at(E, z): return np.where((E > 0) & (Rf(E) > z), Ef(Rf(E) - z), 0)
def dose(m, t, phi_E, Eg, z1=0, z2=30, sscale=1.0):
    Em, Rm = pst(m); Rm = Rm / sscale
    r = np.exp(np.interp(np.log(Eg), np.log(Em), np.log(Rm))) - t
    eo = np.where(r > Rm.min(), np.exp(np.interp(np.log(np.maximum(r, Rm.min())), np.log(Rm), np.log(Em))), 0)
    return trapz(phi_E * (E_at(eo, z1) - E_at(eo, z2)), Eg) / (z2 - z1)
def pdf(sp, Eg): return -np.gradient(SPECTRA[sp](Eg), Eg)
GR = {"kold1956": [30, 60, 120, 170, 250, 400, 1000, 5000], "kold1972": [30, 60, 120, 170, 250, 400, 1000, 5000],
      "king1972": [10, 60, 120, 170, 250, 400, 1000], "webber1956": [10, 60, 120, 170, 250, 400, 2000]}
print("\nC. Sampling: 100 logarithmic bins per group (flat) vs continuous spectrum - ratio of the mean dose 0-30")
cmax = 0
for sp, e in GR.items():
    line = f"   {sp:<11}"
    for m, t in (("Al", 20), ("POLYETHYLENE", 20), ("Pb", 5)):
        Eg = np.geomspace(e[0], e[-1], 40001); exact = dose(m, t, pdf(sp, Eg), Eg)
        approx = np.zeros_like(Eg)
        for lo, hi in zip(e[:-1], e[1:]):
            ed = np.geomspace(lo, hi, 101); w = SPECTRA[sp](ed[:-1]) - SPECTRA[sp](ed[1:])
            idx = np.clip(np.searchsorted(ed, Eg, side="right") - 1, 0, 99); sel = (Eg >= lo) & (Eg < hi)
            approx[sel] = (w / np.diff(ed))[idx[sel]]
        r = dose(m, t, approx, Eg) / exact; cmax = max(cmax, abs(r - 1)); line += f"  {m[:4]} {t}: {r:.4f}"
    print(line)
out.append(["Spectrum sampling (100 bins/group)", "D", f"< {100*cmax:.2f} %", "analytic CSDA"])

print("\nD. Spectrum limits (analytic CSDA, mean dose 0-30, behind Al): missing fraction")
for sp, hi in (("kold1956", 5000), ("kold1972", 5000), ("kold1989", 5000), ("webber1956", 2000), ("king1972", 1000)):
    Eg = np.geomspace(30, 10000, 60001); full = {t: dose("Al", t, pdf(sp, Eg), Eg) for t in (5, 20)}
    Ec = np.geomspace(30, hi, 60001); cut = {t: dose("Al", t, pdf(sp, Ec), Ec) for t in (5, 20)}
    above = SPECTRA[sp](10000)
    print(f"   {sp:<11} > {hi} MeV: {100*(1-cut[5]/full[5]):.2f} % (5 g/cm²), {100*(1-cut[20]/full[20]):.2f} % (20 g/cm²)"
          f"  | Phi(>10 GeV) = {above:.2e} cm^-2 (outside PSTAR)")
    if sp.startswith("kold"):
        out.append([f"Upper spectrum limit ({sp}, 5 GeV)", "D", f"{100*(1-cut[20]/full[20]):.2f} % (20 g/cm2)", "analytic, up to 10 GeV"])
out.append(["Lower spectrum limit (Koldobskiy, 30 MeV)", "D", "0 % for primaries", "E_th >= 50 MeV for all shields"])

print("\nE. Shield stopping power +-2 % (ICRU 49): change of the mean-dose ratio to Al (only the material changes)")
emax = 0
for sp in ("kold1956", "kold1972", "king1972"):
    Eg = np.geomspace(30 if sp.startswith("kold") else 10, 5000 if sp.startswith("kold") else 1000, 40001); ph = pdf(sp, Eg)
    line = f"   {sp:<10}"
    for m in ("POLYETHYLENE", "Pb"):
        for t in (5, 20):
            r0 = dose(m, t, ph, Eg) / dose("Al", t, ph, Eg); r1 = dose(m, t, ph, Eg, sscale=1.02) / dose("Al", t, ph, Eg)
            d = r1 / r0 - 1; emax = max(emax, abs(d)); line += f"  {m[:4]} {t}: {100*d:+.2f} %"
    print(line)
out.append(["Shield stopping power PSTAR/Geant4 (+-2 %)", "D/D_Al", f"<= {100*emax:.1f} %", "analytic CSDA"])

# ---------------- F. Digitisation ----------------
print("\nF. NASA digitisation: depth shift 0.03-0.19 g/cm2 x slope d lnD/dz")
d = np.genfromtxt(f"{B}/spe_dose_depth.csv", delimiter=",", names=True)
sl0 = abs(np.log(d["dose_cGy"][50] / d["dose_cGy"][0]) / (d["depth_g_cm2"][50] - d["depth_g_cm2"][0]))
slm = abs(np.log(d["dose_cGy"][-1] / d["dose_cGy"][100]) / (d["depth_g_cm2"][-1] - d["depth_g_cm2"][100]))
print(f"   slope: {100*sl0:.1f} %/(g/cm2) at the entrance, {100*slm:.1f} %/(g/cm2) at depth -> error <= {100*sl0*0.19:.1f} % at 0, <= {100*slm*0.04:.1f} % elsewhere")
out.append(["Digitisation of the NASA curves", "reference", f"<= {100*sl0*0.19:.1f} % (0), <= {100*slm*0.04:.1f} % (elsewhere)", "+ reading on a log axis ~1-2 %"])

with open("analysis/uncertainty_budget.csv", "w", newline="") as fo:
    w = csv.writer(fo); w.writerow(["source", "quantity", "estimate", "comment"]); w.writerows(out)
print("\nWritten: analysis/uncertainty_budget.csv")
