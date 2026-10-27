"""pb_threshold_dose.py - Skin dose (0-0.3 g/cm2) per source proton near the transmission threshold:
GATE (DoseActor) vs CSDA from the actual phase-space protons (real E_out and angle)
and vs straight-line CSDA (analytic model, as in thresholds_pstar.py)."""
import glob, json, os
import numpy as np, uproot
MEV_PER_G_TO_GY = 1.602176634e-10
def pst(m):
    a = np.array([[float(x) for x in l.split()] for l in open(f"pstar/G4_{m}.txt") if len(l.split()) == 7]); return a[:, 0], a[:, 4]
Ew, Rw = pst("WATER")
R = lambda E: np.exp(np.interp(np.log(np.maximum(E, Ew[0])), np.log(Ew), np.log(Rw)))
Ei = lambda r: np.exp(np.interp(np.log(np.maximum(r, Rw[0])), np.log(Rw), np.log(Ew)))
dep = lambda E, t: E - np.where(R(E) > t, Ei(R(E) - t), 0.0)
SK = 0.3
print(f"{'Mat':<5}{'E':>6}{'GATE':>10}{'CSDA phsp':>11}{'CSDA str.':>11}{'GATE/phsp':>11}{'GATE/str.':>11}   (MeV/g per proton/cm2)")
for d in sorted(glob.glob("spe_runs/pb_threshold/*_E*"), key=lambda p: (p.split("/")[-1][:2], float(p.split("E")[-1]))):
    if not os.path.exists(f"{d}/setup.json"): continue
    s = json.load(open(f"{d}/setup.json")); mat = s["material"].replace("G4_", ""); E0 = 0.5*(s["emin_MeV"]+s["emax_MeV"])
    phi = s["n"] / s["source_area_cm2"]                                    # protons/cm2 in the run
    hdr = {l.split("=")[0].strip(): l.split("=",1)[1].strip() for l in open(f"{d}/dose_{s['material']}_dose.mhd") if "=" in l}
    D = np.fromfile(f"{d}/{hdr['ElementDataFile']}", "<f8")[:3].mean() / MEV_PER_G_TO_GY / phi
    with uproot.open(glob.glob(f"{d}/phsp_*.root")[0]) as f:
        p = f[f.keys()[0]].arrays(["ParticleName", "KineticEnergy", "Direction_Z"], library="np")
    sel = (p["ParticleName"].astype(str) == "proton") & (p["Direction_Z"] > 0)
    E, c = p["KineticEnergy"][sel], np.clip(p["Direction_Z"][sel], 0.02, 1)
    ph = (dep(E, SK / c) / c).sum() / 900.0 / phi / SK                       # per proton/cm2
    Em, Rm = pst(mat); r = np.exp(np.interp(np.log(E0), np.log(Em), np.log(Rm))) - 5.0
    eo = np.exp(np.interp(np.log(r), np.log(Rm), np.log(Em))) if r > Rm.min() else 0.0
    an = dep(np.array(eo), SK) / SK if eo > 0 else 0.0
    q = lambda a, b: f"{a/b:>11.2f}" if b > 0 else f"{'—':>11}"
    print(f"{mat:<5}{E0:>6.0f}{D:>10.2f}{ph:>11.2f}{float(an):>11.2f}{q(D,ph)}{q(D,an)}")
