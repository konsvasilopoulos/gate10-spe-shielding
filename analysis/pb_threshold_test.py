"""pb_threshold_test.py - Transmitted fraction and exit energy near the threshold (GATE vs PSTAR/CSDA)."""
import glob, json, os
import numpy as np
import uproot
def pst(m):
    a = np.array([[float(x) for x in l.split()] for l in open(f"pstar/G4_{m}.txt") if len(l.split()) == 7])
    return a[:, 0], a[:, 4]
print(f"{'Mat':<6}{'E (MeV)':>8}{'transm. GATE':>15}{'CSDA':>7}{'E_out GATE (median)':>22}{'E_out CSDA':>12}{'<1/cos>':>10}")
for d in sorted(glob.glob("spe_runs/pb_threshold/*_E*"), key=lambda p: (p.split("/")[-1].split("_")[0], float(p.split("E")[-1]))):
    if not os.path.exists(f"{d}/setup.json"): continue
    s = json.load(open(f"{d}/setup.json")); mat = s["material"].replace("G4_", "")
    E0 = 0.5 * (s["emin_MeV"] + s["emax_MeV"])
    with uproot.open(glob.glob(f"{d}/phsp_*.root")[0]) as f:
        t = f[f.keys()[0]]
        p = t.arrays(["ParticleName", "KineticEnergy", "Direction_Z"], library="np")
    sel = (p["ParticleName"].astype(str) == "proton") & (p["Direction_Z"] > 0)
    ntar = s["n"] * (30 / 80) ** 2                    # source protons aimed at the 30x30 cm2 plane
    frac = sel.sum() / ntar
    E, R = pst(mat); r = np.exp(np.interp(np.log(E0), np.log(E), np.log(R))) - 5.0
    eo = np.exp(np.interp(np.log(r), np.log(R), np.log(E))) if r > R.min() else 0.0
    print(f"{mat:<6}{E0:>8.1f}{frac:>15.3f}{(1.0 if r > 0 else 0.0):>7.0f}"
          f"{(np.median(p['KineticEnergy'][sel]) if sel.any() else 0):>22.2f}{eo:>12.2f}"
          f"{(np.mean(1/np.clip(p['Direction_Z'][sel],0.02,1)) if sel.any() else 0):>10.3f}")
