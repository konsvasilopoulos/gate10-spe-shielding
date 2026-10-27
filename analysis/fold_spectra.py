"""fold_spectra.py - Dose for ANY SPE spectrum without new simulations (coarse energy groups).
Each material run stores the dose per energy group g (spe_dose_depth.csv). Response per unit group fluence:
    R_g = D_g / Phi_g(source)      (with the shape of the source spectrum inside the group)
For a target spectrum f:  D(f) = sum_g Phi_g(f) * R_g * c_g,   c_g = A_g(f) / A_g(source)
where A_g = analytic (CSDA) dose per proton of the group with the shape of each spectrum (in-group shape correction).
Sources: the Koldobskiy 1956 and 1972 runs (same groups 30-5000 MeV). Quantities: skin D (0-0.3), D at 5 g/cm2 (BFO), mean D 0-30.
Usage: python3 analysis/fold_spectra.py [crossval|predict1989|check1989|scan]"""
import os, sys, csv, json
import numpy as np
sys.path.insert(0, os.getcwd()); sys.path.insert(0, "analysis")
from spectra import SPECTRA
from gle_table import T2, F as FGLE, DATES
trapz = getattr(np, "trapezoid", None) or np.trapz
G = [30, 60, 120, 170, 250, 400, 1000, 5000]
WIN = {"skin": (0, 0.3), "BFO": (4, 6), "mean": (0, 30)}
MATS = ["POLYETHYLENE", "PLEXIGLASS", "Al", "Ti", "Fe", "Cu", "Pb"]
SRC = {"kold1956": "spe_runs/materials_kold1956", "kold1972": "spe_runs/materials_kold1972"}

def pst(m):
    a = np.array([[float(x) for x in l.split()] for l in open(f"pstar/G4_{m}.txt") if len(l.split()) == 7]); return a[:, 0], a[:, 4]
Ew, Rw = pst("WATER")
Rf = lambda E: np.exp(np.interp(np.log(np.maximum(E, Ew[0])), np.log(Ew), np.log(Rw)))
Ef = lambda r: np.exp(np.interp(np.log(np.maximum(r, Rw[0])), np.log(Rw), np.log(Ew)))
def E_at(E, z): return np.where((E > 0) & (Rf(E) > z), Ef(Rf(E) - z), 0.0)
PST = {m: pst(m) for m in MATS}
def A_group(m, t, Ffun, lo, hi, win):
    """analytic dose (MeV/g) per proton of the group, with the shape given by the integral spectrum Ffun."""
    Em, Rm = PST[m]; Eg = np.geomspace(lo, hi, 3001); phi = -np.gradient(Ffun(Eg), Eg)
    phi = np.clip(phi, 0, None); n = trapz(phi, Eg)
    if n <= 0: return 0.0
    r = np.exp(np.interp(np.log(Eg), np.log(Em), np.log(Rm))) - t
    eo = np.where(r > Rm.min(), np.exp(np.interp(np.log(np.maximum(r, Rm.min())), np.log(Rm), np.log(Em))), 0)
    z1, z2 = win
    return trapz(phi * (E_at(eo, z1) - E_at(eo, z2)), Eg) / (z2 - z1) / n

def load_resp(src):
    """R[(m,t)][win][g] = dose (cGy) per unit group fluence (cm^-2) with the shape of the source spectrum."""
    out = {}
    for m in MATS:
        for t in (5, 10, 20):
            f = f"{SRC[src]}/{m}_{t}/spe_dose_depth.csv"
            if not os.path.exists(f): continue
            d = np.genfromtxt(f, delimiter=",", names=True); z = d["depth_g_cm2"]
            out[(m, t)] = {}
            for w, (lo, hi) in WIN.items():
                sel = (z >= lo) & (z <= hi) if w != "skin" else (z < 0.3)
                rg = []
                for a, b in zip(G[:-1], G[1:]):
                    col = f"{a}{b}MeV_cGy"
                    if col not in d.dtype.names: col = [c for c in d.dtype.names if c.startswith(f"{a}") and c.endswith("MeV_cGy")][0]
                    phi = SPECTRA[src](a) - SPECTRA[src](b)
                    rg.append(d[col][sel].mean() / phi)
                out[(m, t)][w] = np.array(rg)
    return out
RESP = {s: load_resp(s) for s in SRC}

def predict(Ftarget, src, m, t, w, correct=True):
    R = RESP[src][(m, t)][w]; tot = 0.0
    for i, (a, b) in enumerate(zip(G[:-1], G[1:])):
        phi = Ftarget(a) - Ftarget(b); c = 1.0
        if correct:
            As = A_group(m, t, SPECTRA[src], a, b, WIN[w]); At = A_group(m, t, Ftarget, a, b, WIN[w])
            if As > 1e-6 * max(At, 1e-30) and As > 0: c = At / As
        tot += phi * R[i] * c
    return tot
def actual(spec, m, t, w):
    d = np.genfromtxt(f"spe_runs/materials_{spec}/{m}_{t}/spe_dose_depth.csv", delimiter=",", names=True); z = d["depth_g_cm2"]
    lo, hi = WIN[w]; sel = (z >= lo) & (z <= hi) if w != "skin" else (z < 0.3)
    return d["dose_cGy"][sel].mean()

mode = sys.argv[1] if len(sys.argv) > 1 else "crossval"
if mode == "crossval":
    print("Cross-validation: prediction of one spectrum from the runs of the other (prediction / actual GATE)")
    rows = []
    for src, tgt in (("kold1956", "kold1972"), ("kold1972", "kold1956")):
        for corr in (False, True):
            dev = []
            for m in MATS:
                for t in (5, 10, 20):
                    for w in WIN:
                        p = predict(SPECTRA[tgt], src, m, t, w, corr); a = actual(tgt, m, t, w)
                        dev.append(p / a); rows.append([src, tgt, corr, m, t, w, p, a, p / a])
            dev = np.array(dev)
            print(f"  {src} -> {tgt}, shape correction {'yes' if corr else 'no'}: median {np.median(dev):.3f}, "
                  f"range {dev.min():.3f}-{dev.max():.3f}, |deviation| <= 5 %: {100*np.mean(abs(dev-1) <= 0.05):.0f} %")
    # ratios to Al
    print("\n  Ratios to Al (mean D), prediction / actual, with correction:")
    for src, tgt in (("kold1956", "kold1972"), ("kold1972", "kold1956")):
        line = f"  {src}→{tgt}:"
        for m in ("POLYETHYLENE", "Pb"):
            for t in (5, 20):
                pr = predict(SPECTRA[tgt], src, m, t, "mean") / predict(SPECTRA[tgt], src, "Al", t, "mean")
                ac = actual(tgt, m, t, "mean") / actual(tgt, "Al", t, "mean")
                line += f"  {m[:4]} {t}: {pr:.3f}/{ac:.3f}"
        print(line)
    with open("analysis/fold_crossval.csv", "w", newline="") as fo:
        wr = csv.writer(fo); wr.writerow(["source", "target", "shape_correction", "material", "areal", "window", "pred_cGy", "actual_cGy", "ratio"]); wr.writerows(rows)
elif mode == "predict1989":
    import datetime
    rows = []; print("BLIND prediction for Koldobskiy 1989 (GLE 42-45), before the run finished - " + datetime.datetime.now().isoformat(timespec="minutes"))
    for m in MATS:
        for t in (5, 10, 20):
            for w in WIN:
                p56 = predict(SPECTRA["kold1989"], "kold1956", m, t, w); p72 = predict(SPECTRA["kold1989"], "kold1972", m, t, w)
                rows.append([m, t, w, p56, p72, 0.5 * (p56 + p72)])
    with open("analysis/fold_predict_1989.csv", "w", newline="") as fo:
        wr = csv.writer(fo); wr.writerow(["material", "areal", "window", "pred_from_1956_cGy", "pred_from_1972_cGy", "pred_mean_cGy"]); wr.writerows(rows)
    P = {(r[0], r[1], r[2]): r[5] for r in rows}
    for t in (5, 10, 20):
        print(f"  {t:>2} g/cm2: Al mean D {P[('Al',t,'mean')]:.2f} cGy | PE/Al {P[('POLYETHYLENE',t,'mean')]/P[('Al',t,'mean')]:.3f} | "
              f"Pb/Al {P[('Pb',t,'mean')]/P[('Al',t,'mean')]:.3f} · BFO Al {1.5*10*P[('Al',t,'BFO')]:.0f} mGy-Eq")
    print("Written: analysis/fold_predict_1989.csv (do NOT modify until the direct run is available)")
elif mode == "check1989":
    P = {(r["material"], int(r["areal"]), r["window"]): float(r["pred_mean_cGy"]) for r in csv.DictReader(open("analysis/fold_predict_1989.csv"))}
    dev = []; print("Prediction / GATE for 1989:")
    for (m, t, w), p in P.items():
        a = actual("kold1989", m, t, w); dev.append(p / a)
    dev = np.array(dev); print(f"  all: median {np.median(dev):.3f}, range {dev.min():.3f}-{dev.max():.3f}")
    for t in (5, 10, 20):
        for m in ("POLYETHYLENE", "Pb"):
            pr = P[(m, t, "mean")] / P[("Al", t, "mean")]; ac = actual("kold1989", m, t, "mean") / actual("kold1989", "Al", t, "mean")
            print(f"  {m:<13}{t:>3} g/cm2: ratio to Al, prediction {pr:.3f} | GATE {ac:.3f}")
elif mode == "scan":
    rows = []
    for gle in sorted(T2):
        Ff = lambda E, g=gle: FGLE(E, g)
        f30, f200 = float(Ff(30)), float(Ff(200))
        if f30 <= 0 or f200 <= 0: continue
        hard = f200 / f30
        r = {}
        for m in ("POLYETHYLENE", "Al", "Pb"):
            for t in (5, 10, 20):
                r[(m, t)] = 0.5 * (predict(Ff, "kold1956", m, t, "mean") + predict(Ff, "kold1972", m, t, "mean"))
                r[(m, t, "BFO")] = 0.5 * (predict(Ff, "kold1956", m, t, "BFO") + predict(Ff, "kold1972", m, t, "BFO"))
        rows.append([gle, DATES.get(gle, ""), f30, f200, hard] + [r[(m, t)] / r[("Al", t)] for m in ("POLYETHYLENE", "Pb") for t in (5, 10, 20)]
                    + [r[("Al", t)] for t in (5, 10, 20)] + [1.5 * 10 * r[(m, 20, "BFO")] for m in ("POLYETHYLENE", "Al")])
    hdr = ["GLE", "date", "F30_cm2", "F200_cm2", "hardness_F200_over_F30", "PE_Al_5", "PE_Al_10", "PE_Al_20", "Pb_Al_5", "Pb_Al_10", "Pb_Al_20",
           "Al_meanD_5_cGy", "Al_meanD_10_cGy", "Al_meanD_20_cGy", "BFO_PE20_mGyEq", "BFO_Al20_mGyEq"]
    with open("analysis/fold_gle_scan.csv", "w", newline="") as fo:
        wr = csv.writer(fo); wr.writerow(hdr); wr.writerows(rows)
    h = np.array([r[4] for r in rows]); pe20 = np.array([r[7] for r in rows]); pb20 = np.array([r[10] for r in rows])
    print(f"{len(rows)} GLEs. Correlation of log(hardness) with PE/Al(20): r = {np.corrcoef(np.log(h), pe20)[0,1]:.3f}, with Pb/Al(20): r = {np.corrcoef(np.log(h), pb20)[0,1]:.3f}")
    for r in sorted(rows, key=lambda x: -x[3])[:12]:
        print(f"  GLE {r[0]:>2} {r[1]:<11} Phi>200 {r[3]:.2e}  hardness {r[4]:.3f}  PE/Al(20) {r[7]:.2f}  Pb/Al(20) {r[10]:.2f}  BFO Al20 {r[15]:.0f} mGy-Eq")
