"""response_matrix.py - Response matrix from set 4 (spe_runs/response_matrix) and its validation.
  build    : R_D(E, z), R_H(E, z) = dose / dose equivalent per unit fluence (Gy cm2, Sv cm2) for each material, thickness,
             narrow energy window (25 windows 30-5000 MeV) and depth (1 mm) -> analysis/response_matrix.npz
  validate : folding with the 5 directly simulated spectra (Koldobskiy 1956/1972/1989, Webber, King) and comparison
             with GATE for D and H at the skin (0-0.3), BFO (4-6) and mean 0-30 g/cm2.
  fold <GLE number or spectrum name> : D and H for any spectrum.
  scan [i0 i1] : all GLEs of Koldobskiy et al. (2021) -> analysis/rm_gle_scan.csv
  build_thick : adds the 30 and 40 g/cm2 shields of run_set5_thick.sh
NOTE: build needs the raw runs of run_set4_response.sh (~2-3 h); validate/fold/scan only need response_matrix.npz.
Folding:  D(z) = sum_b Phi_b * R_b(z) * c_b,  c_b = in-window shape correction (analytic CSDA; small, windows x1.23)."""
import os, sys, glob, json, csv
import numpy as np
sys.path.insert(0, os.getcwd()); sys.path.insert(0, "analysis")
from spectra import SPECTRA
trapz = getattr(np, "trapezoid", None) or np.trapz
BASE = "spe_runs/response_matrix"; OUT = "analysis/response_matrix.npz"
MATS = ["POLYETHYLENE", "PLEXIGLASS", "Al", "Ti", "Fe", "Cu", "Pb"]; AREALS = [5, 10, 20]
WIN = {"skin": (0, 0.3), "BFO": (4, 6), "mean": (0, 30)}

def read_mhd(path):
    hdr = {l.split("=")[0].strip(): l.split("=", 1)[1].strip() for l in open(path) if "=" in l}
    dims = [int(x) for x in hdr["DimSize"].split()]
    a = np.fromfile(os.path.join(os.path.dirname(path), hdr["ElementDataFile"]), "<f8")
    if sorted(dims)[:2] != [1, 1]: a = a.reshape(dims[2], dims[1], dims[0]).sum(axis=(1, 2))
    return a
def Q(L):
    L = np.asarray(L, float); q = np.ones_like(L)
    m = (L >= 10) & (L <= 100); q[m] = 0.32 * L[m] - 2.2
    m = L > 100; q[m] = 300 / np.sqrt(L[m]); return q
def run_DH(d):
    s = json.load(open(f"{d}/setup.json")); mat = s["material"]
    dose = read_mhd(f"{d}/dose_{mat}_dose.mhd"); edep = read_mhd(f"{d}/dose_{mat}_edep.mhd")
    extra = np.zeros_like(dose)
    for npth in glob.glob(f"{d}/let_*numerator*.mhd"):
        num = read_mhd(npth); den = read_mhd(npth.replace("numerator", "denominator"))
        with np.errstate(divide="ignore", invalid="ignore"): L = np.where(den > 0, num / den, 0)
        extra += den * (Q(L) - 1)
    with np.errstate(divide="ignore", invalid="ignore"): H = dose * (1 + np.where(edep > 0, extra / edep, 0))
    per_flu = s["source_area_cm2"] / s["n"]            # Gy -> Gy cm2 per unit fluence
    return s["emin_MeV"], s["emax_MeV"], dose * per_flu, H * per_flu

def build_thick(base="spe_runs/response_matrix_thick", areals=(30, 40)):
    """Adds the 30 and 40 g/cm2 shields (run_set5_thick.sh) to response_matrix.npz, keeping 5/10/20."""
    M = dict(np.load(OUT)); n0 = len(M)
    for m in MATS:
        for t in areals:
            dirs = sorted(glob.glob(f"{base}/{m}_{t}/g*"), key=lambda p: float(p.split("/g")[-1].split("-")[0]))
            dirs = [d for d in dirs if os.path.exists(f"{d}/setup.json")]
            if len(dirs) != len(M["edges"]) - 1: print(f"  [!] {m} {t}: {len(dirs)} windows"); continue
            rows = [run_DH(d) for d in dirs]
            e = np.array([r[0] for r in rows] + [rows[-1][1]])
            assert np.allclose(e, M["edges"], rtol=1e-3), "different energy windows"
            M[f"D_{m}_{t}"] = np.array([r[2] for r in rows]); M[f"H_{m}_{t}"] = np.array([r[3] for r in rows])
            print(f"  {m:<13}{t:>3} g/cm2: {len(rows)} windows")
    np.savez(OUT, **M); print(f"Written: {OUT} ({n0} -> {len(M)} arrays)")

def build():
    R = {}; edges = None
    for m in MATS:
        for t in AREALS:
            dirs = sorted(glob.glob(f"{BASE}/{m}_{t}/g*"), key=lambda p: float(p.split("/g")[-1].split("-")[0]))
            dirs = [d for d in dirs if os.path.exists(f"{d}/setup.json")]
            if not dirs: continue
            rows = [run_DH(d) for d in dirs]
            e = np.array([r[0] for r in rows] + [rows[-1][1]])
            if edges is None: edges = e
            R[f"D_{m}_{t}"] = np.array([r[2] for r in rows]); R[f"H_{m}_{t}"] = np.array([r[3] for r in rows])
            print(f"  {m:<13}{t:>3} g/cm2: {len(rows)} windows")
    nz = next(iter(R.values())).shape[1]
    np.savez(OUT, edges=edges, z=(np.arange(nz) + 0.5) / 10, **R)
    print(f"Written: {OUT}")

# ---- analytic in-window shape correction ----
def pst(m):
    a = np.array([[float(x) for x in l.split()] for l in open(f"pstar/G4_{m}.txt") if len(l.split()) == 7]); return a[:, 0], a[:, 4]
Ew, Rw = pst("WATER")
Rf = lambda E: np.exp(np.interp(np.log(np.maximum(E, Ew[0])), np.log(Ew), np.log(Rw)))
Ef = lambda r: np.exp(np.interp(np.log(np.maximum(r, Rw[0])), np.log(Rw), np.log(Ew)))
def E_at(E, z): return np.where((E > 0) & (Rf(E) > z), Ef(Rf(E) - z), 0.0)
_P = {}
def A_bin(m, t, Ffun, lo, hi, win):
    if m not in _P: _P[m] = pst(m)
    Em, Rm = _P[m]; Eg = np.geomspace(lo, hi, 801); phi = np.clip(-np.gradient(Ffun(Eg), Eg), 0, None); n = trapz(phi, Eg)
    if n <= 0: return 0.0
    r = np.exp(np.interp(np.log(Eg), np.log(Em), np.log(Rm))) - t
    eo = np.where(r > Rm.min(), np.exp(np.interp(np.log(np.maximum(r, Rm.min())), np.log(Rm), np.log(Em))), 0)
    return trapz(phi * (E_at(eo, win[0]) - E_at(eo, win[1])), Eg) / (win[1] - win[0]) / n

_M = {}
def fold(Ffun, m, t, q="D", win="mean", correct=True):
    if not _M: _M.update(dict(np.load(OUT)))
    M = _M; e = M["edges"]; z = M["z"]; R = M[f"{q}_{m}_{t}"]
    lo, hi = WIN[win]; sel = (z < 0.3) if win == "skin" else ((z >= lo) & (z <= hi))
    tot = 0.0
    for b in range(len(e) - 1):
        phi = float(Ffun(e[b]) - Ffun(e[b + 1])); c = 1.0
        if correct:
            As = A_bin(m, t, SPECTRA["logflat"], e[b], e[b + 1], WIN[win]); At = A_bin(m, t, Ffun, e[b], e[b + 1], WIN[win])
            if As > 0: c = At / As
        tot += phi * R[b][sel].mean() * c
    return tot * 100            # cGy / cSv

def validate():
    DIR = {"kold1956": "spe_runs/materials_kold1956", "kold1972": "spe_runs/materials_kold1972", "kold1989": "spe_runs/materials_kold1989",
           "webber1956": "spe_runs/materials", "king1972": "spe_runs/materials_king1972"}
    rows = []
    for sp, b in DIR.items():
        dev = {"D": [], "H": []}
        for m in MATS:
            for t in AREALS:
                f = f"{b}/{m}_{t}/dose_eq_depth.csv"
                if not os.path.exists(f): continue
                a = np.genfromtxt(f, delimiter=",", names=True); z = a["depth_g_cm2"]
                for win, (lo, hi) in WIN.items():
                    sel = (z < 0.3) if win == "skin" else ((z >= lo) & (z <= hi))
                    for q, col in (("D", "D_cGy"), ("H", "H_cSv")):
                        p = fold(SPECTRA[sp], m, t, q, win); g = a[col][sel].mean()
                        dev[q].append(p / g); rows.append([sp, m, t, win, q, p, g, p / g])
        for q in ("D", "H"):
            v = np.array(dev[q]); print(f"  {sp:<11} {q}: median {np.median(v):.3f}, range {v.min():.3f}-{v.max():.3f}, |dev| <= 2 %: {100*np.mean(abs(v-1)<=0.02):.0f} %")
    with open("analysis/response_matrix_validation.csv", "w", newline="") as fo:
        w = csv.writer(fo); w.writerow(["spectrum", "material", "areal", "window", "quantity", "folded", "GATE", "ratio"]); w.writerows(rows)
    print("Note: Webber/King start at 10 MeV; 10-30 MeV protons do not cross >= 5 g/cm2 (E_th >= 50 MeV).")

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "build"
    if mode == "build": build()
    elif mode == "validate": validate()
    elif mode == "build_thick": build_thick()
    elif mode == "fold":
        arg = sys.argv[2]
        if arg.isdigit():
            from gle_table import F; Ffun = lambda E: F(E, int(arg))
        else: Ffun = SPECTRA[arg]
        for t in AREALS:
            al = {q: fold(Ffun, "Al", t, q) for q in "DH"}
            print(f"{t:>2} g/cm²: Al D {al['D']:.2f} cGy, H {al['H']:.2f} cSv · " + " · ".join(
                f"{m[:4]} {fold(Ffun, m, t, 'D') / al['D']:.2f}/{fold(Ffun, m, t, 'H') / al['H']:.2f}" for m in MATS if m != "Al"))


def scan(i0=0, i1=None):
    """All GLEs (Koldobskiy et al. 2021) folded with the response matrix: D and H. In chunks: scan i0 i1 (appends to the csv)."""
    from gle_table import T2, F, DATES
    rows = []
    for gle in sorted(T2)[i0:i1]:
        Ff = lambda E, g=gle: F(E, g)
        r = [gle, DATES.get(gle, ""), float(Ff(30)), float(Ff(200))]
        for t in AREALS:
            al = {q: fold(Ff, "Al", t, q) for q in "DH"}
            r += [al["D"], al["H"]] + [fold(Ff, m, t, q) / al[q] for m in ("POLYETHYLENE", "Pb") for q in "DH"]
        r += [1.5 * 10 * fold(Ff, m, 20, "D", "BFO") for m in ("POLYETHYLENE", "Al")]
        r += [10 * fold(Ff, m, 20, "H", "BFO") for m in ("POLYETHYLENE", "Al")]
        rows.append(r)
    hdr = ["GLE", "date", "F30_cm2", "F200_cm2"]
    for t in AREALS:
        hdr += [f"Al_D_{t}_cGy", f"Al_H_{t}_cSv", f"PE_Al_D_{t}", f"PE_Al_H_{t}", f"Pb_Al_D_{t}", f"Pb_Al_H_{t}"]
    hdr += ["BFO_PE20_mGyEq", "BFO_Al20_mGyEq", "BFO_H_PE20_mSv", "BFO_H_Al20_mSv"]
    new = not os.path.exists("analysis/rm_gle_scan.csv") or i0 == 0
    with open("analysis/rm_gle_scan.csv", "w" if new else "a", newline="") as fo:
        w = csv.writer(fo)
        if new: w.writerow(hdr)
        w.writerows(rows)
    print(f"Written: analysis/rm_gle_scan.csv ({len(rows)} GLEs)")

if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "scan":
    scan(int(sys.argv[2]) if len(sys.argv) > 2 else 0, int(sys.argv[3]) if len(sys.argv) > 3 else None)
