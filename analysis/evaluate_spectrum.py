"""evaluate_spectrum.py - Shielding of any solar-particle-event proton spectrum without new simulations.

Folds a user spectrum with the GATE response matrix (analysis/response_matrix.npz: nine materials, 5-40 g/cm2,
25 energy windows 30-5000 MeV) and reports, for each material and areal density:
  - BFO dose (4-6 g/cm2 of water) in mGy and mGy-Eq (RBE 1.5), and the mean dose over 0-30 g/cm2 of water;
  - ratios to Al at equal areal density (mean and BFO dose, mean dose equivalent);
  - the threshold-rule estimate (E_th,m / E_th,Al)^(-gamma), with gamma the local slope of the integral spectrum
    at E* = sqrt(E_th,m E_th,Al) and E_th the CSDA range threshold from PSTAR;
  - the minimum areal density for the NASA 30-day BFO limit of 250 mGy-Eq (log interpolation in 5-40 g/cm2).
Geometry of the matrix: broad parallel beam at normal incidence on a slab, water phantom behind it (see the article).

Spectrum input (one of):
  --band J0 g1 g2 R0      Band function in rigidity, integral fluence J(>R) in cm^-2, R and R0 in GV
                          (the form of Tylka & Dietrich 2009 / Raukunen et al. 2018; repeat --band to sum components)
  --gle N                 Koldobskiy et al. (2021) reconstruction of GLE N (analysis/gle_table.py)
  --spectrum NAME         built-in: webber1956, king1972, kold1956, kold1972, kold1989
  --table FILE            CSV or whitespace file with columns E (MeV) and integral fluence J(>E) (cm^-2);
                          add --differential if the second column is dN/dE (cm^-2 MeV^-1)
Options: --materials Al POLYETHYLENE ... (default: all nine), --areals 5 10 20 30 40, --csv OUT.csv
The spectrum must cover 30 MeV; above its last tabulated energy the fluence is taken as zero.

Python use:  from evaluate_spectrum import evaluate;  res = evaluate(F)   with F(E_MeV) -> J(>E) in cm^-2
Usage examples (from the repository root):
  python analysis/evaluate_spectrum.py --spectrum kold1972
  python analysis/evaluate_spectrum.py --band 1.75e8 1.76 5.04 0.566
  python analysis/evaluate_spectrum.py --table my_event.csv --materials Al POLYETHYLENE WATER --csv my_event_out.csv
"""
import os, sys, csv, argparse
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
os.chdir(ROOT)
import response_matrix as rm

MATS = ["POLYETHYLENE", "LITHIUM_HYDRIDE", "WATER", "PLEXIGLASS", "Al", "Ti", "Fe", "Cu", "Pb"]
AREALS = [5, 10, 20, 30, 40]
LIMIT, RBE, MP2 = 250.0, 1.5, 1876.544


# ---------------------------------------------------------------- spectra
def band(J0, g1, g2, R0):
    """Integral fluence J(>E) (cm^-2) of a Band function in rigidity (GV); returns a function of E (MeV)."""
    Rb = (g2 - g1) * R0
    def F(E):
        R = np.sqrt(np.asarray(E, float) * (np.asarray(E, float) + MP2)) / 1000.0
        return np.where(R <= Rb, J0 * R ** (-g1) * np.exp(-R / R0), J0 * R ** (-g2) * Rb ** (g2 - g1) * np.exp(g1 - g2))
    return F


def tabulated(E, y, differential=False):
    """Integral fluence from a table (log-log interpolation); a differential table is integrated first."""
    E = np.asarray(E, float); y = np.asarray(y, float); o = np.argsort(E); E, y = E[o], y[o]
    if differential:
        Eg = np.geomspace(E[0], E[-1], 4000); yg = np.exp(np.interp(np.log(Eg), np.log(E), np.log(np.maximum(y, 1e-300))))
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (yg[1:] + yg[:-1]) * np.diff(Eg))]); E, y = Eg, cum[-1] - cum
    y = np.maximum(y, 1e-300)
    def F(x):
        x = np.asarray(x, float)
        v = np.exp(np.interp(np.log(np.clip(x, E[0], E[-1])), np.log(E), np.log(y)))
        return np.where(x > E[-1], 0.0, v)
    F.emin, F.emax = E[0], E[-1]
    return F


def read_table(path):
    rows = []
    for line in open(path):
        p = line.replace(",", " ").replace(";", " ").split()
        try: rows.append((float(p[0]), float(p[1])))
        except (ValueError, IndexError): continue                      # header or comment
    a = np.array(rows); return a[:, 0], a[:, 1]


# ---------------------------------------------------------------- physics helpers
def e_th(m, t):
    """CSDA threshold (MeV): energy whose PSTAR range equals the areal density t (g/cm2)."""
    a = np.array([[float(x) for x in l.split()] for l in open(f"pstar/G4_{m}.txt") if len(l.split()) == 7])
    return float(np.exp(np.interp(np.log(t), np.log(a[:, 4]), np.log(a[:, 0]))))


def slope(F, E):
    """gamma(E) = -dlnJ/dlnE of the integral spectrum."""
    return float(-(np.log(F(E * 1.05)) - np.log(F(E / 1.05))) / np.log(1.05 ** 2))


def t_limit(t, b):
    """Areal density (g/cm2) where the BFO dose b(t) (mGy-Eq) crosses LIMIT; log interpolation."""
    t, b = np.asarray(t, float), np.asarray(b, float)
    if b[0] <= LIMIT: return f"<={t[0]:g}"
    if b[-1] > LIMIT: return f">{t[-1]:g}"
    i = int(np.where(b <= LIMIT)[0][0])
    return f"{np.interp(np.log(LIMIT), [np.log(b[i]), np.log(b[i - 1])], [t[i], t[i - 1]]):.1f}"


# ---------------------------------------------------------------- main evaluation
def evaluate(F, mats=MATS, areals=AREALS):
    """Returns (rows, limits): rows = list of dicts per material and areal density; limits = {material: areal density}."""
    if float(F(30.0)) <= 0: raise ValueError("the spectrum must have protons above 30 MeV")
    mats = list(dict.fromkeys(["Al"] + list(mats)))
    rows, limits, cache = [], {}, {}
    for m in mats:
        for t in areals:
            cache[(m, t)] = {q: rm.fold(F, m, t, q[0], q[1]) for q in (("D", "BFO"), ("D", "mean"), ("H", "mean"))}
    for m in mats:
        bfo = []
        for t in areals:
            v, a = cache[(m, t)], cache[("Al", t)]
            ea, em = e_th("Al", t), e_th(m, t); es = np.sqrt(ea * em); g = slope(F, es)
            d_bfo = 10 * v[("D", "BFO")]                                   # cGy -> mGy
            rows.append(dict(material=m, areal_g_cm2=t, BFO_D_mGy=d_bfo, BFO_mGyEq=RBE * d_bfo, mean_D_mGy=10 * v[("D", "mean")],
                             ratio_mean_D=v[("D", "mean")] / a[("D", "mean")], ratio_BFO_D=v[("D", "BFO")] / a[("D", "BFO")],
                             ratio_mean_H=v[("H", "mean")] / a[("H", "mean")], E_th_MeV=em, E_star_MeV=es, gamma=g,
                             rule_ratio=(em / ea) ** (-g)))
            bfo.append(RBE * d_bfo)
        limits[m] = t_limit(areals, bfo) if len(areals) > 1 else "n/a (give several areal densities)"
    return rows, limits


def report(F, rows, limits, label):
    print(f"Spectrum: {label}")
    print(f"  J(>30 MeV) = {float(F(30)):.3g} cm^-2,  J(>100) = {float(F(100)):.3g},  J(>200) = {float(F(200)):.3g},  "
          f"gamma at 163 MeV = {slope(F, 163.0):.2f}")
    print(f"\n{'material':<16}{'g/cm2':>6}{'BFO mGy':>11}{'BFO mGy-Eq':>12}{'mean/Al':>9}{'BFO/Al':>8}{'H mean/Al':>10}"
          f"{'E_th':>7}{'gamma':>7}{'rule':>7}")
    for r in rows:
        print(f"{r['material']:<16}{r['areal_g_cm2']:>6}{r['BFO_D_mGy']:>11.4g}{r['BFO_mGyEq']:>12.4g}{r['ratio_mean_D']:>9.3f}"
              f"{r['ratio_BFO_D']:>8.3f}{r['ratio_mean_H']:>10.3f}{r['E_th_MeV']:>7.0f}{r['gamma']:>7.2f}{r['rule_ratio']:>7.3f}")
    print(f"\nMinimum areal density (g/cm2) for BFO < {LIMIT:g} mGy-Eq: " + ", ".join(f"{m} {v}" for m, v in limits.items()))
    print("Ratios are to Al at the same areal density. The rule estimates the mean-dose ratio; it was tested at 5-20 g/cm2 and is\n"
          "within about 5 % for the light materials at 10-20 g/cm2, less accurate at 5 g/cm2 and for Pb (see the article).")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--band", nargs=4, type=float, action="append", metavar=("J0", "G1", "G2", "R0"))
    g.add_argument("--gle", type=int); g.add_argument("--spectrum"); g.add_argument("--table")
    ap.add_argument("--differential", action="store_true")
    ap.add_argument("--materials", nargs="+", default=MATS); ap.add_argument("--areals", nargs="+", type=float, default=AREALS)
    ap.add_argument("--csv")
    a = ap.parse_args(argv)
    if a.band:
        comps = [band(*p) for p in a.band]; F = lambda E: sum(c(E) for c in comps); label = f"Band function, {len(comps)} component(s)"
    elif a.gle is not None:
        from gle_table import F as Fg, valid
        if not valid(a.gle): sys.exit(f"GLE {a.gle}: no valid Koldobskiy et al. (2021) parameters")
        F = lambda E: Fg(E, a.gle); label = f"GLE {a.gle}, Koldobskiy et al. (2021)"
    elif a.spectrum:
        from spectra import SPECTRA
        F = SPECTRA[a.spectrum]; label = a.spectrum
    else:
        E, y = read_table(a.table); F = tabulated(E, y, a.differential); label = a.table
        if F.emin > 30: sys.exit(f"the table starts at {F.emin:g} MeV; it must reach down to 30 MeV")
        if F.emax < 1000: print(f"Note: the table ends at {F.emax:g} MeV; the fluence above is taken as zero.")
    areals = [int(t) if float(t).is_integer() else t for t in a.areals]
    bad = [t for t in areals if t not in AREALS]
    if bad: sys.exit(f"areal densities must be among {AREALS} (response matrix); got {bad}")
    rows, limits = evaluate(F, a.materials, areals)
    report(F, rows, limits, label)
    if a.csv:
        with open(a.csv, "w", newline="") as fo:
            w = csv.DictWriter(fo, fieldnames=list(rows[0]) + ["min_areal_for_limit"]); w.writeheader()
            for r in rows: w.writerow({**{k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()}, "min_areal_for_limit": limits[r["material"]]})
        print(f"Written: {a.csv}")


if __name__ == "__main__":
    main()
