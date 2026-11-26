"""gle_rule_td.py - The threshold rule of gle_rule.py tested on a second, independent set of GLE reconstructions:
the Band-function fits of Tylka & Dietrich (2009), extended through 2012, for 53 GLEs of 1956-2012
(Raukunen et al. 2018, Table 2; analysis/data/tylka_dietrich_band.csv). Where a GLE was analysed as separate
GLE and energetic-storm-particle (ESP) components, the two are summed into the event-integrated spectrum.
The rule needs only the slope of the spectrum near the thresholds, so the normalisation convention of J0 does not
enter the ratios.
For each event: folded mean-dose ratio to Al (PE, LiH, water, Pb at 5/10/20 g/cm2), the rule value and gamma(E*).
Also compares PE/Al with the Koldobskiy et al. (2021) reconstruction of the same GLE (analysis/fold_gle_scan.csv).
Writes analysis/gle_rule_td.csv.  Usage (from the repository root):  python analysis/gle_rule_td.py"""
import os, sys, csv
import numpy as np
sys.path.insert(0, os.getcwd()); sys.path.insert(0, "analysis")
import response_matrix as rm

MP2 = 1876.544
ETH = {(r["material"], int(float(r["areal_g_cm2"]))): float(r["E_th_MeV"]) for r in csv.DictReader(open("analysis/thresholds_pstar.csv"))}
MATS = ("POLYETHYLENE", "LITHIUM_HYDRIDE", "WATER", "Pb")


def band(E, J0, g1, g2, R0):
    """Integral fluence J(>E) (cm^-2) of the Band function in rigidity (GV); E in MeV."""
    R = np.sqrt(np.asarray(E, float) * (np.asarray(E, float) + MP2)) / 1000.0
    Rb = (g2 - g1) * R0
    lo = J0 * R ** (-g1) * np.exp(-R / R0)
    hi = J0 * R ** (-g2) * Rb ** (g2 - g1) * np.exp(g1 - g2)
    return np.where(R <= Rb, lo, hi)


def load():
    ev = {}
    for r in csv.DictReader(l for l in open("analysis/data/tylka_dietrich_band.csv") if not l.startswith("#")):
        ev.setdefault(int(r["GLE"]), []).append(tuple(float(r[k]) for k in ("J0", "g1", "g2", "R0")))
    return {g: (lambda E, p=p: sum(band(E, *q) for q in p)) for g, p in sorted(ev.items())}


def gam(F, E): return float(-(np.log(F(E * 1.05)) - np.log(F(E / 1.05))) / np.log(1.05 ** 2))


def main():
    ev = load()
    k21 = {int(r["GLE"]): float(r["PE_Al_20"]) for r in csv.DictReader(open("analysis/fold_gle_scan.csv"))}
    rows = [["GLE", "material", "areal", "gamma_at_Estar", "rule_ratio", "fold_ratio", "K21_fold_ratio"]]
    for g, F in ev.items():
        for m in MATS:
            for t in (5, 10, 20):
                ea, em = ETH[("Al", t)], ETH[(m, t)]; ga = gam(F, np.sqrt(ea * em))
                fo = rm.fold(F, m, t, "D", "mean") / rm.fold(F, "Al", t, "D", "mean")
                kk = k21.get(g, "") if (m == "POLYETHYLENE" and t == 20) else ""
                rows.append([g, m, t, round(ga, 4), round((em / ea) ** (-ga), 4), round(fo, 4), kk])
    with open("analysis/gle_rule_td.csv", "w", newline="") as fo:
        csv.writer(fo).writerows(rows)
    R = rows[1:]
    print(f"Tylka & Dietrich fits: {len(ev)} GLEs")
    for m in MATS:
        for t in (5, 10, 20):
            s = [r for r in R if r[1] == m and r[2] == t]
            d = np.array([r[5] / r[4] for r in s]); x = [r[3] for r in s]; y = np.log([r[5] for r in s])
            print(f"  {m:<16}{t:>3} g/cm2: ratio {min(r[5] for r in s):.2f}-{max(r[5] for r in s):.2f}, "
                  f"gamma {min(x):.2f}-{max(x):.2f}, r(gamma, ln ratio) = {np.corrcoef(x, y)[0, 1]:.3f}, "
                  f"fold/rule median {np.median(d):.3f} ({d.min():.3f}-{d.max():.3f})")
    c = [(r[0], r[5], r[6]) for r in R if r[1] == "POLYETHYLENE" and r[2] == 20 and r[6] != ""]
    dd = np.array([a / b for _, a, b in c])
    print(f"  PE/Al at 20 g/cm2, T&D / K21 for the same GLE ({len(c)} GLEs): median {np.median(dd):.3f}, range {dd.min():.3f}-{dd.max():.3f}")
    for g in (5, 24, 25, 43, 44, 45):
        z = [x for x in c if x[0] == g]
        if z: print(f"     GLE {g}: T&D {z[0][1]:.3f}, K21 {z[0][2]:.3f}")
    print("Written: analysis/gle_rule_td.csv")


if __name__ == "__main__":
    main()
