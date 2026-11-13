"""oltaris_larc.py - The OLTARIS historical spectra "Feb 1956 (LaRC)" and "Aug 1972 (LaRC)" as alternative
reconstructions of the February 1956 and August 1972 events.
Source (OLTARIS team, personal communication, Nov 2026): Wilson & Denn 1976, NASA TN D-8290, with differential
fluences (cm-2 MeV-1, E in MeV)
  Feb 1956 (LaRC): f(E) = 1.5e9 exp(-(E-10)/25)/25 + 3e8 exp(-(E-100)/320)/320
  Aug 1972 (LaRC): f(E) = 2.2e7 exp(-(E-100)/30)
(larc_differential below; tests/test_oltaris.py checks that the saved boundary fluxes follow them).

OLTARIS runs (TARIS 5.01, 12 Nov 2026): Al slab 5 + 5 + 10 g/cm2 followed by water 0.3 + 4.7 + 5 + 20 g/cm2, as in
oltaris_compare.py, with "Save external differential flux" = Yes -> oltaris/Al_LaRC56/, oltaris/Al_LaRC72/.
The boundary proton spectrum (li_flux.dat) is folded with the GATE response matrix and compared with
  (a) HZETRN for the same Al slab (D and H at all interfaces),
  (b) the spectra of the article (Webber and Koldobskiy GLE 5 for 1956; King and Koldobskiy GLE 24+25 for 1972):
      integral fluence, slope gamma near the thresholds, ratios to Al, threshold rule, BFO limit.
Usage (from the repository root):  python analysis/oltaris_larc.py
  -> analysis/oltaris_larc.csv (one row per spectrum), analysis/oltaris_larc_hzetrn.csv (GATE vs HZETRN, Al slab)
"""
import os, sys, csv
import numpy as np
sys.path.insert(0, os.getcwd()); sys.path.insert(0, "analysis")
import response_matrix as rm
import nasa_limits_thick as nl
from spectra import SPECTRA
from oltaris_compare import read_table, POINTS, BAD, gate_point

ETH = {(r["material"], int(float(r["areal_g_cm2"]))): float(r["E_th_MeV"])
       for r in csv.DictReader(open("analysis/thresholds_pstar.csv"))}


def larc_differential(tag, E):
    """Differential proton fluence (cm-2 MeV-1) of the OLTARIS LaRC spectra (Wilson & Denn 1976)."""
    E = np.asarray(E, float)
    if tag == "LaRC56":
        return 1.5e9 * np.exp(-(E - 10) / 25) / 25 + 3e8 * np.exp(-(E - 100) / 320) / 320
    if tag == "LaRC72":
        return 2.2e7 * np.exp(-(E - 100) / 30)
    raise ValueError(tag)


def flux_integral(tag):
    E, phi, parts = read_table(f"oltaris/Al_{tag}/li_flux.dat")
    p = phi[parts.index("proton")]; ok = p > 0
    E, p = E[ok], p[ok]
    Eg = np.geomspace(E[0], E[-1], 20001); pg = np.exp(np.interp(np.log(Eg), np.log(E), np.log(p)))
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (pg[1:] + pg[:-1]) * np.diff(Eg))])
    f = lambda x: np.interp(np.log(np.clip(np.asarray(x, float), E[0], E[-1])), np.log(Eg), cum[-1] - cum)
    return f, E[-1]


def gam(F, E): return float(-(np.log(F(E * 1.05)) - np.log(F(E / 1.05))) / np.log(1.05 ** 2))


def main():
    hz = [["spectrum", "point", "GATE_D_mGy", "HZETRN_D_mGy", "D_ratio", "GATE_H_mSv", "HZETRN_H_mSv", "H_ratio"]]
    rows = [["spectrum", "Phi_gt30", "Phi_gt100", "Phi_gt300", "Phi_gt1000", "gamma_150MeV", "gamma_Estar_PE5", "gamma_Estar_PE10",
             "gamma_Estar_PE20", "PE_Al_mean_D_5", "PE_Al_mean_D_10", "PE_Al_mean_D_20", "Pb_Al_mean_D_20", "rule_PE_Al_20",
             "PE_Al_BFO_D_20", "BFO_Al20_mGyEq", "min_areal_Al", "min_areal_PE", "min_areal_Pb"]]
    L56, e56 = flux_integral("LaRC56"); L72, e72 = flux_integral("LaRC72")
    groups = {"1956": [("LaRC", L56), ("Webber", SPECTRA["webber1956"]), ("K21", SPECTRA["kold1956"])],
              "1972": [("LaRC", L72), ("King", SPECTRA["king1972"]), ("K21", SPECTRA["kold1972"])]}
    print(f"Highest energy in li_flux.dat: 1956 LaRC {e56:.0f} MeV, 1972 LaRC {e72:.0f} MeV")

    print("\n(a) GATE (response matrix, LaRC boundary spectrum) / HZETRN, Al slab")
    for ev, tag, F in (("1956", "LaRC56", L56), ("1972", "LaRC72", L72)):
        _, oD = read_table(f"oltaris/Al_{tag}/li_dose_database_mGy.dat"); _, oH = read_table(f"oltaris/Al_{tag}/li_doseq_database_mSv.dat")
        for lab, t, win, i in POINTS:
            gD, gH = gate_point(F, "Al", t, "D", win), gate_point(F, "Al", t, "H", win)
            print(f"  {ev} {lab:<20} D {gD:9.2f} / {oD[i]:9.2f} = {gD / oD[i]:.3f}   H {gH:9.2f} / {oH[i]:9.2f} = {gH / oH[i]:.3f}")
            hz.append([f"{ev} LaRC", lab, round(gD, 4), round(oD[i], 4), round(gD / oD[i], 3), round(gH, 4), round(oH[i], 4), round(gH / oH[i], 3)])

    print("\n(b) Spectra through the GATE response matrix")
    for ev, lst in groups.items():
        print(f"  --- {ev}")
        print(f"  {'':<8}" + "".join(f"{h:>10}" for h in ("F>30", "F>100", "F>300", "F>1000", "g(150)", "gPE5", "gPE10", "gPE20")))
        spec = {}
        for name, F in lst:
            g = [gam(F, np.sqrt(ETH[("Al", t)] * ETH[("POLYETHYLENE", t)])) for t in (5, 10, 20)]
            vals = [float(F(e)) for e in (30, 100, 300, 1000)]
            print(f"  {name:<8}" + "".join(f"{v:>10.3g}" for v in vals) + f"{gam(F, 150):>10.2f}" + "".join(f"{x:>10.2f}" for x in g))
            spec[name] = [f"{v:.4g}" for v in vals] + [round(gam(F, 150), 3)] + [round(x, 3) for x in g]
        print(f"  {'':<8}" + "".join(f"{h:>13}" for h in ("PE/Al 5", "PE/Al 10", "PE/Al 20", "Pb/Al 20", "rule PE 20", "PE/Al BFO20", "Al20 BFO", "min Al", "min PE", "min Pb")))
        for name, F in lst:
            r = [rm.fold(F, "POLYETHYLENE", t, "D", "mean") / rm.fold(F, "Al", t, "D", "mean") for t in (5, 10, 20)]
            pb = rm.fold(F, "Pb", 20, "D", "mean") / rm.fold(F, "Al", 20, "D", "mean")
            ea, ep = ETH[("Al", 20)], ETH[("POLYETHYLENE", 20)]; rule = (ep / ea) ** (-gam(F, np.sqrt(ea * ep)))
            bfo = rm.fold(F, "POLYETHYLENE", 20, "D", "BFO") / rm.fold(F, "Al", 20, "D", "BFO")
            bal = 1.5 * 10 * rm.fold(F, "Al", 20, "D", "BFO")
            mins = [nl.tmin([1.5 * 10 * rm.fold(F, m, t, "D", "BFO") for t in nl.T]) for m in ("Al", "POLYETHYLENE", "Pb")]
            print(f"  {name:<8}" + "".join(f"{v:>13.3f}" for v in r + [pb, rule, bfo]) + f"{bal:>13.1f}" + "".join(f"{m:>13}" for m in mins))
            rows.append([f"{ev} {name}"] + spec[name] + [round(v, 3) for v in r + [pb, rule, bfo]] + [round(bal, 1)] + mins)
    with open("analysis/oltaris_larc.csv", "w", newline="") as fo:
        csv.writer(fo).writerows(rows)
    with open("analysis/oltaris_larc_hzetrn.csv", "w", newline="") as fo:
        csv.writer(fo).writerows(hz)
    print("\nWritten: analysis/oltaris_larc.csv, analysis/oltaris_larc_hzetrn.csv")


if __name__ == "__main__":
    main()
