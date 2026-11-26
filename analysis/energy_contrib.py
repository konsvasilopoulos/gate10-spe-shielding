"""energy_contrib.py - Which primary proton energies give the BFO dose, and where the difference between materials
comes from (cf. Blunier et al. 2026, Space Weather 24, e2025SW004662, Fig. 5).
From the response matrix: contribution of each 25 energy windows (30-5000 MeV) to the BFO dose (4-6 g/cm2 of
water) and to the mean dose over 0-30 g/cm2, D_b = Phi_b R_b c_b, for the five directly simulated spectra,
Al/PE/Pb at 5, 10 and 20 g/cm2.
Reported per case: primary energies below which 5 %, 50 % and 95 % of the dose is delivered (log-interpolated
within a window), the share of the dose from primaries below the CSDA threshold E_th (dose carried by secondaries),
and for PE the energy range that holds 90 % of the dose difference D_Al - D_PE.
Writes analysis/energy_contrib.csv.  Usage (from the repository root):  python analysis/energy_contrib.py"""
import os, sys, csv
import numpy as np
sys.path.insert(0, os.getcwd()); sys.path.insert(0, "analysis")
import response_matrix as rm
from spectra import SPECTRA

ETH = {(r["material"], int(float(r["areal_g_cm2"]))): float(r["E_th_MeV"]) for r in csv.DictReader(open("analysis/thresholds_pstar.csv"))}
SPEC = [("K21 Feb 1956", "kold1956"), ("Webber Feb 1956", "webber1956"), ("K21 Aug 1972", "kold1972"),
        ("King Aug 1972", "king1972"), ("K21 Sep-Oct 1989", "kold1989")]


def bins(F, m, t, win="BFO"):
    if not rm._M: rm._M.update(dict(np.load(rm.OUT)))
    M = rm._M; e = M["edges"]; z = M["z"]; R = M[f"D_{m}_{t}"]
    lo, hi = rm.WIN[win]; sel = (z >= lo) & (z <= hi)
    out = []
    for b in range(len(e) - 1):
        phi = float(F(e[b]) - F(e[b + 1]))
        As = rm.A_bin(m, t, SPECTRA["logflat"], e[b], e[b + 1], rm.WIN[win]); At = rm.A_bin(m, t, F, e[b], e[b + 1], rm.WIN[win])
        c = At / As if As > 0 else 1.0
        out.append(phi * R[b][sel].mean() * c)
    return e, np.array(out)


def pct(e, d, p):
    """Primary energy below which a fraction p of sum(d) is delivered (log-linear within a window)."""
    c = np.concatenate([[0.0], np.cumsum(d)]) / d.sum()
    i = int(np.searchsorted(c, p)); i = min(max(i, 1), len(e) - 1)
    f = (p - c[i - 1]) / (c[i] - c[i - 1]) if c[i] > c[i - 1] else 0.0
    return float(np.exp(np.log(e[i - 1]) + f * (np.log(e[i]) - np.log(e[i - 1]))))


def below(e, d, E):
    """Share of sum(d) from primaries below E (log-linear split of the window containing E)."""
    s = 0.0
    for b in range(len(d)):
        if e[b + 1] <= E: s += d[b]
        elif e[b] < E: s += d[b] * np.log(E / e[b]) / np.log(e[b + 1] / e[b])
    return s / d.sum()


def main():
    rows = [["window", "spectrum", "material", "areal", "E5_MeV", "E50_MeV", "E95_MeV", "E_th_MeV", "share_below_Eth", "diff_E5_MeV", "diff_E95_MeV"]]
    for lab, k in SPEC:
        F = SPECTRA[k]
        for win in ("BFO", "mean"):
          for t in (5, 10, 20):
            e, dAl = bins(F, "Al", t, win)
            for m in ("Al", "POLYETHYLENE", "Pb"):
                e, d = (e, dAl) if m == "Al" else bins(F, m, t, win)
                r = [win, lab, m, t] + [round(pct(e, d, p), 1) for p in (0.05, 0.5, 0.95)] + [round(ETH[(m, t)], 1), round(below(e, d, ETH[(m, t)]), 4)]
                if m == "POLYETHYLENE":
                    diff = np.clip(dAl - d, 0, None)
                    r += [round(pct(e, diff, 0.05), 1), round(pct(e, diff, 0.95), 1)]
                else:
                    r += ["", ""]
                rows.append(r)
                print("  " + "  ".join(str(x) for x in r))
    with open("analysis/energy_contrib.csv", "w", newline="") as fo:
        csv.writer(fo).writerows(rows)
    print("Written: analysis/energy_contrib.csv")


if __name__ == "__main__":
    main()
