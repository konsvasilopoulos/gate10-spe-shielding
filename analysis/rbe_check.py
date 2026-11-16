"""rbe_check.py - Effect of particle-specific RBE values on the BFO dose, from the particle-resolved HZETRN (OLTARIS)
doses at 5 g/cm2 of water behind 20 g/cm2 of shield.

The article converts the BFO dose to mGy-Eq with RBE = 1.5 (protons above 2 MeV) for the whole dose.
NCRP Report 132 (2000), used for the NASA short-term limits, gives 1.5 for protons and 2.5 for helium and heavier
ions; electrons, photons, muons and pions are given RBE = 1 here. Deuterons and tritons are counted with the protons.
HZETRN deposits the dose of neutrons through the charged particles they produce (neutron column = 0), as GATE does,
so this check does not cover a separate neutron RBE.
For the August 1972 reconstruction the two OLTARIS projects (GLE 24 and GLE 25) are summed after the J0 rescaling
of oltaris_compare.py.
Usage (from the repository root):  python analysis/rbe_check.py  -> analysis/rbe_check_hzetrn.csv
"""
import os, sys, csv
import numpy as np
sys.path.insert(0, os.getcwd()); sys.path.insert(0, "analysis")
from oltaris_compare import read_table, SCALE, EXTRA

BFO = 5                                   # index of 5 g/cm2 of water behind 20 g/cm2 in the OLTARIS depth grid
RBE = {"H": 1.5, "He": 2.5, "Z>=3": 2.5, "pi/mu": 1.0, "e/gamma": 1.0}
CASES = [("1972 King", ["King72"]), ("1956 Webber", ["Webber56"]), ("1972 K21", ["GLE24", "GLE25"]),
         ("Oct 1989 Tylka", ["Tylka89"]), ("1956 LaRC", ["LaRC56"]), ("1972 LaRC", ["LaRC72"])]
GROUP = {"proton": "H", "deuteron": "H", "triton": "H", "helion": "He", "alpha": "He",
         "mu-": "pi/mu", "mu+": "pi/mu", "pi-": "pi/mu", "pi+": "pi/mu",
         "electron": "e/gamma", "positron": "e/gamma", "photon": "e/gamma", "neutron": "n"}


def shares(mat, tags):
    tot = None
    for t in tags:
        d = f"oltaris/{mat}_{t}"
        if not os.path.isdir(d): return None
        _, v, parts = read_table(f"{d}/li_dose_part_database_mGy.dat")
        x = v[:, BFO] / (SCALE.get(t, 1.0) * EXTRA.get((mat, t), 1.0))
        tot = x if tot is None else tot + x
    g = {k: 0.0 for k in list(RBE) + ["n"]}
    for p, x in zip(parts, tot / tot.sum()):
        g[GROUP.get(p, "Z>=3")] += x
    return g


def main():
    rows = [["spectrum", "shield", "H", "He", "Z>=3", "pi/mu", "e/gamma", "neutron", "BFO_change_pct"]]
    for lab, tags in CASES:
        for mat in ("Al", "PE", "Pb"):
            g = shares(mat, tags)
            if g is None: continue
            dch = 100 * (sum(RBE[k] * g[k] for k in RBE) / 1.5 - 1)
            rows.append([lab, mat] + [round(g[k], 4) for k in ("H", "He", "Z>=3", "pi/mu", "e/gamma", "n")] + [round(dch, 2)])
            print(f"  {lab:<15}{mat:<3} " + "  ".join(f"{k} {g[k]:.4f}" for k in ("H", "He", "Z>=3", "pi/mu", "e/gamma", "n")) + f"   BFO {dch:+.1f} %")
    with open("analysis/rbe_check_hzetrn.csv", "w", newline="") as fo:
        csv.writer(fo).writerows(rows)
    print("Written: analysis/rbe_check_hzetrn.csv")


if __name__ == "__main__":
    main()
