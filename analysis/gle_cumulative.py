"""gle_cumulative.py - Does the size of an event predict the benefit of a light shield, and what is the benefit
summed over all reconstructed GLEs?
For the 58 GLEs of Koldobskiy et al. (2021) (analysis/gle_table.py) and, as a check, the 59 Band fits of
Tylka & Dietrich (Raukunen et al. 2018; analysis/data/tylka_dietrich_band.csv), behind 20 g/cm2 (and 10 g/cm2):
  - correlation between the size of the event (log of the Al BFO dose) and the PE/Al ratio of the BFO dose;
  - the ratio of the BFO doses summed over all events (PE, LiH, water, Pb to Al): the benefit a shield of that
    material would have given over the whole GLE record;
  - the share of the summed Al BFO dose carried by the five largest events.
Writes analysis/gle_cumulative.csv.  Usage (from the repository root):  python analysis/gle_cumulative.py"""
import os, sys, csv
import numpy as np
sys.path.insert(0, os.getcwd()); sys.path.insert(0, "analysis")
import response_matrix as rm
from gle_table import F as Fk, T2, valid
from gle_rule_td import load as load_td

MATS = ["Al", "POLYETHYLENE", "LITHIUM_HYDRIDE", "WATER", "Pb"]


def doses(events, t):
    return {g: {m: 10 * rm.fold(F, m, t, "D", "BFO") for m in MATS} for g, F in events.items()}


def main():
    sets = {"Koldobskiy2021": {g: (lambda E, g=g: Fk(E, g)) for g in sorted(T2) if valid(g)},
            "TylkaDietrich": load_td()}
    rows = [["set", "areal", "n_events", "r_logAl_vs_PE_Al", "PE_Al_median", "PE_Al_min", "PE_Al_max",
             "sum_PE_Al", "sum_LiH_Al", "sum_water_Al", "sum_Pb_Al", "top5_share", "largest_GLE", "largest_PE_Al"]]
    for name, ev in sets.items():
        for t in (10, 20):
            d = doses(ev, t); G = list(d)
            al = np.array([d[g]["Al"] for g in G]); pe = np.array([d[g]["POLYETHYLENE"] for g in G]) / al
            r = np.corrcoef(np.log(al), pe)[0, 1]
            s = {m: sum(d[g][m] for g in G) / al.sum() for m in MATS}
            o = np.argsort(-al); top5 = al[o[:5]].sum() / al.sum()
            rows.append([name, t, len(G), round(r, 3), round(float(np.median(pe)), 3), round(pe.min(), 3), round(pe.max(), 3),
                         round(s["POLYETHYLENE"], 3), round(s["LITHIUM_HYDRIDE"], 3), round(s["WATER"], 3), round(s["Pb"], 3),
                         round(top5, 3), G[o[0]], round(pe[o[0]], 3)])
            print(f"  {name:<15}{t:>3} g/cm2, {len(G)} GLEs: r(log D_Al, PE/Al) = {r:+.2f}; PE/Al per event {pe.min():.2f}-{pe.max():.2f} "
                  f"(median {np.median(pe):.2f}); summed: PE {s['POLYETHYLENE']:.2f}, LiH {s['LITHIUM_HYDRIDE']:.2f}, water {s['WATER']:.2f}, "
                  f"Pb {s['Pb']:.2f}; five largest events {100 * top5:.0f} % of the Al dose; largest GLE {G[o[0]]} (PE/Al {pe[o[0]]:.2f})")
    with open("analysis/gle_cumulative.csv", "w", newline="") as fo:
        csv.writer(fo).writerows(rows)
    print("Written: analysis/gle_cumulative.csv")


if __name__ == "__main__":
    main()
