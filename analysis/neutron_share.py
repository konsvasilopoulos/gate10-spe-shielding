"""neutron_share.py - Share of the dose due to neutrons, from run_neutron_share.sh.
For each spectrum and shield at 20 g/cm2: D(full) and D(noN) (all hadronic processes of neutrons inactivated) at the
BFO depth (4-6 g/cm2 of water) and averaged over 0-30 g/cm2. Neutron share f_n = 1 - D(noN)/D(full); it includes
everything that neutrons produce (recoil protons, fragments, capture and inelastic photons).
Effect on the BFO dose in mGy-Eq if this part is weighted with a neutron RBE instead of 1.5:
dBFO = f_n (RBE_n / 1.5 - 1), for RBE_n = 3.5 (5-50 MeV) and 6.0 (1-5 MeV), NCRP Report 132.
The quoted uncertainty treats the two runs as independent and overestimates it (the runs are correlated).
The simulation outputs are part of the Zenodo dataset; without them the script stops with a message.
Usage (from the repository root):  python analysis/neutron_share.py  -> analysis/neutron_share.csv
"""
import os, sys, json, glob, csv
import numpy as np

BASE = "spe_runs/neutron_share"


def read_mhd(path):
    hdr = dict(l.split(" = ", 1) for l in open(path).read().splitlines() if " = " in l)
    dt = {"MET_DOUBLE": np.float64, "MET_FLOAT": np.float32}[hdr["ElementType"].strip()]
    return np.fromfile(os.path.join(os.path.dirname(path), hdr["ElementDataFile"].strip()), dtype=dt).ravel()


def event_dose(d):
    """Event dose (cGy) per depth voxel and its variance, summed over the energy groups."""
    tot = var = None; ng = 0
    for g in sorted(glob.glob(f"{d}/g*/setup.json")):
        s = json.load(open(g)); gd = os.path.dirname(g)
        D = read_mhd(f"{gd}/dose_{s['material']}_dose.mhd"); u = read_mhd(f"{gd}/dose_{s['material']}_dose_uncertainty.mhd")
        x = D * s["group_fluence_cm2"] * s["source_area_cm2"] / s["n"] * 100.0           # Gy -> event cGy
        tot = x if tot is None else tot + x
        var = (x * u) ** 2 if var is None else var + (x * u) ** 2
        ng += 1
    return tot, var, ng


def window(tot, var, z0, z1, dz=0.1):
    i0, i1 = int(round(z0 / dz)), int(round(z1 / dz))
    return tot[i0:i1].mean(), np.sqrt(var[i0:i1].sum()) / (i1 - i0)


def main():
    if not os.path.isdir(BASE):
        sys.exit(f"{BASE} not found (simulation outputs are in the Zenodo dataset); analysis/neutron_share.csv holds the results.")
    rows = [["spectrum", "material", "window", "D_full_cGy", "D_noN_cGy", "neutron_share", "share_sigma_uncorrelated",
             "dBFO_RBE3.5_pct", "dBFO_RBE6_pct", "groups_full", "groups_noN"]]
    for sp in ("kold1972", "kold1956"):
        for m in ("Al", "POLYETHYLENE", "Pb"):
            f, n = f"{BASE}/{sp}/{m}_20_full", f"{BASE}/{sp}/{m}_20_noN"
            if not (os.path.isdir(f) and os.path.isdir(n)): continue
            tf, vf, gf = event_dose(f); tn, vn, gn = event_dose(n)
            for lab, (z0, z1) in (("BFO", (4.0, 6.0)), ("mean", (0.0, 30.0))):
                Df, sf = window(tf, vf, z0, z1); Dn, sn = window(tn, vn, z0, z1)
                fn = 1 - Dn / Df; sig = (Dn / Df) * np.hypot(sf / Df, sn / Dn)
                d35, d6 = 100 * fn * (3.5 / 1.5 - 1), 100 * fn * (6.0 / 1.5 - 1)
                rows.append([sp, m, lab, round(Df, 4), round(Dn, 4), round(fn, 4), round(sig, 4), round(d35, 2), round(d6, 2), gf, gn])
                print(f"  {sp:<10}{m:<14}{lab:<6} D {Df:8.3f} / {Dn:8.3f} cGy  f_n {100 * fn:5.2f} +- {100 * sig:4.2f} %"
                      f"  dBFO {d35:+5.1f} % (RBE 3.5)  {d6:+5.1f} % (RBE 6)")
    with open("analysis/neutron_share.csv", "w", newline="") as fo:
        csv.writer(fo).writerows(rows)
    print("Written: analysis/neutron_share.csv")


if __name__ == "__main__":
    main()
