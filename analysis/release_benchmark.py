"""release_benchmark.py - GATE 10.0.0 (Geant4 11.2.1) and GATE 10.1.1 (Geant4 11.4.2) against the external
references used in the paper. The set of comparisons was fixed before the 10.1.1 runs (3 Nov 2026):
  1) NASA benchmark (Heinbockel et al. 2009): Webber 1956, 20 g/cm2 Al, D and H in water vs FLUKA, HETC-HEDS, HZETRN
  2) HZETRN (OLTARIS): King 1972, Webber 1956 and Koldobskiy 1972 (GLE 24+25), Al / PE / Pb at 20 g/cm2,
     D and H at the water entrance and at 0.3, 5 (BFO), 10 and 30 g/cm2 of water.
Both releases are compared from DIRECT simulations (same settings, seed 1234): analysis/release_data/<version>/.
(The paper's Table 2 uses the response matrix of 10.0.0; the direct 10.0.0 runs agree with it within 1.5 % in D.)
10.1.1 runs: run_release_check.sh (Koldobskiy 1972) and run_release_benchmark.sh (NASA benchmark, King, Webber).
Usage:  python analysis/release_benchmark.py  -> analysis/release_benchmark_nasa.csv, analysis/release_benchmark_hzetrn.csv
"""
import os, sys, csv
import numpy as np
sys.path.insert(0, os.getcwd()); sys.path.insert(0, "analysis")
import oltaris_compare as oc

R = "analysis/release_data"
VERS = ["v10.0.0", "v10.1.1"]
MATS = {"Al": "Al", "POLYETHYLENE": "PE", "Pb": "Pb"}
SPEC = {"king1972": ["King72"], "webber1956": ["Webber56"], "kold1972": ["GLE24", "GLE25"]}
# water points behind 20 g/cm2: (label, GATE depth window g/cm2, index in the OLTARIS depth array)
PTS = [("entrance", (0, 0.1), 3), ("0.3", (0.2, 0.4), 4), ("5 (BFO)", (4.5, 5.5), 5), ("10", (9.5, 10.5), 6), ("30", (29.5, 30), 7)]


def gate_point(v, sp, m, win):
    a = np.genfromtxt(f"{R}/{v}/{sp}/{m}_20_dose_eq_depth.csv", delimiter=",", names=True)
    z = a["depth_g_cm2"]; s = (z >= win[0]) & (z <= win[1])
    return 10 * a["D_cGy"][s].mean(), 10 * a["H_cSv"][s].mean()          # mGy, mSv


def nasa():
    rows = [["version", "quantity", "depth_g_cm2", "GATE", "FLUKA", "HETC", "HZETRN", "GATE_over_HZETRN",
             "within_code_spread", "deviation_from_spread_pct"]]
    summary = {}
    for v in VERS:
        d = list(csv.DictReader(open(f"{R}/{v}/nasa_benchmark/compare_nasa.csv")))
        h = list(csv.DictReader(open(f"{R}/{v}/nasa_benchmark/dose_eq_compare_nasa.csv")))
        for q, data, g, refs in (("D (cGy)", d, "GATE_cGy", ("FLUKA_cGy", "HETC_cGy", "HZETRN_cGy")),
                                 ("H (cSv)", h, "H_cSv", ("H_FLUKA", "H_HETC_corr", "H_HZETRN"))):
            n_in = 0; devs = []
            for r in data:
                x = float(r[g]); ref = [float(r[k]) for k in refs]; lo, hi = min(ref), max(ref)
                dev = 0.0 if lo <= x <= hi else 100 * ((x - lo) / lo if x < lo else (x - hi) / hi)
                n_in += dev == 0.0; devs.append(dev)
                rows.append([v, q, r["depth_g_cm2"], round(x, 4)] + [round(y, 4) for y in ref] +
                            [round(x / ref[2], 3), dev == 0.0, round(dev, 1)])
            summary[(v, q)] = (n_in, len(data), min(devs), max(devs))
    with open("analysis/release_benchmark_nasa.csv", "w", newline="") as fo:
        csv.writer(fo).writerows(rows)
    print("NASA benchmark (Webber 1956, 20 g/cm2 Al): depths within the FLUKA/HETC-HEDS/HZETRN spread")
    for (v, q), (n, tot, lo, hi) in summary.items():
        print(f"  {v}  {q:<8} {n}/{tot} depths inside; deviation from the spread {lo:+.1f} to {hi:+.1f} %")
    for q, k in (("D (cGy)", 7), ("H (cSv)", 7)):
        rr = [r for r in rows[1:] if r[1] == q]
        print(f"  GATE/HZETRN {q}: " + "  ".join(f"{v}: {min(r[k] for r in rr if r[0] == v):.3f}-{max(r[k] for r in rr if r[0] == v):.3f}" for v in VERS))


def hzetrn():
    rows = [["version", "spectrum", "material", "point", "GATE_D_mGy", "HZETRN_D_mGy", "D_ratio",
             "GATE_H_mSv", "HZETRN_H_mSv", "H_ratio"]]
    res = {}
    for v in VERS:
        for sp, tags in SPEC.items():
            for m, l in MATS.items():
                oD = sum(oc.oltaris(l, t)[0] for t in tags); oH = sum(oc.oltaris(l, t)[1] for t in tags)
                for lab, win, i in PTS:
                    gD, gH = gate_point(v, sp, m, win); ok = oD[i] < oc.BAD
                    res[(v, sp, l, lab)] = (gD, oD[i] if ok else np.nan, gH, oH[i] if ok else np.nan)
                    rows.append([v, sp, l, lab, round(gD, 4), round(oD[i], 4) if ok else "", round(gD / oD[i], 3) if ok else "",
                                 round(gH, 4), round(oH[i], 4) if ok else "", round(gH / oH[i], 3) if ok else ""])
    with open("analysis/release_benchmark_hzetrn.csv", "w", newline="") as fo:
        csv.writer(fo).writerows(rows)
    print("\nHZETRN (OLTARIS), 20 g/cm2: GATE/HZETRN at the BFO depth (D | H) and range over all water points")
    for sp in SPEC:
        for l in MATS.values():
            line = f"  {sp:<11}{l:<3}"
            for v in VERS:
                b = res[(v, sp, l, "5 (BFO)")]
                allD = [res[(v, sp, l, p[0])][0] / res[(v, sp, l, p[0])][1] for p in PTS]
                allH = [res[(v, sp, l, p[0])][2] / res[(v, sp, l, p[0])][3] for p in PTS]
                line += (f"  {v}: BFO {b[0] / b[1]:.3f} | {b[2] / b[3]:.3f}, all D {np.nanmin(allD):.2f}-{np.nanmax(allD):.2f}"
                         f" H {np.nanmin(allH):.2f}-{np.nanmax(allH):.2f}")
            print(line)
    print("\nRatios to Al at the BFO depth (D | H): v10.0.0, v10.1.1, HZETRN")
    for sp in SPEC:
        for l in ("PE", "Pb"):
            out = []
            for v in VERS:
                x, a = res[(v, sp, l, "5 (BFO)")], res[(v, sp, "Al", "5 (BFO)")]
                out.append(f"{x[0] / a[0]:.3f} | {x[2] / a[2]:.3f}")
            x, a = res[(VERS[0], sp, l, "5 (BFO)")], res[(VERS[0], sp, "Al", "5 (BFO)")]
            out.append(f"{x[1] / a[1]:.3f} | {x[3] / a[3]:.3f}")
            print(f"  {sp:<11}{l}/Al   " + "   ".join(out))


if __name__ == "__main__":
    nasa()
    hzetrn()
    print("\nWritten: analysis/release_benchmark_nasa.csv, analysis/release_benchmark_hzetrn.csv")
