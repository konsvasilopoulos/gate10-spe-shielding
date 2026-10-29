"""oltaris_compare.py - Comparison of GATE 10 with HZETRN (OLTARIS, TARIS 5.01) for the same slab geometry.

OLTARIS slabs (outermost layer first): shield of 5 + 5 + 10 g/cm2 (Al 2.699, PE 0.94, Pb 11.35 g/cm3)
followed by water 0.3 + 4.7 + 5 + 20 g/cm2. Dose (mGy) and dose equivalent (ICRP 60, mSv) are reported at
every interface: behind 5, 10 and 20 g/cm2 of shield (= water entrance), and at 0.3, 5, 10 and 30 g/cm2 of water.
The OLTARIS tables are in oltaris/<material>_<spectrum>/ (li_dose_database_mGy.dat, li_doseq_database_mSv.dat;
li_flux.dat for the Tylka spectrum).

Spectra (OLTARIS input):
  King72   user-defined, exponential in energy: J0 = 7.9e9, E0 = 26.5 MeV, E1 = 30 MeV
  Webber56 user-defined, exponential in rigidity (MV): J0 = 1e9, R0 = 100 MV, R1 = 239.1 MV
  GLE24/25 user-defined Band function = Koldobskiy et al. (2021) with R2 = infinity. OLTARIS defines J0 with the
           rigidity in GV; the runs were made with J0 given for R in MV, so the results are divided by 1000^gamma1
           (exact, the dose is linear in J0). Equivalent OLTARIS inputs: GLE24 J0 = 1.492e9, R0 = 152 MV,
           g1 = 2.02, g2 = 11.96; GLE25 J0 = 9.920e9, R0 = 39 MV, g1 = 1.53, g2 = 5.07.
           PE_GLE25 was run with J0 = 3.8595 instead of 3.8595e14 (typing error); its results are scaled by 1e14
           (checked: the dose at the slab surface differs from the Al and Pb runs by exactly this factor).
  Tylka89  historical "Sum of October 1989 Tylka Band fits" (protons, 0.01-2500 MeV); the boundary spectrum is
           read from li_flux.dat and folded with the GATE response matrix.
GATE values: response matrix (analysis/response_matrix.npz) folded with the same spectra, averaged over
0-0.1 g/cm2 (entrance), 0.2-0.4, 4.5-5.5, 9.5-10.5 and 29.5-30 g/cm2 of water.
Known OLTARIS problem: with the user-defined Lead material the dose at the Pb/water interface and 0.3 g/cm2 behind
it is unphysical (~1e88-1e95 mGy) for every spectrum; these two points are excluded.

Usage (from the repository root):  python analysis/oltaris_compare.py
  -> analysis/oltaris_compare.csv (all points), analysis/oltaris_ratios.csv (PE/Al, Pb/Al),
     analysis/oltaris_tylka.csv (Tylka vs Koldobskiy, October 1989)
"""
import os, sys, csv
import numpy as np
sys.path.insert(0, os.getcwd()); sys.path.insert(0, "analysis")
import response_matrix as rm
from spectra import SPECTRA, koldobskiy_gle_integral, koldobskiy_sum

ODIR = "oltaris"
MATS = {"Al": "Al", "POLYETHYLENE": "PE", "Pb": "Pb"}
# (label, shield areal density, GATE window in water, index in the OLTARIS depth array)
POINTS = [("behind 5 g/cm2", 5, (0, 0.1), 1), ("behind 10 g/cm2", 10, (0, 0.1), 2), ("behind 20 g/cm2", 20, (0, 0.1), 3),
          ("20 + water 0.3", 20, (0.2, 0.4), 4), ("20 + water 5 (BFO)", 20, (4.5, 5.5), 5),
          ("20 + water 10", 20, (9.5, 10.5), 6), ("20 + water 30", 20, (29.5, 30), 7)]
SCALE = {"GLE24": 1000 ** 2.02, "GLE25": 1000 ** 1.53}
EXTRA = {("PE", "GLE25"): 1e-14}
BAD = 1e30


def read_table(path):
    """OLTARIS ASCII table -> (depth or energy axis, values[, particle names])."""
    tok = " ".join(l for l in open(path) if not l.startswith("C") and l.strip()).split()
    n = int(tok[1]); axis = np.array(tok[2:2 + n], float); rest = tok[2 + n:]
    if rest[0] == "1" and len(rest) > 2 + int(rest[1]) and not _isnum(rest[2]):
        npart = int(rest[1]); parts = rest[2:2 + npart]; vals = np.array(rest[2 + npart + 4:], float).reshape(npart, n)
        return axis, vals, parts
    return axis, np.array(rest[2:], float)


def _isnum(s):
    try: float(s); return True
    except ValueError: return False


def oltaris(mat, tag):
    d = f"{ODIR}/{mat}_{tag}"
    _, D = read_table(f"{d}/li_dose_database_mGy.dat"); _, H = read_table(f"{d}/li_doseq_database_mSv.dat")
    s = SCALE.get(tag, 1.0) * EXTRA.get((mat, tag), 1.0)
    return D / s, H / s


def tylka_integral():
    E, phi, parts = read_table(f"{ODIR}/Al_Tylka89/li_flux.dat")
    p = phi[parts.index("proton")]
    Eg = np.geomspace(E[0], E[-1], 20001); pg = np.exp(np.interp(np.log(Eg), np.log(E), np.log(p)))
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (pg[1:] + pg[:-1]) * np.diff(Eg))])
    return lambda x: np.interp(np.log(np.clip(np.asarray(x, float), E[0], E[-1])), np.log(Eg), cum[-1] - cum)


def gate_point(F, m, t, q, win):
    rm.WIN["_pt"] = win
    return 10 * rm.fold(F, m, t, q, "_pt")          # mGy / mSv


SPEC = {"King72": (SPECTRA["king1972"], ["King72"]),
        "Webber56": (SPECTRA["webber1956"], ["Webber56"]),
        "Kold72": (SPECTRA["kold1972"], ["GLE24", "GLE25"]),        # Koldobskiy Aug 1972 = GLE 24 + 25
        "Tylka89": (None, ["Tylka89"])}


def main():
    SPEC["Tylka89"] = (tylka_integral(), ["Tylka89"])
    rows, res = [], {}
    for sp, (F, tags) in SPEC.items():
        for m, l in MATS.items():
            oD = sum(oltaris(l, t)[0] for t in tags); oH = sum(oltaris(l, t)[1] for t in tags)
            for lab, t, win, i in POINTS:
                gD, gH = gate_point(F, m, t, "D", win), gate_point(F, m, t, "H", win)
                ok = oD[i] < BAD
                res[(sp, l, lab)] = (gD, oD[i] if ok else np.nan, gH, oH[i] if ok else np.nan)
                rows.append([sp, l, lab, round(gD, 4), round(oD[i], 4) if ok else "", round(gD / oD[i], 3) if ok else "",
                             round(gH, 4), round(oH[i], 4) if ok else "", round(gH / oH[i], 3) if ok else "",
                             "" if ok else "OLTARIS numerical error at the Pb/water interface"])
    with open("analysis/oltaris_compare.csv", "w", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(["spectrum", "material", "point", "GATE_D_mGy", "OLTARIS_D_mGy", "D_ratio", "GATE_H_mSv",
                    "OLTARIS_H_mSv", "H_ratio", "note"]); w.writerows(rows)
    rr = []
    print(f"{'':<9}{'point':<20}{'PE/Al D GATE':>13}{'OLT':>7}{'Pb/Al D GATE':>14}{'OLT':>7}{'PE/Al H GATE':>14}{'OLT':>7}{'Pb/Al H GATE':>14}{'OLT':>7}")
    for sp in SPEC:
        for lab, *_ in POINTS:
            al = res[(sp, "Al", lab)]; out = []
            for k in (0, 2):
                for l in ("PE", "Pb"):
                    x = res[(sp, l, lab)]
                    out += [x[k] / al[k], x[k + 1] / al[k + 1]]
            # order: PE D gate, PE D olt, Pb D gate, Pb D olt, PE H gate, PE H olt, Pb H gate, Pb H olt
            rr.append([sp, lab] + [round(v, 3) if np.isfinite(v) else "" for v in out])
            print(f"{sp:<9}{lab:<20}" + "".join(f"{v:>13.3f}{'':>1}" if j % 2 == 0 else f"{v:>7.3f}" for j, v in enumerate(out)))
    with open("analysis/oltaris_ratios.csv", "w", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(["spectrum", "point", "PE_Al_D_GATE", "PE_Al_D_OLTARIS", "Pb_Al_D_GATE", "Pb_Al_D_OLTARIS",
                    "PE_Al_H_GATE", "PE_Al_H_OLTARIS", "Pb_Al_H_GATE", "Pb_Al_H_OLTARIS"]); w.writerows(rr)
    print("\nGATE/OLTARIS, absolute dose (range over all valid points):")
    for sp in SPEC:
        print(f"  {sp:<9}" + "  ".join(
            f"{l}: D {np.nanmin([res[(sp, l, p[0])][0] / res[(sp, l, p[0])][1] for p in POINTS]):.2f}-"
            f"{np.nanmax([res[(sp, l, p[0])][0] / res[(sp, l, p[0])][1] for p in POINTS]):.2f}, H "
            f"{np.nanmin([res[(sp, l, p[0])][2] / res[(sp, l, p[0])][3] for p in POINTS]):.2f}-"
            f"{np.nanmax([res[(sp, l, p[0])][2] / res[(sp, l, p[0])][3] for p in POINTS]):.2f}" for l in ("Al", "PE", "Pb")))

    # Tylka vs Koldobskiy (October 1989 = GLE 43-45) through the GATE response matrix, 5-40 g/cm2
    import nasa_limits_thick as nl
    Ft = SPEC["Tylka89"][0]; Ko = koldobskiy_sum([43, 44, 45])
    trows = [["quantity", "Tylka_OLTARIS", "Koldobskiy_GLE43-45", "ratio"]]
    for e in (30, 60, 100, 200, 500, 1000, 2000):
        a, b = float(Ft(e)), float(Ko(e)); trows.append([f"Phi(>{e} MeV) cm^-2", f"{a:.4g}", f"{b:.4g}", round(a / b, 3)])
    for m in ("Al", "POLYETHYLENE", "Pb"):
        ba = [1.5 * 10 * rm.fold(Ft, m, t, "D", "BFO") for t in nl.T]; bb = [1.5 * 10 * rm.fold(Ko, m, t, "D", "BFO") for t in nl.T]
        trows.append([f"BFO behind 20 g/cm2 {m} (mGy-Eq)", round(ba[2], 1), round(bb[2], 1), round(ba[2] / bb[2], 3)])
        trows.append([f"min areal density {m} for 250 mGy-Eq (g/cm2)", nl.tmin(ba), nl.tmin(bb), ""])
    for m in ("POLYETHYLENE", "Pb"):
        for win in ("BFO", "mean"):
            a = rm.fold(Ft, m, 20, "D", win) / rm.fold(Ft, "Al", 20, "D", win)
            b = rm.fold(Ko, m, 20, "D", win) / rm.fold(Ko, "Al", 20, "D", win)
            trows.append([f"{m}/Al D {win}, 20 g/cm2", round(a, 3), round(b, 3), round(a / b, 3)])
    with open("analysis/oltaris_tylka.csv", "w", newline="") as fo:
        csv.writer(fo).writerows(trows)
    print("\nTylka (OLTARIS) vs Koldobskiy GLE 43-45, October 1989:")
    for r in trows[1:]: print(f"  {r[0]:<48}{str(r[1]):>10}{str(r[2]):>10}{str(r[3]):>8}")
    print("Written: analysis/oltaris_compare.csv, analysis/oltaris_ratios.csv, analysis/oltaris_tylka.csv")


if __name__ == "__main__":
    main()
