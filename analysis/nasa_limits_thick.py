"""nasa_limits_thick.py - Minimum areal density for BFO < 250 mGy-Eq (NASA 30-day limit) from the response matrix
at 5, 10, 20, 30 and 40 g/cm2 (logarithmic interpolation), for the 5 spectra and 9 materials. BFO = D(4-6 g/cm2) x RBE 1.5.
Needs response_matrix.npz with the 30/40 g/cm2 shields (response_matrix.py build_thick).
Writes analysis/nasa_limits_thick.csv."""
import sys, csv
import numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "analysis")
import response_matrix as rm
from spectra import SPECTRA
T = [5, 10, 20, 30, 40]; LIM = 250.0
SP = ["kold1956", "webber1956", "kold1972", "king1972", "kold1989"]
def tmin(b):
    b = np.array(b); y = np.log(b); t = np.array(T, float)
    if b[0] <= LIM: return "<=5"
    if b[-1] > LIM: return ">40"
    i = np.where(b <= LIM)[0][0]
    return f"{np.interp(np.log(LIM), [y[i], y[i-1]], [t[i], t[i-1]]):.1f}"


def main():
    rows = []
    print(f"{'':<14}" + "".join(f"{s:>12}" for s in SP))
    res = {}
    for m in rm.MATS:
        line = f"{m:<14}"
        for sp in SP:
            b = [1.5 * 10 * rm.fold(SPECTRA[sp], m, t, "D", "BFO") for t in T]
            tm = tmin(b); res[(m, sp)] = (b, tm); line += f"{tm:>12}"
            rows.append([sp, m] + [round(x, 1) for x in b] + [tm])
        print(line)
    print("\nBFO (mGy-Eq) for Al and PE at 20/30/40 g/cm2:")
    for sp in SP:
        for m in ("Al", "POLYETHYLENE"):
            b = res[(m, sp)][0]; print(f"  {sp:<11}{m:<13} " + " / ".join(f"{x:.0f}" for x in b[2:]))
    with open("analysis/nasa_limits_thick.csv", "w", newline="") as fo:
        w = csv.writer(fo); w.writerow(["spectrum", "material"] + [f"BFO_{t}_mGyEq" for t in T] + ["min_areal_g_cm2"]); w.writerows(rows)
    print("Written: analysis/nasa_limits_thick.csv")


if __name__ == "__main__":
    main()
