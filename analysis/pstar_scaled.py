"""pstar_scaled.py - PSTAR-format stopping-power and range table for a material that PSTAR does not list
(here lithium hydride), by scaling the PSTAR table of a chemically similar reference material with the ratio of
Bethe stopping powers:

    S_m(E) = S_ref,PSTAR(E) * S_Bethe(E; Z/A_m, I_m) / S_Bethe(E; Z/A_ref, I_ref)

The Bethe formula is used without shell, Barkas, Bloch and density-effect corrections; these mostly cancel in the
ratio of two low-Z materials. The ratio is evaluated at E >= 1 MeV and held constant below (the range below 1 MeV is
about 0.002 g/cm2 and does not matter for shields of 5-40 g/cm2). The CSDA range is the integral of 1/S.
The method is checked by building water from polyethylene and comparing with the PSTAR table of water.

Lithium hydride (G4_LITHIUM_HYDRIDE): density 0.82 g/cm3, Z/A = 4/7.949 = 0.50321, I = 36.5 eV (Geant4/NIST).
Usage (from the repository root):  python analysis/pstar_scaled.py   -> pstar/G4_LITHIUM_HYDRIDE.txt
"""
import numpy as np

ME, MP, K = 0.51099895, 938.27209, 0.307075     # MeV, MeV, MeV cm2/mol

MAT = {  # Z/A, I (eV)
    "POLYETHYLENE": (0.57034, 57.4),
    "WATER": (0.55509, 78.0),          # Geant4 value; PSTAR uses 75 eV
    "WATER_PSTAR": (0.55509, 75.0),
    "LITHIUM_HYDRIDE": (4.0 / (6.941 + 2 * 0 + 1.00794), 36.5),
}


def bethe(E, za, I_eV):
    g = 1 + E / MP; b2 = 1 - 1 / g ** 2
    tmax = 2 * ME * b2 * g ** 2 / (1 + 2 * g * ME / MP + (ME / MP) ** 2)
    return K * za / b2 * (0.5 * np.log(2 * ME * b2 * g ** 2 * tmax / (I_eV * 1e-6) ** 2) - b2)


def read(m):
    a = np.array([[float(x) for x in l.split()] for l in open(f"pstar/G4_{m}.txt") if len(l.split()) == 7])
    return a


def scaled(ref, target, ref_bethe=None):
    a = read(ref); E = a[:, 0]
    rz, ri = MAT[ref_bethe or ref]; tz, ti = MAT[target]
    Ec = np.maximum(E, 1.0)
    ratio = bethe(Ec, tz, ti) / bethe(Ec, rz, ri)
    S_el, S_nuc = a[:, 1] * ratio, a[:, 2] * (tz / rz)
    S = S_el + S_nuc
    # CSDA range: reference range at the first point scaled, then trapezoidal integration of 1/S in E
    R = np.empty_like(E); R[0] = a[0, 4] * (a[0, 3] / S[0])
    R[1:] = R[0] + np.cumsum(0.5 * (1 / S[1:] + 1 / S[:-1]) * np.diff(E))
    proj = R * a[:, 6]                     # projected range with the detour factor of the reference
    return np.column_stack([E, S_el, S_nuc, S, R, proj, a[:, 6]])


def write(m, t):
    with open(f"pstar/G4_{m}.txt", "w") as fo:
        for r in t:
            fo.write("  " + "  ".join(f"{x:.3E}" for x in r[:6]) + f"  {r[6]:.4f}\n")


def main():
    # check: water from polyethylene (with the PSTAR I-value of water) against the PSTAR table of water
    w = scaled("POLYETHYLENE", "WATER_PSTAR"); p = read("WATER")
    for e in (50, 100, 150, 200, 300):
        i = np.argmin(abs(p[:, 0] - e))
        print(f"  water from PE, {p[i, 0]:5.0f} MeV: S {w[i, 3] / p[i, 3]:.4f}, CSDA range {w[i, 4] / p[i, 4]:.4f} (scaled / PSTAR)")
    t = scaled("POLYETHYLENE", "LITHIUM_HYDRIDE")
    write("LITHIUM_HYDRIDE", t)
    pe = read("POLYETHYLENE")
    for e in (50, 100, 150, 200):
        i = np.argmin(abs(t[:, 0] - e))
        print(f"  LiH {t[i, 0]:5.0f} MeV: S/rho {t[i, 3]:.3f} MeV cm2/g (PE {pe[i, 3]:.3f}), CSDA range {t[i, 4]:.2f} g/cm2 (PE {pe[i, 4]:.2f})")
    print("Written: pstar/G4_LITHIUM_HYDRIDE.txt")


if __name__ == "__main__":
    main()
