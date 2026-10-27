"""nuclear_attenuation.py - Validation of the hadronic part: probability that primary protons (150 MeV)
undergo a nuclear interaction in the shield, GATE vs measured reaction cross sections sigma_R(E).
Reference: empirical formula fitted to measurements (Carlson 1996 compilation and others),
arXiv:2602.21192 (2026): sigma = sigma_h(A)*X(E)*Y_coul(Z,E). p-p has no inelastic channel below 290 MeV.
P_int = 1 - exp(-sum_i n_i int sigma_i(E(x)) dx), with E(x) from CSDA (PSTAR) inside the slab.
GATE: 1 - (primary survival behind the material / survival in air), from runs_origin/hadronic_origin.csv."""
import csv
import numpy as np
NA = 6.02214076e23
def sig_h(A): return 47.4 * A**0.676 * (1 - 0.045 * np.cos(0.42 * np.sqrt(A)) + 1.8e-5 * A**1.63)
def X(E): return 1 - 0.3 * np.exp(-E / 3742) * np.sin(68 * E**-0.83) + 0.5 * np.exp(-E / 134)
def Y(Z, E):
    k = sum(d / (1 + np.exp(t * (zt - Z))) for d, t, zt in ((-0.3, 3, 19), (-0.6, 0.8, 28), (-0.4, 10, 70)))
    return 1 + k * np.exp(-E / 24)
def sigma_mb(Z, A, E): return sig_h(A) * X(E) * Y(Z, E)
EL = {"H": (1, 1.008), "C": (6, 12.011), "O": (8, 15.999), "Al": (13, 26.982), "Ti": (22, 47.867), "Fe": (26, 55.845), "Cu": (29, 63.546), "Pb": (82, 207.2)}
MAT = {  # material: (composition by mass, areal density g/cm2)
    "G4_POLYETHYLENE": ({"H": 0.1437, "C": 0.8563}, 2 * 0.94), "G4_PLEXIGLASS": ({"H": 0.0805, "C": 0.5999, "O": 0.3196}, 4 * 1.19),
    "G4_Al": ({"Al": 1}, 2 * 2.699), "G4_Ti": ({"Ti": 1}, 2 * 4.54), "G4_Fe": ({"Fe": 1}, 1.1 * 7.874),
    "G4_Cu": ({"Cu": 1}, 1 * 8.96), "G4_Pb": ({"Pb": 1}, 0.4 * 11.35)}
def pst(m):
    a = np.array([[float(x) for x in l.split()] for l in open(f"pstar/{m}.txt") if len(l.split()) == 7]); return a[:, 0], a[:, 4]
g = {r["material"]: float(r["primary_survival"]) for r in csv.DictReader(open("runs_origin/hadronic_origin.csv"))}
nt = 250000 * 0.5625
print(f"{'Material':<16}{'g/cm2':>7}{'E_out':>7}{'sigR(150) mb':>13}{'P GATE %':>12}{'P expect %':>11}{'GATE/exp':>12}")
rows = []
for m, (comp, t) in MAT.items():
    E, R = pst(m); R0 = np.exp(np.interp(np.log(150), np.log(E), np.log(R)))
    x = np.linspace(0, t, 2001); Ex = np.exp(np.interp(np.log(R0 - x), np.log(R), np.log(E)))
    mu = sum(w * NA / EL[el][1] * (0 if el == "H" else sigma_mb(*EL[el], Ex) * 1e-27) for el, w in comp.items())   # cm²/g
    P = 1 - np.exp(-(getattr(np, "trapezoid", None) or np.trapz)(mu, x))
    pg = 1 - g[m] / g["G4_AIR"]; se = np.sqrt(pg * (1 - pg) / nt)
    s150 = " / ".join(f"{sigma_mb(*EL[el], 150.0):.0f}" for el in comp if el != "H")
    print(f"{m:<16}{t:>7.2f}{Ex[-1]:>7.1f}{s150:>13}{100*pg:>8.2f}±{100*se:.2f}{100*P:>11.2f}{pg/P:>12.3f}")
    rows.append([m, t, Ex[-1], s150, pg, se, P, pg / P])
with open("analysis/nuclear_attenuation.csv", "w", newline="") as fo:
    w = csv.writer(fo); w.writerow(["material", "areal_g_cm2", "E_out_MeV", "sigmaR_150MeV_mb", "P_int_GATE", "P_int_GATE_se", "P_int_expected", "ratio"]); w.writerows(rows)
r = np.array([x[7] for x in rows]); print(f"\nGATE/expected: mean {r.mean():.3f}, range {r.min():.3f}-{r.max():.3f}")
print("The formula uncertainty is ~5 % (deviation from data); GATE also counts elastic scattering out of the plane (small).")
