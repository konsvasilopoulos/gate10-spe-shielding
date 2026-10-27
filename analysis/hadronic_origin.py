"""
hadronic_origin.py - Origin of the hadronic contribution to the dose in the first voxel (0-1 mm of water),
150 MeV, validation geometry.

For each material:
  - hadronic contribution = Edep(QGSP_BIC_EMZ, this run) - Edep(EM-only), MeV per incident proton
  - charged secondaries (p, d, t, 3He, alpha, ions) entering the water, by creator process
    (elastic / inelastic) and creation place (shield or air), and an estimate of the energy
    they deposit in the first mm (CSDA, PSTAR water, with range scaling for d, t, 3He, alpha).
Usage:  python analysis/hadronic_origin.py
"""
import os, sys, glob, json
import numpy as np
import uproot

MATS = ["G4_AIR", "G4_POLYETHYLENE", "G4_PLEXIGLASS", "G4_Al", "G4_Ti", "G4_Fe", "G4_Cu", "G4_Pb"]
EM_ONLY = {"G4_AIR": 0.5409, "G4_POLYETHYLENE": 0.5725, "G4_PLEXIGLASS": 0.6307, "G4_Al": 0.6152,
           "G4_Ti": 0.6748, "G4_Fe": 0.6609, "G4_Cu": 0.6637, "G4_Pb": 0.5771}   # MeV per incident proton, cut 0.1 mm (EM-only)
FRAC = (30.0 / 40.0) ** 2
VOXEL = 0.1                                   # g/cm2 (1 mm of water)

pw = np.array([[float(x) for x in l.split()] for l in open("pstar/G4_WATER.txt") if len(l.split()) == 7])
Ew, Rw = pw[:, 0], pw[:, 4]
R_p = lambda E: np.exp(np.interp(np.log(np.maximum(E, Ew[0])), np.log(Ew), np.log(Rw)))
E_p = lambda R: np.exp(np.interp(np.log(np.maximum(R, Rw[0])), np.log(Rw), np.log(Ew)))
# range scaling at equal velocity: R_x(E) = (m_x/m_p)/z_x^2 * R_p(E*m_p/m_x)
SCALE = {"proton": (1, 1), "deuteron": (2, 1), "triton": (3, 1), "He3": (3, 2), "alpha": (4, 2)}


def dep_first_voxel(name, E, cosz):
    t = VOXEL / np.clip(cosz, 0.05, 1)             # path length in the voxel
    if name not in SCALE:                          # heavier ions: very short range -> all energy deposited
        return E
    A, z = SCALE[name]
    k = A / z ** 2
    R = k * R_p(E / A)
    Eleft = np.where(R > t, A * E_p((R - t) / k), 0.0)
    return E - Eleft


def read_mhd_first(path):
    hdr = {l.split("=")[0].strip(): l.split("=", 1)[1].strip() for l in open(path) if "=" in l}
    dt = np.dtype({"MET_DOUBLE": "f8", "MET_FLOAT": "f4"}[hdr["ElementType"]]).newbyteorder("<")
    return float(np.fromfile(os.path.join(os.path.dirname(path), hdr["ElementDataFile"]), dtype=dt)[0])


print(f"{'Material':<16}{'hadr. D':>9}{'n sec/p':>10}{'elast%':>8}{'+sec shld':>12}{'(elast)':>9}"
      f"{'p surviv.':>11}{'-p loss':>11}{'=water':>8}")
print(" " * 16 + "(MeV/p)" + " " * 18 + "(MeV/p, 1st mm)" + " " * 10 + "(MeV/p)   (D - sec + loss)")
rows = []
for m in MATS:
    d = f"runs_origin/output_{m}"
    if not os.path.exists(f"{d}/setup.json"):
        print(f"{m:<16} - missing"); continue
    s = json.load(open(f"{d}/setup.json")); nt = s["n"] * FRAC
    edep = read_mhd_first(glob.glob(f"{d}/dose_{m}_edep.mhd")[0]) / nt
    had = edep - EM_ONLY[m]
    with uproot.open(glob.glob(f"{d}/phsp_*.root")[0]) as f:
        t = f[f.keys()[0]]
        p = t.arrays(["ParticleName", "KineticEnergy", "Direction_Z", "TrackCreatorProcess",
                      "TrackVertexPosition_Z"], library="np")
    name = p["ParticleName"].astype(str); proc = p["TrackCreatorProcess"].astype(str)
    fwd = p["Direction_Z"] > 0
    primary = np.isin(proc, ["none", "primary", ""])
    charged_sec = fwd & ~primary & ~np.isin(name, ["e-", "e+", "gamma", "neutron", "nu_e", "anti_nu_e"])
    elastic = np.char.find(np.char.lower(proc.astype("U")), "elastic") >= 0
    elastic &= np.char.find(np.char.lower(proc.astype("U")), "inelastic") < 0
    in_shield = p["TrackVertexPosition_Z"] < 0          # the shield is at z < 0
    dep = np.zeros(len(name))
    for nm in np.unique(name[charged_sec]):
        sel = charged_sec & (name == nm)
        dep[sel] = dep_first_voxel(nm, p["KineticEnergy"][sel], p["Direction_Z"][sel])
    n_sec = charged_sec.sum() / nt
    f_el = 100 * (charged_sec & elastic).sum() / max(charged_sec.sum(), 1)
    f_sh = 100 * (charged_sec & in_shield).sum() / max(charged_sec.sum(), 1)
    e_sh = dep[charged_sec & in_shield].sum() / nt
    e_el = dep[charged_sec & in_shield & elastic].sum() / nt
    e_in = e_sh - e_el
    # primary protons reaching the water: nuclear reactions in the shield remove them
    surv = (fwd & primary & (name == "proton")).sum() / nt
    loss = (1.0 - surv) * EM_ONLY[m]              # dose lost because primaries are missing (EM-only p/p ~ 1)
    water = had - e_sh + loss
    by = {}
    for nm in ["proton", "deuteron", "triton", "He3", "alpha"]:
        by[nm] = dep[charged_sec & in_shield & (name == nm)].sum() / nt
    by["ions Z>=3"] = dep[charged_sec & in_shield & ~np.isin(name, list(SCALE))].sum() / nt
    rows.append((m, had, n_sec, f_el, f_sh, e_sh, e_el, e_in, surv, loss, water) + tuple(by.values()))
    print(f"{m:<16}{had:>9.4f}{n_sec:>10.4f}{f_el:>8.1f}{e_sh:>12.4f}{e_el:>9.4f}{surv:>11.4f}{loss:>11.4f}{water:>8.4f}")
    print(" " * 18 + "shield secondaries by species: " + " | ".join(f"{k} {v:.4f}" for k, v in by.items()))
    if m == "G4_AIR":
        print("   creator processes (charged secondaries, forward):",
              dict(zip(*np.unique(proc[charged_sec], return_counts=True))))

print("\nBalance: hadronic D = (+) secondaries from the shield + (+) reactions in the water - (-) primaries lost in the shield.")
print("The '=water' column should be about the same for all materials (~ the air value) if the balance closes.")
import csv
with open("runs_origin/hadronic_origin.csv", "w", newline="") as fo:
    w = csv.writer(fo)
    w.writerow(["material", "hadronic_increment_MeV_per_p", "charged_secondaries_per_p", "elastic_pct",
                "from_shield_pct", "Edep_shield_secondaries_MeV_per_p", "Edep_elastic", "Edep_inelastic",
                "primary_survival", "primary_loss_MeV_per_p", "water_term", "Edep_p", "Edep_d", "Edep_t",
                "Edep_He3", "Edep_alpha", "Edep_ionsZ3"])
    w.writerows(rows)
