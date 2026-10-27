"""check_spectrum_sampling.py - Check of the opengate spectrum_histogram sampling (~1 min).
Proton source in vacuum, 10-60 MeV group of King 1972 with 100 logarithmic bins:
  (a) old way: weights = fluence inside the bin   (b) corrected: weights = fluence / bin width
Compares the fraction of protons above 49.7 MeV (threshold of 5 g/cm2 Pb) with the exact value.
Also: protons that crossed in the old Pb 5 g/cm2 run (10-60 group) against the two predictions."""
import os, sys, glob, json
import numpy as np, uproot
sys.path.insert(0, os.getcwd())
import opengate as gate
from spectra import spectrum_bins, SPECTRA

J = SPECTRA["king1972"]; TH = 49.7
true = (J(TH) - J(60)) / (J(10) - J(60))
edges, w, _ = spectrum_bins("king1972", 10, 60, 100)

def run(tag, weights):
    MeV, cm, m = gate.g4_units.MeV, gate.g4_units.cm, gate.g4_units.m
    sim = gate.Simulation(); sim.output_dir = f"tests_archive/sampling_{tag}"; sim.random_seed = 7
    sim.world.size = [1 * m] * 3; sim.world.material = "G4_Galactic"
    pl = sim.add_volume("Box", name="plane"); pl.size = [10 * cm, 10 * cm, 1 * cm]; pl.material = "G4_Galactic"
    src = sim.add_source("GenericSource", "p"); src.particle = "proton"; src.n = 200000
    src.energy.type = "spectrum_histogram"; src.energy.spectrum_energy_bin_edges = list(edges * MeV)
    src.energy.spectrum_weights = list(weights / weights.sum()); src.energy.spectrum_histogram_interpolation = None
    src.position.type = "point"; src.position.translation = [0, 0, -20 * cm]
    src.direction.type = "momentum"; src.direction.momentum = [0, 0, 1]
    ph = sim.add_actor("PhaseSpaceActor", "ph"); ph.attached_to = "plane"; ph.attributes = ["KineticEnergy"]
    ph.output_filename = "ph.root"
    sim.run(start_new_process=True)
    with uproot.open(f"{sim.output_dir}/ph.root") as f:
        E = f[f.keys()[0]].arrays(["KineticEnergy"], library="np")["KineticEnergy"]
    return np.mean(E > TH), np.mean(E)

def main():
    Eg = np.linspace(10, 60, 200001); pdf = -np.gradient(J(Eg), Eg)
    print(f"Exact fraction >{TH} MeV: {100*true:.2f} %   | exact mean E: {np.sum(pdf*Eg)/np.sum(pdf):.2f} MeV")
    for tag, ww in (("old_flux_in_bin", w), ("fixed_density", w / np.diff(edges))):
        fr, me = run(tag, ww)
        print(f"  {tag:<16}: >{TH} MeV {100*fr:6.2f} %  (ratio to exact {fr/true:.2f}) | mean E {me:.2f} MeV")

    d = "spe_runs/materials_king1972/Pb_5/g10-60"
    if glob.glob(f"{d}/phsp_*.root"):
        s = json.load(open(f"{d}/setup.json"))
        with uproot.open(glob.glob(f"{d}/phsp_*.root")[0]) as f:
            p = f[f.keys()[0]].arrays(["ParticleName", "Direction_Z", "KineticEnergy"], library="np")
        n = np.sum((p["ParticleName"].astype(str) == "proton") & (p["Direction_Z"] > 0) & (p["KineticEnergy"] > 1)) / (s["n"] * (30 / 80) ** 2)
        tilt = np.sum((pdf * Eg)[Eg > TH]) / np.sum(pdf * Eg)
        print(f"\nOld run Pb 5 g/cm2, group 10-60: protons into water per source proton {n:.4f}"
              f"  | prediction, exact spectrum {0.97*true:.4f} | prediction, tilted {0.97*tilt:.4f}")


if __name__ == "__main__":
    main()
