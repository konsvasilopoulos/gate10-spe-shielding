import os
import sys
import json
import argparse
import numpy as np
import opengate as gate
from spectra import spectrum_bins

# =============================================================================
# run_spectrum.py - SPE proton spectrum on a shield slab followed by a water phantom.
#
# Same validated geometry as run_material.py, with three differences:
#   1. The source energy follows the event spectrum (e.g. Webber 1956) within one
#      energy group [--emin, --emax]. Each group is a separate run, so that more
#      events can be spent where they are needed (high energies).
#   2. The shield thickness is given as areal density in g/cm2 (--areal).
#   3. The water phantom thickness is variable (--phantom-cm, default 30).
#
# Normalisation (written to setup.json): each run represents the event fluence
# inside the group, Phi_group (cm^-2). Event dose = run dose x Phi_group x A_source / N.
#
# Example (NASA TP-2009-215560 benchmark: 20 g/cm2 Al + 30 g/cm2 water):
#   python run_spectrum.py G4_Al --areal 20 --emin 150 --emax 2000 --n 250000 \
#          --out spe_runs/webber1956_Al20/g150-2000
# All groups and materials: bash run_materials_spe.sh
# =============================================================================

DENSITY = {  # g/cm3, Geant4/NIST values
    "G4_AIR": 0.00120479, "G4_POLYETHYLENE": 0.94, "G4_PLEXIGLASS": 1.19, "G4_Al": 2.699,
    "G4_Ti": 4.54, "G4_Fe": 7.874, "G4_Cu": 8.96, "G4_Pb": 11.35, "G4_WATER": 1.0,
}

SOURCE_XY, SLAB_XY, DETECTOR_XY = 40, 60, 30   # cm (defaults; SOURCE_XY and SLAB_XY can be changed with --source-xy/--slab-xy)
PLANE_THICKNESS = 0.1                          # cm (1 mm)
SPLIT = {"proton": "proton", "neutron": "neutron", "photon_electron": ["gamma", "e-", "e+"]}

# --- Dose equivalent -----------------------------------------------------------------------
# One LETActor (dose-averaged, LET in water) per (particle, kinetic-energy bin), writing the
# numerator sum(Edep*L) and the denominator sum(Edep). Then H = D + sum_bin D_bin*(Q(Lbar_bin) - 1).
# Bins are dense around the knees of the ICRP-60 Q(L) (L = 10 and 100 keV/um); between them
# Q is linear in L, so the dose-averaged LET of the bin gives the exact result.
# Above the last edge of each particle L < 10 keV/um, so Q = 1 (no actor needed).
LET_BINS = {  # MeV, kinetic energy of the particle (not per nucleon)
    "proton":   [0, 1, 3, 4, 5, 8],
    "deuteron": [0, 2, 6, 8, 10, 16],
    "triton":   [0, 3, 9, 12, 15, 24],
    "He3":      [0, 1, 2, 3, 5, 15, 40, 60, 80, 150],
    "alpha":    [0, 0.5, 1, 2, 4, 6, 20, 50, 80, 110, 200],
}
# Everything except e+-, gamma, n, p, d, t, He3, alpha (heavier fragments and recoils) -> one actor
LET_REST_EXCLUDE = ["e-", "e+", "gamma", "neutron", "proton", "deuteron", "triton", "He3", "alpha"]


def _let_actor(sim, name, nz, dz, cm, mm):
    a = sim.add_actor("LETActor", name)
    a.attached_to = "Water_Phantom"
    a.size = [1, 1, nz]
    a.spacing = [DETECTOR_XY * cm, DETECTOR_XY * cm, dz * mm]
    a.averaging_method = "dose_average"
    a.score_in = "G4_WATER"
    a.output_filename = f"{name}.mhd"
    a.user_output.let.set_write_to_disk(True, item=0)   # numerator   sum Edep*L
    a.user_output.let.set_write_to_disk(True, item=1)   # denominator sum Edep
    return a


def add_let_actors(sim, mat, nz, dz, cm, mm, MeV):
    for part, edges in LET_BINS.items():
        for lo, hi in zip(edges[:-1], edges[1:]):
            tag = f"let_{part}_{lo:g}-{hi:g}MeV"
            a = _let_actor(sim, tag, nz, dz, cm, mm)
            fp = sim.add_filter("ParticleFilter", f"fp_{tag}")
            fp.particle = part
            fe = sim.add_filter("KineticEnergyFilter", f"fe_{tag}")
            fe.energy_min = lo * MeV
            fe.energy_max = hi * MeV
            a.filters.append(fp)
            a.filters.append(fe)
            a.filters_boolean_operator = "and"
    a = _let_actor(sim, "let_rest_all", nz, dz, cm, mm)
    for nm in LET_REST_EXCLUDE:
        f = sim.add_filter("ParticleFilter", f"fx_rest_{nm}")
        f.particle = nm
        f.policy = "reject"
        a.filters.append(f)
    a.filters_boolean_operator = "and"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("material")
    ap.add_argument("--areal", type=float, required=True, help="shield areal density (g/cm2)")
    ap.add_argument("--spectrum", default="webber1956")
    ap.add_argument("--emin", type=float, required=True, help="lower edge of the energy group (MeV)")
    ap.add_argument("--emax", type=float, required=True, help="upper edge of the energy group (MeV)")
    ap.add_argument("--nbins", type=int, default=100, help="number of spectrum bins inside the group")
    ap.add_argument("--phantom-cm", type=float, default=30.0)
    ap.add_argument("--dz", type=float, default=1.0, help="voxel size along z (mm)")
    ap.add_argument("--n", type=int, default=250000)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--cut", type=float, default=0.1, help="production cut (mm)")
    ap.add_argument("--physics", default="QGSP_BIC_EMZ")
    ap.add_argument("--split", action="store_true")
    ap.add_argument("--dose-eq", action="store_true",
                    help="LET actors per particle and energy, for dose equivalent (analyze_dose_eq.py)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--no-phsp", action="store_true", help="no PhaseSpaceActor (e.g. for the response matrix)")
    ap.add_argument("--particle", default="proton", choices=["proton", "alpha"],
                    help="source particle. For alpha, --emin/--emax and the spectrum are per nucleon (MeV/n); the group "
                         "fluence is the proton fluence at the same MeV/n (He/H = 1, scaled afterwards)")
    ap.add_argument("--theta-deg", type=float, default=0.0,
                    help="incidence angle to the z axis (degrees). The source is shifted and widened along x so that the "
                         "central 30x30 cm2 column receives a uniform fluence (isotropic-incidence test)")
    ap.add_argument("--source-xy", type=float, default=SOURCE_XY, help="source side length (cm)")
    ap.add_argument("--slab-xy", type=float, default=SLAB_XY, help="shield slab side length (cm)")
    ap.add_argument("--water-xy", type=float, default=DETECTOR_XY,
                    help="water side length (cm). If > 30, the water is widened and dose is scored in the central 30x30 cm2")
    args = ap.parse_args()

    mat = args.material
    layers = None
    if ":" in mat:
        # Layered shield, e.g. "G4_Al:10+G4_POLYETHYLENE:10" (the first layer faces the source)
        layers = []
        for part in mat.split("+"):
            mname, ar = part.split(":")
            if mname not in DENSITY:
                sys.exit(f"Unknown material {mname}. Add it to DENSITY.")
            layers.append((mname, float(ar)))
        if abs(sum(a for _, a in layers) - args.areal) > 1e-6:
            sys.exit("The layer areal densities must add up to --areal.")
        short = {"G4_POLYETHYLENE": "PE", "G4_PLEXIGLASS": "PMMA"}
        mat = "_".join(f"{short.get(m, m.replace('G4_', ''))}{a:g}" for m, a in layers)
        thick = sum(a / DENSITY[m] for m, a in layers)
    else:
        if mat not in DENSITY:
            sys.exit(f"Unknown material {mat}. Add it to DENSITY.")
        thick = args.areal / DENSITY[mat]                      # cm
    nz = int(round(args.phantom_cm * 10 / args.dz))
    if abs(nz * args.dz - args.phantom_cm * 10) > 1e-6:
        sys.exit("--dz does not divide the phantom thickness exactly.")

    edges, weights, group_fluence = spectrum_bins(args.spectrum, args.emin, args.emax, args.nbins)

    m, cm, mm, MeV = gate.g4_units.m, gate.g4_units.cm, gate.g4_units.mm, gate.g4_units.MeV

    sim = gate.Simulation()
    sim.output_dir = args.out
    sim.random_seed = args.seed

    world = sim.world
    tan_t = float(np.tan(np.radians(args.theta_deg)))
    src_x = args.source_xy + (thick + args.phantom_cm) * tan_t          # cm, source side along x
    src_x0 = -(20 + thick + args.phantom_cm / 2) * tan_t               # cm, source centre along x
    half = max(150.0, abs(src_x0) + src_x / 2 + 20, args.slab_xy / 2 + 20, args.water_xy / 2 + 20)
    world.size = [2 * half * cm, 2 * half * cm, 3 * m]                  # 3 m cube for normal incidence
    world.material = "G4_AIR"

    if layers is None:
        shield = sim.add_volume("Box", name="Shielding_Slab")
        shield.size = [args.slab_xy * cm, args.slab_xy * cm, thick * cm]
        shield.translation = [0, 0, -thick * cm / 2]
        shield.material = mat
    else:
        z0 = -thick
        for i, (mname, ar) in enumerate(layers):
            tl = ar / DENSITY[mname]
            sl = sim.add_volume("Box", name=f"Shielding_Slab_{i}")
            sl.size = [args.slab_xy * cm, args.slab_xy * cm, tl * cm]
            sl.translation = [0, 0, (z0 + tl / 2) * cm]
            sl.material = mname
            z0 += tl

    wide = args.water_xy > DETECTOR_XY + 1e-9
    plane_z = PLANE_THICKNESS * cm / 2
    phantom_z = PLANE_THICKNESS * cm + args.phantom_cm * cm / 2
    if wide:
        # Wide water (and wide 1 mm air gap) with the central 30x30 cm2 column as a daughter
        # volume, so that the actors and the normalisation stay the same.
        gap = sim.add_volume("Box", name="Air_Gap_Wide")
        gap.size = [args.water_xy * cm, args.water_xy * cm, PLANE_THICKNESS * cm]
        gap.translation = [0, 0, plane_z]
        gap.material = "G4_AIR"
        water = sim.add_volume("Box", name="Water_Wide")
        water.size = [args.water_xy * cm, args.water_xy * cm, args.phantom_cm * cm]
        water.translation = [0, 0, phantom_z]
        water.material = "G4_WATER"

    plane = sim.add_volume("Box", name="Scoring_Plane")
    plane.size = [DETECTOR_XY * cm, DETECTOR_XY * cm, PLANE_THICKNESS * cm]
    plane.translation = [0, 0, 0] if wide else [0, 0, plane_z]
    plane.material = "G4_AIR"
    if wide:
        plane.mother = "Air_Gap_Wide"

    phantom = sim.add_volume("Box", name="Water_Phantom")
    phantom.size = [DETECTOR_XY * cm, DETECTOR_XY * cm, args.phantom_cm * cm]
    phantom.translation = [0, 0, 0] if wide else [0, 0, phantom_z]
    phantom.material = "G4_WATER"
    if wide:
        phantom.mother = "Water_Wide"

    sim.physics_manager.physics_list_name = args.physics
    for p in ["gamma", "electron", "positron", "proton"]:
        setattr(sim.physics_manager.global_production_cuts, p, args.cut * mm)

    source = sim.add_source("GenericSource", name="SPE_Protons")
    source.particle = args.particle
    A_nuc = {"proton": 1, "alpha": 4}[args.particle]          # total kinetic energy = A x (MeV/n)
    source.energy.type = "spectrum_histogram"                      # verified with opengate 10.0.0
    source.energy.spectrum_energy_bin_edges = list(edges * A_nuc * MeV)    # in Geant4 units (total energy)
    # NOTE: opengate (GateGenericSource.cpp) treats spectrum_weights as a DENSITY (per MeV) and
    # multiplies them by the bin width. Our weights are the fluence INSIDE each bin, so we divide by
    # the bin width. (Without this, with logarithmic bins the spectrum inside each group is tilted
    # by a factor proportional to E - fixed on 22/10/2026.)
    dens = weights / np.diff(edges)
    source.energy.spectrum_weights = list(dens / dens.sum())
    source.energy.spectrum_histogram_interpolation = None         # flat inside each bin (as in opengate test010)
    source.position.type = "box"
    source.position.size = [src_x * cm, args.source_xy * cm, 0]
    source.position.translation = [src_x0 * cm, 0, -thick * cm - 20 * cm]
    source.direction.type = "momentum"
    source.direction.momentum = [float(np.sin(np.radians(args.theta_deg))), 0, float(np.cos(np.radians(args.theta_deg)))]
    source.n = args.n

    sim.add_actor("SimulationStatisticsActor", "Stats")

    dose = sim.add_actor("DoseActor", "Dose")
    dose.attached_to = "Water_Phantom"
    dose.size = [1, 1, nz]
    dose.spacing = [DETECTOR_XY * cm, DETECTOR_XY * cm, args.dz * mm]
    dose.edep.active = True
    dose.dose.active = True
    dose.dose_uncertainty.active = True
    dose.output_filename = f"dose_{mat}.mhd"

    if not args.no_phsp:
        phsp = sim.add_actor("PhaseSpaceActor", "PhaseSpace")
        phsp.attached_to = "Scoring_Plane"
        phsp.attributes = ["EventID", "ParticleName", "PDGCode", "KineticEnergy", "Position", "Direction"]
        phsp.output_filename = f"phsp_{mat}.root"

    if args.split:
        for tag, particles in SPLIT.items():
            names = particles if isinstance(particles, list) else [particles]
            a = sim.add_actor("DoseActor", f"Dose_{tag}")
            a.attached_to = "Water_Phantom"
            a.size = [1, 1, nz]
            a.spacing = [DETECTOR_XY * cm, DETECTOR_XY * cm, args.dz * mm]
            a.dose.active = True
            a.output_filename = f"dose_{mat}_{tag}.mhd"
            # ParticleFilter takes ONE particle name. For a group (gamma, e-, e+) use one filter
            # per particle, combined with "or" (a list matches nothing).
            for nm in names:
                f = sim.add_filter("ParticleFilter", f"f_{tag}_{nm}")
                f.particle = nm
                a.filters.append(f)
            a.filters_boolean_operator = "or"

    if args.dose_eq:
        add_let_actors(sim, mat, nz, args.dz, cm, mm, MeV)

    sim.run()

    with open(os.path.join(args.out, "setup.json"), "w") as fo:
        json.dump({
            "material": mat, "layers": layers, "areal_g_cm2": args.areal, "thickness_cm": thick,
            "spectrum": args.spectrum, "emin_MeV": args.emin, "emax_MeV": args.emax, "particle": args.particle,
            "energy_per_nucleon": args.particle != "proton",
            "group_fluence_cm2": group_fluence,           # event fluence inside the group (cm^-2)
            "source_area_cm2": src_x * args.source_xy, "n": args.n, "theta_deg": args.theta_deg,
            "source_xy_cm": args.source_xy, "slab_xy_cm": args.slab_xy, "water_xy_cm": args.water_xy,
            "phantom_cm": args.phantom_cm, "dz_mm": args.dz,
            "physics": args.physics, "weights_as_density": True, "seed": args.seed, "cut_mm": args.cut, "split": args.split,
            "dose_eq": args.dose_eq,
        }, fo, indent=1)
    print(f"--> Done: {mat} {args.areal} g/cm2 ({thick:.3f} cm), "
          f"group {args.emin}-{args.emax} MeV, Phi = {group_fluence:.3e} cm^-2")


if __name__ == "__main__":
    main()
