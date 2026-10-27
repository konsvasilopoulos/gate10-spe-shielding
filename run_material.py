import os
import sys
import json
import argparse
import opengate as gate

# Usage:
#   python run_material.py G4_Al                  -> total dose + phase space + setup.json (cut 0.1 mm)
#   python run_material.py G4_Al --split          -> + dose per particle type
#   python run_material.py G4_Cu --seed 4321 --physics G4EmStandardPhysics_option4
#   python run_material.py G4_Cu --cut 0.7        -> old production cut, for comparison
#   python run_material.py G4_Ti --dz 0.2         -> 0.2 mm voxels along z (1000 voxels)
#
# ALL slab thicknesses are defined ONLY here.
THICKNESS_CM = {
    "G4_AIR": 2,
    "G4_POLYETHYLENE": 2,
    "G4_PLEXIGLASS": 4,
    "G4_Al": 2,
    "G4_Ti": 2,
    "G4_Fe": 1.1,
    "G4_Cu": 1,
    "G4_Pb": 0.4,
}

SOURCE_XY, SLAB_XY, DETECTOR_XY = 40, 60, 30   # cm
PLANE_THICKNESS = 0.1                          # cm (1 mm)
PHANTOM_THICKNESS = 20                         # cm
E0_MEV = 150
DZ_MM = 1.0                                    # default voxel size along z (mm); change with --dz

# Filter syntax for the per-particle dose split
SPLIT = {"proton": "proton", "neutron": "neutron", "photon_electron": ["gamma", "e-", "e+"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("material")
    ap.add_argument("--split", action="store_true")
    ap.add_argument("--phsp-origin", action="store_true",
                    help="phase space with creator process/position of each particle (hadronic-origin study)")
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--n", type=int, default=250000)
    ap.add_argument("--cut", type=float, default=0.1,
                    help="production cut in mm (default 0.1)")
    ap.add_argument("--physics", default="QGSP_BIC_EMZ")
    ap.add_argument("--dz", type=float, default=DZ_MM, help="voxel size along z in mm (default 1)")
    ap.add_argument("--source-xy", type=float, default=SOURCE_XY, help="source side length (cm)")
    ap.add_argument("--slab-xy", type=float, default=SLAB_XY, help="slab side length (cm)")
    ap.add_argument("--water-xy", type=float, default=DETECTOR_XY,
                    help="water side length (cm); if > 30, dose is scored in the central 30x30 cm2")
    ap.add_argument("--out", default=None, help="output folder (default: output_<material>)")
    args = ap.parse_args()

    mat = args.material
    if mat not in THICKNESS_CM:
        sys.exit(f"Unknown material {mat}. Add it to THICKNESS_CM.")
    thick = THICKNESS_CM[mat]
    out = args.out or f"output_{mat}"
    NZ = int(round(PHANTOM_THICKNESS * 10 / args.dz))
    if abs(NZ * args.dz - PHANTOM_THICKNESS * 10) > 1e-6:
        sys.exit(f"--dz {args.dz} mm does not divide the {PHANTOM_THICKNESS} cm phantom exactly.")

    m = gate.g4_units.m
    cm = gate.g4_units.cm
    mm = gate.g4_units.mm
    MeV = gate.g4_units.MeV

    sim = gate.Simulation()
    sim.output_dir = out
    sim.random_seed = args.seed

    world = sim.world
    world.size = [2 * m, 2 * m, 2 * m]
    world.material = "G4_AIR"

    shield = sim.add_volume("Box", name="Shielding_Slab")
    shield.size = [args.slab_xy * cm, args.slab_xy * cm, thick * cm]
    shield.translation = [0, 0, -thick * cm / 2]
    shield.material = mat

    wide = args.water_xy > DETECTOR_XY + 1e-9
    plane_z = PLANE_THICKNESS * cm / 2
    phantom_z = PLANE_THICKNESS * cm + PHANTOM_THICKNESS * cm / 2
    if wide:
        # wide water: the central 30x30 cm2 column is a daughter volume (same actors and normalisation)
        gap = sim.add_volume("Box", name="Air_Gap_Wide")
        gap.size = [args.water_xy * cm, args.water_xy * cm, PLANE_THICKNESS * cm]
        gap.translation = [0, 0, plane_z]
        gap.material = "G4_AIR"
        water = sim.add_volume("Box", name="Water_Wide")
        water.size = [args.water_xy * cm, args.water_xy * cm, PHANTOM_THICKNESS * cm]
        water.translation = [0, 0, phantom_z]
        water.material = "G4_WATER"

    plane = sim.add_volume("Box", name="Scoring_Plane")
    plane.size = [DETECTOR_XY * cm, DETECTOR_XY * cm, PLANE_THICKNESS * cm]
    plane.translation = [0, 0, 0] if wide else [0, 0, plane_z]
    plane.material = "G4_AIR"
    if wide:
        plane.mother = "Air_Gap_Wide"

    phantom = sim.add_volume("Box", name="Water_Phantom")
    phantom.size = [DETECTOR_XY * cm, DETECTOR_XY * cm, PHANTOM_THICKNESS * cm]
    phantom.translation = [0, 0, 0] if wide else [0, 0, phantom_z]
    phantom.material = "G4_WATER"
    if wide:
        phantom.mother = "Water_Wide"

    sim.physics_manager.physics_list_name = args.physics

    # production cut (same for all relevant particles)
    for p in ["gamma", "electron", "positron", "proton"]:
        setattr(sim.physics_manager.global_production_cuts, p, args.cut * mm)

    source = sim.add_source("GenericSource", name="Proton_Beam")
    source.particle = "proton"
    source.energy.mono = E0_MEV * MeV
    source.position.type = "box"
    source.position.size = [args.source_xy * cm, args.source_xy * cm, 0]
    source.position.translation = [0, 0, -thick * cm - 20 * cm]
    source.direction.type = "momentum"
    source.direction.momentum = [0, 0, 1]
    source.n = args.n

    sim.add_actor("SimulationStatisticsActor", "Stats")

    dose = sim.add_actor("DoseActor", "Dose")
    dose.attached_to = "Water_Phantom"
    dose.size = [1, 1, NZ]
    dose.spacing = [DETECTOR_XY * cm, DETECTOR_XY * cm, args.dz * mm]
    dose.edep.active = True
    dose.dose.active = True
    dose.dose_uncertainty.active = True
    dose.output_filename = f"dose_{mat}.mhd"

    phsp = sim.add_actor("PhaseSpaceActor", "PhaseSpace")
    phsp.attached_to = "Scoring_Plane"
    phsp.attributes = ["EventID", "ParticleName", "PDGCode", "KineticEnergy", "Position", "Direction"]
    if args.phsp_origin:
        phsp.attributes += ["TrackCreatorProcess", "TrackCreatorModelName", "TrackVertexPosition",
                            "TrackVertexKineticEnergy", "ParentParticleName"]
    phsp.output_filename = f"phsp_{mat}.root"

    if args.split:
        for tag, particles in SPLIT.items():
            names = particles if isinstance(particles, list) else [particles]
            a = sim.add_actor("DoseActor", f"Dose_{tag}")
            a.attached_to = "Water_Phantom"
            a.size = [1, 1, NZ]
            a.spacing = [DETECTOR_XY * cm, DETECTOR_XY * cm, args.dz * mm]
            a.dose.active = True
            a.output_filename = f"dose_{mat}_{tag}.mhd"
            # one ParticleFilter per particle, combined with "or" (the filter takes a single name)
            for nm in names:
                f = sim.add_filter("ParticleFilter", f"f_{tag}_{nm}")
                f.particle = nm
                a.filters.append(f)
            a.filters_boolean_operator = "or"

    sim.run()

    with open(os.path.join(out, "setup.json"), "w") as fo:
        json.dump({"material": mat, "thickness_cm": thick, "n": args.n,
                   "energy_MeV": E0_MEV, "physics": args.physics, "seed": args.seed,
                   "cut_mm": args.cut, "dz_mm": args.dz, "source_xy_cm": args.source_xy,
                   "slab_xy_cm": args.slab_xy, "water_xy_cm": args.water_xy}, fo)
    print(f"--> Done: {mat} ({out}), cut {args.cut} mm, voxel {args.dz} mm")


if __name__ == "__main__":
    main()
