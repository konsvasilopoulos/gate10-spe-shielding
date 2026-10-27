#!/usr/bin/env bash
# run_validation_150MeV.sh - 150 MeV monoenergetic validation (stopping power, range shift, nuclear attenuation).
#   Air (reference) + 7 materials, three independent seeds, QGSP_BIC_EMZ, cut 0.1 mm, 250 000 protons each.
#   Output: runs/s<seed>/output_<material>/
#   Then:   python validation/compare_pstar.py   and   python validation/seed_stats.py runs
# Options (environment variables):
#   SEEDS="1234 2345 3456"   PHYSICS=G4EmStandardPhysics_option4 (EM-only entrance-dose check)
#   EXTRA_ARGS="--dz 0.2"    (voxel test)   or   EXTRA_ARGS="--source-xy 80 --slab-xy 120 --water-xy 120" (wide geometry)
#   BASE=runs_vox02 / runs_wide
# Usage:  bash run_validation_150MeV.sh      (~1 h on 4 cores)
set -u
# Keeps the Mac awake while the script runs (also with the lid closed); stops by itself at the end.
if command -v caffeinate >/dev/null 2>&1; then caffeinate -i -w $$ & fi

SEEDS="${SEEDS:-1234 2345 3456}"
PHYSICS="${PHYSICS:-QGSP_BIC_EMZ}"
BASE="${BASE:-runs}"
EXTRA_ARGS="${EXTRA_ARGS:-}"
MATERIALS="G4_AIR G4_POLYETHYLENE G4_PLEXIGLASS G4_Al G4_Ti G4_Fe G4_Cu G4_Pb"

for s in ${SEEDS}; do
    for m in ${MATERIALS}; do
        out="${BASE}/s${s}/output_${m}"
        if [ -f "${out}/setup.json" ]; then echo "[skip] ${out}"; continue; fi
        mkdir -p "${BASE}/s${s}"
        python run_material.py "${m}" --seed "${s}" --physics "${PHYSICS}" --out "${out}" ${EXTRA_ARGS} \
            > "${BASE}/s${s}/${m}.log" 2>&1 && echo "[ok]   ${out}" || echo "[FAIL] ${out} -> see ${BASE}/s${s}/${m}.log"
    done
done
echo "Finished: $(date)"
