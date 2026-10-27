#!/usr/bin/env bash
# run_quick_tests.sh - two short tests (~30-35 min in total; keeps the Mac awake):
#   A. Lateral geometry: source 120x120, slab and water 200x200 cm (instead of 80/120/120),
#      Al, PE, Pb at 20 g/cm2, Koldobskiy 1956 (hardest spectrum, most secondaries)
#   B. Layered shields of 20 g/cm2, Koldobskiy 1972: Al 10 + PE 10 and PE 10 + Al 10
#      (the first layer faces the source), compared with Al 20 and PE 20
# Usage:  bash run_quick_tests.sh 2>&1 | tee spe_runs/quick_tests.log      then: python analysis/quick_tests.py
set -u
if command -v caffeinate >/dev/null 2>&1; then caffeinate -i -w $$ & fi
echo "===== A. Lateral geometry | $(date) ====="
SPECTRUM=kold1956 AREALS=20 MATERIALS="G4_Al G4_POLYETHYLENE G4_Pb" \
    WIDE="--source-xy 120 --slab-xy 200 --water-xy 200" BASE_OVERRIDE=spe_runs/lateral_kold1956 bash run_materials_spe.sh
echo "===== B. Layered shields | $(date) ====="
SPECTRUM=kold1972 AREALS=20 MATERIALS="G4_Al:10+G4_POLYETHYLENE:10 G4_POLYETHYLENE:10+G4_Al:10" \
    BASE_OVERRIDE=spe_runs/layered_kold1972 bash run_materials_spe.sh
echo "===== Finished: $(date) ====="
python analysis/quick_tests.py
