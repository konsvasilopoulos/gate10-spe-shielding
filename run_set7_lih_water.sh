#!/usr/bin/env bash
# =============================================================================
# run_set7_lih_water.sh - SET 7 (November 2026): lithium hydride and water as shield materials.
#   0. Short test (~2 min) that the new materials are built: LiH, water, Al + 2 g/cm2 cortical bone.
#   1. Direct runs of G4_LITHIUM_HYDRIDE (0.82 g/cm3) and G4_WATER at 5/10/20 g/cm2 for kold1972, kold1956,
#      kold1989, king1972 and webber1956, with the settings of sets 1-3 and into the same folders
#      (spe_runs/materials*), so that compare_materials_spe.py gives the nine-material tables.
#   The borated-polyethylene test of the article (5 % boron, 20 g/cm2, K21 1956 and 1972) needs a user-defined
#   material with its own ionisation potential, which opengate 10.0.0 does not support; it was run with
#   opengate 10.1.1 and a custom-material library that is not part of this repository (composition in README).
#   GATE 10.0.0, JOBS=4, ~3 h on 4 cores of an Apple M5 (measured: 210 min including the borated-PE runs).
# Usage:  bash run_set7_lih_water.sh 2>&1 | tee spe_runs/set7.log
#   then: python compare_materials_spe.py spe_runs/materials_<spectrum>   (Webber: spe_runs/materials)
# If interrupted, run the same command again; it resumes where it stopped.
# =============================================================================
set -u
cd "$(dirname "$0")"
# Keeps the Mac awake while the script runs (also with the lid closed); stops by itself at the end.
if command -v caffeinate >/dev/null 2>&1; then caffeinate -i -w $$ & fi
T0=$(date +%s); echo "Set 7 started: $(date)"

echo "===== 0. Material test (2000 protons, 150-170 MeV) ====="
mkdir -p spe_runs/set7_test
for m in G4_LITHIUM_HYDRIDE G4_WATER "G4_Al:20+G4_BONE_CORTICAL_ICRP:2"; do
    ar=20; [[ "$m" == *BONE* ]] && ar=22
    tag=$(echo "$m" | sed 's/G4_//g; s/[:+]/_/g')
    python run_spectrum.py "$m" --areal $ar --spectrum kold1972 --emin 150 --emax 170 --n 2000 --seed 1234 --no-phsp \
        --source-xy 80 --slab-xy 120 --water-xy 120 --out spe_runs/set7_test/$tag > spe_runs/set7_test/$tag.log 2>&1
    if [ -f spe_runs/set7_test/$tag/setup.json ]; then echo "[ok] $m"; else echo "[FAIL] $m -> spe_runs/set7_test/$tag.log"; exit 1; fi
done

for sp in kold1972 kold1956 kold1989 king1972 webber1956; do
    t=$(date +%s); echo "===== 1. LiH + water, ${sp} - $(date '+%H:%M') ====="
    SPECTRUM=$sp MATERIALS="G4_LITHIUM_HYDRIDE G4_WATER" bash run_materials_spe.sh
    echo "      time ${sp}: $(( ($(date +%s) - t) / 60 )) min"
done
echo "===== Set 7 finished: $(date) - $(( ($(date +%s) - T0) / 60 )) min ====="
