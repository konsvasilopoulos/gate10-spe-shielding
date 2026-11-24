#!/usr/bin/env bash
# =============================================================================
# run_set8_lih_water_rm.sh - SET 8 (November 2026), after set 7.
#   1. Response matrix for lithium hydride and water: 5/10/20 g/cm2 (spe_runs/response_matrix) and 30/40 g/cm2
#      (spe_runs/response_matrix_thick), the same 25 windows and events as sets 4 and 5.
#      Then: python analysis/response_matrix.py add   (adds the two materials to analysis/response_matrix.npz)
#   2. SPE alpha particles for LiH and water, kold1972/1956/1989, 5/10/20 g/cm2, NSCALE 0.1 (as set 6).
#      Then: python analysis/alpha_contribution.py
#   3. Bone layer: 2 g/cm2 of G4_BONE_CORTICAL_ICRP between the shield (Al, PE, Pb at 20 g/cm2) and the water,
#      kold1972 and kold1956 -> spe_runs/bone_<spectrum>/. Then: python compare_materials_spe.py spe_runs/bone_<spectrum>
#      and python analysis/bone_bpe.py (ratios to Al with and without the bone layer).
#   GATE 10.0.0, JOBS=4, measured 199 min on 4 cores of an Apple M5.
# Usage:  bash run_set8_lih_water_rm.sh 2>&1 | tee spe_runs/set8.log
# If interrupted, run the same command again; it resumes where it stopped.
# =============================================================================
set -u
cd "$(dirname "$0")"
# Keeps the Mac awake while the script runs (also with the lid closed); stops by itself at the end.
if command -v caffeinate >/dev/null 2>&1; then caffeinate -i -w $$ & fi
T0=$(date +%s); echo "Set 8 started: $(date)"
W="30:36.8:80000 36.8:45.2:80000 45.2:55.4:80000 55.4:68:80000 68:83.5:80000 83.5:102.4:80000 102.4:125.7:80000 125.7:154.2:80000 154.2:189.2:80000 189.2:232.2:80000 232.2:284.9:80000 284.9:349.6:80000 349.6:429:80000 429:526.4:80000 526.4:646:80000 646:792.7:80000 792.7:972.7:80000 972.7:1193.6:80000 1193.6:1464.6:80000 1464.6:1797.2:80000 1797.2:2205.3:80000 2205.3:2706.1:80000 2706.1:3320.6:80000 3320.6:4074.7:80000 4074.7:5000:80000"

t=$(date +%s); echo "===== 1a. Response matrix LiH + water, 5/10/20 g/cm2 - $(date '+%H:%M') ====="
SPECTRUM=logflat NSCALE=1 NO_FLUENCE=1 EXTRA_ARGS="--no-phsp" BASE_OVERRIDE=spe_runs/response_matrix \
GROUPS_N_OVERRIDE="$W" MATERIALS="G4_LITHIUM_HYDRIDE G4_WATER" bash run_materials_spe.sh
echo "      time: $(( ($(date +%s) - t) / 60 )) min"

t=$(date +%s); echo "===== 1b. Response matrix LiH + water, 30/40 g/cm2 - $(date '+%H:%M') ====="
SPECTRUM=logflat NSCALE=1 AREALS="30 40" NO_FLUENCE=1 EXTRA_ARGS="--no-phsp" BASE_OVERRIDE=spe_runs/response_matrix_thick \
GROUPS_N_OVERRIDE="$W" MATERIALS="G4_LITHIUM_HYDRIDE G4_WATER" bash run_materials_spe.sh
echo "      time: $(( ($(date +%s) - t) / 60 )) min"

for sp in kold1972 kold1956 kold1989; do
  t=$(date +%s); echo "===== 2. alpha, LiH + water, ${sp} - $(date '+%H:%M') ====="
  SPECTRUM=${sp} NSCALE=0.1 MATERIALS="G4_LITHIUM_HYDRIDE G4_WATER" NO_FLUENCE=1 EXTRA_ARGS="--particle alpha --no-phsp" \
  BASE_OVERRIDE=spe_runs/alpha_${sp} bash run_materials_spe.sh
  echo "      time: $(( ($(date +%s) - t) / 60 )) min"
done

for sp in kold1972 kold1956; do
  t=$(date +%s); echo "===== 3. Bone layer 2 g/cm2, ${sp} - $(date '+%H:%M') ====="
  SPECTRUM=${sp} AREALS="22" NO_FLUENCE=1 EXTRA_ARGS="--no-phsp" BASE_OVERRIDE=spe_runs/bone_${sp} \
  MATERIALS="G4_Al:20+G4_BONE_CORTICAL_ICRP:2 G4_POLYETHYLENE:20+G4_BONE_CORTICAL_ICRP:2 G4_Pb:20+G4_BONE_CORTICAL_ICRP:2" \
  bash run_materials_spe.sh
  echo "      time: $(( ($(date +%s) - t) / 60 )) min"
done
echo "===== Set 8 finished: $(date) - $(( ($(date +%s) - T0) / 60 )) min ====="
