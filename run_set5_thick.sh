#!/usr/bin/env bash
# run_set5_thick.sh - SET 5: shields of 30 and 40 g/cm2, so that the NASA-limit table has no "> 20" entries.
#   A. Direct check: Koldobskiy 1989, Al and PE at 30 g/cm2 (7 groups, D and H) -> test of the folding above 20 g/cm2.
#   B. Response matrix: 7 materials (LiH and water: run_set8_lih_water_rm.sh) x 30/40 g/cm2 x 25 narrow windows (30 MeV - 5 GeV), 80 000 protons each.
#   Then: python3 analysis/response_matrix.py build_thick   and   python3 analysis/nasa_limits_thick.py
# Usage:  bash run_set5_thick.sh 2>&1 | tee spe_runs/set5.log     (~1-2.5 h, keeps the Mac awake)
# If interrupted, run the same command again; it resumes where it stopped.
set -u
# Keeps the Mac awake while the script runs (also with the lid closed); stops by itself at the end.
if command -v caffeinate >/dev/null 2>&1; then caffeinate -i -w $$ & fi
echo "===== Set 5A: direct check kold1989, Al/PE 30 g/cm2 | $(date) ====="
SPECTRUM=kold1989 AREALS="30" MATERIALS="G4_Al G4_POLYETHYLENE" NO_FLUENCE=1 EXTRA_ARGS="--no-phsp" \
BASE_OVERRIDE=spe_runs/materials_kold1989_thick bash run_materials_spe.sh
echo "===== Set 5B: response matrix 30/40 g/cm2 | $(date) ====="
SPECTRUM=logflat NSCALE=1 AREALS="30 40" NO_FLUENCE=1 EXTRA_ARGS="--no-phsp" BASE_OVERRIDE=spe_runs/response_matrix_thick \
GROUPS_N_OVERRIDE="30:36.8:80000 36.8:45.2:80000 45.2:55.4:80000 55.4:68:80000 68:83.5:80000 83.5:102.4:80000 102.4:125.7:80000 125.7:154.2:80000 154.2:189.2:80000 189.2:232.2:80000 232.2:284.9:80000 284.9:349.6:80000 349.6:429:80000 429:526.4:80000 526.4:646:80000 646:792.7:80000 792.7:972.7:80000 972.7:1193.6:80000 1193.6:1464.6:80000 1464.6:1797.2:80000 1797.2:2205.3:80000 2205.3:2706.1:80000 2706.1:3320.6:80000 3320.6:4074.7:80000 4074.7:5000:80000" \
    bash run_materials_spe.sh
echo "===== End of set 5: $(date) ====="
