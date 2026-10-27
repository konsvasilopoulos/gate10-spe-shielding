#!/usr/bin/env bash
# run_set4_response.sh - SET 4: response matrix.
#   7 materials x 5/10/20 g/cm2 x 25 narrow logarithmic energy windows (30 MeV - 5 GeV), 80 000 protons each,
#   with D, H (LET actors) and the per-particle split; no phase space.
#   Then: python3 analysis/response_matrix.py build   and   python3 analysis/response_matrix.py validate
# Usage:  bash run_set4_response.sh 2>&1 | tee spe_runs/set4.log     (~2.5-3 h, keeps the Mac awake)
# If interrupted, run the same command again; it resumes where it stopped.
set -u
if command -v caffeinate >/dev/null 2>&1; then caffeinate -i -w $$ & fi
echo "===== Response matrix | $(date) ====="
SPECTRUM=logflat NSCALE=1 NO_FLUENCE=1 EXTRA_ARGS="--no-phsp" BASE_OVERRIDE=spe_runs/response_matrix \
GROUPS_N_OVERRIDE="30:36.8:80000 36.8:45.2:80000 45.2:55.4:80000 55.4:68:80000 68:83.5:80000 83.5:102.4:80000 102.4:125.7:80000 125.7:154.2:80000 154.2:189.2:80000 189.2:232.2:80000 232.2:284.9:80000 284.9:349.6:80000 349.6:429:80000 429:526.4:80000 526.4:646:80000 646:792.7:80000 792.7:972.7:80000 972.7:1193.6:80000 1193.6:1464.6:80000 1464.6:1797.2:80000 1797.2:2205.3:80000 2205.3:2706.1:80000 2706.1:3320.6:80000 3320.6:4074.7:80000 4074.7:5000:80000" \
    bash run_materials_spe.sh
echo "===== End of set 4: $(date) ====="
