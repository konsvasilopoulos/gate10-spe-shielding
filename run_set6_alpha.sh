#!/usr/bin/env bash
# run_set6_alpha.sh - SET 6: SPE alpha particles.
#   Koldobskiy 1956, 1972 and 1989 spectra taken per nucleon, with He/H = 1 at the same MeV/n
#   (results scale linearly with the actual He/H ratio, typically 0.02-0.06 in large gradual events).
#   Al, PE, Pb x 5/10/20 g/cm2, the same 7 groups (30-5000 MeV/n), D and H (LET actors), wide geometry.
#   NSCALE = 0.1 (alphas cost more per particle; their contribution is small).
#   Output: spe_runs/alpha_<spectrum>/<material>_<g/cm2>/...      Then: python3 analysis/alpha_contribution.py
# Usage:  bash run_set6_alpha.sh 2>&1 | tee spe_runs/set6.log     (~2 h estimated, keeps the Mac awake)
# If interrupted, run the same command again; it resumes where it stopped.
set -u
# Keeps the Mac awake while the script runs (also with the lid closed); stops by itself at the end.
if command -v caffeinate >/dev/null 2>&1; then caffeinate -i -w $$ & fi
for sp in kold1972 kold1956 kold1989; do
  echo "===== Set 6: alpha, ${sp} | $(date) ====="
  SPECTRUM=${sp} NSCALE=0.1 MATERIALS="G4_Al G4_POLYETHYLENE G4_Pb" NO_FLUENCE=1 EXTRA_ARGS="--particle alpha --no-phsp" \
  BASE_OVERRIDE=spe_runs/alpha_${sp} bash run_materials_spe.sh
done
echo "===== End of set 6: $(date) ====="
python3 analysis/alpha_contribution.py
