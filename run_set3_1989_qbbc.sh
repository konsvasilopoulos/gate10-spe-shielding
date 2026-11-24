#!/usr/bin/env bash
# run_set3_1989_qbbc.sh - SET 3:
#   1. Materials with the September-October 1989 SPE series (Koldobskiy et al. 2021, GLE 42-45), 7 materials x 5/10/20 g/cm2 (LiH and water: run_set7_lih_water.sh)
#      (NASA design-reference event for storm shelters, Townsend et al. 2018)
#   2. Physics-list test of the RATIOS to Al: QBBC instead of QGSP_BIC_EMZ, PE/Al/Pb at 20 g/cm2, Koldobskiy 1972
# Usage:  bash run_set3_1989_qbbc.sh 2>&1 | tee spe_runs/set3.log     (~3.5-4 h, keeps the Mac awake)
# If interrupted, run the same command again; it resumes where it stopped.
set -u
if command -v caffeinate >/dev/null 2>&1; then caffeinate -i -w $$ & fi
echo "===== Materials kold1989 | $(date) ====="
SPECTRUM=kold1989 bash run_materials_spe.sh
echo "===== QBBC, kold1972, 20 g/cm² · $(date) ====="
PHYSICS=QBBC SPECTRUM=kold1972 AREALS=20 MATERIALS="G4_POLYETHYLENE G4_Al G4_Pb" \
    BASE_OVERRIDE=spe_runs/materials_QBBC_kold1972 bash run_materials_spe.sh
echo "===== End of set 3: $(date) ====="
