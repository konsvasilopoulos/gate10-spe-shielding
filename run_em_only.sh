#!/usr/bin/env bash
# run_em_only.sh - Material comparison at 20 g/cm2 with electromagnetic physics ONLY (G4EmStandardPhysics_option4),
# for the Koldobskiy spectra (set SPECTRA="..." for others). Splits the PE advantage into stopping and nuclear parts.
# Usage:  bash run_em_only.sh      (then: python analysis/decompose_em_nuclear.py)
set -u
if command -v caffeinate >/dev/null 2>&1; then caffeinate -i -w $$ & fi
for sp in ${SPECTRA:-kold1956 kold1972}; do
    PHYSICS=G4EmStandardPhysics_option4 SPECTRUM=$sp AREALS=20 \
        BASE_OVERRIDE=spe_runs/materials_EMonly_${sp} bash run_materials_spe.sh
done
python analysis/decompose_em_nuclear.py | tee spe_runs/decompose_em_nuclear.txt
