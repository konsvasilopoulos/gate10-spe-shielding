#!/usr/bin/env bash
# run_hadronic_origin.sh - Origin of the hadronic contribution to the entrance dose (150 MeV validation).
# 8 runs (air + 7 materials), QGSP_BIC_EMZ, cut 0.1 mm, validation geometry, phase space with the creator
# process and position of each particle. Then: python analysis/hadronic_origin.py
# Usage:  bash run_hadronic_origin.sh     (~5-10 min)
set -u
if command -v caffeinate >/dev/null 2>&1; then caffeinate -i -w $$ & fi
JOBS="${JOBS:-4}"
mkdir -p runs_origin/logs
run_one() {
    m="$1"; out="runs_origin/output_${m}"
    [ -f "${out}/setup.json" ] && { echo "[skip] $m"; return 0; }
    python run_material.py "$m" --seed 1234 --cut 0.1 --physics QGSP_BIC_EMZ --n 250000 --split --phsp-origin \
        --out "${out}" > "runs_origin/logs/${m}.log" 2>&1
    [ -f "${out}/setup.json" ] && echo "[ok]   $m  $(date '+%H:%M')" || echo "[FAIL] $m -> runs_origin/logs/${m}.log"
}
export -f run_one
for m in G4_AIR G4_POLYETHYLENE G4_PLEXIGLASS G4_Al G4_Ti G4_Fe G4_Cu G4_Pb; do echo "$m"; done \
    | xargs -P "${JOBS}" -L 1 bash -c 'run_one "$@"' _
python analysis/hadronic_origin.py | tee runs_origin/hadronic_origin_result.txt
