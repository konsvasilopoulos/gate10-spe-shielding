#!/usr/bin/env bash
# run_pb_threshold_test.sh - Proton transmission near the threshold of 5 g/cm2 Pb (and Al for comparison):
# narrow 1 MeV energy windows; transmitted fraction and exit energy compared with PSTAR (CSDA).
# Usage:  bash run_pb_threshold_test.sh   (~10 min)   then: python analysis/pb_threshold_test.py
set -u
if command -v caffeinate >/dev/null 2>&1; then caffeinate -i -w $$ & fi
JOBS="${JOBS:-4}"
mkdir -p spe_runs/pb_threshold/logs
run_one() {
    mat="$1"; e="$2"; lo=$(python -c "print($e-0.5)"); hi=$(python -c "print($e+0.5)")
    out="spe_runs/pb_threshold/${mat#G4_}_E${e}"
    [ -f "${out}/setup.json" ] && return 0
    python run_spectrum.py "$mat" --areal 5 --spectrum king1972 --emin "$lo" --emax "$hi" --nbins 5 --n 40000 --seed 1234 \
        --source-xy 80 --slab-xy 120 --water-xy 120 --phantom-cm 5 --out "$out" > "spe_runs/pb_threshold/logs/${mat#G4_}_${e}.log" 2>&1
    [ -f "${out}/setup.json" ] && echo "[ok] $mat $e MeV" || echo "[FAIL] $mat $e"
}
export -f run_one
{ for e in 46 48 50 52 54 56 58 60; do echo "G4_Pb $e"; done; for e in 64 66 68 70 72 74; do echo "G4_Al $e"; done; } \
    | xargs -P "$JOBS" -L 1 bash -c 'run_one "$@"' _
python analysis/pb_threshold_test.py | tee spe_runs/pb_threshold/pb_threshold_result.txt
