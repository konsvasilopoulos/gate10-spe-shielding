#!/usr/bin/env bash
# =============================================================================
# run_dose_eq.sh - Webber 1956, 20 g/cm2 Al + 30 g/cm2 water, wide geometry, with LET actors
# for dose equivalent (ICRP-60). Output: spe_runs/webber1956_Al20_wide_H
# Usage:  JOBS=4 bash run_dose_eq.sh       (runs analyze_dose_eq.py at the end)
# Events: NSCALE multiplies the N per group (default 1).
# =============================================================================
set -u
# Keeps the Mac awake while the script runs (also with the lid closed); stops by itself at the end.
if command -v caffeinate >/dev/null 2>&1; then caffeinate -i -w $$ & fi
JOBS="${JOBS:-4}"   # parallel runs (number of fast CPU cores)
NSCALE="${NSCALE:-1}"
BASE="spe_runs/webber1956_Al20_wide_H"
WIDE="--source-xy 80 --slab-xy 120 --water-xy 120"
GROUPS_N="10:60:800000 60:120:800000 120:170:1000000 170:250:1000000 250:400:1000000 400:2000:1000000"
mkdir -p "${BASE}/logs"
export BASE WIDE NSCALE

run_group() {
    IFS=: read -r emin emax n <<< "$1"
    n=$(python -c "print(int(${n} * ${NSCALE}))")
    out="${BASE}/g${emin}-${emax}"; log="${BASE}/logs/g${emin}-${emax}.log"
    if [ -f "${out}/setup.json" ]; then echo "[skip] ${emin}–${emax} MeV"; return 0; fi
    echo "[start] ${emin}–${emax} MeV (N=${n})  $(date '+%H:%M:%S')"; t0=$(date +%s)
    python run_spectrum.py G4_Al --areal 20 --emin "${emin}" --emax "${emax}" --n "${n}" --seed 1234 \
        --split --dose-eq ${WIDE} --out "${out}" > "${log}" 2>&1
    if [ -f "${out}/setup.json" ]; then echo "[ok]    ${emin}–${emax} MeV  ($(( $(date +%s) - t0 )) s)"
    else echo "[FAIL]  ${emin}-${emax} -> see ${log}"; echo "${emin}-${emax}" >> "${BASE}/failed.txt"; fi
}
export -f run_group
rm -f "${BASE}/failed.txt"
for g in ${GROUPS_N}; do echo "$g"; done | xargs -P "${JOBS}" -L 1 bash -c 'run_group "$@"' _
if [ -f "${BASE}/failed.txt" ]; then echo "Failed:"; cat "${BASE}/failed.txt"; fi
python analyze_dose_eq.py "${BASE}" 2>&1 | tee "${BASE}/dose_eq_result.txt"
