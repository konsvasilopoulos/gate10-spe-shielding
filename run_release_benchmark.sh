#!/usr/bin/env bash
# run_release_benchmark.sh - All comparisons with external references repeated with GATE 10.1.1 (Geant4 11.4.2).
# The set was fixed before the runs: (1) NASA benchmark, Webber 1956, 20 g/cm2 Al (same settings as run_dose_eq.sh)
#   -> spe_runs/gate101_webber1956_Al20_wide_H ; (2) HZETRN cases: King 1972 and Webber 1956, Al / PE / Pb at 20 g/cm2
#   (same settings as set 2) -> spe_runs/gate101_king1972, spe_runs/gate101_webber1956
#   (Koldobskiy 1972: run_release_check.sh). Analysis: python analysis/release_benchmark.py
# Uses the virtual environment with opengate 10.1.1 (VENV, default .venv-gate101). About 2-2.5 h with JOBS=4.
set -u
cd "$(dirname "$0")"
if command -v caffeinate >/dev/null 2>&1; then caffeinate -i -w $$ & fi
source "${VENV:-.venv-gate101}/bin/activate"
python -c "from importlib.metadata import version; print('GATE', version('opengate'))"
JOBS="${JOBS:-4}"

echo "=== 1/3  NASA benchmark (Webber 1956, 20 g/cm2 Al)  ($(date '+%H:%M'))"
BASE="spe_runs/gate101_webber1956_Al20_wide_H"
WIDE="--source-xy 80 --slab-xy 120 --water-xy 120"
GROUPS_N="400:2000:1000000 250:400:1000000 170:250:1000000 120:170:1000000 60:120:800000 10:60:800000"
mkdir -p "${BASE}/logs"
export BASE WIDE
run_group() {
    IFS=: read -r emin emax n <<< "$1"
    out="${BASE}/g${emin}-${emax}"; log="${BASE}/logs/g${emin}-${emax}.log"
    if [ -f "${out}/setup.json" ]; then echo "[skip] ${emin}-${emax} MeV"; return 0; fi
    python run_spectrum.py G4_Al --areal 20 --spectrum webber1956 --emin "${emin}" --emax "${emax}" --n "${n}" --seed 1234 \
        --split --dose-eq ${WIDE} --out "${out}" > "${log}" 2>&1 \
        && echo "[ok]   ${emin}-${emax} MeV" || echo "[FAIL] ${emin}-${emax} -> ${log}"
}
export -f run_group
for g in ${GROUPS_N}; do echo "$g"; done | xargs -P "${JOBS}" -L 1 bash -c 'run_group "$@"' _
python analyze_spectrum.py "${BASE}" > "${BASE}/analyze_spectrum.txt" 2>&1
python analyze_dose_eq.py "${BASE}" > "${BASE}/dose_eq_result.txt" 2>&1
python compare_nasa.py "${BASE}" 2>&1 | tail -30

for sp in king1972 webber1956; do
    echo "=== Materials ${sp}: Al / PE / Pb, 20 g/cm2  ($(date '+%H:%M'))"
    SPECTRUM=${sp} MATERIALS="G4_Al G4_POLYETHYLENE G4_Pb" AREALS="20" \
    BASE_OVERRIDE=spe_runs/gate101_${sp} EXTRA_ARGS="--no-phsp" NO_FLUENCE=1 JOBS="${JOBS}" \
        bash run_materials_spe.sh
done
echo "Finished: $(date)"
deactivate
