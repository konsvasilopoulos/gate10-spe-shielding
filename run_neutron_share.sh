#!/usr/bin/env bash
# =============================================================================
# run_neutron_share.sh - Share of the BFO dose due to neutrons (neutron RBE check, Supplementary Section S7).
#   Koldobskiy 1972 and 1956, Al / PE / Pb at 20 g/cm2, two runs with the SAME seed and Geant4 in MT mode
#   (one thread, one seed per event -> correlated runs):
#     full : standard physics
#     noN  : all hadronic processes of neutrons inactivated (neutrons leave without depositing energy)
#   Neutron share = 1 - D(noN) / D(full). Dose only (no LET actors, no phase space).
#   Output: spe_runs/neutron_share/<spectrum>/<material>_20_<full|noN>/g<emin>-<emax>/
#   A one-minute test first checks that the inactivation works (process dump in the log).
#   Time: ~45 min with JOBS=4 (Apple M5). Rerun to resume; groups with setup.json are skipped.
# Usage:  bash run_neutron_share.sh        (then: python analysis/neutron_share.py)
# =============================================================================
set -u
cd "$(dirname "$0")"
if command -v caffeinate >/dev/null 2>&1; then caffeinate -i -w $$ & fi
JOBS="${JOBS:-4}"
NSCALE="${NSCALE:-0.25}"
BASE="spe_runs/neutron_share"
WIDE="--source-xy 80 --slab-xy 120 --water-xy 120"
GROUPS_N="1000:5000:600000 400:1000:1000000 250:400:1000000 170:250:1000000 120:170:1000000 60:120:800000 30:60:800000"
mkdir -p "${BASE}/logs"

echo "=== Test: 150-170 MeV, Al 20 g/cm2, 2000 protons, with and without neutron interactions"
for mode in full noN; do
    flag=""; [ "${mode}" = "noN" ] && flag="--no-neutron-hadronic"
    python run_spectrum.py G4_Al --areal 20 --spectrum kold1972 --emin 150 --emax 170 --n 2000 --seed 1234 --mt --no-phsp \
        ${WIDE} ${flag} --out "${BASE}/test_${mode}" > "${BASE}/logs/test_${mode}.log" 2>&1
    [ -f "${BASE}/test_${mode}/setup.json" ] && echo "[ok] test ${mode}" || { echo "[FAIL] test ${mode} -> ${BASE}/logs/test_${mode}.log"; exit 1; }
done
grep "process\[" "${BASE}/logs/test_noN.log" | head -8

export WIDE BASE NSCALE
run_one() {   # $1 spectrum, $2 material, $3 mode, $4 emin:emax:n
    sp="$1"; mat="$2"; mode="$3"; IFS=: read -r emin emax n <<< "$4"
    n=$(python -c "print(int(${n} * ${NSCALE}))")
    tag=$(echo "${mat}" | sed 's/G4_//')
    out="${BASE}/${sp}/${tag}_20_${mode}/g${emin}-${emax}"; log="${BASE}/logs/${sp}_${tag}_${mode}_g${emin}-${emax}.log"
    if [ -f "${out}/setup.json" ]; then echo "[skip] ${sp} ${tag} ${mode} ${emin}-${emax}"; return 0; fi
    flag=""; [ "${mode}" = "noN" ] && flag="--no-neutron-hadronic"
    t0=$(date +%s)
    python run_spectrum.py "${mat}" --areal 20 --spectrum "${sp}" --emin "${emin}" --emax "${emax}" --n "${n}" --seed 1234 \
        --mt --no-phsp ${WIDE} ${flag} --out "${out}" > "${log}" 2>&1
    if [ -f "${out}/setup.json" ]; then echo "[ok]    ${sp} ${tag} ${mode} ${emin}-${emax} MeV ($(( $(date +%s) - t0 )) s)"
    else echo "[FAIL]  ${sp} ${tag} ${mode} ${emin}-${emax} -> ${log}"; fi
}
export -f run_one
for sp in kold1972 kold1956; do
  for mat in G4_Al G4_POLYETHYLENE G4_Pb; do
    for mode in full noN; do
      for g in ${GROUPS_N}; do echo "${sp} ${mat} ${mode} ${g}"; done
    done
  done
done | xargs -P "${JOBS}" -L 1 bash -c 'run_one "$@"' _
python analysis/neutron_share.py
