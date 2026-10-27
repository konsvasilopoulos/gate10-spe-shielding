#!/usr/bin/env bash
# =============================================================================
# run_materials_spe.sh - Material comparison for one SPE spectrum (SPECTRUM=..., default webber1956), wide geometry,
# 7 materials x 3 areal densities, with dose per particle type and dose equivalent.
#   Output: spe_runs/materials[_<spectrum>]/<material>_<g/cm2>/g<emin>-<emax>/
#   Then:   python compare_materials_spe.py <output folder>
#
# Usage:  SPECTRUM=kold1956 JOBS=4 bash run_materials_spe.sh
#   JOBS   : parallel runs (about the number of CPU cores - 1)
#   NSCALE : event multiplier (default 0.5 -> half of the benchmark statistics,
#            ~1 % statistical uncertainty on the dose). Cost ~22 CPU-min per (material, thickness).
# If interrupted, run the same command again: it resumes (runs with a setup.json are skipped).
# =============================================================================
set -u
# Keeps the Mac awake while the script runs (also with the lid closed); stops by itself at the end.
if command -v caffeinate >/dev/null 2>&1; then caffeinate -i -w $$ & fi
JOBS="${JOBS:-4}"   # parallel runs (number of fast CPU cores)
SPECTRUM="${SPECTRUM:-webber1956}"
NSCALE="${NSCALE:-0.5}"
PHYSICS="${PHYSICS:-QGSP_BIC_EMZ}"   # e.g. G4EmStandardPhysics_option4 for electromagnetic physics only
MATERIALS="${MATERIALS:-G4_POLYETHYLENE G4_PLEXIGLASS G4_Al G4_Ti G4_Fe G4_Cu G4_Pb}"
AREALS="${AREALS:-5 10 20}"
WIDE="${WIDE:---source-xy 80 --slab-xy 120 --water-xy 120}"
if [ "${SPECTRUM}" = "webber1956" ]; then
    BASE="spe_runs/materials"
    GROUPS_N="400:2000:1000000 250:400:1000000 170:250:1000000 120:170:1000000 10:60:800000 60:120:800000"
elif [ "${SPECTRUM}" = "king1972" ]; then
    BASE="spe_runs/materials_${SPECTRUM}"
    # the 1972 King spectrum is soft: the 400-1000 MeV group has negligible fluence, so few events
    GROUPS_N="250:400:1000000 400:1000:200000 170:250:1000000 120:170:1000000 10:60:800000 60:120:800000"
else
    # Koldobskiy et al. 2021 reconstructions (kold1956, kold1972, kold1989): valid for E > 30 MeV,
    # hard spectra -> groups up to 5 GeV
    BASE="spe_runs/materials_${SPECTRUM}"
    GROUPS_N="1000:5000:600000 400:1000:1000000 250:400:1000000 170:250:1000000 120:170:1000000 60:120:800000 30:60:800000"
fi
BASE="${BASE_OVERRIDE:-${BASE}}"
GROUPS_N="${GROUPS_N_OVERRIDE:-${GROUPS_N}}"   # e.g. narrow windows for the response matrix
EXTRA_ARGS="${EXTRA_ARGS:-}"   # extra arguments for run_spectrum.py (e.g. --no-phsp)
mkdir -p "${BASE}/logs"
export WIDE NSCALE SPECTRUM BASE PHYSICS EXTRA_ARGS

run_one() {   # $1 material, $2 g/cm2, $3 emin:emax:n
    mat="$1"; ar="$2"; IFS=: read -r emin emax n <<< "$3"
    n=$(python -c "print(int(${n} * ${NSCALE}))")
    tag=$(echo "${mat}" | sed 's/G4_POLYETHYLENE/PE/g; s/G4_PLEXIGLASS/PMMA/g; s/G4_//g; s/:\([0-9.]*\)/\1/g; s/+/_/g')
    [ "${tag}" = "PE" ] && [ "${mat}" = "G4_POLYETHYLENE" ] && tag="POLYETHYLENE"
    [ "${tag}" = "PMMA" ] && [ "${mat}" = "G4_PLEXIGLASS" ] && tag="PLEXIGLASS"
    out="${BASE}/${tag}_${ar}/g${emin}-${emax}"
    log="${BASE}/logs/${tag}_${ar}_g${emin}-${emax}.log"
    if [ -f "${out}/setup.json" ]; then echo "[skip] ${mat} ${ar} ${emin}–${emax}"; return 0; fi
    t0=$(date +%s)
    python run_spectrum.py "${mat}" --areal "${ar}" --emin "${emin}" --emax "${emax}" --n "${n}" --seed 1234 --spectrum "${SPECTRUM}" --physics "${PHYSICS}" \
        --split --dose-eq ${WIDE} ${EXTRA_ARGS} --out "${out}" > "${log}" 2>&1
    if [ -f "${out}/setup.json" ]; then echo "[ok]    ${mat} ${ar} g/cm² ${emin}–${emax} MeV ($(( $(date +%s) - t0 )) s)  $(date '+%H:%M')"
    else echo "[FAIL]  ${mat} ${ar} ${emin}–${emax} -> ${log}"; echo "${mat} ${ar} ${emin}-${emax}" >> ${BASE}/failed.txt; fi
}
export -f run_one
rm -f ${BASE}/failed.txt
total=0
{
  # slow groups (high energies) first, for better load balancing
  for g in ${GROUPS_N}; do
    for ar in ${AREALS}; do for m in ${MATERIALS}; do echo "$m $ar $g"; done; done
  done
} | xargs -P "${JOBS}" -L 1 bash -c 'run_one "$@"' _
echo "Finished: $(date)"
if [ -f ${BASE}/failed.txt ]; then echo "Failed:"; cat ${BASE}/failed.txt; fi
python compare_materials_spe.py "${BASE}" | tee "${BASE}/materials_result.txt"
# plane fluence from the phase-space files, then phsp -> _to_delete/
if [ -z "${NO_FLUENCE:-}" ]; then bash extract_fluence.sh "${BASE}"; fi
