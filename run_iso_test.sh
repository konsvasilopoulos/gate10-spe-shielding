#!/usr/bin/env bash
# run_iso_test.sh - Isotropic-incidence test:
#   Al, PE, Pb at 20 g/cm2, Koldobskiy 1972 and 1956, energy groups >= 120 MeV (only these cross 20 g/cm2).
#   Angles: 0 deg (reference, same geometry) and the 4 equal-probability angles of an isotropic field
#   crossing a plane (cosine law, p(mu) = 2 mu): mu = sqrt((k+0.5)/4) -> 69.3, 52.2, 37.8, 20.7 deg.
#   Isotropic = mean of the 4 angles. D only (no LET actors, no phase space) for speed.
#   Slab and water 600x600 cm so that the oblique paths fit.
#   Output: spe_runs/iso_<spectrum>/<material>_20_t<angle>/g<emin>-<emax>/     Then: python3 analysis/iso_test.py
# Usage:  bash run_iso_test.sh 2>&1 | tee spe_runs/iso_test.log      (~1-1.5 h, keeps the Mac awake)
# If interrupted, run the same command again; it resumes where it stopped.
# Keeps the Mac awake while the script runs (also with the lid closed); stops by itself at the end.
if command -v caffeinate >/dev/null 2>&1; then caffeinate -i -w $$ & fi
JOBS="${JOBS:-4}"
NSCALE="${NSCALE:-0.3}"
SPECTRA="${SPECTRA:-kold1972 kold1956}"
ANGLES="0 20.7 37.8 52.2 69.3"
GROUPS_N="1000:5000:600000 400:1000:1000000 250:400:1000000 170:250:1000000 120:170:1000000"
export NSCALE

run_one() {   # $1 spectrum, $2 material, $3 angle, $4 emin:emax:n
    sp="$1"; mat="$2"; th="$3"; IFS=: read -r emin emax n <<< "$4"
    n=$(python -c "print(int(${n} * ${NSCALE}))")
    tag=$(echo "${mat}" | sed 's/G4_//')
    out="spe_runs/iso_${sp}/${tag}_20_t${th}/g${emin}-${emax}"
    log="spe_runs/iso_${sp}/logs/${tag}_t${th}_g${emin}-${emax}.log"
    mkdir -p "spe_runs/iso_${sp}/logs"
    if [ -f "${out}/setup.json" ]; then echo "[skip] ${sp} ${tag} ${th} deg ${emin}-${emax}"; return 0; fi
    t0=$(date +%s)
    python run_spectrum.py "${mat}" --areal 20 --emin "${emin}" --emax "${emax}" --n "${n}" --seed 1234 --spectrum "${sp}" \
        --source-xy 80 --slab-xy 600 --water-xy 600 --theta-deg "${th}" --no-phsp --out "${out}" > "${log}" 2>&1
    if [ -f "${out}/setup.json" ]; then echo "[ok]    ${sp} ${tag} ${th} deg ${emin}-${emax} MeV ($(( $(date +%s) - t0 )) s)  $(date '+%H:%M')"
    else echo "[FAIL]  ${sp} ${tag} ${th} deg ${emin}-${emax} -> ${log}"; fi
}
export -f run_one
echo "===== Isotropic-incidence test | $(date) ====="
for sp in ${SPECTRA}; do
  for g in ${GROUPS_N}; do for th in ${ANGLES}; do for m in G4_Al G4_POLYETHYLENE G4_Pb; do echo "$sp $m $th $g"; done; done; done
done | xargs -P "${JOBS}" -L 1 bash -c 'run_one "$@"' _
echo "===== Finished: $(date) ====="
python3 analysis/iso_test.py
