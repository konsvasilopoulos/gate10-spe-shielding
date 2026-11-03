#!/usr/bin/env bash
# run_release_check.sh - Check of the latest release, GATE 10.1.1 (Geant4 11.4.2), against GATE 10.0.0 (Geant4 11.2.1).
#   1) Koldobskiy 1972 SPE, Al / PE / Pb at 20 g/cm2, same settings as set 1 (wide geometry, --split --dose-eq,
#      NSCALE 0.5, seed 1234), without phase space            -> spe_runs/gate101_kold1972/
#   2) 150 MeV validation (PSTAR): air + 7 materials, seed 1234 -> runs_gate101/s1234/
#   Uses a separate virtual environment with opengate 10.1.1 (Python >= 3.10), e.g.
#      python3.13 -m venv .venv-gate101 && .venv-gate101/bin/pip install opengate==10.1.1 uproot matplotlib
#   Then: python validation/seed_stats.py runs_gate101 ; results in analysis/release_check.csv
#   About 1-1.5 h with JOBS=4; resumes if interrupted.
set -u
cd "$(dirname "$0")"
if command -v caffeinate >/dev/null 2>&1; then caffeinate -i -w $$ & fi
VENV="${VENV:-.venv-gate101}"
source "${VENV}/bin/activate"
python -c "from importlib.metadata import version; print('GATE', version('opengate'))"
JOBS="${JOBS:-4}"

echo "=== 1/2  Koldobskiy 1972, Al / PE / Pb, 20 g/cm2  ($(date '+%H:%M'))"
SPECTRUM=kold1972 MATERIALS="G4_Al G4_POLYETHYLENE G4_Pb" AREALS="20" \
BASE_OVERRIDE=spe_runs/gate101_kold1972 EXTRA_ARGS="--no-phsp" NO_FLUENCE=1 JOBS="${JOBS}" \
    bash run_materials_spe.sh

echo "=== 2/2  150 MeV validation  ($(date '+%H:%M'))"
BASE=runs_gate101/s1234
mkdir -p "${BASE}"
val_one() {
    m="$1"; out="${BASE}/output_${m}"
    if [ -f "${out}/setup.json" ]; then echo "[skip] ${out}"; return 0; fi
    python run_material.py "${m}" --seed 1234 --physics QGSP_BIC_EMZ --out "${out}" > "${BASE}/${m}.log" 2>&1 \
        && echo "[ok]   ${m}  $(date '+%H:%M')" || echo "[FAIL] ${m} -> ${BASE}/${m}.log"
}
export -f val_one; export BASE
echo "G4_AIR G4_POLYETHYLENE G4_PLEXIGLASS G4_Al G4_Ti G4_Fe G4_Cu G4_Pb" | tr ' ' '\n' | xargs -P "${JOBS}" -I{} bash -c 'val_one "$@"' _ {}
echo "Finished: $(date)"
deactivate
