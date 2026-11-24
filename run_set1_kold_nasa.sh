#!/usr/bin/env bash
# run_set1_kold_nasa.sh - SET 1 (after the spectrum-sampling fix of 22/10/2026):
#   1. Materials with Koldobskiy 1972 and 1956 (7 materials x 5/10/20 g/cm2; LiH and water: run_set7_lih_water.sh)
#   2. NASA benchmark: Webber 1956, Al 20 g/cm2 + water, D and H
# Earlier results are NOT deleted: they are moved to spe_runs/_old_tilted/.
# Usage:  bash run_set1_kold_nasa.sh 2>&1 | tee spe_runs/set1.log     (~7-10 h, keeps the Mac awake)
# If interrupted, run the same command again; it resumes where it stopped.
set -u
if command -v caffeinate >/dev/null 2>&1; then caffeinate -i -w $$ & fi
mkdir -p spe_runs/_old_tilted
for d in webber1956_Al20_wide webber1956_Al20_wide_H; do
    if [ -d "spe_runs/$d" ] && ! grep -q weights_as_density "spe_runs/$d"/*/setup.json 2>/dev/null; then
        mv "spe_runs/$d" "spe_runs/_old_tilted/$d" && echo "[moved] $d -> spe_runs/_old_tilted/"
    fi
done
for sp in kold1972 kold1956; do
    echo "===== Materials ${sp} | $(date) ====="
    SPECTRUM=$sp bash run_materials_spe.sh
done
echo "===== Benchmark NASA · $(date) ====="
bash run_dose_eq.sh
python analyze_spectrum.py spe_runs/webber1956_Al20_wide_H > spe_runs/webber1956_Al20_wide_H/analyze_spectrum.txt 2>&1
python compare_nasa.py spe_runs/webber1956_Al20_wide_H 2>&1 | tail -30
echo "===== End of set 1: $(date) ====="
