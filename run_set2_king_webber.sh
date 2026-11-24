#!/usr/bin/env bash
# run_set2_king_webber.sh - SET 2 (after the spectrum-sampling fix of 22/10/2026):
#   Materials with King 1972 and Webber 1956 (7 materials x 5/10/20 g/cm2; LiH and water: run_set7_lih_water.sh)
# Earlier results are NOT deleted: they are moved to spe_runs/_old_tilted/.
# Usage:  bash run_set2_king_webber.sh 2>&1 | tee spe_runs/set2.log    (~8-10 h, keeps the Mac awake)
# If interrupted, run the same command again; it resumes where it stopped.
set -u
if command -v caffeinate >/dev/null 2>&1; then caffeinate -i -w $$ & fi
mkdir -p spe_runs/_old_tilted
for d in materials materials_king1972 materials_king1972_20hiN materials_EMonly_king1972 materials_EMonly_webber1956 pb_skin_test; do
    if [ -d "spe_runs/$d" ] && ! grep -q weights_as_density "spe_runs/$d"/*/*/setup.json 2>/dev/null; then
        mv "spe_runs/$d" "spe_runs/_old_tilted/$d" && echo "[moved] $d -> spe_runs/_old_tilted/"
    fi
done
for sp in king1972 webber1956; do
    echo "===== Materials ${sp} | $(date) ====="
    SPECTRUM=$sp bash run_materials_spe.sh
done
echo "===== End of set 2: $(date) ====="
