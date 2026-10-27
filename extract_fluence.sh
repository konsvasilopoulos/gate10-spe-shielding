#!/usr/bin/env bash
# extract_fluence.sh - For each configuration in a folder (e.g. spe_runs/materials) runs analyze_spectrum.py,
# which writes spe_fluence_plane.csv (fluence of all species at the water entrance), spe_dose_depth.csv
# and spe_summary.png. Once the fluence is written, the phsp_*.root files are moved to _to_delete/ (large).
# Usage:  bash extract_fluence.sh spe_runs/materials
set -u
BASE="${1:-spe_runs/materials}"
for d in "${BASE}"/*_*/; do
    d="${d%/}"; [ "$(basename "$d")" = "logs" ] && continue
    ls "$d"/g*/phsp_*.root >/dev/null 2>&1 || { echo "[--] $(basename "$d"): no phase-space files"; continue; }
    rm -f "$d/spe_fluence_plane.csv"
    python analyze_spectrum.py "$d" > "$d/analyze_spectrum.txt" 2>&1
    if [ -s "$d/spe_fluence_plane.csv" ]; then
        for f in "$d"/g*/phsp_*.root; do mkdir -p "_to_delete/$(dirname "$f")"; mv "$f" "_to_delete/$f"; done
        echo "[ok] $(basename "$d"): fluence saved, phsp -> _to_delete/"
    else
        echo "[!!] $(basename "$d"): fluence not written (see $d/analyze_spectrum.txt) - phsp kept"
    fi
done
