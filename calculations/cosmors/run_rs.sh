#!/bin/bash
# Roadmap step 3.4: the openCOSMORS step for every input in rs/ (scripts/cosmors_henry.py prepare).
# Each run takes well under a second, so this runs on the Mahameru login node, not as a job:
#   rsync -a rs/ mzakariah@login2.hpc.brin.go.id:cosmors/rs/
#   ssh mzakariah@login2.hpc.brin.go.id 'bash -s' < run_rs.sh
#   rsync -a mzakariah@login2.hpc.brin.go.id:cosmors/rs/'*_out.json' rs/
# then: python scripts/cosmors_henry.py results
BIN=/mgpfs/apps/nuclear/apps/orca/6.1.1/openCOSMORS
cd ~/cosmors/rs || exit 1
n=0
for j in *.json; do
  case "$j" in *_out.json) continue;; esac
  "$BIN" "$j" > "${j%.json}.log" 2>&1 || echo "[FAIL] $j"
  n=$((n + 1))
done
echo "$n openCOSMORS runs; warnings:"
grep -L '"warnings": *\[\]' ./*_out.json 2>/dev/null || echo "  none"
