#!/usr/bin/env bash
# Run ORCA inputs from calculations/ions/ one after another on this machine.
#
# Each job runs in its own folder under $ORCA_WORK (default ~/orca_runs), so ORCA's
# temporary files stay out of the repository. When a job finishes, its .out file is
# copied back next to the input. Jobs run on one core (no MPI needed).
#
# Usage, from the repository root:
#   scripts/run_orca_local.sh PF6 HSO4 C3F7CO2          # ion file names without .inp
#   ORCA=~/orca_6_1_1/orca scripts/run_orca_local.sh DCA
set -u
# macOS default stack (8 MB) makes ORCA 6.1.1 crash in the analytic r2SCAN frequency
# step (segfault in orca_prop, XC Hessian). The macOS maximum (64 MB) fixes it.
ulimit -s 65520 2>/dev/null || ulimit -s hard
ORCA="${ORCA:-$HOME/orca_6_1_1/orca}"
WORK="${ORCA_WORK:-$HOME/orca_runs}"
REPO="$(cd "$(dirname "$0")/.." && pwd)"
IONS="$REPO/calculations/ions"

for name in "$@"; do
  inp="$IONS/$name.inp"
  if [[ ! -f "$inp" ]]; then echo "[skip] no input $inp"; continue; fi
  dir="$WORK/$name"
  mkdir -p "$dir"
  cp "$inp" "$dir/"
  start=$(date +%s)
  echo "[run ] $name  ($(date '+%H:%M:%S'))"
  (cd "$dir" && "$ORCA" "$name.inp" > "$name.out" 2> "$name.err")
  secs=$(( $(date +%s) - start ))
  if grep -q "ORCA TERMINATED NORMALLY" "$dir/$name.out"; then
    cp "$dir/$name.out" "$IONS/$name.out"
    echo "[done] $name  ${secs}s"
  else
    echo "[FAIL] $name  ${secs}s  see $dir/$name.out and $dir/$name.err"
  fi
done
