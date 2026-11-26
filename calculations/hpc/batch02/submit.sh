#!/bin/bash
# Batch 02: analytic r2SCAN-3c frequency checks for [N1,8,8,8]+, [P6,6,6,14]+ and
# [oleyl-MIM]+ (80, 101 and 65 atoms), whose laptop optimisations converged but whose
# frequency step crashed. The three serial ORCA runs in jobs.txt run side by side, one core
# each. Batch 01 (53 atoms) took 1 h; the 101-atom run may take most of a day, hence 23 h.
#
# Usage on the login node, in this folder:
#   sbatch submit.sh
#
#SBATCH --job-name=batch02
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=28G
#SBATCH --time=23:00:00
#SBATCH --output=%x_%j.log

module purge
module load nuclear/orca/6.1.1
ORCA="$(command -v orca)"
echo "host=$(hostname) orca=$ORCA start=$(date)"

run_one() {
  local name="$1"
  local work="$HOME/_scratch/batch02_${SLURM_JOB_ID}/$name"
  mkdir -p "$work"
  cp "$SLURM_SUBMIT_DIR/$name.inp" "$work/"
  ( cd "$work" && "$ORCA" "$name.inp" > "$name.out" )
  cp "$work/$name.out" "$SLURM_SUBMIT_DIR/"
  if grep -q "ORCA TERMINATED NORMALLY" "$work/$name.out"; then
    echo "[done] $name $(date '+%H:%M')"
  else
    echo "[FAIL] $name $(date '+%H:%M')"
  fi
}

while read -r name; do
  [ -n "$name" ] && run_one "$name" &
done < "$SLURM_SUBMIT_DIR/jobs.txt"
wait
echo "all finished $(date)"
