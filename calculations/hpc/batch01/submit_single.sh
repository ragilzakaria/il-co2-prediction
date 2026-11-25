#!/bin/bash
# Batch 01 as ONE SLURM job (for accounts whose submit limit rejects job arrays):
# the five serial ORCA runs listed in jobs.txt run side by side, one core each.
# See submit.sh for what each run is.
#
# Usage on the login node, in this folder:
#   sbatch submit_single.sh
#
#SBATCH --job-name=batch01
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=30G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.log

module purge
module load nuclear/orca/6.1.1
ORCA="$(command -v orca)"
echo "host=$(hostname) orca=$ORCA start=$(date)"

run_one() {
  local name="$1"
  local work="$HOME/_scratch/batch01_${SLURM_JOB_ID}/$name"
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
