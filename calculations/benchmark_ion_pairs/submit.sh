#!/bin/bash
# Ion-pair benchmark on Mahameru: 528 single-point ORCA runs (see scripts/make_ionpair_benchmark.py),
# run as many serial ORCA processes side by side inside ONE job (the student account allows
# 2 running / 4 submitted jobs). Small inputs and r2SCAN-3c run first, so an input problem shows
# up within minutes in the log.
#
# Usage on the login node, in this folder (it must contain inputs/):
#   sbatch submit.sh
#
#SBATCH --job-name=ipbench
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=110G
#SBATCH --time=23:00:00
#SBATCH --output=%x_%j.log

module purge
module load nuclear/orca/6.1.1
export ORCA="$(command -v orca)"
export SUBMIT="$SLURM_SUBMIT_DIR"
export WORK="$HOME/_scratch/ipbench_${SLURM_JOB_ID}"
mkdir -p "$WORK" "$SUBMIT/outputs"
echo "host=$(hostname) cores=$SLURM_CPUS_PER_TASK orca=$ORCA start=$(date)"

run_one() {
  local name="$1"
  local dir="$WORK/$name"
  [ -s "$SUBMIT/outputs/$name.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$SUBMIT/outputs/$name.out" && return 0
  mkdir -p "$dir" && cp "$SUBMIT/inputs/$name.inp" "$dir/" && cd "$dir" || return 1
  "$ORCA" "$name.inp" > "$name.out" 2> "$name.err"
  cp "$name.out" "$SUBMIT/outputs/"
  if grep -q "ORCA TERMINATED NORMALLY" "$name.out"; then echo "[done] $name"; rm -rf "$dir"; else echo "[FAIL] $name"; fi
}
export -f run_one

# Order: r2SCAN-3c before wB97X-D3, smaller files first.
cd "$SUBMIT/inputs" || exit 1
{ ls -S -r *__r2scan3c__*.inp; ls -S -r *__wb97xd3__*.inp; } | sed 's/\.inp$//' \
  | xargs -P "$SLURM_CPUS_PER_TASK" -I{} bash -c 'run_one "$@"' _ {}

echo "finished $(date)"
echo "normal terminations: $(grep -l 'ORCA TERMINATED NORMALLY' "$SUBMIT"/outputs/*.out | wc -l) of $(ls "$SUBMIT"/inputs/*.inp | wc -l)"
