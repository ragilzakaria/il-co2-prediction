#!/bin/bash
# Ion-pair benchmark, one method per job: a smaller request than submit.sh, so the
# scheduler can start it in a gap (backfill) instead of waiting for 32 free cores.
# Finished outputs in outputs/ are skipped, so a job that hits its time limit can
# simply be submitted again and continues where it stopped.
#
# Usage on the login node, in this folder (it must contain inputs/):
#   sbatch --job-name=ipb_r2scan --cpus-per-task=8  --mem=28G --time=04:00:00 submit_part.sh r2scan3c
#   sbatch --job-name=ipb_wb97x  --cpus-per-task=16 --mem=56G --time=12:00:00 submit_part.sh wb97xd3
#
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=%x_%j.log

METHOD="${1:?give the method: r2scan3c or wb97xd3}"
module purge
module load nuclear/orca/6.1.1
export ORCA="$(command -v orca)"
export SUBMIT="$SLURM_SUBMIT_DIR"
export WORK="$HOME/_scratch/ipbench_${SLURM_JOB_ID}"
mkdir -p "$WORK" "$SUBMIT/outputs"
echo "method=$METHOD host=$(hostname) cores=$SLURM_CPUS_PER_TASK orca=$ORCA start=$(date)"

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

# Smaller files first, so an input problem shows up early in the log.
cd "$SUBMIT/inputs" || exit 1
ls -S -r *__${METHOD}__*.inp | sed 's/\.inp$//' \
  | xargs -P "$SLURM_CPUS_PER_TASK" -I{} bash -c 'run_one "$@"' _ {}

echo "finished $(date)"
echo "normal terminations ($METHOD): $(grep -l 'ORCA TERMINATED NORMALLY' "$SUBMIT"/outputs/*__${METHOD}__*.out 2>/dev/null | wc -l) of $(ls "$SUBMIT"/inputs/*__${METHOD}__*.inp | wc -l)"
