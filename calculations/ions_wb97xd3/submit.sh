#!/bin/bash
# Roadmap steps 2.2 and 2.3: wB97X-D3/def2-TZVPD single points of all 83 ions on their r2SCAN-3c
# geometries (inputs from scripts/make_ion_singlepoints.py), with a parallel-ORCA timing test first.
#
#   Part 1, parallel test (alone on the allocation, one run at a time):
#     [Tf2N], [BMIM] and [P4,4,4,4] on 1, 4 and 16 cores.
#   Part 2, all other ions: serial ORCA runs side by side, largest first.
#
# Finished outputs in outputs/ are skipped, so if the job hits its time limit, submit it again.
#
# Usage on the login node, in this folder (it must contain inputs/):
#   sbatch submit.sh
#
#SBATCH --job-name=ionsp
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=23:00:00
#SBATCH --output=%x_%j.log

module purge
module load nuclear/orca/6.1.1
module load openmpi4/4.1.4          # needed for %pal (parallel ORCA)
export ORCA="$(command -v orca)"    # parallel ORCA must be called with its full path
export SUBMIT="$SLURM_SUBMIT_DIR"
export WORK="$HOME/_scratch/ionsp_${SLURM_JOB_ID}"
mkdir -p "$WORK" "$SUBMIT/outputs"
echo "host=$(hostname) tasks=$SLURM_NTASKS orca=$ORCA mpirun=$(command -v mpirun) start=$(date)"

run_one() {
  local name="$1"
  local dir="$WORK/$name"
  [ -s "$SUBMIT/outputs/$name.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$SUBMIT/outputs/$name.out" && return 0
  mkdir -p "$dir" && cp "$SUBMIT/inputs/$name.inp" "$dir/" && cd "$dir" || return 1
  local t0=$(date +%s)
  "$ORCA" "$name.inp" > "$name.out" 2> "$name.err"
  local t=$(( $(date +%s) - t0 ))
  cp "$name.out" "$SUBMIT/outputs/"
  [ -s "$name.err" ] && cp "$name.err" "$SUBMIT/outputs/"
  if grep -q "ORCA TERMINATED NORMALLY" "$name.out"; then echo "[done] $name ${t}s"; rm -rf "$dir"; else echo "[FAIL] $name ${t}s"; fi
}
export -f run_one

echo "== part 1: parallel test $(date)"
for ion in Tf2N BMIM "P4,4,4,4"; do
  for suffix in "" __pal4 __pal16; do
    run_one "${ion}${suffix}"
  done
done

echo "== part 2: all ions, serial, $SLURM_NTASKS at a time $(date)"
cd "$SUBMIT/inputs" || exit 1
ls -S *.inp | grep -v "__pal" | sed 's/\.inp$//' \
  | xargs -P "$SLURM_NTASKS" -I{} bash -c 'run_one "$@"' _ {}

echo "finished $(date)"
echo "normal terminations: $(grep -l 'ORCA TERMINATED NORMALLY' "$SUBMIT"/outputs/*.out | wc -l) of $(ls "$SUBMIT"/inputs/*.inp | wc -l)"
