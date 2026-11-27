#!/bin/bash
# Roadmap step 3.5: ion pairs (110) and anion-CO2 complexes (23). Starts from
# scripts/make_ion_pairs.py (starts/, systems.tsv); steps between stages in pipeline.py.
#
#   Stage A: every start (about 22 per system) optimised with GFN2-xTB (ORCA native xTB, Fermi smearing 1000 K)
#   Stage B: the lowest 3 distinct structures within 3 kcal/mol re-optimised with r2SCAN-3c
#   Stage C: wB97X-D3/def2-TZVPD single points of the lowest complex and of both ions at the
#            complex geometry -> interaction energy (results.tsv)
# Structures in which a proton moved between the ions, or an ion broke up, are dropped at each
# stage and reported in the log. Serial ORCA runs, 16 side by side, largest first.
# Finished runs are skipped, so if the job hits its time limit, submit it again.
#
# Usage on the login node, in this folder (it must contain starts/, systems.tsv, pipeline.py):
#   sbatch submit.sh
#
#SBATCH --job-name=ionpairs
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=23:00:00
#SBATCH --output=%x_%j.log

module purge
module load nuclear/orca/6.1.1
export ORCA="$(command -v orca)"
export SUBMIT="$SLURM_SUBMIT_DIR"
export WORK="$HOME/_scratch/ionpairs_${SLURM_JOB_ID}"
mkdir -p "$WORK"
echo "host=$(hostname) tasks=$SLURM_NTASKS orca=$ORCA python=$(python3 --version 2>&1) start=$(date)"

ok() { [ -s "$1" ] && grep -q "ORCA TERMINATED NORMALLY" "$1"; }
export -f ok

# run_one <stage> <name>: run ORCA in scratch, copy .out and .xyz to <stage>/outputs
run_one() {
  local stage="$1" name="$2" dir="$WORK/$1/$2"
  ok "$SUBMIT/$stage/outputs/$name.out" && return 0
  rm -rf "$dir"; mkdir -p "$dir" && cp "$SUBMIT/$stage/inputs/$name.inp" "$dir/" && cd "$dir" || return 1
  local t0=$(date +%s)
  "$ORCA" "$name.inp" > "$name.out" 2> "$name.err"
  local t=$(( $(date +%s) - t0 ))
  cp "$name.out" "$SUBMIT/$stage/outputs/"
  [ -s "$name.xyz" ] && cp "$name.xyz" "$SUBMIT/$stage/outputs/"
  if ok "$name.out"; then cd "$WORK" && rm -rf "$dir"; else echo "[FAIL] $stage $name ${t}s"; fi
}
export -f run_one

run_stage() {
  local stage="$1"
  mkdir -p "$SUBMIT/$stage/outputs"
  echo "== stage $stage: $(ls "$SUBMIT/$stage/inputs/"*.inp | wc -l) runs, $(date)"
  cd "$SUBMIT/$stage/inputs" || exit 1
  ls -S *.inp | sed 's/\.inp$//' | xargs -P "$SLURM_NTASKS" -I{} bash -c 'run_one "$1" "$2"' _ "$stage" {}
  echo "   done: $(grep -l 'ORCA TERMINATED NORMALLY' "$SUBMIT/$stage/outputs/"*.out 2>/dev/null | wc -l), $(date)"
  cd "$SUBMIT" || exit 1
}

cd "$SUBMIT" || exit 1
[ -d A/inputs ] || python3 pipeline.py prepare-a
run_stage A
python3 pipeline.py select-a
run_stage B
python3 pipeline.py select-b
run_stage C
python3 pipeline.py results
echo "finished $(date)"
