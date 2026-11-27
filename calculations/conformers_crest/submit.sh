#!/bin/bash
# Roadmap step 3.1: CREST conformer search (iMTD-GC, GFN2-xTB, gas phase) for all 83 computed ions.
# Inputs from scripts/make_conformer_search.py: inputs/<ion>.xyz and jobs.txt (ion and CREST flags,
# largest ion first, so the long searches start at once). CREST 3.0.2, official static build, in
# ~/software/crest-3.0.2 (downloaded from github.com/crest-lab/crest, release v3.0.2).
#
# 4 ions at a time, 4 threads each. Kept per ion in outputs/:
#   <ion>.crest.out          CREST log
#   <ion>.best.xyz           lowest conformer (crest_best.xyz)
#   <ion>.conformers.xyz     all conformers within 6 kcal/mol (crest_conformers.xyz)
# Finished ions are skipped, so if the job hits its time limit, submit it again. Partition
# medium-small (3 days): the first job (576149) on short (23 h) finished only 4 of 83 ions, and an
# ion cut off by the time limit starts again from zero.
#
# Usage on the login node, in this folder (it must contain inputs/ and jobs.txt):
#   sbatch submit.sh
#
#SBATCH --job-name=crest
#SBATCH --partition=medium-small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --output=%x_%j.log

export CREST="$HOME/software/crest-3.0.2/crest/crest"
export THREADS=4
export SUBMIT="$SLURM_SUBMIT_DIR"
export WORK="$HOME/_scratch/crest_${SLURM_JOB_ID}"
export OMP_STACKSIZE=1G
ulimit -s unlimited
mkdir -p "$WORK" "$SUBMIT/outputs"
echo "host=$(hostname) cpus=$SLURM_CPUS_PER_TASK crest=$("$CREST" --version 2>&1 | grep -m1 -o 'Version [0-9.]*') start=$(date)"

run_one() {
  local name="$1"; shift
  local flags="$*"
  local dir="$WORK/$name"
  [ -s "$SUBMIT/outputs/$name.best.xyz" ] && grep -q "CREST terminated normally" "$SUBMIT/outputs/$name.crest.out" && return 0
  rm -rf "$dir"; mkdir -p "$dir" && cp "$SUBMIT/inputs/$name.xyz" "$dir/start.xyz" && cd "$dir" || return 1
  local t0=$(date +%s)
  "$CREST" start.xyz $flags --T "$THREADS" > crest.out 2>&1
  local t=$(( $(date +%s) - t0 ))
  cp crest.out "$SUBMIT/outputs/$name.crest.out"
  [ -s crest_best.xyz ] && cp crest_best.xyz "$SUBMIT/outputs/$name.best.xyz"
  [ -s crest_conformers.xyz ] && cp crest_conformers.xyz "$SUBMIT/outputs/$name.conformers.xyz"
  if grep -q "CREST terminated normally" crest.out && [ -s crest_best.xyz ]; then
    echo "[done] $name ${t}s conformers=$(grep -c '^ *-\?[0-9]*\.[0-9]*$' crest_conformers.xyz 2>/dev/null) flags=$flags"
    cd "$WORK" && rm -rf "$dir"
  else
    echo "[FAIL] $name ${t}s flags=$flags"
  fi
}
export -f run_one

# jobs.txt: "<ion>\t<flags>", largest first
tr '\t' ' ' < "$SUBMIT/jobs.txt" \
  | xargs -P $(( SLURM_CPUS_PER_TASK / THREADS )) -L 1 bash -c 'run_one "$@"' _

echo "finished $(date)"
echo "ions done: $(ls "$SUBMIT"/outputs/*.best.xyz 2>/dev/null | wc -l) of $(wc -l < "$SUBMIT/jobs.txt")"
