#!/bin/bash
# Roadmap step 3.2: r2SCAN-3c re-optimisation of the lowest CREST conformers of every ion, then
# frequencies of the lowest one (inputs from scripts/make_conformer_reopt.py).
#
#   Stage 1: every candidate (<ion>__c0 = current structure, c1..c5 = CREST conformers) is
#            optimised; serial ORCA runs, 16 side by side, largest first.
#   Stage 2: per ion, the lowest converged candidate is kept (written to selected.tsv) and its
#            frequencies are computed on 4 cores, 4 ions at a time.
#
# Kept in outputs/: <ion>__c<k>.out and .xyz (optimisations), <ion>.best.xyz, <ion>.freq.out.
# Finished runs are skipped, so if the job hits its time limit, submit it again.
#
# Usage on the login node, in this folder (it must contain inputs/, freq_templates/ and ions.txt):
#   sbatch submit.sh
#
#SBATCH --job-name=reopt
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=23:00:00
#SBATCH --output=%x_%j.log

module purge
module load nuclear/orca/6.1.1
module load openmpi4/4.1.4          # stage 2 runs ORCA in parallel (%pal)
export ORCA="$(command -v orca)"    # parallel ORCA must be called with its full path
export SUBMIT="$SLURM_SUBMIT_DIR"
export WORK="$HOME/_scratch/reopt_${SLURM_JOB_ID}"
mkdir -p "$WORK" "$SUBMIT/outputs"
echo "host=$(hostname) tasks=$SLURM_NTASKS orca=$ORCA start=$(date)"

ok() { [ -s "$1" ] && grep -q "ORCA TERMINATED NORMALLY" "$1"; }
export -f ok

# run_one <name> <input file>: run ORCA in a scratch folder, copy .out and .xyz back
run_one() {
  local name="$1" inp="$2" dir="$WORK/$1"
  ok "$SUBMIT/outputs/$name.out" && return 0
  rm -rf "$dir"; mkdir -p "$dir" && cp "$inp" "$dir/$name.inp" || return 1
  [ -n "$3" ] && cp "$3" "$dir/"
  cd "$dir" || return 1
  local t0=$(date +%s)
  "$ORCA" "$name.inp" > "$name.out" 2> "$name.err"
  local t=$(( $(date +%s) - t0 ))
  cp "$name.out" "$SUBMIT/outputs/"
  [ -s "$name.xyz" ] && cp "$name.xyz" "$SUBMIT/outputs/"
  [ -s "$name.err" ] && cp "$name.err" "$SUBMIT/outputs/"
  if ok "$name.out"; then echo "[done] $name ${t}s"; cd "$WORK" && rm -rf "$dir"; else echo "[FAIL] $name ${t}s"; fi
}
export -f run_one

echo "== stage 1: optimisations $(date)"
cd "$SUBMIT/inputs" || exit 1
ls -S *.inp | sed 's/\.inp$//' \
  | xargs -P "$SLURM_NTASKS" -I{} bash -c 'run_one "$1" "$SUBMIT/inputs/$1.inp"' _ {}

echo "== selection $(date)"
cd "$SUBMIT/outputs" || exit 1
printf "ion\tbest\tE_best_Eh\tE_c0_Eh\tbest_minus_c0_kcal\tn_converged\tn_candidates\n" > "$SUBMIT/selected.tsv"
while IFS=$'\t' read -r ion charge natoms; do
  best=""; ebest=""; e0=""; nconv=0; ncand=0
  for f in "$ion"__c*.out; do
    [ -e "$f" ] || continue
    ncand=$((ncand + 1))
    ok "$f" && grep -q "THE OPTIMIZATION HAS CONVERGED" "$f" || continue
    nconv=$((nconv + 1))
    e=$(grep "FINAL SINGLE POINT ENERGY" "$f" | tail -1 | awk '{print $NF}')
    tag=${f%.out}; tag=${tag##*__}
    [ "$tag" = c0 ] && e0=$e
    if [ -z "$ebest" ] || awk -v a="$e" -v b="$ebest" 'BEGIN{exit !(a < b)}'; then best=$tag; ebest=$e; fi
  done
  if [ -z "$best" ]; then echo "[NOSEL] $ion: no converged optimisation"; continue; fi
  cp "${ion}__${best}.xyz" "$ion.best.xyz"
  d=$(awk -v a="$ebest" -v b="$e0" 'BEGIN{if (b == "") print "NA"; else printf "%.2f", (a - b) * 627.509}')
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$ion" "$best" "$ebest" "$e0" "$d" "$nconv" "$ncand" >> "$SUBMIT/selected.tsv"
done < "$SUBMIT/ions.txt"
column -t "$SUBMIT/selected.tsv"

echo "== stage 2: frequencies of the selected conformers $(date)"
cut -f1 "$SUBMIT/ions.txt" | while read -r ion; do [ -s "$SUBMIT/outputs/$ion.best.xyz" ] && echo "$ion"; done \
  | xargs -P 4 -I{} bash -c 'run_one "$1.freq" "$SUBMIT/freq_templates/$1.freq.inp" "$SUBMIT/outputs/$1.best.xyz"' _ {}

echo "finished $(date)"
echo "optimisations done: $(grep -l 'THE OPTIMIZATION HAS CONVERGED' "$SUBMIT"/outputs/*__c*.out 2>/dev/null | wc -l) of $(ls "$SUBMIT"/inputs/*.inp | wc -l)"
echo "frequencies done:   $(grep -l 'ORCA TERMINATED NORMALLY' "$SUBMIT"/outputs/*.freq.out 2>/dev/null | wc -l) of $(wc -l < "$SUBMIT/ions.txt")"
