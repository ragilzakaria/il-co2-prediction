#!/bin/bash
# Batch 01 on Mahameru: the ORCA jobs left over from the laptop reruns (25 Nov 2026).
#
#   P4,4,4,4.freq  C12PhSO3.freq  HMMIM.freq   analytic frequencies on geometries already
#                                              optimised on the laptop (macOS ORCA crashes
#                                              in this step for ions above ~30 atoms)
#   DBPO4                                      full r2SCAN-3c Opt Freq (rebuilt input, 0953722)
#   C3F7CO2                                    TightOpt Freq DefGrid3 (removes a -7 cm-1 mode)
#
# One serial ORCA 6.1.1 job per line of jobs.txt, run as a SLURM job array in scratch.
# Each <name>.inp gives <name>.out, copied back to this folder.
#
# Usage on the login node, in this folder:
#   sbatch submit.sh
#
#SBATCH --job-name=batch01
#SBATCH --partition=short
#SBATCH --array=1-5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=20:00:00
#SBATCH --output=%x_%A_%a.log

NAME="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SLURM_SUBMIT_DIR/jobs.txt")"

module purge
module load nuclear/orca/6.1.1
ORCA="$(command -v orca)"

WORK="$HOME/_scratch/batch01_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$WORK"
cp "$SLURM_SUBMIT_DIR/$NAME.inp" "$WORK/"
cd "$WORK" || exit 1

echo "job=$NAME host=$(hostname) orca=$ORCA start=$(date)"
"$ORCA" "$NAME.inp" > "$NAME.out"
echo "exit=$? end=$(date)"

cp "$NAME.out" "$SLURM_SUBMIT_DIR/"
grep -q "ORCA TERMINATED NORMALLY" "$NAME.out" && echo "status: terminated normally" || echo "status: FAILED"
grep -A14 "VIBRATIONAL FREQUENCIES" "$NAME.out" | grep "cm\*\*-1" | sed -n 7,9p
