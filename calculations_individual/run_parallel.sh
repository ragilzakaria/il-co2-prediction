#!/bin/bash

# --- CONFIGURATION ---
MAX_JOBS=8
ORCA_PATH="/Users/ragilzakaria/Library/orca_6_0_0/orca"
MPI_PATH=$(brew --prefix open-mpi)
GCC_PATH=$(brew --prefix gcc)
LIB_PATH="$MPI_PATH/lib:$GCC_PATH/lib/gcc/current"
# --- END OF CONFIGURATION ---


# FIX: The loop now only runs on the specific list of failed jobs.
# NOTE: Make sure these filenames exactly match what is in your folder.
for file in \
  DMIM.inp \
  Cl_m.inp \
  HMMIM.inp \
  "N1,8,8,8.inp" \
  oleyl_mMIM.inp \
  "P1,4,4,4.inp" \
  "P4,4,4,4.inp" \
  "P6,6,6,14.inp" \
; do
  # Run the ORCA command in the background
  echo "==> Re-running calculation for $file..."
  DYLD_LIBRARY_PATH="$LIB_PATH" "$ORCA_PATH" "$file" > "${file%.inp}.out" &

  # Check the number of running jobs and wait if the max is reached
  while (( $(jobs -p | wc -l) >= MAX_JOBS )); do
    sleep 10
  done
done

# Wait for all remaining background jobs to finish
wait

echo "--- RERUN COMPLETE ---"