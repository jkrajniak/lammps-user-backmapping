#!/usr/bin/env bash
# Full Tier B protocol in MPI for all three systems.
# Records: completion status, wall time, final thermo, any ERROR.
# Outputs to /home/azureuser/sc/mpi-full/.
set -u

LMP=/home/azureuser/sc/lammps/build-mpi/lmp
OUT=/home/azureuser/sc/mpi-full
mkdir -p "$OUT"
SUMMARY="$OUT/full_mpi_summary.txt"
: > "$SUMMARY"

run_case() {
  local name="$1" dir="$2" inp="$3" np="$4" tag="$5"
  local log="$OUT/${name}_${tag}.lammps"
  echo "=== $name np=$np ($tag) started $(date +%T) ===" | tee -a "$SUMMARY"
  if [[ "$np" -eq 1 ]]; then
    (cd "$dir" && "$LMP" -in "$inp" -log "$log" >"$OUT/${name}_${tag}.out" 2>&1)
  else
    (cd "$dir" && mpirun -np "$np" "$LMP" -in "$inp" -log "$log" >"$OUT/${name}_${tag}.out" 2>&1)
  fi
  local rc=$?
  local wall errors final_t final_e
  wall=$(grep "Total wall time" "$log" 2>/dev/null | sed 's/Total wall time://')
  errors=$(grep -c "^ERROR" "$log" 2>/dev/null)
  final_t=$(grep -E "^ +[0-9]+ " "$log" 2>/dev/null | tail -1 | awk '{print $2}')
  final_e=$(grep -E "^ +[0-9]+ " "$log" 2>/dev/null | tail -1 | awk '{print $5}')
  echo "  rc=$rc wall=[$wall] errors=$errors final_temp=$final_t final_etotal=$final_e" | tee -a "$SUMMARY"
}

# Dodecane full Tier B (in.dodecane): 1, 4, 8 ranks
run_case dodecane /home/azureuser/sc/pe-mpi-test in.dodecane 1 serial
run_case dodecane /home/azureuser/sc/pe-mpi-test in.dodecane 4 np4
run_case dodecane /home/azureuser/sc/pe-mpi-test in.dodecane 8 np8

# PE full Tier B (in.pe_robust): 1, 4 ranks
run_case pe /home/azureuser/sc/pe-mpi in.pe_robust 1 serial
run_case pe /home/azureuser/sc/pe-mpi in.pe_robust 4 np4

# rim135 full Tier B (in.rim135): 1, 4 ranks
run_case rim135 /home/azureuser/sc/rim135 in.rim135 1 serial
run_case rim135 /home/azureuser/sc/rim135 in.rim135 4 np4

echo "=== ALL DONE $(date +%T) ===" | tee -a "$SUMMARY"
