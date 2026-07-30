#!/usr/bin/env bash
# Speedup benchmark: dodecane + rim135 at 1,2,4,8 ranks.
# Records wall time of the timed ramp block; writes speedup.csv.
set -euo pipefail

LMP="${LMP:-/path/to/lammps/build-mpi/lmp}"
OUT="${OUT:-./mpi-speedup}"
mkdir -p "$OUT"
CSV="$OUT/speedup.csv"
echo "system,ranks,wall_s,atoms" > "$CSV"

run_one() {
  local sys="$1" dir="$2" inp="$3" atoms="$4" np="$5"
  local log="$OUT/${sys}_${np}rank.lammps"
  if [[ "$np" -eq 1 ]]; then
    (cd "$dir" && "$LMP" -in "$inp" -log "$log" >/dev/null 2>&1)
  else
    (cd "$dir" && mpirun -np "$np" "$LMP" -in "$inp" -log "$log" >/dev/null 2>&1)
  fi
  # Loop time of the timed 500-step run block (last "Loop time" in log)
  local t
  t=$(grep "Loop time" "$log" | tail -1 | awk "{print \$3}")
  echo "$sys,$np,$t,$atoms" | tee -a "$CSV"
}

DOD_DIR="${DOD_DIR:-/path/to/pe-mpi-test}"
RIM_DIR="${RIM_DIR:-/path/to/rim135}"

for np in 1 2 4 8; do
  echo "--- dodecane np=$np ---"
  run_one dodecane "$DOD_DIR" in.dodecane_bench 9000 "$np"
done

for np in 1 2 4 8; do
  echo "--- rim135 np=$np ---"
  run_one rim135 "$RIM_DIR" in.rim135_bench 13653 "$np"
done

echo "Done. CSV at $CSV"
cat "$CSV"
