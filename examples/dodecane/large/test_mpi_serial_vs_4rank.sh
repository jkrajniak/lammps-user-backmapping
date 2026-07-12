#!/usr/bin/env bash
# MPI correctness test: serial vs N-rank dodecane backmap.
# Runs the short in.dodecane_mpi input serially and with mpirun -np NRANK,
# then compares the final data files (positions + per-atom lambda).
#
# Usage: ./test_mpi_serial_vs_4rank.sh [NRANK]
# Requires a real-MPI LAMMPS build (not MPI STUBS). Run on the VM.
set -euo pipefail

NRANK="${1:-4}"
DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR"

LMP="${LMP:-lmp}"
MPIRUN="${MPIRUN:-mpirun}"

echo "=== Serial run ==="
"$LMP" -in in.dodecane_mpi -log log.mpi_serial.lammps
mv -f dodecane_hybrid_mpi.data dodecane_hybrid_serial.data

echo "=== ${NRANK}-rank run ==="
"$MPIRUN" -np "$NRANK" "$LMP" -in in.dodecane_mpi -log log.mpi_${NRANK}rank.lammps
mv -f dodecane_hybrid_mpi.data dodecane_hybrid_${NRANK}rank.data

echo "=== Compare (positions within 1e-5 Å, lambda within 1e-10) ==="
uv run compare_mpi_data.py \
  --serial dodecane_hybrid_serial.data \
  --parallel dodecane_hybrid_${NRANK}rank.data \
  --pos-tol 1e-5 --lambda-tol 1e-10 \
  --report mpi_comparison_report.txt

echo "Done. See mpi_comparison_report.txt"
