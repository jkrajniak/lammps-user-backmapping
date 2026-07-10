#!/usr/bin/env bash
# PE Tier C on Azure VM: reference + backmapped AT RDF production (1 ns each).
# Compare step: run compare_rdf_blocks.py locally after rsync (needs numpy).
set -euo pipefail

DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR"

LMP="${LMP:-/home/azureuser/sc/lammps/build-omp/lmp}"
OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export OMP_NUM_THREADS
LMP_OMP=(-sf omp -pk omp "$OMP_NUM_THREADS")
PYTHON="${PYTHON:-python3}"

echo "=== Extract AT from hybrid ==="
"$PYTHON" extract_at_frame.py pe_hybrid.data pe_at.data --cg-types 1 2

echo "=== Build independent AT reference ==="
if "$PYTHON" -c "import numpy" 2>/dev/null; then
  "$PYTHON" build_at_reference.py --gro pe_single.gro --output pe_at_ref.data --n-mol 10
else
  echo "SKIP build_at_reference (no numpy) — using pre-synced pe_at_ref.data"
fi

echo "=== Reference RDF production (LMP=$LMP OMP=$OMP_NUM_THREADS) ==="
"$LMP" "${LMP_OMP[@]}" -in in.pe_at_ref -log log.pe_at_ref.lammps 2>&1 | tee log.pe_at_ref.out

echo "=== Backmapped AT RDF production ==="
"$LMP" "${LMP_OMP[@]}" -in in.pe_at -log log.pe_at.lammps 2>&1 | tee log.pe_at.out

echo "=== Tier C LAMMPS complete ==="
echo "Rsync rdf_reference_long.dat and rdf_backmap_long.dat to local, then:"
echo "  uv run compare_rdf_blocks.py --backmap rdf_backmap_long.dat --reference rdf_reference_long.dat ..."
