#!/usr/bin/env bash
# Tier C pipeline: hybrid AT production (frozen CG) → RDF validation vs paper ref
set -euo pipefail
cd "$(dirname "$0")"
# Serial build/lmp for backmap correctness; use build-omp for faster AT-only prod later.
LMP="${LMP:-/home/azureuser/sc/lammps/build-omp/lmp}"
LMP_SUFFIX_OMP="${LMP_SUFFIX_OMP:--sf omp -pk omp ${OMP_NUM_THREADS:-8}}"
OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export OMP_NUM_THREADS
PYTHON="${PYTHON:-python3}"
REF_DIR="${REF_DIR:-ref}"

mkdir -p "$REF_DIR"
for f in ref_C_C.xvg ref_C_N.xvg ref_O_H.xvg; do
  if [[ -f "$f" && ! -f "$REF_DIR/$f" ]]; then
    mv "$f" "$REF_DIR/"
  fi
done

echo "=== Hybrid AT production from melamine_hybrid.data (500 ps) ==="
echo "LMP=$LMP OMP_NUM_THREADS=$OMP_NUM_THREADS"
"$LMP" $LMP_SUFFIX_OMP -in in.melamine_hybrid_tierc -log log.melamine.tierc.lammps 2>&1 | tee log.melamine.tierc.out

echo "=== RDF validation ==="
if "$PYTHON" -c "import numpy, matplotlib" 2>/dev/null; then
  "$PYTHON" compare_melamine_structure.py \
    --dump dump.at_prod \
    --final-data melamine_hybrid_final.data \
    --ref-dir "$REF_DIR" \
    --mode hybrid \
    --min-step 50000 \
    --plot rdf_comparison.png \
    --report structural_validation_report.txt
else
  echo "SKIP: numpy/matplotlib not available — run compare_melamine_structure.py locally after rsync"
fi

echo "=== Tier C complete ==="
