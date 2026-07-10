#!/usr/bin/env bash
# PE Tier C: extract AT, build reference, run RDF production, compare.
set -euo pipefail

LMP="${LMP:-/Users/jakubkrajniak/Work/Science/lammps/build/lmp}"
DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR"

echo "=== Extract AT from hybrid ==="
uv run extract_at_frame.py pe_hybrid.data pe_at.data --cg-types 1 2

echo "=== Build independent AT reference ==="
uv run build_at_reference.py --gro pe_single.gro --output pe_at_ref.data --n-mol 10

echo "=== Reference RDF production ==="
"$LMP" -in in.pe_at_ref -log log.pe_at_ref.lammps

echo "=== Backmapped AT RDF production ==="
"$LMP" -in in.pe_at -log log.pe_at.lammps

echo "=== Compare RDFs ==="
uv run compare_rdf_blocks.py \
  --backmap rdf_backmap_long.dat \
  --reference rdf_reference_long.dat \
  --plot rdf_comparison_long.png \
  --report rdf_comparison_long.txt

echo "Done. See rdf_comparison_long.txt"
