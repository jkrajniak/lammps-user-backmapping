#!/usr/bin/env bash
# PE Tier C extended: 5 ns RDF production from equilibrated finals (both legs).
set -euo pipefail
cd "$(dirname "$0")"

LMP="${LMP:-/path/to/lammps/build-omp/lmp}"
OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export OMP_NUM_THREADS
LMP_OMP=(-sf omp -pk omp "$OMP_NUM_THREADS")

DIH_LINE="1 0.717018 -2.217713 2.905285 3.135783 -0.731287 -6.271589"

fix_dihedral_coeffs() {
  local f="$1"
  python3 - "$f" "$DIH_LINE" <<'PY'
import re, sys
path, line = sys.argv[1], sys.argv[2]
text = open(path).read()
pat = r"(Dihedral Coeffs #[^\n]*\n)\n+(Atoms)"
if not re.search(pat, text):
    sys.exit(0)
text = re.sub(pat, rf"\1\n{line}\n\n\2", text, count=1)
open(path, "w").write(text)
print(f"Fixed empty Dihedral Coeffs in {path}")
PY
}

fix_dihedral_coeffs pe_at_ref_final.data
fix_dihedral_coeffs pe_at_final.data

echo "=== Reference 5 ns production ==="
"$LMP" "${LMP_OMP[@]}" -in in.pe_at_ref_prod5ns -log log.pe_at_ref_prod5ns.lammps 2>&1 | tee log.pe_at_ref_prod5ns.out

echo "=== Backmapped 5 ns production ==="
"$LMP" "${LMP_OMP[@]}" -in in.pe_at_prod5ns -log log.pe_at_prod5ns.lammps 2>&1 | tee log.pe_at_prod5ns.out

echo "=== Done. Compare with rdf_*_5ns.dat ==="
