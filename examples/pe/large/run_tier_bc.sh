#!/usr/bin/env bash
# PE Tier B (backmapping) and Tier C (AT RDFs vs an independent reference).
#   LMP=/path/to/lmp NP=8 bash run_tier_bc.sh backmap     # Tier B, then Tier C of the backmapped melt
#   LMP=/path/to/lmp NP=8 bash run_tier_bc.sh reference   # independent AT reference
# Every force field comes from backmap-prep; the scripts here are protocol only.
set -euo pipefail
cd "$(dirname "$0")"
LMP="${LMP:?set LMP to a LAMMPS binary with the BACKMAP package}"
NP="${NP:-8}"
PREP="${PREP:-uv run backmap-prep}"
lmp() { mpirun -np "$NP" "$LMP" "$@"; }

$PREP build settings.yaml
L=$(awk '/xlo xhi/ {print $2 - $1; exit}' pe.data)

case "${1:?backmap or reference}" in
  backmap)
    lmp -in in.pe -log log.pe.lammps
    $PREP at-system settings.yaml --from pe_hybrid.data
    lmp -in in.pe_at -log log.pe_at.lammps ${PROD_STEPS:+-var prod_steps $PROD_STEPS}
    ;;
  reference)
    $PREP at-system settings.yaml --reference 75 --seed 42
    lmp -in in.pe_at_ref -log log.pe_at_ref.lammps -var L "$L" ${PROD_STEPS:+-var prod_steps $PROD_STEPS}
    ;;
  *) echo "unknown stage $1" >&2; exit 1 ;;
esac
