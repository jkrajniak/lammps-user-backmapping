#!/usr/bin/env bash
# PE (all atom, published pe) Tier B (backmapping) and Tier C (AT RDFs vs an independent reference).
#   LMP=/path/to/lmp NP=8 bash run_tier_bc.sh backmap     # Tier B, then Tier C of the backmapped melt
#   LMP=/path/to/lmp NP=8 bash run_tier_bc.sh reference   # independent AT reference
# Every force field comes from backmap-prep; the scripts here are protocol only.
set -euo pipefail
cd "$(dirname "$0")"
LMP="${LMP:?set LMP to a LAMMPS binary with the BACKMAP package}"
NP="${NP:-8}"
PREP="${PREP:-uv run backmap-prep}"
lmp() { mpirun -np "$NP" "$LMP" "$@"; }
# Spot or preemptible machines: RESUMABLE=/path/to/scripts/review/resumable_lmp.sh makes the long
# atomistic stage restart from its newest checkpoint, and every finished step leaves a .done.<step>
# marker, so calling this script again after an interruption continues where it stopped.
RESUMABLE="${RESUMABLE:-}"
once() { local m=".done.$1"; shift; [ -e "$m" ] || { "$@" && touch "$m"; }; }
lmp_at() { if [ -n "$RESUMABLE" ]; then LMP="$LMP" NP="$NP" "$RESUMABLE" "$@"; else shift; lmp "$@"; fi; }

once build $PREP build settings.yaml
L=$(awk '/xlo xhi/ {print $2 - $1; exit}' pe_aa.data)

case "${1:?backmap or reference}" in
  backmap)
    once hybrid lmp -in in.pe_aa -log log.pe_aa.lammps
    once atsystem $PREP at-system settings.yaml --from pe_aa_hybrid.data
    lmp_at at -in in.pe_aa_at -log log.pe_aa_at.lammps ${PROD_STEPS:+-var prod_steps $PROD_STEPS}
    ;;
  reference)
    once atsystem_ref $PREP at-system settings.yaml --reference 75 --seed 42
    # not restartable mid-run (high-temperature compression and cooling stages); finished runs are skipped
    once reference lmp -in in.pe_aa_at_ref -log log.pe_aa_at_ref.lammps -var L "$L" ${PROD_STEPS:+-var prod_steps $PROD_STEPS}
    ;;
  *) echo "unknown stage $1" >&2; exit 1 ;;
esac
