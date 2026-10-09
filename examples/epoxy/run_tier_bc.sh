#!/usr/bin/env bash
# RIM135 (epoxy) Tier B (backmapping) and Tier C (AT continuation for the structure analysis).
#   LMP=/path/to/lmp NP=8 bash run_tier_bc.sh backmap
# Tier C continues the backmapped AT system with the protocol of the published
# GROMACS reference (JCC 2018 data, rim135_small/backmapping/aa/template_eq): em,
# annealing NPT 1000 -> 298 K (5 ns), NVT 1 ns, NPT 5 ns, NVT 5 ns, RDF 2.5-5 ns.
# There is no reference stage here. Every force field comes from backmap-prep; the inputs
# here are protocol only (AT production: ../common/in.at_reference_protocol).
set -euo pipefail
cd "$(dirname "$0")"
LMP="${LMP:?set LMP to a LAMMPS binary with the BACKMAP package}"
NP="${NP:-8}"
PREP="${PREP:-uv run backmap-prep}"
AT_IN="${AT_IN:-../common/in.at_reference_protocol}"
lmp() { mpirun -np "$NP" "$LMP" "$@"; }
# Spot or preemptible machines: RESUMABLE=/path/to/scripts/review/resumable_lmp.sh makes the long
# atomistic stage restart from its newest checkpoint, and every finished step leaves a .done.<step>
# marker, so calling this script again after an interruption continues where it stopped.
RESUMABLE="${RESUMABLE:-}"
once() { local m=".done.$1"; shift; [ -e "$m" ] || { "$@" && touch "$m"; }; }
lmp_at() { if [ -n "$RESUMABLE" ]; then LMP="$LMP" NP="$NP" "$RESUMABLE" "$@"; else shift; lmp "$@"; fi; }
prepare_at() {
  $PREP at-system settings.v2.yaml --from rim135_hybrid.data &&
    grep "^pair_coeff" rim135.at.ff.lmp > rim135.at.pair_coeffs
}

case "${1:?backmap}" in
  backmap)
    once build $PREP build settings.v2.yaml
    once hybrid lmp -in in.rim135 -log log.rim135.lammps
    once atsystem prepare_at
    lmp_at at -in "$AT_IN" -log log.rim135_at.lammps -var prefix rim135 -var anneal 1 ${AT_VARS:-}
    ;;
  *) echo "unknown stage $1" >&2; exit 1 ;;
esac
