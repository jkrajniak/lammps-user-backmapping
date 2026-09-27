#!/usr/bin/env bash
# Melamine network Tier B (backmapping) and Tier C (AT continuation for the structure analysis).
#   LMP=/path/to/lmp NP=8 bash run_tier_bc.sh backmap
# Tier C continues the backmapped AT system with the protocol of the published
# GROMACS reference (bakery examples/network_backmapping/mf/backmapping/template_eq):
# em, NPT 5 ns, NVT 10 ns, RDF 5-10 ns.
# There is no reference stage here. Every force field comes from backmap-prep; the inputs
# here are protocol only (AT production: ../../common/in.at_reference_protocol).
set -euo pipefail
cd "$(dirname "$0")"
LMP="${LMP:?set LMP to a LAMMPS binary with the BACKMAP package}"
NP="${NP:-8}"
PREP="${PREP:-uv run backmap-prep}"
AT_IN="${AT_IN:-../../common/in.at_reference_protocol}"
lmp() { mpirun -np "$NP" "$LMP" "$@"; }

case "${1:?backmap}" in
  backmap)
    $PREP build settings.yaml
    lmp -in in.melamine_network -log log.melamine_network.lammps
    $PREP at-system settings.yaml --from melamine_network_hybrid.data
    grep "^pair_coeff" melamine_network.at.ff.lmp > melamine_network.at.pair_coeffs
    lmp -in "$AT_IN" -log log.melamine_network_at.lammps -var prefix melamine_network -var nvt_steps 5000000 -var window_steps 5000000 ${AT_VARS:-}
    ;;
  *) echo "unknown stage $1" >&2; exit 1 ;;
esac
