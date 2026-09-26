#!/usr/bin/env bash
# MARTINI POPC Tier B: backmap in LAMMPS, then continue the AT system in GROMACS
# with the settings of the independent Slipids reference.
#   LMP=/path/to/lmp GMX=/path/to/gmx NP=8 bash run_tier_b.sh
# Keeps dump.backmap (lambda per atom every 1000 steps) for the transition figure.
set -euo pipefail
cd "$(dirname "$0")"
LMP="${LMP:?set LMP to a LAMMPS binary with the BACKMAP, MOLECULE, EXTRA-MOLECULE packages}"
GMX="${GMX:?set GMX to a GROMACS binary}"
NP="${NP:-8}"
PREP="${PREP:-uv run backmap-prep}"
PY="${PY:-uv run python}"

bash fetch_sources.sh
$PY make_cg_topology.py
$PREP build settings.yaml
mpirun -np "$NP" "$LMP" -in in.popc -log log.popc.lammps
$PY to_gromacs.py popc_hybrid.data at_backmapped.gro

mkdir -p at_continue && cd at_continue
for f in ../gromacs/*.mdp ../gromacs/topol_at.top ../gromacs/tip3p.itp ../POPC.itp; do cp "$f" .; done
ln -sfn ../amber99sb-ildn_slipids.ff amber99sb-ildn_slipids.ff
"$GMX" grompp -f em.mdp -c ../at_backmapped.gro -p topol_at.top -o em.tpr
"$GMX" mdrun -deffnm em -nt "$NP"
"$GMX" grompp -f eq.mdp -c em.gro -p topol_at.top -o eq.tpr
"$GMX" mdrun -deffnm eq -nt "$NP"
"$GMX" grompp -f prod.mdp -c eq.gro -t eq.cpt -p topol_at.top -o prod.tpr
"$GMX" mdrun -deffnm prod -nt "$NP"
