#!/usr/bin/env bash
# Review item B2 (R1-m3): strong scaling of the lambda-ramp block with the reduced communication cutoff.
# One system, ranks RANKS (default 1 2 4 8), REPS repetitions each. The measured block is the hybrid at
# lambda = 0.5 with the ramp active and the AT atoms integrated (the same case as the hybrid case of b7_cost.sh):
# WARM warm-up steps, then STEPS steps whose loop time is read from the LAMMPS log.
# Inputs from a finished Tier B run directory: <prefix>_hybrid.data with <prefix>.ff.lmp / <prefix>.backmap.lmp.
# Run on a quiet machine: ranks are bound to cores, other jobs distort the result. b2_summary.py collects the logs.
#   LMP=/path/to/lmp bash b2_scaling.sh <run-dir> <prefix> <T> <outdir>
set -euo pipefail
LMP="${LMP:?set LMP}"
RANKS="${RANKS:-1 2 4 8}"
REPS="${REPS:-3}"
STEPS="${STEPS:-500}"
WARM="${WARM:-200}"
DT_HYB="${DT_HYB:-0.1}"   # the generated protocol's ramp timestep
BIND="${BIND:---bind-to core --map-by core}"
src="$1"; prefix="$2"; T="$3"; out="$4"
# Safe to call again after an interruption: a finished case (log with "Total wall time") is skipped,
# an interrupted case runs again from its start.
mkdir -p "$out"
for f in "$prefix"_hybrid.data "$prefix.ff.lmp" "$prefix.backmap.lmp" pairs.dat; do
  [ -e "$src/$f" ] && cp "$src/$f" "$out/"
done
cp "$src"/table_*.table "$out/" 2>/dev/null || true
cd "$out"

# A data file written by write_data carries coefficient sections; they need the styles defined before read_data,
# and the force field include sets all coefficients anyway, so they are dropped.
awk '/^[A-Za-z]/ {skip = ($0 ~ /Coeffs/)} !skip' "${prefix}_hybrid.data" > "${prefix}_hybrid_b2.data"
sed -E "s/^(fix bm all backmap .*) lambda0 [^ ]+(.*)$/\1 lambda0 0.5\2/" "$prefix.backmap.lmp" > backmap_b2.lmp
cat > in.b2 <<EOT
units real
atom_style full
boundary p p p
read_data ${prefix}_hybrid_b2.data
include ${prefix}.ff.lmp
include backmap_b2.lmp
fix_modify bm active yes
thermo 1000
velocity at_atoms create $T 4928459 mom yes rot yes dist gaussian loop geom
timestep $DT_HYB
fix integrate at_atoms nvt temp $T $T $(python3 -c "print(100 * $DT_HYB)")
run $WARM
run $STEPS
EOT

for np in $RANKS; do
  for rep in $(seq 1 "$REPS"); do
    log="log.b2_np${np}_r${rep}"
    if grep -qs "^Total wall time" "$log"; then echo "b2 np=$np rep=$rep already finished"; continue; fi
    # shellcheck disable=SC2086
    mpirun -np "$np" $BIND "$LMP" -in in.b2 -log "$log" -screen none
  done
done
grep -h "comm_modify" "$prefix.ff.lmp" "$prefix.backmap.lmp" > b2_comm_cutoff.txt || true
echo B2_DONE
