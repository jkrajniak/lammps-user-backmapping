#!/usr/bin/env bash
# Review item B7 (R1-4): cost of pure CG, the hybrid during the lambda ramp,
# and pure AT for one system, same machine and rank count.
# Inputs from a finished Tier B run directory: <prefix>_hybrid.data with
# <prefix>.ff.lmp / <prefix>.backmap.lmp, <prefix>_at.data with <prefix>.at.ff.lmp.
# Pure CG: the hybrid frame with the AT atoms deleted and the CG terms of the
# hybrid force field as stock styles (cg_part_of_hybrid.py).
# Each case: STEPS NVT steps after a short warm-up; LAMMPS "Performance" line
# and loop time are collected by the caller from log.b7_<case>.
#   LMP=/path/to/lmp NP=8 bash b7_cost.sh <run-dir> <prefix> <T> <outdir>
set -euo pipefail
LMP="${LMP:?set LMP}"
NP="${NP:-8}"
HERE="$(cd "$(dirname "$0")" && pwd)"
STEPS="${STEPS:-2000}"
WARM="${WARM:-200}"
DT_CG="${DT_CG:-5.0}"
DT_HYB="${DT_HYB:-0.1}"   # the generated protocol's ramp timestep
DT_AT="${DT_AT:-1.0}"
src="$1"; prefix="$2"; T="$3"; out="$4"
# Safe to call again: a finished measurement (b7_summary.txt) is not repeated; an interrupted one runs
# all three cases again, because timings of a partly interrupted run are not comparable.
if [ -s "$out/b7_summary.txt" ]; then echo "b7 already finished: $out/b7_summary.txt"; exit 0; fi
mkdir -p "$out"
cp -r "$src/." "$out/"
cd "$out"

timed() { # name data setup-lines dt group
  cat > "in.b7_$1" <<EOT
units real
atom_style full
boundary p p p
read_data $2
$3
thermo 1000
velocity $5 create $T 4928459 mom yes rot yes dist gaussian loop geom
timestep $4
fix integrate $5 nvt temp $T $T $(python3 -c "print(100 * $4)")
run $WARM
run $STEPS
EOT
  mpirun -np "$NP" "$LMP" -in "in.b7_$1" -log "log.b7_$1" -screen none
}

# Pure CG: hybrid frame without its AT atoms, CG terms as stock styles
python3 "$HERE/cg_part_of_hybrid.py" "$prefix.ff.lmp" > cg_ff.lmp
cg_types=$(grep -oE "cg_type( [0-9]+)+" "$prefix.backmap.lmp" | cut -d" " -f2-)
timed cg "${prefix}_hybrid.data" "group cgb type $cg_types
group atd subtract all cgb
delete_atoms group atd bond yes mol no
include cg_ff.lmp" "$DT_CG" all
# Hybrid during the ramp: lambda from 0.5, ramp active, AT atoms integrated
sed -E "s/^(fix bm all backmap .*) lambda0 [^ ]+(.*)$/\1 lambda0 0.5\2/" "$prefix.backmap.lmp" > backmap_b7.lmp
timed hybrid "${prefix}_hybrid.data" "include $prefix.ff.lmp
include backmap_b7.lmp
fix_modify bm active yes" "$DT_HYB" at_atoms
# Pure AT
timed at "${prefix}_at.data" "include $prefix.at.ff.lmp" "$DT_AT" all

for c in cg hybrid at; do
  # "Loop time of T on P procs for S steps with N atoms" of the timed run
  atoms=$(grep "^Loop time" "log.b7_$c" | tail -1 | awk '{print $(NF-1)}')
  loop=$(grep "^Loop time" "log.b7_$c" | tail -1 | awk '{print $4}')
  perf=$(grep "^Performance" "log.b7_$c" | tail -1)
  echo "$prefix $c atoms=$atoms ranks=$NP steps=$STEPS loop_s=$loop s_per_step=$(python3 -c "print($loop / $STEPS)") | $perf"
done | tee b7_summary.txt
