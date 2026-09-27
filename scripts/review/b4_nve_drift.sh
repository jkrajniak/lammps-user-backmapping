#!/usr/bin/env bash
# Review item B4 (R1-1): energy conservation of the hybrid at fixed lambda.
# From a backmapped hybrid frame (<prefix>_hybrid.data, written by the generated
# protocol) and its AT-only counterpart (<prefix>_at.data from `at-system`):
# for lambda in LAMBDAS (ramp inactive, no CapForce) and the plain AT system,
# thermalize (NVT, 0.5 fs), then run NVE at each timestep in DTS and write the
# total energy; b4_drift.py fits the drift per atom per ns. In the hybrid only
# the AT atoms are integrated: fix backmap places each bead on its fragment COM
# and zeroes its velocity, so the conserved energy is PE + KE(AT).
#   LMP=/path/to/lmp NP=4 bash b4_nve_drift.sh <run-dir-with-outputs> <prefix> <T> <outdir>
set -euo pipefail
LMP="${LMP:?set LMP}"
NP="${NP:-4}"
LAMBDAS="${LAMBDAS:-0.0 0.5 1.0}"
DTS="${DTS:-1.0 0.5}"
NVE_PS="${NVE_PS:-100}"
THERM_PS="${THERM_PS:-20}"
src="$1"; prefix="$2"; T="$3"; out="$4"
[ -e "$out" ] && { echo "$out exists, not overwriting" >&2; exit 1; }
mkdir -p "$out"
for f in "$prefix"_hybrid.data "$prefix".ff.lmp "$prefix".backmap.lmp "$prefix"_at.data "$prefix".at.ff.lmp pairs.dat; do
  [ -e "$src/$f" ] && cp "$src/$f" "$out/"
done
cp "$src"/table_*.table "$out/" 2>/dev/null || true
cd "$out"

run_case() { # name data include-lines dt group
  local name="$1" data="$2" includes="$3" dt="$4" grp="$5"
  local therm=$(python3 -c "print(int(round($THERM_PS * 1000 / 0.5)))")
  local nve=$(python3 -c "print(int(round($NVE_PS * 1000 / $dt)))")
  local every=$(python3 -c "print(max(1, int(round(100 / $dt))))")
  cat > "in.b4_$name" <<EOT
units real
atom_style full
boundary p p p
read_data $data
$includes
compute tgrp $grp temp
thermo_style custom step time c_tgrp pe ke etotal
thermo_modify format float %.12g norm no
thermo $every
reset_timestep 0
velocity $grp create $T 4928459 mom yes rot yes dist gaussian
timestep 0.5
fix therm $grp nvt temp $T $T 50.0
fix_modify therm temp tgrp
run $therm
unfix therm
reset_timestep 0
timestep $dt
fix integrate $grp nve
run $nve
EOT
  mpirun -np "$NP" "$LMP" -in "in.b4_$name" -log "log.b4_$name" -screen none
}

for dt in $DTS; do
  for lam in $LAMBDAS; do
    sed -E "s/^(fix bm all backmap .*) lambda0 [^ ]+(.*)$/\1 lambda0 ${lam}\2/" "$prefix.backmap.lmp" > "backmap_l${lam}.lmp"
    run_case "hyb_l${lam}_dt${dt}" "${prefix}_hybrid.data" "include $prefix.ff.lmp
include backmap_l${lam}.lmp" "$dt" at_atoms
  done
  run_case "at_dt${dt}" "${prefix}_at.data" "include $prefix.at.ff.lmp" "$dt" all
done
echo B4_DONE
