#!/usr/bin/env bash
# Review item B1 (R2-3): static MPI parity of the hybrid force field.
# For each system: read the generated initial hybrid frame, set lambda to
# LAMBDA (default 0.5, so CG and AT forces both act), `run 0`, and dump per-atom
# forces, bead COMs and bead CG forces (fix backmap peratom full) and the
# energies, for several rank counts, processor grids and a tiled RCB
# decomposition. compare_mpi_parity.py compares every run with the 1-rank run.
#   LMP=/path/to/lmp OUT=/path/to/outdir bash b1_mpi_parity.sh <example-dir> <settings> <prefix> [...]
# e.g. bash b1_mpi_parity.sh examples/dodecane/large settings.yaml dodecane
set -euo pipefail
LMP="${LMP:?set LMP to a LAMMPS binary with the BACKMAP package}"
OUT="${OUT:?set OUT to an output directory}"
PREP="${PREP:-uv run backmap-prep}"
LAMBDA="${LAMBDA:-0.5}"
HERE="$(cd "$(dirname "$0")" && pwd)"
# name:ranks:processors-or-tiled
CONFIGS="${CONFIGS:-r1:1:1*1*1 r2x:2:2*1*1 r2z:2:1*1*2 r4x:4:4*1*1 r4z:4:1*1*4 r4xy:4:2*2*1 r4t:4:tiled r8:8:2*2*2 r8x:8:8*1*1 r8t:8:tiled}"

while [ $# -ge 3 ]; do
  src="$1"; settings="$2"; prefix="$3"; shift 3
  work="$OUT/$prefix"
  [ -e "$work" ] && { echo "$work exists, not overwriting" >&2; exit 1; }
  mkdir -p "$OUT" && cp -r "$src" "$work"
  (cd "$work" && $PREP build "$settings" > build.log 2>&1)
  # lambda fixed at LAMBDA, per-bead output
  sed -i -E "s/^(fix bm all backmap .*) lambda0 [^ ]+(.*)$/\1 lambda0 ${LAMBDA}\2 peratom full/" "$work/$prefix.backmap.lmp"
  grep -q "peratom full" "$work/$prefix.backmap.lmp"
  for cfg in $CONFIGS; do
    IFS=: read -r name np grid <<< "$cfg"
    if [ "$grid" = tiled ]; then
      decomp="comm_style tiled"; balance="balance 1.0 rcb"
    else
      decomp="processors ${grid//\*/ }"; balance=""
    fi
    cat > "$work/in.b1_$name" <<EOT
units real
atom_style full
boundary p p p
$decomp
read_data $prefix.data
include $prefix.ff.lmp
include $prefix.backmap.lmp
$balance
thermo_style custom step pe evdwl ecoul ebond eangle edihed eimp
thermo_modify format float %.17g
dump d all custom 1 forces_$name.dump id type x y z fx fy fz f_bm[1] f_bm[2] f_bm[3] f_bm[4] f_bm[5] f_bm[6] f_bm[7]
dump_modify d sort id format float %.17g
run 0
print "RESULT pe=\$(pe:%.17g) evdwl=\$(evdwl:%.17g) ecoul=\$(ecoul:%.17g) ebond=\$(ebond:%.17g) eangle=\$(eangle:%.17g) edihed=\$(edihed:%.17g) eimp=\$(eimp:%.17g)"
EOT
    (cd "$work" && mpirun -np "$np" "$LMP" -in "in.b1_$name" -log "log.b1_$name" -screen none)
  done
  uv run --no-project --with numpy python "$HERE/compare_mpi_parity.py" "$work" --ref r1 | tee "$work/parity.txt"
done
