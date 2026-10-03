#!/usr/bin/env bash
# Review item B8 (R1-3): how the quality of the CG model carries into the
# backmapped structure. Dodecane melt (500 molecules, 298 K), three CG models
# that differ only in the nonbonded bead-bead tables:
#   ibi  published IBI tables (control)
#   lj   12-6 LJ fitted to the IBI first well (b8_cg_tables.py)
#   wca  repulsion only (b8_cg_tables.py)
#   dbi  direct Boltzmann inversion of the bead RDF of the independent atomistic reference
#        (b8_dbi_tables.py; needs DBI_DATA = the reference's final data file, DBI_TRAJ = its DCD)
# Per model: CG NVT at the published density (bead RDF over the second half;
# the mean CG pressure is logged), the hybrid built from the model's own
# equilibrated frame, the generated backmapping protocol, and the AT
# continuation (in.dodecane_at) with a dense trajectory over its first NPT
# stage and PROD_STEPS of production. No CG NPT: the IBI tables carry no
# pressure correction, so the CG melt expands under NPT by construction and
# the CG density is not a property of the model.
#   LMP=/path/to/lmp NP=8 OUT=/path/to/outdir bash b8_degraded_cg.sh <examples-dir> [variants...]
# Safe to call again after an interruption (see below).
# SMOKE=1 stops each variant after the hybrid rebuild (CG_STEPS small for a quick check).
set -euo pipefail
LMP="${LMP:?set LMP}"
OUT="${OUT:?set OUT}"
NP="${NP:-8}"
PREP="${PREP:-uv run backmap-prep}"
PY="${PY:-uv run --no-project --with numpy python}"
PY_DBI="${PY_DBI:-uv run --no-project --with numpy --with MDAnalysis python}"
CG_STEPS="${CG_STEPS:-500000}"   # 1 ns at 2 fs
PROD_STEPS="${PROD_STEPS:-1000000}"
HERE="$(cd "$(dirname "$0")" && pwd)"
EXAMPLES="$(cd "${1:?examples directory}" && pwd)"
shift
VARIANTS="${*:-ibi lj wca}"
export LMP NP
lmp() { mpirun -np "$NP" "$LMP" "$@"; }

# Spot or preemptible machines: every stage leaves a .done.<stage> marker in the model directory
# and the atomistic stage is resumable (RESUMABLE, default: resumable_lmp.sh next to this script), so
# calling this script again after an interruption continues where it stopped.
RESUMABLE="${RESUMABLE:-$HERE/resumable_lmp.sh}"
once() { local m=".done.$1"; shift; [ -e "$m" ] || { "$@" && touch "$m"; }; }

tables() { # variant
  case "$1" in
    ibi) ;;
    dbi) $PY_DBI "$HERE/b8_dbi_tables.py" --dir . --data "${DBI_DATA:?set DBI_DATA}" --traj "${DBI_TRAJ:?set DBI_TRAJ}" ;;
    *) $PY "$HERE/b8_cg_tables.py" --dir . --variant "$1" ;;
  esac
}
# CG model = the CG part of the generated hybrid force field
cg_model() {
  $PREP build settings.yaml > build_published_frame.log 2>&1 &&
    python3 "$HERE/cg_part_of_hybrid.py" dodecane.ff.lmp > cg_ff.lmp
}
cg_nvt() {
  local cg_types ta tb
  cg_types=$(grep -oE "cg_type( [0-9]+)+" dodecane.backmap.lmp | cut -d" " -f2-)
  read -r ta tb <<< "$cg_types"
  cat > in.b8_cg_nvt <<EOF
units real
atom_style full
boundary p p p
read_data dodecane.data
group cgb type $cg_types
group atd subtract all cgb
delete_atoms group atd bond yes mol no
include cg_ff.lmp
velocity all create 298.0 48279 mom yes rot yes dist gaussian loop geom
timestep 2.0
thermo 5000
thermo_style custom step temp pe press density
fix integ all nvt temp 298.0 298.0 100.0
run $((CG_STEPS / 2))
compute rdf_cg all rdf 150 $ta $ta $ta $tb $tb $tb
fix rdf_out all ave/time 100 $((CG_STEPS / 200)) $CG_STEPS c_rdf_cg[*] file rdf_cg.dat mode vector
variable pcg equal press
fix pcg_out all ave/time 100 $((CG_STEPS / 200)) $CG_STEPS v_pcg file pressure_cg.dat
run $((CG_STEPS / 2))
write_data cg_equil_nvt.data
EOF
  lmp -in in.b8_cg_nvt -log log.b8_cg_nvt -screen none
}
model_frame() { # hybrid from the model's own equilibrated CG frame
  cp cg_conf.gro cg_conf_published.gro &&
    python3 "$HERE/b8_data_to_gro.py" cg_equil_nvt.data cg_conf.gro &&
    $PREP build settings.yaml > build.log 2>&1
}
hybrid() { lmp -in in.dodecane -log log.dodecane.lammps -screen none; }
at_system() { $PREP at-system settings.yaml --from dodecane_hybrid.data > at_system.log 2>&1; }
at_input() {
  # dense trajectory over the first AT stage (NPT, 100 ps): relaxation after backmapping
  awk '/^fix +integrate all npt/ {print; print "dump early all dcd 200 traj_early.s${seg}.dcd"; print "dump_modify early unwrap yes"; next}
       /^unfix +integrate/ && !done {print "undump early"; done=1} {print}' in.dodecane_at > in.b8_at
  grep -q "dump early" in.b8_at && grep -q "undump early" in.b8_at
}

for v in $VARIANTS; do
  W="$OUT/$v"
  mkdir -p "$OUT"
  [ -e "$W" ] || cp -r "$EXAMPLES/dodecane/large" "$W"
  cd "$W"
  if [ -e ".done.all" ]; then cd - > /dev/null; continue; fi
  once tables tables "$v"
  once cgmodel cg_model
  once cgnvt cg_nvt
  once frame model_frame
  if [ -n "${SMOKE:-}" ]; then echo "$v smoke stop after the rebuild" >> "$OUT/b8.log"; cd - > /dev/null; continue; fi
  once hybrid hybrid
  once atsystem at_system
  once atinput at_input
  "$RESUMABLE" at -in in.b8_at -log log.b8_at.lammps -var prod_steps "$PROD_STEPS"
  touch .done.all
  echo "$v done $(date -u)" >> "$OUT/b8.log"
  cd - > /dev/null
done
