#!/usr/bin/env bash
# Review item B3 (R1-1) control: potential energy of the backmapped and of the independent reference
# system at the SAME density. The final atomistic configuration of a run is rescaled uniformly (atom
# coordinates and box) to RHO g/cm3 and run in NVT: EQUIL steps equilibration, then PROD steps
# production at 1 fs (thermo every 5000 steps). Resumable: restart files every CKPT steps, run under
# resumable_lmp.sh; safe to call again after an interruption.
#   LMP=... NP=2 bash b3_matched_density.sh <final.data> <at.ff.lmp> <outdir> <T> <RHO> [EQUIL PROD]
# Energy terms of the production: b3_energy_terms.py <log_backmap> <log_reference> <n_atoms> --from-step EQUIL
# with log.matched.s*.lammps of the two output directories.
set -euo pipefail
LMP="${LMP:?}"; NP="${NP:-2}"
HERE="$(cd "$(dirname "$0")" && pwd)"
RESUMABLE="${RESUMABLE:-$HERE/resumable_lmp.sh}"
data="$1"; ff="$2"; out="$3"; T="$4"; RHO="$5"; EQ="${6:-100000}"; PR="${7:-300000}"
CKPT="${CKPT:-50000}"
mkdir -p "$out"

prepare() {
  cp "$data" "$out/start.data"
  cp "$ff" "$out/at.ff.lmp"
  cp "$(dirname "$data")"/*.pairs14.dat "$out/" 2> /dev/null || true
  # the force-field file defines every style and coefficient: drop the coefficient sections of write_data output
  python3 - "$out/start.data" <<'PY2'
import re, sys
from pathlib import Path
f = Path(sys.argv[1])
out, skip = [], False
for ln in f.read_text().splitlines():
    if re.match(r"^(Pair|PairIJ|Bond|Angle|Dihedral|Improper) Coeffs", ln):
        skip = True
        continue
    if skip and re.match(r"^[A-Za-z]", ln):
        skip = False
    if not skip:
        out.append(ln)
f.write_text("\n".join(out) + "\n")
PY2
  # scale factor from the mass density of the file
  python3 - "$out/start.data" "$RHO" > "$out/scale.txt" <<'PY3'
import sys
t = open(sys.argv[1]).read().splitlines()
rho = float(sys.argv[2])
lo, mass = [], {}
for l in t:
    p = l.split()
    if len(p) >= 4 and p[2] in ("xlo", "ylo", "zlo"):
        lo.append((float(p[0]), float(p[1])))
i = next(k for k, l in enumerate(t) if l.startswith("Masses")) + 2
while t[i].strip():
    p = t[i].split(); mass[int(p[0])] = float(p[1]); i += 1
j = next(k for k, l in enumerate(t) if l.startswith("Atoms")) + 2
M = 0.0
while t[j].strip():
    p = t[j].split(); M += mass[int(p[2])]; j += 1
vol = 1.0
for a, b in lo:
    vol *= b - a
cur = M / 6.02214076e23 / (vol * 1e-24)
print(f"{cur:.6f} {(cur / rho) ** (1 / 3):.8f}")
PY3
  read -r cur f < "$out/scale.txt"
  cat > "$out/in.matched" <<EOT
variable resume index 0
variable resume_file index none
variable seg index 1
variable ckpt index $CKPT
units real
atom_style full
boundary p p p
if "\${resume} == 0" then "read_data start.data" else "read_restart \${resume_file}"
include at.ff.lmp
thermo 5000
thermo_style custom step temp pe ke etotal ebond eangle edihed evdwl ecoul press density vol
thermo_modify flush yes
restart \${ckpt} rst.s\${seg}.a rst.s\${seg}.b
if "\${resume} == 1" then "jump SELF resume"
change_box all x scale $f y scale $f z scale $f remap
velocity all create $T 4928459 mom yes rot yes dist gaussian
label resume
timestep 1.0
fix integrate all nvt temp $T $T 100.0
if "\$(step) >= $EQ" then "jump SELF prod"
run $EQ upto
label prod
run $((EQ + PR)) upto
EOT
  echo "start density $cur -> $RHO, scale factor $f" | tee "$out/README.txt"
}
[ -e "$out/.done.prep" ] || { prepare && touch "$out/.done.prep"; }

cd "$out"
export LMP NP
"$RESUMABLE" matched -in in.matched -log log.matched.lammps
