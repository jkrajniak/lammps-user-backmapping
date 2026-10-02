#!/usr/bin/env bash
# Review item B3 (R1-1) control: potential energy of the backmapped and the independent reference
# system at the SAME density. Each final AT configuration is rescaled (atom coordinates, uniform) to
# RHO g/cm3 and run in NVT (EQUIL steps equilibration, PROD steps production at 1 fs, thermo every 5000).
#   LMP=... NP=2 bash b3_matched_density.sh <final.data> <at.ff.lmp> <outdir> <T> <RHO> [EQUIL PROD]
set -euo pipefail
LMP="${LMP:?}"; NP="${NP:-2}"
data="$1"; ff="$2"; out="$3"; T="$4"; RHO="$5"; EQ="${6:-100000}"; PR="${7:-300000}"
mkdir -p "$out"; cp "$data" "$out/start.data"; cp "$ff" "$out/at.ff.lmp"
cp "$(dirname "$data")"/*.pairs14.dat "$out/" 2>/dev/null || true
# the force-field file defines every style and coefficient: drop the coefficient sections of write_data output
python3 - "$out/start.data" <<'PY2'
import re, sys
from pathlib import Path
f = Path(sys.argv[1]); lines = f.read_text().splitlines()
out, skip = [], False
for ln in lines:
    if re.match(r"^(Pair|PairIJ|Bond|Angle|Dihedral|Improper) Coeffs", ln):
        skip = True; continue
    if skip and re.match(r"^[A-Za-z]", ln):
        skip = False
    if not skip:
        out.append(ln)
f.write_text("\n".join(out) + "\n")
PY2
# current density and box from the data file and the log of the run that wrote it are not needed: use mass and volume
python3 - "$out/start.data" "$RHO" > "$out/scale.txt" <<'PY'
import sys
t=open(sys.argv[1]).read().splitlines(); rho=float(sys.argv[2])
lo=[];mass={}
for i,l in enumerate(t):
    p=l.split()
    if len(p)>=4 and p[2] in("xlo","ylo","zlo"): lo.append((float(p[0]),float(p[1])))
i=next(k for k,l in enumerate(t) if l.startswith("Masses"))+2
while t[i].strip():
    p=t[i].split(); mass[int(p[0])]=float(p[1]); i+=1
j=next(k for k,l in enumerate(t) if l.startswith("Atoms"))+2
M=0.0
while t[j].strip():
    p=t[j].split(); M+=mass[int(p[2])]; j+=1
vol=1.0
for a,b in lo: vol*=b-a
cur=M/6.02214076e23/(vol*1e-24)
L=[b-a for a,b in lo]
f=(cur/rho)**(1/3)
print(f"{cur:.6f} {f:.8f} "+" ".join(f"{x*f:.6f}" for x in L))
PY
read -r cur f lx ly lz < "$out/scale.txt"
cat > "$out/in.matched" <<EOT
units real
atom_style full
boundary p p p
read_data start.data
include at.ff.lmp
change_box all x scale $f y scale $f z scale $f remap
thermo 5000
thermo_style custom step temp pe ke etotal ebond eangle edihed evdwl ecoul press density vol
velocity all create $T 4928459 mom yes rot yes dist gaussian
timestep 1.0
fix integrate all nvt temp $T $T 100.0
run $EQ
reset_timestep 0
run $PR
EOT
echo "start density $cur -> $RHO, scale factor $f" | tee "$out/README.txt"
cd "$out" && mpirun -np "$NP" "$LMP" -in in.matched -log log.matched -screen none
