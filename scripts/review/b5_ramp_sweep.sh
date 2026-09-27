#!/usr/bin/env bash
# Review item B5 (R1-2): ramp-time sweep. One job = one example backmapped with
# simulation.alpha = ALPHA and simulation.rng_seed = SEED (the generator sets
# every stage length and seed from these), then the example's own Tier B/C
# script (run_tier_bc.sh backmap) with a shortened AT production (PROD_STEPS).
# Metrics come from the logs and the production outputs (b5_metrics.py).
#   LMP=/path/to/lmp PREP="uv run backmap-prep" NP=4 PROD_STEPS=500000 \
#     bash b5_ramp_sweep.sh <example-dir> <settings> <alpha> <seed> <outdir>
set -euo pipefail
: "${LMP:?set LMP}"
export LMP NP="${NP:-4}" PREP="${PREP:-uv run backmap-prep}" PROD_STEPS="${PROD_STEPS:-500000}"
src="$1"; settings="$2"; alpha="$3"; seed="$4"; out="$5"
[ -e "$out" ] && { echo "$out exists, not overwriting" >&2; exit 1; }
cp -r "$src" "$out"
python3 - "$src" "$out/$settings" "$alpha" "$seed" <<'PY'
import re, sys
from pathlib import Path
src, settings, alpha, seed = Path(sys.argv[1]).resolve(), Path(sys.argv[2]), sys.argv[3], sys.argv[4]
text = settings.read_text()
# relative settings paths must still resolve from the copy
text = re.sub(r"^(\s*(?:data_dir|tables_dir|forcefield_dir|bakery_xml):\s*)([^\s#]+)",
              lambda m: m.group(1) + (m.group(2) if m.group(2).startswith("/") else str((src / m.group(2)).resolve())),
              text, flags=re.MULTILINE)
# floats with a decimal point (YAML 1.1 reads 3e-5 as a string)
for key, value in (("alpha", f"{float(alpha):.6e}"), ("rng_seed", str(int(seed)))):
    pattern = rf"^(\s+){key}:[ \t]*[^\s#]+"
    if re.search(pattern, text, re.MULTILINE):
        text = re.sub(pattern, rf"\g<1>{key}: {value}", text, count=1, flags=re.MULTILINE)
    else:
        text = re.sub(r"^(simulation:\s*\n)", rf"\g<1>  {key}: {value}\n", text, count=1, flags=re.MULTILINE)
settings.write_text(text)
PY
grep -nE "^\s+(alpha|rng_seed):" "$out/$settings"
cd "$out"
echo "alpha=$alpha seed=$seed start $(date -u)" > b5.info
# a failed run is a result (B5 counts failures): record it, do not abort
if bash run_tier_bc.sh backmap > run.out 2>&1; then rc=0; else rc=$?; fi
echo "rc=$rc end $(date -u)" >> b5.info
