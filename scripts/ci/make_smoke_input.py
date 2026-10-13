"""Shorten a generated backmapping input for the release smoke test.

The `robust` protocol written by backmap-prep has fixed stage lengths (relaxation, ramp, staged
NVT, production). This script divides every ``run N`` by ``--scale``, multiplies the lambda rate
``alpha`` of the backmap include by the same factor (the ramp then still reaches lambda = 1), prints
``EQ done`` / ``RAMP done`` / ``PRODUCTION done`` after the relaxation, ramp and last run, and
samples the thermo output every 100 steps.

    make_smoke_input.py in.pet pet.backmap.lmp --scale 10 --out in.pet.smoke --out-backmap pet.backmap.smoke.lmp
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

RUN = re.compile(r"^run\s+(\d+)\s*$")
ALPHA = re.compile(r"(\balpha\s+)([0-9.eE+-]+)")
MIN_STEPS = 50


def shorten(text: str, backmap_name: str, smoke_backmap_name: str, scale: int) -> str:
    """Return the shortened input."""
    lines = text.splitlines()
    run_idx = [i for i, ln in enumerate(lines) if RUN.match(ln)]
    if len(run_idx) < 3:
        raise ValueError(f"expected at least three run commands, found {len(run_idx)}")
    out: list[str] = []
    for i, line in enumerate(lines):
        m = RUN.match(line)
        if m:
            steps = max(MIN_STEPS, int(m.group(1)) // scale)
            out.append(f"run {steps}")
            if i == run_idx[0]:
                out.append('print "EQ done"')
            elif i == run_idx[1]:
                out.append('print "RAMP done"')
            elif i == run_idx[-1]:
                out.append('print "PRODUCTION done"')
            continue
        if line.startswith("include ") and line.split()[1] == backmap_name:
            line = f"include {smoke_backmap_name}"
        elif re.match(r"^thermo\s+\d+\s*$", line):
            line = "thermo 100"
        out.append(line)
    return "\n".join(out) + "\n"


def scale_alpha(text: str, scale: int) -> str:
    if not ALPHA.search(text):
        raise ValueError("no alpha in the backmap include")
    return ALPHA.sub(lambda m: f"{m.group(1)}{float(m.group(2)) * scale:.10g}", text, count=1)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("input")
    ap.add_argument("backmap")
    ap.add_argument("--scale", type=int, default=10, help="divide every run length by this factor")
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-backmap", required=True)
    args = ap.parse_args()
    if args.scale < 1:
        raise SystemExit("--scale must be at least 1")
    src, bm = Path(args.input), Path(args.backmap)
    text = shorten(src.read_text(), bm.name, Path(args.out_backmap).name, args.scale)
    Path(args.out).write_text(text)
    Path(args.out_backmap).write_text(scale_alpha(bm.read_text(), args.scale))
    print(f"scale={args.scale}")


if __name__ == "__main__":
    main()
