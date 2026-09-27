"""Compare the B1 static MPI parity runs of one system with the 1-rank run.

Reads ``forces_<name>.dump`` and ``log.b1_<name>`` written by b1_mpi_parity.sh
in ``workdir`` and reports, per run, the largest deviation from the reference
run of:

- per-atom total force (AT atoms after the CG force redistribution, CG beads),
  relative to the largest force magnitude in the system;
- bead COM positions (fix backmap peratom full, columns 2-4), absolute, in A;
- bead CG forces before redistribution (columns 5-7), relative as above;
- each energy term, relative to its magnitude.

Deviations near 1e-12 relative are summation-order round-off.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

ENERGIES = ("pe", "evdwl", "ecoul", "ebond", "eangle", "edihed", "eimp")


def read_dump(path: Path) -> tuple[list[str], np.ndarray]:
    lines = path.read_text().splitlines()
    head = next(i for i, ln in enumerate(lines) if ln.startswith("ITEM: ATOMS"))
    cols = lines[head].split()[2:]
    data = np.loadtxt(lines[head + 1 :], ndmin=2)
    return cols, data[np.argsort(data[:, 0])]


def read_energies(path: Path) -> dict[str, float]:
    match = re.search(r"^RESULT (.*)$", path.read_text(), re.MULTILINE)
    if match is None:
        raise ValueError(f"no RESULT line in {path}")
    return {k: float(v) for k, v in (kv.split("=") for kv in match.group(1).split())}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("workdir", type=Path)
    parser.add_argument("--ref", default="r1")
    parser.add_argument("--tol", type=float, default=1e-10, help="relative force tolerance")
    args = parser.parse_args()

    cols, ref = read_dump(args.workdir / f"forces_{args.ref}.dump")
    idx = {c: i for i, c in enumerate(cols)}
    f_cols = [idx["fx"], idx["fy"], idx["fz"]]
    com_cols = [idx["f_bm[2]"], idx["f_bm[3]"], idx["f_bm[4]"]]
    fcg_cols = [idx["f_bm[5]"], idx["f_bm[6]"], idx["f_bm[7]"]]
    fix_line = next(
        ln
        for bm in args.workdir.glob("*.backmap.lmp")
        for ln in bm.read_text().splitlines()
        if ln.startswith("fix bm all backmap")
    )
    cg_types = {int(t) for t in re.search(r"cg_type ((?:\d+ ?)+)", fix_line).group(1).split()}
    is_bead = np.isin(ref[:, idx["type"]].astype(int), sorted(cg_types))
    f_scale = np.max(np.linalg.norm(ref[:, f_cols], axis=1))
    fcg_scale = np.max(np.linalg.norm(ref[is_bead][:, fcg_cols], axis=1)) if is_bead.any() else 1.0
    e_ref = read_energies(args.workdir / f"log.b1_{args.ref}")

    print(
        f"system {args.workdir.name}: {len(ref)} atoms, {int(is_bead.sum())} beads, "
        f"max |f| {f_scale:.6g} kcal/mol/A, lambda {ref[0, idx['f_bm[1]']]:g}"
    )
    print(f"{'run':6s} {'force rel':>10s} {'COM abs A':>10s} {'CG f rel':>10s} {'energy rel':>10s}")
    worst = 0.0
    for dump in sorted(args.workdir.glob("forces_*.dump")):
        name = dump.stem.removeprefix("forces_")
        if name == args.ref:
            continue
        _, run = read_dump(dump)
        if not np.array_equal(run[:, 0], ref[:, 0]):
            print(f"{name}: atom ids differ from the reference run", file=sys.stderr)
            return 1
        d_f = np.max(np.abs(run[:, f_cols] - ref[:, f_cols])) / f_scale
        d_com = (
            np.max(np.abs(run[is_bead][:, com_cols] - ref[is_bead][:, com_cols]))
            if is_bead.any()
            else 0.0
        )
        d_fcg = (
            np.max(np.abs(run[is_bead][:, fcg_cols] - ref[is_bead][:, fcg_cols])) / fcg_scale
            if is_bead.any()
            else 0.0
        )
        e_run = read_energies(args.workdir / f"log.b1_{name}")
        d_e = max(
            abs(e_run[k] - e_ref[k]) / max(abs(e_ref[k]), 1e-300)
            for k in ENERGIES
            if e_ref[k] != 0.0
        )
        worst = max(worst, d_f, d_fcg)
        print(f"{name:6s} {d_f:10.2e} {d_com:10.2e} {d_fcg:10.2e} {d_e:10.2e}")
    verdict = "PASS" if worst <= args.tol else "FAIL"
    print(f"{verdict}: largest relative force deviation {worst:.2e} (tolerance {args.tol:.0e})")
    return 0 if verdict == "PASS" else 2


if __name__ == "__main__":
    sys.exit(main())
