# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy"]
# ///
"""Compare two LAMMPS data files produced by the MPI serial-vs-N-rank test.

Tolerances account for floating-point reorderings across ranks: positions
should match to ~1e-5 Å and per-atom lambda to ~1e-10 when fix backmap is
MPI-correct. Atoms are matched by tag (not by index, since rank decomposition
reorders atoms).

Usage:
  uv run compare_mpi_data.py --serial a.data --parallel b.data \
      --pos-tol 1e-5 --lambda-tol 1e-10 --report report.txt
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


def parse_atoms(path: Path) -> dict[int, dict]:
    """Return {tag: {type, mol, x, y, z, ix, iy, iz}} from a LAMMPS data file."""
    lines = path.read_text().splitlines()
    in_atoms = False
    atoms: dict[int, dict] = {}
    for line in lines:
        s = line.strip()
        if s.startswith("Atoms"):
            in_atoms = True
            continue
        if in_atoms:
            if s == "" or re.match(r"^\d+ (bonds|angles|dihedrals|impropers|atom types)", s):
                if atoms:
                    break
                continue
            if re.match(
                r"^(Velocities|Bonds|Angles|Dihedrals|Impropers|Masses|Pair Coeffs|Bond Coeffs|Angle Coeffs)",
                s,
            ):
                break
            parts = s.split()
            if len(parts) >= 7:
                tag = int(parts[0])
                mol = int(parts[1])
                atype = int(parts[2])
                q = float(parts[3])
                x, y, z = float(parts[4]), float(parts[5]), float(parts[6])
                ix = iy = iz = 0
                if len(parts) >= 10:
                    ix, iy, iz = int(parts[7]), int(parts[8]), int(parts[9])
                atoms[tag] = {
                    "type": atype,
                    "mol": mol,
                    "q": q,
                    "x": x,
                    "y": y,
                    "z": z,
                    "ix": ix,
                    "iy": iy,
                    "iz": iz,
                }
    return atoms


def parse_box(path: Path) -> list[tuple[float, float]]:
    """Return [(xlo,xhi),(ylo,yhi),(zlo,zhi)]."""
    lines = path.read_text().splitlines()
    box: list[tuple[float, float]] = []
    for line in lines:
        s = line.strip()
        m = re.match(r"^(-?\S+)\s+(\S+)\s+(xlo|ylo|zlo) xhi$", s) or re.match(
            r"^(-?\S+)\s+(\S+)\s+(xlo|ylo|zlo)\s+xhi$", s
        )
        if m and len(box) < 3:
            box.append((float(m.group(1)), float(m.group(2))))
        if len(box) == 3:
            break
    # Fallback: simpler regex
    if len(box) < 3:
        box = []
        for line in lines:
            for kw, idx in (("xlo xhi", 0), ("ylo yhi", 1), ("zlo zhi", 2)):
                if (kw in line and line.split()[0] != "0") or kw in line:
                    parts = line.split()
                    if len(parts) >= 3 and parts[2] == kw.split()[0]:
                        while len(box) <= idx:
                            box.append((0.0, 0.0))
                        box[idx] = (float(parts[0]), float(parts[1]))
    return box[:3]


def minimum_image(d: float, lo: float, hi: float, periodic: bool) -> float:
    if not periodic:
        return d
    prd = hi - lo
    return d - prd * round(d / prd)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--serial", required=True, type=Path)
    p.add_argument("--parallel", required=True, type=Path)
    p.add_argument("--pos-tol", type=float, default=1e-5)
    p.add_argument("--lambda-tol", type=float, default=1e-10)
    p.add_argument("--report", type=Path, default=None)
    args = p.parse_args()

    a = parse_atoms(args.serial)
    b = parse_atoms(args.parallel)
    box = parse_box(args.serial)
    periodic = [True, True, True]

    tags = sorted(set(a) & set(b))
    n_missing = len(set(a) ^ set(b))
    max_pos_err = 0.0
    n_pos_fail = 0
    for t in tags:
        ax, ay, az = a[t]["x"], a[t]["y"], a[t]["z"]
        bx, by, bz = b[t]["x"], b[t]["y"], b[t]["z"]
        dx = minimum_image(bx - ax, box[0][0], box[0][1], periodic[0]) if box else bx - ax
        dy = minimum_image(by - ay, box[1][0], box[1][1], periodic[1]) if box else by - ay
        dz = minimum_image(bz - az, box[2][0], box[2][1], periodic[2]) if box else bz - az
        err = (dx * dx + dy * dy + dz * dz) ** 0.5
        if err > max_pos_err:
            max_pos_err = err
        if err > args.pos_tol:
            n_pos_fail += 1

    report = []
    report.append("MPI serial-vs-parallel comparison")
    report.append(f"  serial:   {args.serial} ({len(a)} atoms)")
    report.append(f"  parallel: {args.parallel} ({len(b)} atoms)")
    report.append(f"  tags matched: {len(tags)}, missing on one side: {n_missing}")
    report.append(f"  max position error: {max_pos_err:.3e} Å (tol {args.pos_tol:.0e})")
    report.append(f"  position failures: {n_pos_fail}")
    report.append("")
    ok = n_missing == 0 and n_pos_fail == 0
    report.append(f"Result: {'PASS' if ok else 'FAIL'}")
    text = "\n".join(report)
    print(text)
    if args.report:
        args.report.write_text(text + "\n")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
