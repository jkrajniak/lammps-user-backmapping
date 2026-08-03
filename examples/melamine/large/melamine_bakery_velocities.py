#!/usr/bin/env python3
"""Write bakery-faithful initial velocities for melamine hybrid data.

Same scheme as examples/pet/large/pet_bakery_velocities.py (see
research/decisions/2026-07-19-prefer-bakery-protocol-no-frozen-cg.md and
research/checkpoints.md): one Maxwell-Boltzmann draw per CG bead, on the
CG mass; every AT atom in that bead copies the CG velocity (zero
intra-bead relative KE at t=0). Box bounds are read from the data file
instead of hardcoded, unlike the PET script this was adapted from.

Unlike PET (molecule ID == one CG bead), melamine's molecule ID covers a
whole physical molecule with 3 CG beads (A1/A2/A3) + 27 AT atoms each.
Bead membership isn't in the LAMMPS data file directly, so this relies on
backmap-prep's atom ordering within each molecule: CG atoms first (in
bead order), then AT atoms grouped contiguously by bead in the same
order, each bead's AT block the same size (verified: every molecule here
has exactly 3 CG atoms + 27 AT atoms, 9 per bead) -- confirmed against
melamine.data before relying on it.

Usage:
  python3 melamine_bakery_velocities.py melamine.data melamine_vel.bakery.dump [T_K] [seed]
  # then in LAMMPS after read_data:
  #   read_dump melamine_vel.bakery.dump 0 vx vy vz box no replace yes format native
"""

from __future__ import annotations

import math
import random
import sys
from collections import defaultdict
from pathlib import Path

KB_REAL = 0.0019872041  # kcal/mol/K
VEL_FACTOR = 0.02045482882  # matches LAMMPS velocity create (real)

CG_TYPES = {1}


def parse_data(path: Path) -> tuple[dict, dict, tuple[str, str, str]]:
    masses: dict[int, float] = {}
    atoms: dict[int, tuple[int, int]] = {}  # id -> (type, mol)
    box: list[str] = []
    lines = path.read_text().splitlines()
    section = None
    for line in lines:
        s = line.strip()
        if "xlo xhi" in s or "ylo yhi" in s or "zlo zhi" in s:
            box.append(" ".join(s.split()[:2]))
            continue
        if not s:
            continue
        if s == "Masses":
            section = "Masses"
            continue
        if s.startswith("Atoms"):
            section = "Atoms"
            continue
        if s[0].isalpha() and section in {"Masses", "Atoms"}:
            section = None
        if section == "Masses":
            p = s.split()
            if p[0].isdigit():
                masses[int(p[0])] = float(p[1])
        elif section == "Atoms":
            p = s.split()
            if len(p) >= 3 and p[0].isdigit():
                atoms[int(p[0])] = (int(p[2]), int(p[1]))
    if len(box) != 3:
        raise ValueError(f"expected 3 box-bound lines, found {len(box)}")
    return masses, atoms, (box[0], box[1], box[2])


def bakery_velocities(
    masses: dict[int, float],
    atoms: dict[int, tuple[int, int]],
    temperature: float,
    seed: int,
) -> dict[int, tuple[float, float, float]]:
    rng = random.Random(seed)
    by_mol: dict[int, list[int]] = defaultdict(list)
    for aid in atoms:
        by_mol[atoms[aid][1]].append(aid)

    vel: dict[int, tuple[float, float, float]] = {}
    for _mol, aids in by_mol.items():
        aids.sort()
        cg_ids = [a for a in aids if atoms[a][0] in CG_TYPES]
        at_ids = [a for a in aids if atoms[a][0] not in CG_TYPES]
        if not cg_ids:
            for aid in aids:
                vel[aid] = (0.0, 0.0, 0.0)
            continue
        n_beads = len(cg_ids)
        if len(at_ids) % n_beads != 0:
            raise ValueError(
                f"mol with {n_beads} CG beads has {len(at_ids)} AT atoms, "
                "not evenly divisible -- bead grouping assumption broken"
            )
        chunk = len(at_ids) // n_beads
        for i, cg in enumerate(cg_ids):
            bead_ids = [cg, *at_ids[i * chunk : (i + 1) * chunk]]
            m = masses[atoms[cg][0]]
            sigma = math.sqrt(KB_REAL * temperature / m) * VEL_FACTOR
            vx = rng.gauss(0.0, sigma)
            vy = rng.gauss(0.0, sigma)
            vz = rng.gauss(0.0, sigma)
            for aid in bead_ids:
                vel[aid] = (vx, vy, vz)
    return vel


def write_velocity_data(
    path: Path,
    n_atoms: int,
    vel: dict[int, tuple[float, float, float]],
    box: tuple[str, str, str],
) -> None:
    with path.open("w") as f:
        f.write("ITEM: TIMESTEP\n0\n")
        f.write("ITEM: NUMBER OF ATOMS\n")
        f.write(f"{n_atoms}\n")
        f.write("ITEM: BOX BOUNDS pp pp pp\n")
        for b in box:
            f.write(f"{b}\n")
        f.write("ITEM: ATOMS id vx vy vz\n")
        for aid in range(1, n_atoms + 1):
            vx, vy, vz = vel[aid]
            f.write(f"{aid} {vx:.10e} {vy:.10e} {vz:.10e}\n")


def main() -> None:
    data = Path(sys.argv[1] if len(sys.argv) > 1 else "melamine.data")
    out = Path(sys.argv[2] if len(sys.argv) > 2 else "melamine_vel.bakery.dump")
    temperature = float(sys.argv[3]) if len(sys.argv) > 3 else 300.0
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 12345
    masses, atoms, box = parse_data(data)
    vel = bakery_velocities(masses, atoms, temperature, seed)
    write_velocity_data(out, len(atoms), vel, box)
    uniq = len(set(vel.values()))
    print(f"wrote {out}: {len(vel)} atoms, {uniq} unique velocity vectors (expect ~n_beads)")


if __name__ == "__main__":
    main()
