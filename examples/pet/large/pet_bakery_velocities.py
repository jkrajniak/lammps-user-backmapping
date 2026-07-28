#!/usr/bin/env python3
"""Write bakery-faithful initial velocities for PET hybrid data.

Bakery start_backmapping.py: Maxwell-Boltzmann on CG masses; every AT atom
in the bead copies that CG velocity (zero intra-bead relative KE).

LAMMPS pet.data uses molecule ID == CG bead (1 CG + N AT per mol).

Usage:
  python3 pet_bakery_velocities.py pet.data pet_vel.bakery.data [T_K] [seed]
  # then in LAMMPS after read_data:
  #   read_dump pet_vel.bakery.data 0 x y z vx vy vz box yes replace yes
  # or merge Velocities section into a new data file.
"""

from __future__ import annotations

import math
import random
import sys
from collections import defaultdict
from pathlib import Path

# real units: mass g/mol, T K, v Angstrom/fs
# v_rms component: sqrt(kB*T/m) with kB=0.0019872041 kcal/mol/K
# and conversion: energy -> velocity via m*v^2 ~ kT
# LAMMPS real: v in Angstrom/fs; m in g/mol;
# factor: v = sqrt(T/m) * 0.000911337 (approx from LAMMPS velocity create)
# Use same as LAMMPS: sqrt(k_B T / m) with k_B in real units.
KB_REAL = 0.0019872041  # kcal/mol/K
# Convert kcal/mol -> (g/mol)*(Ang/fs)^2 : 1 kcal/mol = 4184 J/mol;
# 1 J = 1 kg m^2/s^2; result: v(Ang/fs) = sqrt(2*E/m) * factor
# Standard LAMMPS: v = sqrt(kT/m) * sqrt(4.184e-4) roughly.
# Use: sigma = sqrt(KB_REAL * T / mass) * 0.020454828 (from LAMMPS src)
VEL_FACTOR = 0.02045482882  # matches LAMMPS velocity create (real)


CG_TYPES = {1, 6, 7, 11, 14, 16}


def parse_data(path: Path) -> tuple[dict, dict, list[str]]:
    masses: dict[int, float] = {}
    atoms: dict[int, tuple[int, int]] = {}  # id -> (type, mol)
    lines = path.read_text().splitlines()
    section = None
    for line in lines:
        s = line.strip()
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
    return masses, atoms, lines


def bakery_velocities(
    masses: dict[int, float],
    atoms: dict[int, tuple[int, int]],
    temperature: float,
    seed: int,
) -> dict[int, tuple[float, float, float]]:
    rng = random.Random(seed)
    by_mol: dict[int, list[int]] = defaultdict(list)
    mol_cg: dict[int, int] = {}
    for aid, (atype, mol) in atoms.items():
        by_mol[mol].append(aid)
        if atype in CG_TYPES:
            mol_cg[mol] = aid

    vel: dict[int, tuple[float, float, float]] = {}
    for mol, aids in by_mol.items():
        cg = mol_cg.get(mol)
        if cg is None:
            # no CG — zero
            for aid in aids:
                vel[aid] = (0.0, 0.0, 0.0)
            continue
        m = masses[atoms[cg][0]]
        sigma = math.sqrt(KB_REAL * temperature / m) * VEL_FACTOR
        vx = rng.gauss(0.0, sigma)
        vy = rng.gauss(0.0, sigma)
        vz = rng.gauss(0.0, sigma)
        for aid in aids:
            vel[aid] = (vx, vy, vz)
    return vel


def write_velocity_data(
    path: Path, n_atoms: int, vel: dict[int, tuple[float, float, float]]
) -> None:
    with path.open("w") as f:
        f.write("ITEM: TIMESTEP\n0\n")
        f.write("ITEM: NUMBER OF ATOMS\n")
        f.write(f"{n_atoms}\n")
        f.write("ITEM: BOX BOUNDS pp pp pp\n")
        f.write("0.0 71.3297000\n0.0 71.3297000\n0.0 71.3297000\n")
        f.write("ITEM: ATOMS id vx vy vz\n")
        for aid in range(1, n_atoms + 1):
            vx, vy, vz = vel[aid]
            f.write(f"{aid} {vx:.10e} {vy:.10e} {vz:.10e}\n")


def main() -> None:
    data = Path(sys.argv[1] if len(sys.argv) > 1 else "pet.data")
    out = Path(sys.argv[2] if len(sys.argv) > 2 else "pet_vel.bakery.dump")
    temperature = float(sys.argv[3]) if len(sys.argv) > 3 else 298.0
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 48279
    masses, atoms, _ = parse_data(data)
    vel = bakery_velocities(masses, atoms, temperature, seed)
    write_velocity_data(out, len(atoms), vel)
    # quick check: unique velocities == n beads
    uniq = len(set(vel.values()))
    print(f"wrote {out}: {len(vel)} atoms, {uniq} unique velocity vectors (expect ~n_beads)")


if __name__ == "__main__":
    main()
