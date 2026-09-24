#!/usr/bin/env python3
"""Rewrite LAMMPS image flags so every bond is the minimum-image bond.

Coordinates stay in the primary cell. Flags are assigned by walking each
molecule's bond tree. A bond whose flags differ by one box vector while the
atoms sit a bond length apart in the cell is what made the 4-rank polyethylene
run diverge.
"""

from __future__ import annotations

import argparse
from collections import defaultdict, deque
from pathlib import Path


def _min_image(delta: float, length: float) -> float:
    if length == 0:
        return delta
    return delta - length * round(delta / length)


def _parse(
    text: str,
) -> tuple[
    list[str],
    dict[int, tuple[float, float, float]],
    list[tuple[int, int]],
    tuple[float, float, float],
]:
    box = [0.0, 0.0, 0.0]
    atoms: dict[int, tuple[float, float, float]] = {}
    bonds: list[tuple[int, int]] = []
    section: str | None = None
    lines = text.splitlines(keepends=True)
    for line in lines:
        stripped = line.strip()
        if stripped.endswith("xlo xhi"):
            box[0] = float(stripped.split()[1]) - float(stripped.split()[0])
        elif stripped.endswith("ylo yhi"):
            box[1] = float(stripped.split()[1]) - float(stripped.split()[0])
        elif stripped.endswith("zlo zhi"):
            box[2] = float(stripped.split()[1]) - float(stripped.split()[0])
        if stripped.startswith("Atoms"):
            section = "atoms"
            continue
        if stripped.startswith("Bonds"):
            section = "bonds"
            continue
        if stripped.startswith(("Angles", "Dihedrals", "Velocities", "Masses")):
            section = None
            continue
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if section == "atoms" and len(parts) >= 10:
            atoms[int(parts[0])] = tuple(float(value) for value in parts[4:7])
        elif section == "bonds" and len(parts) >= 4:
            bonds.append((int(parts[2]), int(parts[3])))
    return lines, atoms, bonds, (box[0], box[1], box[2])


def _assign_flags(
    atoms: dict[int, tuple[float, float, float]],
    bonds: list[tuple[int, int]],
    box: tuple[float, float, float],
) -> dict[int, tuple[int, int, int]]:
    adjacency: dict[int, list[int]] = defaultdict(list)
    for atom_i, atom_j in bonds:
        adjacency[atom_i].append(atom_j)
        adjacency[atom_j].append(atom_i)
    flags = dict.fromkeys(atoms, (0, 0, 0))
    unwrapped = {atom_id: atoms[atom_id] for atom_id in atoms}
    visited: set[int] = set()
    for root in sorted(atoms):
        if root in visited:
            continue
        queue: deque[int] = deque([root])
        visited.add(root)
        while queue:
            current = queue.popleft()
            cx, cy, cz = unwrapped[current]
            sx, sy, sz = atoms[current]
            for neighbor in adjacency[current]:
                if neighbor in visited:
                    continue
                visited.add(neighbor)
                queue.append(neighbor)
                nx, ny, nz = atoms[neighbor]
                delta = (
                    _min_image(nx - sx, box[0]),
                    _min_image(ny - sy, box[1]),
                    _min_image(nz - sz, box[2]),
                )
                unwrapped[neighbor] = (cx + delta[0], cy + delta[1], cz + delta[2])
    for atom_id, (x, y, z) in atoms.items():
        ux, uy, uz = unwrapped[atom_id]
        flags[atom_id] = (
            round((ux - x) / box[0]) if box[0] else 0,
            round((uy - y) / box[1]) if box[1] else 0,
            round((uz - z) / box[2]) if box[2] else 0,
        )
    return flags


def _max_flag_bond(
    atoms: dict[int, tuple[float, float, float]],
    bonds: list[tuple[int, int]],
    flags: dict[int, tuple[int, int, int]],
    box: tuple[float, float, float],
) -> float:
    worst = 0.0
    for atom_i, atom_j in bonds:
        xi, yi, zi = atoms[atom_i]
        xj, yj, zj = atoms[atom_j]
        ixi, iyi, izi = flags[atom_i]
        ixj, iyj, izj = flags[atom_j]
        dx = (xj + ixj * box[0]) - (xi + ixi * box[0])
        dy = (yj + iyj * box[1]) - (yi + iyi * box[1])
        dz = (zj + izj * box[2]) - (zi + izi * box[2])
        worst = max(worst, (dx * dx + dy * dy + dz * dz) ** 0.5)
    return worst


def rewrite(path: Path, destination: Path) -> float:
    lines, atoms, bonds, box = _parse(path.read_text())
    flags = _assign_flags(atoms, bonds, box)
    section: str | None = None
    rewritten: list[str] = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("Atoms"):
            section = "atoms"
            rewritten.append(line)
            continue
        if stripped.startswith(("Bonds", "Angles", "Dihedrals", "Velocities", "Masses")):
            section = None
            rewritten.append(line)
            continue
        parts = stripped.split()
        if section == "atoms" and len(parts) >= 10 and parts[0].isdigit():
            ix, iy, iz = flags[int(parts[0])]
            prefix = " ".join(parts[:7])
            rewritten.append(f"{prefix} {ix} {iy} {iz}\n")
            continue
        rewritten.append(line)
    destination.write_text("".join(rewritten))
    return _max_flag_bond(atoms, bonds, flags, box)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    worst = rewrite(args.source, args.destination)
    print(f"longest flag-unwrapped bond: {worst:.3f} A")
    if worst > 5.0:
        raise SystemExit("image flags still leave a bond longer than 5 A")


if __name__ == "__main__":
    main()
