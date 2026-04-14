#!/usr/bin/env python3
"""Extract AT-only atoms from a backmapped hybrid LAMMPS data file.

Reads a LAMMPS data file containing both CG and AT atoms after backmapping,
filters out CG atoms (types matching cg_types), renumbers remaining atoms/bonds/angles,
and writes a clean pure-AT data file for standalone production runs.

Usage:
    uv run extract_at_frame.py dodecane_hybrid.data dodecane_at.data --cg-types 1 2
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class LAMMPSData:
    header: str = ""
    counts: dict[str, int] = field(default_factory=dict)
    box_bounds: list[tuple[float, float]] = field(default_factory=list)
    masses: list[tuple[int, float, str]] = field(default_factory=list)
    atoms: list[list[str]] = field(default_factory=list)
    bonds: list[list[str]] = field(default_factory=list)
    angles: list[list[str]] = field(default_factory=list)
    dihedrals: list[list[str]] = field(default_factory=list)
    impropers: list[list[str]] = field(default_factory=list)


def parse_data_file(path: Path) -> LAMMPSData:
    """Parse a LAMMPS data file into structured data."""
    data = LAMMPSData()
    lines = path.read_text().splitlines()

    # Header
    data.header = lines[0] if lines else ""

    # Counts — skip blank lines after header
    i = 1
    while i < len(lines) and not lines[i].strip():
        i += 1

    count_re = re.compile(r"^(\d+)\s+(\w+(?:\s+\w+)?)$")
    while i < len(lines):
        m = count_re.match(lines[i].strip())
        if not m:
            break
        data.counts[m.group(2)] = int(m.group(1))
        i += 1

    # Skip blank lines before box bounds
    while i < len(lines) and not lines[i].strip():
        i += 1

    # Box bounds
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            break
        if any(kw in line for kw in ("xlo", "ylo", "zlo")):
            parts = line.split()
            data.box_bounds.append((float(parts[0]), float(parts[1])))
            i += 1
        else:
            break

    # Skip blank line
    while i < len(lines) and not lines[i].strip():
        i += 1

    # Sections — include all possible LAMMPS data file sections
    section_map = {
        "Masses": "masses",
        "Atoms": "atoms",
        "Bonds": "bonds",
        "Angles": "angles",
        "Dihedrals": "dihedrals",
        "Impropers": "impropers",
        "Velocities": None,
        "Pair Coeffs": None,
        "Bond Coeffs": None,
        "Angle Coeffs": None,
        "Dihedral Coeffs": None,
        "Improper Coeffs": None,
        "Special Bonds Coeffs": None,
    }

    current_section = None
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue

        # Check for section header
        found_section = False
        for keyword, attr in section_map.items():
            if line.startswith(keyword):
                current_section = attr
                found_section = True
                i += 1
                break
        if not found_section:
            if current_section is not None:
                getattr(data, current_section).append(line.split())
            i += 1

    return data


def extract_at(
    data: LAMMPSData,
    cg_types: set[int],
) -> LAMMPSData:
    """Filter out CG atoms and renumber everything."""
    out = LAMMPSData(header=data.header + " (AT-only, extracted from backmapped system)")
    out.box_bounds = list(data.box_bounds)

    # Identify CG atom IDs
    cg_atom_ids: set[int] = set()
    at_atom_ids: list[int] = []
    old_to_new: dict[int, int] = {}

    for atom_line in data.atoms:
        atom_id = int(atom_line[0])
        atom_type = int(atom_line[2])
        if atom_type in cg_types:
            cg_atom_ids.add(atom_id)
        else:
            at_atom_ids.append(atom_id)

    at_atom_ids.sort()
    for new_id, old_id in enumerate(at_atom_ids, start=1):
        old_to_new[old_id] = new_id

    # Remap masses — keep only AT types
    new_type_map: dict[int, int] = {}
    new_type_id = 1
    for m in data.masses:
        tid = int(m[0])
        mass = float(m[1])
        comment = m[2] if len(m) > 2 else ""
        if tid not in cg_types:
            new_type_map[tid] = new_type_id
            out.masses.append((new_type_id, mass, comment))
            new_type_id += 1

    # Count remaining types
    n_atom_types = len(out.masses)

    # Remap atoms
    for atom_line in data.atoms:
        old_id = int(atom_line[0])
        if old_id in cg_atom_ids:
            continue
        new_line = list(atom_line)
        new_line[0] = str(old_to_new[old_id])
        new_line[2] = str(new_type_map[int(atom_line[2])])
        out.atoms.append(new_line)

    # Remap bonds
    n_bonds = 0
    bond_type_map: dict[int, int] = {}
    new_bond_type_id = 1
    for bond_line in data.bonds:
        bond_type = int(bond_line[1])
        a1 = int(bond_line[2])
        a2 = int(bond_line[3])
        if a1 in cg_atom_ids or a2 in cg_atom_ids:
            continue
        if bond_type not in bond_type_map:
            bond_type_map[bond_type] = new_bond_type_id
            new_bond_type_id += 1
        new_line = [
            str(n_bonds + 1),
            str(bond_type_map[bond_type]),
            str(old_to_new[a1]),
            str(old_to_new[a2]),
        ]
        out.bonds.append(new_line)
        n_bonds += 1

    # Remap angles
    n_angles = 0
    angle_type_map: dict[int, int] = {}
    new_angle_type_id = 1
    for angle_line in data.angles:
        angle_type = int(angle_line[1])
        a1 = int(angle_line[2])
        a2 = int(angle_line[3])
        a3 = int(angle_line[4])
        if a1 in cg_atom_ids or a2 in cg_atom_ids or a3 in cg_atom_ids:
            continue
        if angle_type not in angle_type_map:
            angle_type_map[angle_type] = new_angle_type_id
            new_angle_type_id += 1
        new_line = [
            str(n_angles + 1),
            str(angle_type_map[angle_type]),
            str(old_to_new[a1]),
            str(old_to_new[a2]),
            str(old_to_new[a3]),
        ]
        out.angles.append(new_line)
        n_angles += 1

    # Update counts
    out.counts["atoms"] = len(out.atoms)
    out.counts["bonds"] = n_bonds
    out.counts["angles"] = n_angles
    out.counts["atom types"] = n_atom_types
    out.counts["bond types"] = len(bond_type_map)
    out.counts["angle types"] = len(angle_type_map)
    if "dihedrals" in data.counts:
        out.counts["dihedrals"] = 0
    if "impropers" in data.counts:
        out.counts["impropers"] = 0

    return out


def write_data_file(data: LAMMPSData, path: Path) -> None:
    """Write a LAMMPS data file from structured data."""
    lines = [data.header, ""]

    # Counts
    for key, val in data.counts.items():
        lines.append(f"{val} {key}")
    lines.append("")

    # Box bounds
    labels = ["xlo xhi", "ylo yhi", "zlo zhi"]
    for idx, (lo, hi) in enumerate(data.box_bounds):
        label = labels[idx] if idx < len(labels) else "lo hi"
        lines.append(f"{lo} {hi} {label}")
    lines.append("")

    # Masses
    if data.masses:
        lines.append("Masses")
        lines.append("")
        for tid, mass, comment in data.masses:
            if comment:
                lines.append(f"{tid} {mass} {comment}")
            else:
                lines.append(f"{tid} {mass}")
        lines.append("")

    # Atoms
    if data.atoms:
        lines.append("Atoms # full")
        lines.append("")
        for atom_line in data.atoms:
            lines.append(" ".join(atom_line))
        lines.append("")

    # Bonds
    if data.bonds:
        lines.append("Bonds")
        lines.append("")
        for bond_line in data.bonds:
            lines.append(" ".join(bond_line))
        lines.append("")

    # Angles
    if data.angles:
        lines.append("Angles")
        lines.append("")
        for angle_line in data.angles:
            lines.append(" ".join(angle_line))
        lines.append("")

    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Extract AT-only atoms from a backmapped hybrid LAMMPS data file"
    )
    parser.add_argument("input", type=Path, help="Hybrid LAMMPS data file")
    parser.add_argument("output", type=Path, help="Output AT-only data file")
    parser.add_argument(
        "--cg-types",
        type=int,
        nargs="+",
        default=[1, 2],
        help="CG atom type IDs to remove (default: 1 2)",
    )
    args = parser.parse_args()

    if not args.input.exists():
        print(f"ERROR: {args.input} not found", file=sys.stderr)
        return 1

    cg_types = set(args.cg_types)
    print(f"Reading {args.input} ...")
    data = parse_data_file(args.input)
    print(
        f"  {data.counts.get('atoms', 0)} atoms, {data.counts.get('bonds', 0)} bonds, "
        f"{data.counts.get('angles', 0)} angles"
    )
    print(f"  CG types to remove: {cg_types}")

    at_data = extract_at(data, cg_types)
    print(
        f"  AT-only: {at_data.counts['atoms']} atoms, "
        f"{at_data.counts['bonds']} bonds, {at_data.counts['angles']} angles"
    )

    write_data_file(at_data, args.output)
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
