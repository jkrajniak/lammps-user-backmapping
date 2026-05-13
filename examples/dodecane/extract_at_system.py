# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Extract a pure atomistic system from a hybrid CG+AT LAMMPS data file.

Strips CG atoms, renumbers types/IDs, and writes a standalone AT data file
suitable for a reference all-atom simulation.

Usage:
  uv run extract_at_system.py [--input dodecane_final.data] [--output dodecane_at.data] [--cg-types 1 2]
"""

from __future__ import annotations

import argparse
import contextlib
import sys
from pathlib import Path


def parse_header(lines: list[str]) -> dict[str, int]:
    """Parse LAMMPS data file header for counts and box."""
    header: dict[str, int] = {}
    for line in lines:
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        for key in ("atoms", "bonds", "angles", "atom types", "bond types", "angle types"):
            if s.endswith(key):
                header[key] = int(s.split()[0])
                break
    return header


def extract_at_system(
    input_path: Path,
    output_path: Path,
    cg_types: set[int],
) -> None:
    text = input_path.read_text()
    lines = text.splitlines()

    # Parse box dimensions
    box_lines: list[str] = []
    for line in lines:
        for tag in ("xlo", "ylo", "zlo"):
            if tag in line:
                box_lines.append(line.strip())

    # Parse sections
    section = None
    atoms_raw: list[list[str]] = []
    bonds_raw: list[list[str]] = []
    angles_raw: list[list[str]] = []
    masses_raw: dict[int, str] = {}

    for line in lines:
        s = line.strip()
        if not s:
            continue

        if s.startswith("Atoms"):
            section = "atoms"
            continue
        if s.startswith("Velocities"):
            section = "velocities"
            continue
        if s.startswith("Bonds"):
            section = "bonds"
            continue
        if s.startswith("Angles"):
            section = "angles"
            continue
        if s.startswith("Masses"):
            section = "masses"
            continue
        if s.startswith(("Pair Coeffs", "Bond Coeffs", "Angle Coeffs")):
            section = "skip"
            continue

        if section == "masses":
            parts = s.split()
            if len(parts) >= 2:
                with contextlib.suppress(ValueError):
                    masses_raw[int(parts[0])] = s
        elif section == "atoms":
            parts = s.split()
            if len(parts) >= 7:
                atoms_raw.append(parts)
        elif section == "bonds":
            parts = s.split()
            if len(parts) >= 4:
                bonds_raw.append(parts)
        elif section == "angles":
            parts = s.split()
            if len(parts) >= 5:
                angles_raw.append(parts)

    # Filter AT atoms and build ID mapping
    at_atoms = [a for a in atoms_raw if int(a[2]) not in cg_types]
    at_types = sorted({int(a[2]) for a in at_atoms})
    type_map = {old: new for new, old in enumerate(at_types, 1)}
    old_ids = {int(a[0]) for a in at_atoms}
    id_map = {int(a[0]): new for new, a in enumerate(sorted(at_atoms, key=lambda x: int(x[0])), 1)}

    # Renumber molecule IDs to be contiguous
    mol_ids = sorted({int(a[1]) for a in at_atoms})
    mol_map = {old: new for new, old in enumerate(mol_ids, 1)}

    # Filter bonds: keep only AT-AT bonds, renumber
    at_bonds = []
    for b in bonds_raw:
        i, j = int(b[2]), int(b[3])
        if i in old_ids and j in old_ids:
            at_bonds.append(b)

    # Filter angles: keep only AT-AT-AT angles, renumber
    at_angles = []
    for a in angles_raw:
        i, j, k = int(a[2]), int(a[3]), int(a[4])
        if i in old_ids and j in old_ids and k in old_ids:
            at_angles.append(a)

    n_at_types = len(at_types)

    with output_path.open("w") as f:
        f.write(f"LAMMPS data file — pure AT system extracted from {input_path.name}\n\n")
        f.write(f"{len(at_atoms)} atoms\n")
        f.write(f"{len(at_bonds)} bonds\n")
        f.write(f"{len(at_angles)} angles\n")
        f.write("0 dihedrals\n")
        f.write("0 impropers\n\n")
        f.write(f"{n_at_types} atom types\n")
        f.write("1 bond types\n")
        f.write("1 angle types\n")
        f.write("0 dihedral types\n")
        f.write("0 improper types\n\n")
        for bl in box_lines:
            f.write(f"{bl}\n")
        f.write("\n")

        f.write("Masses\n\n")
        for old_type in at_types:
            new_type = type_map[old_type]
            mass_val = masses_raw[old_type].split()[1]
            comment = (
                masses_raw[old_type].split("#")[1].strip() if "#" in masses_raw[old_type] else ""
            )
            if comment:
                f.write(f"{new_type} {mass_val} # {comment}\n")
            else:
                f.write(f"{new_type} {mass_val}\n")
        f.write("\n")

        f.write("Atoms # full\n\n")
        for a in sorted(at_atoms, key=lambda x: int(x[0])):
            old_id = int(a[0])
            new_id = id_map[old_id]
            new_mol = mol_map[int(a[1])]
            new_type = type_map[int(a[2])]
            charge = a[3]
            coords = " ".join(a[4:7])
            image_flags = " ".join(a[7:10]) if len(a) >= 10 else "0 0 0"
            f.write(f"{new_id} {new_mol} {new_type} {charge} {coords} {image_flags}\n")
        f.write("\n")

        f.write("Bonds\n\n")
        for idx, b in enumerate(at_bonds, 1):
            new_i = id_map[int(b[2])]
            new_j = id_map[int(b[3])]
            f.write(f"{idx} 1 {new_i} {new_j}\n")
        f.write("\n")

        f.write("Angles\n\n")
        for idx, a in enumerate(at_angles, 1):
            new_i = id_map[int(a[2])]
            new_j = id_map[int(a[3])]
            new_k = id_map[int(a[4])]
            f.write(f"{idx} 1 {new_i} {new_j} {new_k}\n")

    n_mols = len(mol_ids)
    print(f"Extracted {len(at_atoms)} AT atoms ({n_at_types} types, {n_mols} molecules)")
    print(f"  {len(at_bonds)} bonds, {len(at_angles)} angles")
    print(f"  Type mapping: {dict(type_map)}")
    print(f"  Written to {output_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Extract pure AT system from hybrid data file")
    parser.add_argument("--input", type=Path, default=Path("dodecane_final.data"))
    parser.add_argument("--output", type=Path, default=Path("dodecane_at.data"))
    parser.add_argument(
        "--cg-types", type=int, nargs="+", default=[1, 2], help="CG atom type IDs to strip"
    )
    args = parser.parse_args()

    if not args.input.exists():
        print(f"ERROR: {args.input} not found", file=sys.stderr)
        return 1

    extract_at_system(args.input, args.output, set(args.cg_types))
    return 0


if __name__ == "__main__":
    sys.exit(main())
