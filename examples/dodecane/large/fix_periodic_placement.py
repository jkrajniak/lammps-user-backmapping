# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy"]
# ///
"""Make molecules whole and ensure AT atoms are near their CG beads.

For molecules spanning periodic boundaries, backmap-prep places AT atoms
near their CG bead in minimum-image space, but the raw coordinates can
be on opposite sides of the box. This script:

1. Makes each molecule "whole" by unwrapping all atoms relative to the
   first CG bead (removing periodic jumps within the molecule).
2. Wraps the entire molecule back into the box using the molecule's COM.

After this, all bonded atoms within a molecule are geometrically close,
and fix_backmap can safely track CG-AT distances.
"""

import argparse
from pathlib import Path

import numpy as np


def fix_molecules(data_path: Path, output_path: Path, cg_types: set[int]) -> None:
    lines = data_path.read_text().splitlines()

    box = np.zeros((3, 2))
    prd = np.zeros(3)
    atoms: list[dict] = []

    # Parse
    section = None
    for idx, line in enumerate(lines):
        s = line.strip()
        if "xlo xhi" in s:
            parts = s.split()
            box[0] = [float(parts[0]), float(parts[1])]
        elif "ylo yhi" in s:
            parts = s.split()
            box[1] = [float(parts[0]), float(parts[1])]
        elif "zlo zhi" in s:
            parts = s.split()
            box[2] = [float(parts[0]), float(parts[1])]

        if s.startswith("Atoms"):
            section = "atoms"
            idx + 1
            continue
        if section == "atoms" and (s.startswith(("Velocities", "Bonds"))):
            section = None
            continue
        if section == "atoms" and s and not s.startswith("#"):
            parts = s.split()
            try:
                atoms.append(
                    {
                        "line_idx": idx,
                        "parts": parts,
                        "id": int(parts[0]),
                        "mol": int(parts[1]),
                        "type": int(parts[2]),
                        "x": np.array([float(parts[4]), float(parts[5]), float(parts[6])]),
                    }
                )
            except (ValueError, IndexError):
                continue

    prd = box[:, 1] - box[:, 0]
    half = prd / 2.0
    print(f"Read {len(atoms)} atoms, box = {prd}, lo = {box[:, 0]}")

    # Group by molecule
    by_mol: dict[int, list[int]] = {}
    for i, a in enumerate(atoms):
        by_mol.setdefault(a["mol"], []).append(i)

    n_fixed = 0
    for _mol_id, indices in by_mol.items():
        # Sort by atom ID within molecule
        indices.sort(key=lambda i: atoms[i]["id"])

        # Find the first CG bead as reference
        ref_idx = None
        for i in indices:
            if atoms[i]["type"] in cg_types:
                ref_idx = i
                break
        if ref_idx is None:
            continue

        ref_pos = atoms[ref_idx]["x"].copy()

        # Make molecule whole: unwrap each atom relative to the previous
        # connected atom (using sequential order as proxy for connectivity)
        prev_pos = ref_pos.copy()
        for i in indices:
            if i == ref_idx:
                continue
            pos = atoms[i]["x"]
            delta = pos - prev_pos
            for d in range(3):
                if delta[d] > half[d]:
                    pos[d] -= prd[d]
                    n_fixed += 1
                elif delta[d] < -half[d]:
                    pos[d] += prd[d]
                    n_fixed += 1
            prev_pos = pos.copy()

        # Compute molecule COM and wrap the whole molecule into the box
        positions = np.array([atoms[i]["x"] for i in indices])
        com = positions.mean(axis=0)
        shift = np.zeros(3)
        for d in range(3):
            while com[d] + shift[d] < box[d, 0]:
                shift[d] += prd[d]
            while com[d] + shift[d] >= box[d, 1]:
                shift[d] -= prd[d]

        if np.any(shift != 0):
            for i in indices:
                atoms[i]["x"] += shift

    print(f"Fixed {n_fixed} periodic jumps within molecules")

    # Write output
    out_lines = lines.copy()
    for a in atoms:
        parts = a["parts"].copy()
        parts[4] = f"{a['x'][0]:.6f}"
        parts[5] = f"{a['x'][1]:.6f}"
        parts[6] = f"{a['x'][2]:.6f}"
        out_lines[a["line_idx"]] = " ".join(parts)

    output_path.write_text("\n".join(out_lines) + "\n")

    # Verify
    n_outside = 0
    max_cg_at_dist = 0.0
    for _mol_id, indices in by_mol.items():
        indices.sort(key=lambda i: atoms[i]["id"])
        cg = [i for i in indices if atoms[i]["type"] in cg_types]
        at = [i for i in indices if atoms[i]["type"] not in cg_types]
        if not cg or not at:
            continue
        apb = len(at) // len(cg)
        for ci, cg_i in enumerate(cg):
            for ai in range(ci * apb, (ci + 1) * apb):
                d = np.linalg.norm(atoms[at[ai]]["x"] - atoms[cg_i]["x"])
                max_cg_at_dist = max(max_cg_at_dist, d)

        for i in indices:
            pos = atoms[i]["x"]
            for d in range(3):
                if pos[d] < box[d, 0] or pos[d] >= box[d, 1]:
                    n_outside += 1
                    break

    print(f"Max CG-AT distance: {max_cg_at_dist:.1f} A")
    print(f"Atoms outside box: {n_outside}")
    print(f"Written to {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--cg-types", type=int, nargs="+", default=[1, 2])
    args = parser.parse_args()
    output = args.output or args.input
    fix_molecules(args.input, output, set(args.cg_types))


if __name__ == "__main__":
    main()
