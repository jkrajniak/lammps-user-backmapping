# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy"]
# ///
"""Fix AT atom positions: shift each AT group so its COM coincides with its parent CG bead.

The original data file places AT atoms from a reference all-atom structure
that doesn't match the CG bead positions. This script translates each
AT atom pair so their mass-weighted COM equals the CG bead position.
"""

from __future__ import annotations

import numpy as np
from pathlib import Path

MASS = {1: 29.062, 2: 28.054, 3: 15.035, 4: 14.027}

# CG_id → [AT_ids]
CG_TO_AT = {
    1: [7, 8],
    2: [9, 10],
    3: [11, 12],
    4: [13, 14],
    5: [15, 16],
    6: [17, 18],
}


def main():
    data_path = Path("dodecane_single_chain.data")
    text = data_path.read_text()
    lines = text.splitlines()

    # Parse atom positions
    atoms: dict[int, dict] = {}
    in_atoms = False
    for idx, line in enumerate(lines):
        if line.strip() == "Atoms # full":
            in_atoms = True
            continue
        if in_atoms and line.strip() == "":
            if atoms:
                break
            continue
        if in_atoms and line.strip():
            parts = line.split()
            aid = int(parts[0])
            atoms[aid] = {
                "id": aid,
                "mol": int(parts[1]),
                "type": int(parts[2]),
                "charge": float(parts[3]),
                "x": float(parts[4]),
                "y": float(parts[5]),
                "z": float(parts[6]),
                "line_idx": idx,
            }

    print("Before fix — AT COM vs CG position:")
    for cg_id, at_ids in CG_TO_AT.items():
        cg = atoms[cg_id]
        cg_pos = np.array([cg["x"], cg["y"], cg["z"]])

        masses = np.array([MASS[atoms[a]["type"]] for a in at_ids])
        positions = np.array(
            [[atoms[a]["x"], atoms[a]["y"], atoms[a]["z"]] for a in at_ids]
        )
        com = np.average(positions, axis=0, weights=masses)

        displacement = cg_pos - com
        print(
            f"  CG {cg_id} ({cg_pos}) ← AT COM ({com}) | shift = {displacement} | |d| = {np.linalg.norm(displacement):.2f} Å"
        )

        # Apply shift to AT atoms
        for a_id in at_ids:
            atoms[a_id]["x"] += displacement[0]
            atoms[a_id]["y"] += displacement[1]
            atoms[a_id]["z"] += displacement[2]

    print("\nAfter fix — verification:")
    for cg_id, at_ids in CG_TO_AT.items():
        cg = atoms[cg_id]
        cg_pos = np.array([cg["x"], cg["y"], cg["z"]])
        masses = np.array([MASS[atoms[a]["type"]] for a in at_ids])
        positions = np.array(
            [[atoms[a]["x"], atoms[a]["y"], atoms[a]["z"]] for a in at_ids]
        )
        com = np.average(positions, axis=0, weights=masses)
        err = np.linalg.norm(cg_pos - com)
        print(f"  CG {cg_id}: |CG - AT_COM| = {err:.6f} Å")

    # Rewrite atom lines
    for aid, a in atoms.items():
        line_idx = a["line_idx"]
        lines[line_idx] = (
            f"{a['id']} {a['mol']} {a['type']} {a['charge']:.6f} {a['x']:.6f} {a['y']:.6f} {a['z']:.6f}"
        )

    data_path.write_text("\n".join(lines) + "\n")
    print(f"\nFixed data written to {data_path}")


if __name__ == "__main__":
    main()
