# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy"]
# ///
"""Generate single-chain dodecane data file with per-bead molecule IDs.

Each CG bead and its 2 AT atoms share a unique molecule ID so that
fix backmap can correctly map CG→AT per bead.

CG bead mapping:
  mol 1: CG 1 (A1,type1) → AT 7(CH3), 8(CH2)
  mol 2: CG 2 (B1,type2) → AT 9(CH2), 10(CH2)
  mol 3: CG 3 (B2,type2) → AT 11(CH2),12(CH2)
  mol 4: CG 4 (B3,type2) → AT 13(CH2),14(CH2)
  mol 5: CG 5 (B4,type2) → AT 15(CH2),16(CH2)
  mol 6: CG 6 (A2,type1) → AT 17(CH2),18(CH3)
"""

from __future__ import annotations
import numpy as np

MASS = {1: 29.062, 2: 28.054, 3: 15.035, 4: 14.027}

CG_BEADS = [
    # (atom_id, type, x, y, z)
    (1, 1, 36.71, 17.95, 12.06),
    (2, 2, 37.18, 16.62, 9.96),
    (3, 2, 38.71, 16.81, 8.24),
    (4, 2, 41.06, 16.40, 7.69),
    (5, 2, 43.05, 17.53, 7.30),
    (6, 1, 45.24, 18.06, 7.78),
]

# Reference AT positions (from dodecane_single.gro, converted to Angstrom)
AT_REF = [
    # (atom_id, type, x, y, z)
    (7, 3, 36.71, 17.95, 12.06),
    (8, 4, 36.32, 16.52, 11.58),
    (9, 4, 35.89, 15.06, 8.24),
    (10, 4, 36.11, 13.62, 7.67),
    (11, 4, 36.83, 13.71, 4.63),
    (12, 4, 35.37, 14.00, 4.86),
    (13, 4, 36.84, 13.75, 3.06),
    (14, 4, 37.26, 14.86, 2.16),
    (15, 4, 38.05, 16.62, 1.12),
    (16, 4, 38.34, 17.86, 0.29),
    (17, 4, 39.65, 19.66, 0.89),
    (18, 3, 40.24, 20.83, 0.02),
]

# CG bead → AT atom indices (0-based into AT_REF)
BEAD_MAP = {
    0: [0, 1],  # CG1 → AT7, AT8
    1: [2, 3],  # CG2 → AT9, AT10
    2: [4, 5],  # CG3 → AT11,AT12
    3: [6, 7],  # CG4 → AT13,AT14
    4: [8, 9],  # CG5 → AT15,AT16
    5: [10, 11],  # CG6 → AT17,AT18
}


def main():
    # Shift each AT pair so its mass-weighted COM equals CG bead position
    at_atoms = list(AT_REF)
    for cg_idx, at_indices in BEAD_MAP.items():
        cg = CG_BEADS[cg_idx]
        cg_pos = np.array(cg[2:5])

        masses = np.array([MASS[at_atoms[i][1]] for i in at_indices])
        positions = np.array([list(at_atoms[i][2:5]) for i in at_indices])
        com = np.average(positions, axis=0, weights=masses)
        shift = cg_pos - com

        for i in at_indices:
            aid, atype, x, y, z = at_atoms[i]
            at_atoms[i] = (aid, atype, x + shift[0], y + shift[1], z + shift[2])

    # Build atoms list with per-bead molecule IDs
    atoms = []
    for cg_idx, cg in enumerate(CG_BEADS):
        mol_id = cg_idx + 1
        aid, atype, x, y, z = cg
        atoms.append((aid, mol_id, atype, 0.0, x, y, z))

    for cg_idx, at_indices in BEAD_MAP.items():
        mol_id = cg_idx + 1
        for i in at_indices:
            aid, atype, x, y, z = at_atoms[i]
            atoms.append((aid, mol_id, atype, 0.0, x, y, z))

    atoms.sort(key=lambda a: a[0])

    lines = []
    lines.append(
        "LAMMPS data file — single dodecane chain for backmapping verification"
    )
    lines.append("")
    lines.append("18 atoms")
    lines.append("16 bonds")
    lines.append("10 angles")
    lines.append("0 dihedrals")
    lines.append("0 impropers")
    lines.append("")
    lines.append("4 atom types")
    lines.append("4 bond types")
    lines.append("1 angle types")
    lines.append("0 dihedral types")
    lines.append("0 improper types")
    lines.append("")
    lines.append("0.0 57.537600 xlo xhi")
    lines.append("0.0 57.537600 ylo yhi")
    lines.append("0.0 57.537600 zlo zhi")
    lines.append("")
    lines.append("Masses")
    lines.append("")
    lines.append("1 29.062000 # A (CG)")
    lines.append("2 28.054000 # B (CG)")
    lines.append("3 15.035000 # CH3")
    lines.append("4 14.027000 # CH2")
    lines.append("")
    lines.append("Atoms # full")
    lines.append("")
    for a in atoms:
        lines.append(
            f"{a[0]} {a[1]} {a[2]} {a[3]:.6f} {a[4]:.6f} {a[5]:.6f} {a[6]:.6f}"
        )

    lines.append("")
    lines.append("Bonds")
    lines.append("")
    bonds = [
        # Intra-bead AT bonds (type 1 = harmonic)
        (1, 1, 7, 8),
        (2, 1, 9, 10),
        (3, 1, 11, 12),
        (4, 1, 13, 14),
        (5, 1, 15, 16),
        (6, 1, 17, 18),
        # CG-CG tabulated bonds (type 2, 3 = backmap/table cg)
        (7, 2, 1, 2),
        (8, 2, 2, 3),
        (9, 2, 3, 4),
        (10, 2, 4, 5),
        (11, 3, 5, 6),
        # AT-AT cross-bead bonds (type 4 = backmap/harmonic at)
        (12, 4, 8, 9),
        (13, 4, 10, 11),
        (14, 4, 12, 13),
        (15, 4, 14, 15),
        (16, 4, 16, 17),
    ]
    for b in bonds:
        lines.append(f"{b[0]} {b[1]} {b[2]} {b[3]}")

    lines.append("")
    lines.append("Angles")
    lines.append("")
    angles = [
        (1, 1, 7, 8, 9),
        (2, 1, 8, 9, 10),
        (3, 1, 9, 10, 11),
        (4, 1, 10, 11, 12),
        (5, 1, 11, 12, 13),
        (6, 1, 12, 13, 14),
        (7, 1, 13, 14, 15),
        (8, 1, 14, 15, 16),
        (9, 1, 15, 16, 17),
        (10, 1, 16, 17, 18),
    ]
    for a in angles:
        lines.append(f"{a[0]} {a[1]} {a[2]} {a[3]} {a[4]}")

    lines.append("")

    with open("dodecane_single_chain.data", "w") as f:
        f.write("\n".join(lines))

    print("Generated dodecane_single_chain.data with per-bead molecule IDs:")
    for cg_idx, cg in enumerate(CG_BEADS):
        at_ids = [AT_REF[i][0] for i in BEAD_MAP[cg_idx]]
        print(f"  mol {cg_idx + 1}: CG {cg[0]} (type {cg[1]}) + AT {at_ids}")


if __name__ == "__main__":
    main()
