# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy"]
# ///
"""Build independent PE melt (10 chains) for reference RDF at 423 K.

Usage:
  uv run build_at_reference.py [--gro pe_single.gro] [--output pe_at_ref.data] [--n-mol 10]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ATOM_TYPES = {
    "CH3": {"type_id": 1, "mass": 15.035, "sigma": 3.748, "epsilon": 0.207266},
    "CH2": {"type_id": 2, "mass": 14.027, "sigma": 3.905, "epsilon": 0.117997},
}

BOND_K = 633.365900
BOND_R0 = 1.530
ANGLE_K = 124.283120
ANGLE_THETA0 = 111.0
DIHEDRAL_RB = (0.717018, -2.217713, 2.905285, 3.135783, -0.731287, -6.271589)

TYPE_SEQUENCE = ["CH3"] + ["CH2"] * 98 + ["CH3"]
BOX_SIZE = 60.1855
N_MOLECULES = 10


def parse_gro(path: Path) -> np.ndarray:
    lines = path.read_text().splitlines()
    n_atoms = int(lines[1].strip())
    coords = []
    for i in range(2, 2 + n_atoms):
        x = float(lines[i][20:28]) * 10.0
        y = float(lines[i][28:36]) * 10.0
        z = float(lines[i][36:44]) * 10.0
        coords.append([x, y, z])
    return np.array(coords)


def random_rotation_matrix(rng: np.random.Generator) -> np.ndarray:
    q = rng.standard_normal(4)
    q /= np.linalg.norm(q)
    w, x, y, z = q
    return np.array(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
            [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
            [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)],
        ]
    )


def build_system(
    template: np.ndarray,
    n_mol: int,
    box: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, float]:
    com = template.mean(axis=0)
    centered = template - com
    extent = np.max(np.ptp(centered, axis=0))
    spacing = extent * 2.5
    n_side = int(np.ceil(n_mol ** (1 / 3)))
    init_box = max(box * 1.5, n_side * spacing * 1.5)

    all_coords = np.zeros((n_mol, len(template), 3))
    idx = 0
    for ix in range(n_side):
        for iy in range(n_side):
            for iz in range(n_side):
                if idx >= n_mol:
                    break
                rot = random_rotation_matrix(rng)
                offset = np.array([ix, iy, iz], dtype=float) * spacing
                all_coords[idx] = (centered @ rot.T) + offset + init_box / 2
                idx += 1
    return all_coords, init_box


def write_lammps_data(path: Path, all_coords: np.ndarray, box: float) -> None:
    n_mol, n_atoms_per_mol, _ = all_coords.shape
    n_atoms = n_mol * n_atoms_per_mol
    n_bonds_per_mol = n_atoms_per_mol - 1
    n_angles_per_mol = n_atoms_per_mol - 2
    n_dihedrals_per_mol = n_atoms_per_mol - 3
    n_bonds = n_mol * n_bonds_per_mol
    n_angles = n_mol * n_angles_per_mol
    n_dihedrals = n_mol * n_dihedrals_per_mol

    with open(path, "w") as f:
        f.write("LAMMPS data file — independent PE AT reference\n\n")
        f.write(f"{n_atoms} atoms\n{n_bonds} bonds\n{n_angles} angles\n{n_dihedrals} dihedrals\n")
        f.write("0 impropers\n\n")
        f.write("2 atom types\n1 bond types\n1 angle types\n1 dihedral types\n0 improper types\n\n")
        f.write(f"0.0 {box:.6f} xlo xhi\n0.0 {box:.6f} ylo yhi\n0.0 {box:.6f} zlo zhi\n\n")

        f.write("Masses\n\n")
        f.write(f"1 {ATOM_TYPES['CH3']['mass']:.6f} # CH3\n")
        f.write(f"2 {ATOM_TYPES['CH2']['mass']:.6f} # CH2\n\n")

        f.write("Pair Coeffs # lj/cut\n\n")
        f.write(f"1 {ATOM_TYPES['CH3']['epsilon']:.6f} {ATOM_TYPES['CH3']['sigma']:.6f}\n")
        f.write(f"2 {ATOM_TYPES['CH2']['epsilon']:.6f} {ATOM_TYPES['CH2']['sigma']:.6f}\n\n")

        f.write("Bond Coeffs # harmonic\n\n")
        f.write(f"1 {BOND_K:.6f} {BOND_R0:.6f}\n\n")

        f.write("Angle Coeffs # harmonic\n\n")
        f.write(f"1 {ANGLE_K:.6f} {ANGLE_THETA0:.4f}\n\n")

        c0, c1, c2, c3, c4, c5 = DIHEDRAL_RB
        f.write("Dihedral Coeffs # ryckaert\n\n")
        f.write(f"1 {c0:.6f} {c1:.6f} {c2:.6f} {c3:.6f} {c4:.6f} {c5:.6f}\n\n")

        f.write("Atoms # full\n\n")
        atom_id = 0
        for mol_idx in range(n_mol):
            for at_idx in range(n_atoms_per_mol):
                atom_id += 1
                type_name = TYPE_SEQUENCE[at_idx]
                type_id = ATOM_TYPES[type_name]["type_id"]
                x, y, z = all_coords[mol_idx, at_idx]
                f.write(f"{atom_id} {mol_idx + 1} {type_id} 0.000000 {x:.6f} {y:.6f} {z:.6f}\n")
        f.write("\n")

        f.write("Bonds\n\n")
        bond_id = 0
        for mol_idx in range(n_mol):
            base = mol_idx * n_atoms_per_mol
            for j in range(n_bonds_per_mol):
                bond_id += 1
                f.write(f"{bond_id} 1 {base + j + 1} {base + j + 2}\n")
        f.write("\n")

        f.write("Angles\n\n")
        angle_id = 0
        for mol_idx in range(n_mol):
            base = mol_idx * n_atoms_per_mol
            for j in range(n_angles_per_mol):
                angle_id += 1
                f.write(f"{angle_id} 1 {base + j + 1} {base + j + 2} {base + j + 3}\n")
        f.write("\n")

        f.write("Dihedrals\n\n")
        dih_id = 0
        for mol_idx in range(n_mol):
            base = mol_idx * n_atoms_per_mol
            for j in range(n_dihedrals_per_mol):
                dih_id += 1
                f.write(f"{dih_id} 1 {base + j + 1} {base + j + 2} {base + j + 3} {base + j + 4}\n")

    print(f"Written {n_atoms} atoms, {n_bonds} bonds, {n_angles} angles, {n_dihedrals} dihedrals")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build independent PE AT reference")
    parser.add_argument("--gro", type=Path, default=Path("pe_single.gro"))
    parser.add_argument("--output", type=Path, default=Path("pe_at_ref.data"))
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n-mol", type=int, default=N_MOLECULES)
    parser.add_argument("--box", type=float, default=BOX_SIZE)
    args = parser.parse_args()

    if not args.gro.exists():
        print(f"ERROR: {args.gro} not found", file=sys.stderr)
        return 1

    rng = np.random.default_rng(args.seed)
    template = parse_gro(args.gro)
    all_coords, init_box = build_system(template, args.n_mol, args.box, rng)
    write_lammps_data(args.output, all_coords, init_box)
    print(f"Initial box: {init_box:.2f} Å (target melt box: {args.box:.2f} Å)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
