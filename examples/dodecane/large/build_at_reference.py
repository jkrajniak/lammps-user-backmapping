# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy"]
# ///
"""Build an independent atomistic dodecane system for reference RDF.

Reads the single-molecule GROMACS coordinates, replicates 10 copies on a
lattice with random orientations, and writes a LAMMPS data file with the
GROMOS united-atom force field. This system has NO connection to the
backmapping procedure — it serves as the independent AT reference.

Usage:
  uv run build_at_reference.py [--gro dodecane_single.gro] [--output dodecane_at_ref.data] [--seed 42]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

# GROMOS united-atom parameters (from topol_aa.top)
ATOM_TYPES = {
    "CH3": {"type_id": 1, "mass": 15.035, "sigma": 3.748, "epsilon": 0.207266},
    "CH2": {"type_id": 2, "mass": 14.027, "sigma": 3.905, "epsilon": 0.117997},
}

BOND_K = 800.000883  # kcal/(mol·Å²) — converted from 334720 kJ/(mol·nm²)
BOND_R0 = 1.530  # Å

ANGLE_K = 126.673180  # kcal/(mol·rad²) — converted from 530 kJ/(mol·rad²)
ANGLE_THETA0 = 111.0  # degrees

# dodecane type sequence: CH3, CH2*10, CH3
TYPE_SEQUENCE = ["CH3"] + ["CH2"] * 10 + ["CH3"]

BOX_SIZE = 57.5376  # Å (same as CG system: 5.75376 nm)
N_MOLECULES = 500


def parse_gro(path: Path) -> np.ndarray:
    """Read atom positions from a GRO file. Returns coords in Å."""
    lines = path.read_text().splitlines()
    n_atoms = int(lines[1].strip())
    coords = []
    for i in range(2, 2 + n_atoms):
        x = float(lines[i][20:28]) * 10.0  # nm → Å
        y = float(lines[i][28:36]) * 10.0
        z = float(lines[i][36:44]) * 10.0
        coords.append([x, y, z])
    return np.array(coords)


def random_rotation_matrix(rng: np.random.Generator) -> np.ndarray:
    """Generate a uniformly distributed random 3D rotation matrix."""
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
    """Place n_mol copies of the template molecule on a lattice with random rotations.

    Uses a larger initial box to avoid overlaps. Returns (coords, actual_box_size).
    """
    com = template.mean(axis=0)
    centered = template - com

    n_per_side = int(np.ceil(n_mol ** (1.0 / 3.0)))
    # Minimum spacing: molecule length (~15 Å) + buffer
    mol_extent = np.max(np.linalg.norm(centered, axis=1)) * 2
    min_spacing = mol_extent + 4.0  # 4 Å buffer
    init_box = max(box, n_per_side * min_spacing)

    dx = init_box / n_per_side
    dy = init_box / n_per_side
    dz = init_box / n_per_side

    positions = []
    for ix in range(n_per_side):
        for iy in range(n_per_side):
            for iz in range(n_per_side):
                if len(positions) >= n_mol:
                    break
                center = np.array(
                    [
                        (ix + 0.5) * dx,
                        (iy + 0.5) * dy,
                        (iz + 0.5) * dz,
                    ]
                )
                rot = random_rotation_matrix(rng)
                mol_coords = centered @ rot.T + center
                positions.append(mol_coords)

    print(f"  Grid: {n_per_side}³ = {n_per_side**3} slots, using {n_mol}")
    print(f"  Initial box: {init_box:.1f} Å (target: {box:.1f} Å)")
    print("  Initial density: will need NPT compression to reach liquid density")

    return np.array(positions), init_box


def write_lammps_data(
    path: Path,
    all_coords: np.ndarray,
    box: float,
) -> None:
    """Write LAMMPS data file for the pure AT system."""
    n_mol, n_atoms_per_mol, _ = all_coords.shape
    n_atoms = n_mol * n_atoms_per_mol
    n_bonds_per_mol = n_atoms_per_mol - 1
    n_angles_per_mol = n_atoms_per_mol - 2
    n_bonds = n_mol * n_bonds_per_mol
    n_angles = n_mol * n_angles_per_mol

    with path.open("w") as f:
        f.write(f"LAMMPS data file — independent AT dodecane reference ({n_mol} molecules)\n\n")
        f.write(f"{n_atoms} atoms\n")
        f.write(f"{n_bonds} bonds\n")
        f.write(f"{n_angles} angles\n")
        f.write("0 dihedrals\n")
        f.write("0 impropers\n\n")
        f.write("2 atom types\n")
        f.write("1 bond types\n")
        f.write("1 angle types\n")
        f.write("0 dihedral types\n")
        f.write("0 improper types\n\n")
        f.write(f"0.0 {box:.6f} xlo xhi\n")
        f.write(f"0.0 {box:.6f} ylo yhi\n")
        f.write(f"0.0 {box:.6f} zlo zhi\n\n")

        f.write("Masses\n\n")
        f.write(f"1 {ATOM_TYPES['CH3']['mass']:.6f} # CH3\n")
        f.write(f"2 {ATOM_TYPES['CH2']['mass']:.6f} # CH2\n\n")

        f.write("Pair Coeffs # lj/cut\n\n")
        f.write(f"1 {ATOM_TYPES['CH3']['epsilon']:.6f} {ATOM_TYPES['CH3']['sigma']:.6f} # CH3\n")
        f.write(f"2 {ATOM_TYPES['CH2']['epsilon']:.6f} {ATOM_TYPES['CH2']['sigma']:.6f} # CH2\n\n")

        f.write("Bond Coeffs # harmonic\n\n")
        f.write(f"1 {BOND_K:.6f} {BOND_R0:.6f}\n\n")

        f.write("Angle Coeffs # harmonic\n\n")
        f.write(f"1 {ANGLE_K:.6f} {ANGLE_THETA0:.4f}\n\n")

        f.write("Atoms # full\n\n")
        atom_id = 0
        for mol_idx in range(n_mol):
            for at_idx in range(n_atoms_per_mol):
                atom_id += 1
                mol_id = mol_idx + 1
                type_name = TYPE_SEQUENCE[at_idx]
                type_id = ATOM_TYPES[type_name]["type_id"]
                x, y, z = all_coords[mol_idx, at_idx]
                # Wrap into box
                x = x % box
                y = y % box
                z = z % box
                f.write(f"{atom_id} {mol_id} {type_id} 0.000000 {x:.6f} {y:.6f} {z:.6f}\n")
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

    print(f"Written {n_atoms} atoms, {n_bonds} bonds, {n_angles} angles to {path}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build independent AT dodecane system for reference simulation"
    )
    parser.add_argument("--gro", type=Path, default=Path("dodecane_single.gro"))
    parser.add_argument("--output", type=Path, default=Path("dodecane_at_ref.data"))
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n-mol", type=int, default=N_MOLECULES)
    args = parser.parse_args()

    if not args.gro.exists():
        print(f"ERROR: {args.gro} not found", file=sys.stderr)
        return 1

    rng = np.random.default_rng(args.seed)
    template = parse_gro(args.gro)

    print(f"Template molecule: {len(template)} atoms from {args.gro}")
    print(f"Placing {args.n_mol} copies (target box: {BOX_SIZE:.2f} Å)")

    all_coords, init_box = build_system(template, args.n_mol, BOX_SIZE, rng)
    write_lammps_data(args.output, all_coords, init_box)

    return 0


if __name__ == "__main__":
    sys.exit(main())
