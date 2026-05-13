# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy"]
# ///
"""Place atomistic (AT) fragments around equilibrated CG beads.

Takes the equilibrated CG data file and inserts the appropriate AT united-atom
fragments around each CG bead:
  - Terminal A bead (type 1): CH3 (type 3) + CH2 (type 4)
  - Middle B bead  (type 2): CH2 (type 4) + CH2 (type 4)

The two AT atoms per bead are placed along the local chain tangent so that
their center of mass coincides with the CG bead position and they are
separated by the AT bond equilibrium distance (1.53 Å).
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

AT_BOND_EQ = 1.53  # Å — AT C-C bond equilibrium distance

MASS = {1: 29.062, 2: 28.054, 3: 15.035, 4: 14.027}

# AT fragment per bead type and position in chain
# (left_type, right_type) where left is towards chain start, right towards end
FRAG_FIRST = (3, 4)  # terminal A at start: CH3(outer) - CH2(inner)
FRAG_MIDDLE = (4, 4)  # middle B:            CH2 - CH2
FRAG_LAST = (4, 3)  # terminal A at end:   CH2(inner) - CH3(outer)


def parse_cg_data(path: Path) -> tuple[dict, list[dict], list[tuple[int, int]], np.ndarray]:
    """Parse the equilibrated CG data file.

    Returns (header_info, atoms, bonds, box) where box = [Lx, Ly, Lz].
    """
    text = path.read_text()
    lines = text.splitlines()

    box = np.zeros(3)
    atoms: list[dict] = []
    bonds: list[tuple[int, int]] = []
    section = None

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        if stripped.startswith("Atoms"):
            section = "atoms"
            continue
        if stripped.startswith("Bonds"):
            section = "bonds"
            continue
        if stripped.startswith("Velocities"):
            section = "vel"
            continue

        if "xlo" in stripped:
            parts = stripped.split()
            box[0] = float(parts[1]) - float(parts[0])
        elif "ylo" in stripped:
            parts = stripped.split()
            box[1] = float(parts[1]) - float(parts[0])
        elif "zlo" in stripped:
            parts = stripped.split()
            box[2] = float(parts[1]) - float(parts[0])
        elif section == "atoms":
            parts = stripped.split()
            if len(parts) >= 7:
                atoms.append(
                    {
                        "id": int(parts[0]),
                        "mol": int(parts[1]),
                        "type": int(parts[2]),
                        "charge": float(parts[3]),
                        "x": float(parts[4]),
                        "y": float(parts[5]),
                        "z": float(parts[6]),
                        "ix": int(parts[7]) if len(parts) > 7 else 0,
                        "iy": int(parts[8]) if len(parts) > 8 else 0,
                        "iz": int(parts[9]) if len(parts) > 9 else 0,
                    }
                )
        elif section == "bonds":
            parts = stripped.split()
            if len(parts) >= 4:
                bonds.append((int(parts[2]), int(parts[3])))
        elif section == "vel":
            pass

    return {}, atoms, bonds, box


def unwrap(pos: np.ndarray, img: np.ndarray, box: np.ndarray) -> np.ndarray:
    return pos + img * box


def wrap(pos: np.ndarray, box: np.ndarray) -> np.ndarray:
    return pos % box


def trace_chain(mol_atoms: list[dict], mol_bonds: list[tuple[int, int]]) -> list[dict]:
    """Order CG beads along the chain from terminal to terminal."""
    id_to_atom = {a["id"]: a for a in mol_atoms}
    adj: dict[int, list[int]] = defaultdict(list)
    for a1, a2 in mol_bonds:
        adj[a1].append(a2)
        adj[a2].append(a1)

    terminals = [aid for aid, nbrs in adj.items() if len(nbrs) == 1]
    if len(terminals) != 2:
        raise ValueError(f"Expected 2 terminal beads, got {len(terminals)}")

    start = min(terminals)
    chain = [start]
    prev = -1
    current = start
    while True:
        nbrs = [n for n in adj[current] if n != prev]
        if not nbrs:
            break
        prev = current
        current = nbrs[0]
        chain.append(current)

    return [id_to_atom[aid] for aid in chain]


def place_fragments(
    chain: list[dict], box: np.ndarray
) -> tuple[list[dict], list[tuple[int, int, int]], list[tuple[int, int, int, int]]]:
    """Place AT atoms along the chain.

    Returns (at_atoms, intra_bonds, inter_bonds) where:
      at_atoms: list of AT atom dicts (id, mol, type, x, y, z)
      intra_bonds: (bond_type, at1_id, at2_id) within each bead
      inter_bonds: (bond_type, at1_id, at2_id) between adjacent beads
    """
    n = len(chain)
    mol_id = chain[0]["mol"]

    # Unwrap chain positions locally: start from bead 0's wrapped position
    # and use minimum-image displacements to reconstruct a contiguous chain.
    raw = np.array([[a["x"], a["y"], a["z"]] for a in chain])
    positions = np.zeros_like(raw)
    positions[0] = raw[0]  # anchor: use wrapped position
    for i in range(1, n):
        delta = raw[i] - raw[i - 1]
        delta -= np.round(delta / box) * box  # minimum image
        positions[i] = positions[i - 1] + delta

    at_atoms: list[dict] = []
    intra_bonds: list[tuple[int, int, int]] = []
    inter_bonds: list[tuple[int, int, int, int]] = []

    for i in range(n):
        if i == 0:
            tangent = positions[1] - positions[0]
        elif i == n - 1:
            tangent = positions[-1] - positions[-2]
        else:
            tangent = positions[i + 1] - positions[i - 1]

        norm = np.linalg.norm(tangent)
        tangent = np.array([1.0, 0.0, 0.0]) if norm < 1e-10 else tangent / norm

        if i == 0:
            frag = FRAG_FIRST
        elif i == n - 1:
            frag = FRAG_LAST
        else:
            frag = FRAG_MIDDLE

        m_left = MASS[frag[0]]
        m_right = MASS[frag[1]]
        m_total = m_left + m_right

        offset_left = m_right / m_total * AT_BOND_EQ
        offset_right = m_left / m_total * AT_BOND_EQ

        pos_left = positions[i] - offset_left * tangent
        pos_right = positions[i] + offset_right * tangent

        # Store unwrapped positions; we wrap the whole molecule later
        at_atoms.append(
            {
                "mol": mol_id,
                "type": frag[0],
                "x": pos_left[0],
                "y": pos_left[1],
                "z": pos_left[2],
                "bead_idx": i,
                "side": "left",
            }
        )
        at_atoms.append(
            {
                "mol": mol_id,
                "type": frag[1],
                "x": pos_right[0],
                "y": pos_right[1],
                "z": pos_right[2],
                "bead_idx": i,
                "side": "right",
            }
        )

    return at_atoms, intra_bonds, inter_bonds


def build_hybrid(cg_path: Path, output_path: Path) -> None:
    _, atoms, bonds, box = parse_cg_data(cg_path)

    # Group atoms by molecule
    mol_atoms: dict[int, list[dict]] = defaultdict(list)
    for a in atoms:
        mol_atoms[a["mol"]].append(a)

    mol_bonds: dict[int, list[tuple[int, int]]] = defaultdict(list)
    id_to_mol = {a["id"]: a["mol"] for a in atoms}
    for a1, a2 in bonds:
        mol_id = id_to_mol[a1]
        mol_bonds[mol_id].append((a1, a2))

    # Sort molecules by ID
    mol_ids = sorted(mol_atoms.keys())
    n_cg = len(atoms)
    n_beads_per_mol = len(mol_atoms[mol_ids[0]])
    n_at_per_mol = n_beads_per_mol * 2
    n_at = n_cg * 2
    n_cg + n_at

    # Trace chains and place AT atoms
    all_cg_ordered: list[dict] = []
    all_at: list[dict] = []

    for mol_id in mol_ids:
        chain = trace_chain(mol_atoms[mol_id], mol_bonds[mol_id])
        at_atoms, _, _ = place_fragments(chain, box)

        # Renumber CG atoms sequentially within molecule
        for bead in chain:
            all_cg_ordered.append(bead)
        for at in at_atoms:
            all_at.append(at)

    # Assign final IDs: CG atoms first per molecule, then AT atoms
    # Ordering: mol1_cg1..cg6, mol1_at1..at12, mol2_cg1..cg6, mol2_at1..at12, ...
    final_atoms: list[dict] = []
    atom_id = 1

    cg_per_mol = n_beads_per_mol
    at_per_mol = n_at_per_mol

    for mol_idx, _mol_id in enumerate(mol_ids):
        cg_start = mol_idx * cg_per_mol
        at_start = mol_idx * at_per_mol

        # Collect all unwrapped positions for this molecule (CG + AT)
        mol_positions: list[np.ndarray] = []

        for j in range(cg_per_mol):
            a = all_cg_ordered[cg_start + j]
            pos = unwrap(
                np.array([a["x"], a["y"], a["z"]]),
                np.array([a["ix"], a["iy"], a["iz"]]),
                box,
            )
            mol_positions.append(pos)

        for j in range(at_per_mol):
            a = all_at[at_start + j]
            mol_positions.append(np.array([a["x"], a["y"], a["z"]]))

        # Wrap the whole molecule: shift so first CG bead is inside the box,
        # then apply the same shift to all atoms (keeps bonds intact).
        anchor = mol_positions[0].copy()
        anchor_wrapped = wrap(anchor, box)
        shift = anchor_wrapped - anchor
        for k in range(len(mol_positions)):
            mol_positions[k] = mol_positions[k] + shift

        # Wrap each atom into [0, L) with correct image flags
        def wrap_with_image(pos: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
            img = np.floor(pos / box).astype(int)
            wrapped = pos - img * box
            return wrapped, img

        # Write CG atoms
        cg_types_list = [all_cg_ordered[cg_start + j]["type"] for j in range(cg_per_mol)]
        for j in range(cg_per_mol):
            w, img = wrap_with_image(mol_positions[j])
            final_atoms.append(
                {
                    "id": atom_id,
                    "mol": mol_idx + 1,
                    "type": cg_types_list[j],
                    "charge": 0.0,
                    "x": w[0],
                    "y": w[1],
                    "z": w[2],
                    "ix": img[0],
                    "iy": img[1],
                    "iz": img[2],
                }
            )
            atom_id += 1

        # Write AT atoms
        at_types_list = [all_at[at_start + j]["type"] for j in range(at_per_mol)]
        for j in range(at_per_mol):
            w, img = wrap_with_image(mol_positions[cg_per_mol + j])
            final_atoms.append(
                {
                    "id": atom_id,
                    "mol": mol_idx + 1,
                    "type": at_types_list[j],
                    "charge": 0.0,
                    "x": w[0],
                    "y": w[1],
                    "z": w[2],
                    "ix": img[0],
                    "iy": img[1],
                    "iz": img[2],
                }
            )
            atom_id += 1

    # Build bonds
    all_bonds: list[dict] = []
    bond_id = 1

    for mol_idx in range(len(mol_ids)):
        base_cg = mol_idx * (cg_per_mol + at_per_mol) + 1
        base_at = base_cg + cg_per_mol

        # CG bonds (between consecutive CG beads)
        for j in range(cg_per_mol - 1):
            btype = 3 if j == cg_per_mol - 2 else 2  # last CG bond is type 3
            all_bonds.append(
                {"id": bond_id, "type": btype, "a1": base_cg + j, "a2": base_cg + j + 1}
            )
            bond_id += 1

        # AT intra-bead bonds (type 1, harmonic between the 2 AT atoms of each bead)
        for j in range(cg_per_mol):
            at1 = base_at + j * 2
            at2 = base_at + j * 2 + 1
            all_bonds.append({"id": bond_id, "type": 1, "a1": at1, "a2": at2})
            bond_id += 1

        # AT inter-bead bonds (type 4, between right AT of bead j and left AT of bead j+1)
        for j in range(cg_per_mol - 1):
            at_right = base_at + j * 2 + 1  # right atom of bead j
            at_left = base_at + (j + 1) * 2  # left atom of bead j+1
            all_bonds.append({"id": bond_id, "type": 4, "a1": at_right, "a2": at_left})
            bond_id += 1

    # Build angles (AT-AT-AT along the chain)
    all_angles: list[dict] = []
    angle_id = 1

    for mol_idx in range(len(mol_ids)):
        base_at = mol_idx * (cg_per_mol + at_per_mol) + 1 + cg_per_mol

        # Consecutive AT atoms form angles
        for j in range(at_per_mol - 2):
            at1 = base_at + j
            at2 = base_at + j + 1
            at3 = base_at + j + 2
            all_angles.append({"id": angle_id, "type": 1, "a1": at1, "a2": at2, "a3": at3})
            angle_id += 1

    # Write hybrid data file
    with open(output_path, "w") as f:
        f.write("LAMMPS data file for backmapping — hybrid CG+AT\n\n")
        f.write(f"{len(final_atoms)} atoms\n")
        f.write(f"{len(all_bonds)} bonds\n")
        f.write(f"{len(all_angles)} angles\n")
        f.write("0 dihedrals\n0 impropers\n\n")
        f.write("4 atom types\n4 bond types\n1 angle types\n")
        f.write("0 dihedral types\n0 improper types\n\n")
        f.write(f"0.0 {box[0]:.6f} xlo xhi\n")
        f.write(f"0.0 {box[1]:.6f} ylo yhi\n")
        f.write(f"0.0 {box[2]:.6f} zlo zhi\n")
        f.write("\nMasses\n\n")
        f.write("1 29.062000 # A (CG)\n")
        f.write("2 28.054000 # B (CG)\n")
        f.write("3 15.035000 # CH3\n")
        f.write("4 14.027000 # CH2\n")
        f.write("\nAtoms # full\n\n")
        for a in final_atoms:
            f.write(
                f"{a['id']} {a['mol']} {a['type']} {a['charge']:.6f} "
                f"{a['x']:.6f} {a['y']:.6f} {a['z']:.6f} "
                f"{a['ix']} {a['iy']} {a['iz']}\n"
            )
        f.write("\nBonds\n\n")
        for b in all_bonds:
            f.write(f"{b['id']} {b['type']} {b['a1']} {b['a2']}\n")
        f.write("\nAngles\n\n")
        for ang in all_angles:
            f.write(f"{ang['id']} {ang['type']} {ang['a1']} {ang['a2']} {ang['a3']}\n")
        f.write("\n")

    print(
        f"Wrote {len(final_atoms)} atoms ({n_cg} CG + {n_at} AT), "
        f"{len(all_bonds)} bonds, {len(all_angles)} angles → {output_path}"
    )


if __name__ == "__main__":
    cg_equil = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("dodecane_cg_equil.data")
    output = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("dodecane_hybrid.data")
    build_hybrid(cg_equil, output)
