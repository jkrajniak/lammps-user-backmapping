# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Extract pure CG atoms/bonds from the hybrid data file to create
a standalone CG data file for equilibration."""

from __future__ import annotations

import sys
from pathlib import Path


def parse_hybrid_data(path: Path) -> dict:
    """Parse the hybrid LAMMPS data file into sections."""
    text = path.read_text()
    lines = text.splitlines()

    header: list[str] = []
    sections: dict[str, list[str]] = {}
    current_section: str | None = None

    for line in lines:
        stripped = line.strip()
        if stripped in ("Masses", "Atoms # full", "Bonds", "Angles"):
            current_section = stripped.split("#")[0].strip()
            sections[current_section] = []
        elif current_section is not None:
            sections[current_section].append(line)
        else:
            header.append(line)

    return {"header": header, **sections}


def cg_type_map(mass_lines: list[str]) -> dict[int, int]:
    """Hybrid type id -> consecutive CG type id, from the ``(CG)`` marks in the Masses comments."""
    hybrid_ids = sorted(
        int(ln.split()[0]) for ln in mass_lines if ln.strip() and "(CG)" in ln.split("#", 1)[-1]
    )
    if not hybrid_ids:
        raise ValueError("no '(CG)' atom types in the Masses section of the hybrid data file")
    return {old: new for new, old in enumerate(hybrid_ids, 1)}


def chain_ids(n_atoms: int, bonds: list[dict]) -> dict[int, int]:
    """Molecule id per CG atom: connected components of the CG bonds, numbered by first atom.

    In the hybrid file the molecule id of a bead is the bead's own index, so it cannot be reused.
    """
    parent = list(range(n_atoms + 1))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for b in bonds:
        parent[find(b["a1"])] = find(b["a2"])
    numbering: dict[int, int] = {}
    return {i: numbering.setdefault(find(i), len(numbering) + 1) for i in range(1, n_atoms + 1)}


def extract_cg(hybrid_path: Path, output_path: Path) -> None:
    data = parse_hybrid_data(hybrid_path)
    type_map = cg_type_map(data["Masses"])
    cg_types = set(type_map)

    atoms_raw = [
        ln.split() for ln in data["Atoms"] if ln.strip() and not ln.strip().startswith("#")
    ]
    atoms = [
        {
            "id": int(t[0]),
            "mol": int(t[1]),
            "type": int(t[2]),
            "charge": float(t[3]),
            "x": float(t[4]),
            "y": float(t[5]),
            "z": float(t[6]),
        }
        for t in atoms_raw
    ]

    cg_atoms = [a for a in atoms if a["type"] in cg_types]
    cg_ids = {a["id"] for a in cg_atoms}

    old_to_new: dict[int, int] = {}
    for i, a in enumerate(cg_atoms, 1):
        old_to_new[a["id"]] = i

    bonds_raw = [
        ln.split() for ln in data["Bonds"] if ln.strip() and not ln.strip().startswith("#")
    ]
    cg_bonds = []
    for b in bonds_raw:
        a1, a2 = int(b[2]), int(b[3])
        if a1 in cg_ids and a2 in cg_ids:
            cg_bonds.append(
                {"id": len(cg_bonds) + 1, "type": 1, "a1": old_to_new[a1], "a2": old_to_new[a2]}
            )

    mol_of = chain_ids(len(cg_atoms), cg_bonds)

    mass_lines = [ln for ln in data["Masses"] if ln.strip()]

    box_lines = [ln for ln in data["header"] if "xlo" in ln or "ylo" in ln or "zlo" in ln]

    with open(output_path, "w") as f:
        f.write("LAMMPS data file — pure CG for equilibration\n\n")
        f.write(f"{len(cg_atoms)} atoms\n")
        f.write(f"{len(cg_bonds)} bonds\n")
        f.write("0 angles\n0 dihedrals\n0 impropers\n\n")
        f.write(f"{len(type_map)} atom types\n1 bond types\n0 angle types\n")
        f.write("0 dihedral types\n0 improper types\n\n")
        for bl in box_lines:
            f.write(bl + "\n")
        f.write("\nMasses\n\n")
        for ml in mass_lines:
            parts = ml.split()
            if len(parts) >= 2 and int(parts[0]) in type_map:
                f.write(" ".join([str(type_map[int(parts[0])]), *parts[1:]]) + "\n")
        f.write("\nAtoms # full\n\n")
        for a in cg_atoms:
            new_id = old_to_new[a["id"]]
            f.write(
                f"{new_id} {mol_of[new_id]} {type_map[a['type']]} {a['charge']:.6f} "
                f"{a['x']:.6f} {a['y']:.6f} {a['z']:.6f}\n"
            )
        f.write("\nBonds\n\n")
        for b in cg_bonds:
            f.write(f"{b['id']} {b['type']} {b['a1']} {b['a2']}\n")
        f.write("\n")

    print(f"Wrote {len(cg_atoms)} CG atoms, {len(cg_bonds)} CG bonds → {output_path}")


if __name__ == "__main__":
    hybrid = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("dodecane.data")
    output = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("dodecane_cg.data")
    extract_cg(hybrid, output)
