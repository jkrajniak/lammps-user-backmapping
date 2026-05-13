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


def extract_cg(hybrid_path: Path, output_path: Path) -> None:
    cg_types = {1, 2}

    data = parse_hybrid_data(hybrid_path)

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

    mass_lines = [ln for ln in data["Masses"] if ln.strip()]

    box_lines = [ln for ln in data["header"] if "xlo" in ln or "ylo" in ln or "zlo" in ln]

    with open(output_path, "w") as f:
        f.write("LAMMPS data file — pure CG for equilibration\n\n")
        f.write(f"{len(cg_atoms)} atoms\n")
        f.write(f"{len(cg_bonds)} bonds\n")
        f.write("0 angles\n0 dihedrals\n0 impropers\n\n")
        f.write("2 atom types\n1 bond types\n0 angle types\n")
        f.write("0 dihedral types\n0 improper types\n\n")
        for bl in box_lines:
            f.write(bl + "\n")
        f.write("\nMasses\n\n")
        for ml in mass_lines:
            parts = ml.split()
            if len(parts) >= 2 and parts[0] in ("1", "2"):
                f.write(ml + "\n")
        f.write("\nAtoms # full\n\n")
        for a in cg_atoms:
            new_id = old_to_new[a["id"]]
            f.write(
                f"{new_id} {a['mol']} {a['type']} {a['charge']:.6f} "
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
