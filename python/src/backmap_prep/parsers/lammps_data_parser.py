"""LAMMPS data file parser for reading equilibrated CG configurations."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from .gro_parser import GroAtom, GroFile
from .top_parser import AtomType, MoleculeType, TopAtom, Topology

if TYPE_CHECKING:
    from pathlib import Path


@dataclass
class LammpsDataAtom:
    atom_id: int
    mol_id: int
    type_id: int
    charge: float
    x: float
    y: float
    z: float
    ix: int = 0
    iy: int = 0
    iz: int = 0


@dataclass
class LammpsDataBond:
    bond_id: int
    type_id: int
    i: int
    j: int


@dataclass
class LammpsDataFile:
    atoms: list[LammpsDataAtom] = field(default_factory=list)
    bonds: list[LammpsDataBond] = field(default_factory=list)
    box: tuple[float, float, float] = (0.0, 0.0, 0.0)
    n_atom_types: int = 0
    n_bond_types: int = 0


def parse_lammps_data(path: Path) -> LammpsDataFile:
    """Parse a LAMMPS data file (atom_style full).

    Reads atom positions (with image flags if present), bonds, and box
    dimensions. Used to read equilibrated CG frames for the rebuild workflow.
    """
    text = path.read_text()
    lines = text.splitlines()

    result = LammpsDataFile()
    section: str | None = None
    box = [0.0, 0.0, 0.0]

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
            section = "velocities"
            continue
        if stripped.startswith("Masses"):
            section = "masses"
            continue
        if stripped.startswith("Angles"):
            section = "angles"
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
        elif "atom types" in stripped:
            result.n_atom_types = int(stripped.split()[0])
        elif "bond types" in stripped:
            result.n_bond_types = int(stripped.split()[0])
        elif section == "atoms":
            parts = stripped.split()
            if len(parts) < 7:
                continue
            result.atoms.append(
                LammpsDataAtom(
                    atom_id=int(parts[0]),
                    mol_id=int(parts[1]),
                    type_id=int(parts[2]),
                    charge=float(parts[3]),
                    x=float(parts[4]),
                    y=float(parts[5]),
                    z=float(parts[6]),
                    ix=int(parts[7]) if len(parts) > 7 else 0,
                    iy=int(parts[8]) if len(parts) > 8 else 0,
                    iz=int(parts[9]) if len(parts) > 9 else 0,
                )
            )
        elif section == "bonds":
            parts = stripped.split()
            if len(parts) < 4:
                continue
            result.bonds.append(
                LammpsDataBond(
                    bond_id=int(parts[0]),
                    type_id=int(parts[1]),
                    i=int(parts[2]),
                    j=int(parts[3]),
                )
            )
        elif section in ("velocities", "masses", "angles"):
            pass

    result.box = (box[0], box[1], box[2])
    return result


_CG_SYSTEM_SECTIONS = frozenset(
    {"Masses", "Atoms", "Bonds", "Angles", "Dihedrals", "Impropers", "Velocities"}
)


def parse_cg_system(path: Path) -> tuple[GroFile, Topology]:
    """Parse a LAMMPS data file as a full CG system source (``cg_system.format: lammps``).

    Reads box bounds, the ``Masses`` section, and the ``Atoms # full`` section,
    and returns them as :class:`GroFile`/:class:`Topology` objects — the same
    types the GROMACS ``.gro``/``.top`` CG path produces — so ``build_system()``
    needs no further branching once this returns.

    Atom order is assumed to be contiguous per molecule instance: the first
    contiguous run of atoms sharing the first atom's molecule-ID is taken as
    one molecule template (this determines ``cg_atom_count`` the same way
    ``build_system()`` indexes CG atoms via ``cg_start = mol_idx * cg_atom_count``).
    ``Bonds``/``Angles``/etc. sections, if present, are tolerated but not read:
    CG-CG bonded connectivity is configured via ``cross_interactions`` in the
    YAML settings for both CG formats, not from ``cg_system`` itself.
    """
    text = path.read_text()

    box = [0.0, 0.0, 0.0]
    masses: dict[int, float] = {}
    # (atom_id, mol_id, type_id, charge, x, y, z)
    atom_rows: list[tuple[int, int, int, float, float, float, float]] = []
    section: str | None = None

    for raw_line in text.splitlines():
        stripped = raw_line.strip()
        if not stripped:
            continue

        header = stripped.split("#", 1)[0].strip()
        if header in _CG_SYSTEM_SECTIONS:
            section = header
            continue

        body = stripped.split("#", 1)[0].strip()
        if not body:
            continue

        if "xlo" in body and "xhi" in body:
            parts = body.split()
            box[0] = float(parts[1]) - float(parts[0])
        elif "ylo" in body and "yhi" in body:
            parts = body.split()
            box[1] = float(parts[1]) - float(parts[0])
        elif "zlo" in body and "zhi" in body:
            parts = body.split()
            box[2] = float(parts[1]) - float(parts[0])
        elif section == "Masses":
            parts = body.split()
            masses[int(parts[0])] = float(parts[1])
        elif section == "Atoms":
            parts = body.split()
            if len(parts) < 7:
                raise ValueError(f"malformed 'Atoms' line in {path}: {raw_line!r}")
            atom_rows.append(
                (
                    int(parts[0]),
                    int(parts[1]),
                    int(parts[2]),
                    float(parts[3]),
                    float(parts[4]),
                    float(parts[5]),
                    float(parts[6]),
                )
            )
        # Bonds/Angles/Dihedrals/Impropers/Velocities: intentionally not read.

    if not masses:
        raise ValueError(
            f"LAMMPS data file {path} has no 'Masses' section "
            "(required for cg_system.format: lammps)"
        )
    if not atom_rows:
        raise ValueError(
            f"LAMMPS data file {path} has no 'Atoms' section "
            "(required for cg_system.format: lammps)"
        )

    atom_rows.sort(key=lambda r: r[0])

    first_mol_id = atom_rows[0][1]
    cg_atom_count = 0
    for row in atom_rows:
        if row[1] != first_mol_id:
            break
        cg_atom_count += 1
    if len(atom_rows) % cg_atom_count != 0:
        raise ValueError(
            f"LAMMPS data file {path}: {len(atom_rows)} atoms is not a multiple of "
            f"{cg_atom_count} (the apparent per-molecule atom count from the first "
            "molecule-ID block) — check that atoms are grouped contiguously per molecule"
        )
    n_molecules = len(atom_rows) // cg_atom_count

    gro_atoms = [
        GroAtom(
            resid=mol_id,
            resname="CG",
            name=str(type_id),
            index=i + 1,
            x=x,
            y=y,
            z=z,
        )
        for i, (_atom_id, mol_id, type_id, _charge, x, y, z) in enumerate(atom_rows)
    ]
    gro = GroFile(title=path.name, atoms=gro_atoms, box=(box[0], box[1], box[2]))

    template_rows = atom_rows[:cg_atom_count]
    template_atoms = [
        TopAtom(
            index=i + 1,
            type=str(type_id),
            resid=mol_id,
            resname="CG",
            name=str(type_id),
            charge_group=mol_id,
            charge=charge,
            mass=masses.get(type_id, 0.0),
        )
        for i, (_atom_id, mol_id, type_id, charge, _x, _y, _z) in enumerate(template_rows)
    ]

    mol_name = "CG"
    atom_types = {
        str(type_id): AtomType(name=str(type_id), mass=mass, charge=0.0, ptype="V")
        for type_id, mass in masses.items()
    }
    top = Topology(
        atom_types=atom_types,
        molecule_types={mol_name: MoleculeType(name=mol_name, nrexcl=0, atoms=template_atoms)},
        molecules=[(mol_name, n_molecules)],
    )

    return gro, top
