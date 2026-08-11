"""LAMMPS data file parser for reading equilibrated CG configurations."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from .gro_parser import GroAtom, GroFile
from .lammps_script_parser import parse_at_fragment_script
from .top_parser import AtomType, MoleculeType, TopAngle, TopAtom, TopBond, TopDihedral, Topology

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


_DATA_FILE_SECTION_HEADERS = frozenset(
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
        if header in _DATA_FILE_SECTION_HEADERS:
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


def parse_at_fragment(data_path: Path, script_path: Path) -> tuple[GroFile, Topology]:
    """Parse a LAMMPS-native AT fragment (``molecules[].source.format: lammps``).

    Unlike :func:`parse_cg_system`, the `data` file's ``Bonds``/``Angles``/
    ``Dihedrals`` sections ARE read here: `builder.py` consumes
    `at_mol.bonds`/`.angles`/`.dihedrals` directly to build the intra-bead
    bond/angle/dihedral coefficients, unlike the CG side where the
    equivalent sections are parsed but never consumed. Coefficients
    themselves come from *script_path* (see
    :func:`parsers.lammps_script_parser.parse_at_fragment_script`) — a
    `data` file alone cannot supply them.

    Each atom's LAMMPS atom ID (as a string) is used for ``TopAtom.name``,
    mirroring how :func:`parse_cg_system` uses the numeric type ID for
    ``TopAtom.type`` — this is what ``beads[].atoms`` entries reference for
    a LAMMPS-native fragment, instead of a symbolic atom name.

    Atom IDs SHALL be contiguous ``1..N`` with no gaps (matching how
    ``at_mol.atoms[bond.i - 1]``-style indexing already assumes GROMACS
    `.top` files number atoms 1..N via the ``nr`` column).

    Coefficient values are used as-is, with no unit conversion (the
    input script must declare ``units real``) and, for dihedrals, no
    GROMACS RB-to-LAMMPS conversion (`dihedral_coeff` is already in this
    package's native `dihedral_style ryckaert` form).
    """
    text = data_path.read_text()

    masses: dict[int, float] = {}
    # (atom_id, mol_id, type_id, charge, x, y, z)
    atom_rows: list[tuple[int, int, int, float, float, float, float]] = []
    bond_rows: list[tuple[int, int, int]] = []  # (type_id, i, j)
    angle_rows: list[tuple[int, int, int, int]] = []  # (type_id, i, j, k)
    dihedral_rows: list[tuple[int, int, int, int, int]] = []  # (type_id, i, j, k, l)
    section: str | None = None

    for raw_line in text.splitlines():
        stripped = raw_line.strip()
        if not stripped:
            continue

        header = stripped.split("#", 1)[0].strip()
        if header in _DATA_FILE_SECTION_HEADERS:
            section = header
            continue

        body = stripped.split("#", 1)[0].strip()
        if not body:
            continue

        if section == "Masses":
            parts = body.split()
            masses[int(parts[0])] = float(parts[1])
        elif section == "Atoms":
            parts = body.split()
            if len(parts) < 7:
                raise ValueError(f"malformed 'Atoms' line in {data_path}: {raw_line!r}")
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
        elif section == "Bonds":
            parts = body.split()
            if len(parts) < 4:
                raise ValueError(f"malformed 'Bonds' line in {data_path}: {raw_line!r}")
            bond_rows.append((int(parts[1]), int(parts[2]), int(parts[3])))
        elif section == "Angles":
            parts = body.split()
            if len(parts) < 5:
                raise ValueError(f"malformed 'Angles' line in {data_path}: {raw_line!r}")
            angle_rows.append((int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4])))
        elif section == "Dihedrals":
            parts = body.split()
            if len(parts) < 6:
                raise ValueError(f"malformed 'Dihedrals' line in {data_path}: {raw_line!r}")
            dihedral_rows.append(
                (int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4]), int(parts[5]))
            )
        # Impropers/Velocities: not read (builder.py has no AT-fragment improper handling).

    if not masses:
        raise ValueError(
            f"LAMMPS data file {data_path} has no 'Masses' section "
            "(required for molecules[].source.format: lammps)"
        )
    if not atom_rows:
        raise ValueError(
            f"LAMMPS data file {data_path} has no 'Atoms' section "
            "(required for molecules[].source.format: lammps)"
        )

    atom_rows.sort(key=lambda r: r[0])
    expected_ids = list(range(1, len(atom_rows) + 1))
    if [row[0] for row in atom_rows] != expected_ids:
        raise ValueError(
            f"LAMMPS data file {data_path}: atom IDs must be contiguous 1..N with no gaps"
        )

    coeffs = parse_at_fragment_script(script_path)

    # Namespaced distinctly from CG numeric type IDs ("1", "2", ...): builder.py
    # keys its shared `type_map` dict by this string across CG and AT atoms in
    # the same build, so an AT type and a CG type both named "1" would
    # otherwise silently collide and merge into one type.
    def _at_type_name(type_id: int) -> str:
        return f"AT{type_id}"

    gro_atoms = [
        GroAtom(resid=mol_id, resname="AT", name=str(atom_id), index=i + 1, x=x, y=y, z=z)
        for i, (atom_id, mol_id, _type_id, _charge, x, y, z) in enumerate(atom_rows)
    ]
    gro = GroFile(title=data_path.name, atoms=gro_atoms, box=(0.0, 0.0, 0.0))

    top_atoms = [
        TopAtom(
            index=atom_id,
            type=_at_type_name(type_id),
            resid=mol_id,
            resname="AT",
            name=str(atom_id),
            charge_group=mol_id,
            charge=charge,
            mass=masses.get(type_id, 0.0),
        )
        for (atom_id, mol_id, type_id, charge, _x, _y, _z) in atom_rows
    ]

    def _missing(kind: str, type_id: int) -> ValueError:
        return ValueError(
            f"LAMMPS data file {data_path} references {kind} type {type_id}, but "
            f"{script_path} has no matching {kind}_coeff (required for "
            "molecules[].source.format: lammps)"
        )

    top_bonds: list[TopBond] = []
    for type_id, i, j in bond_rows:
        if type_id not in coeffs.bond:
            raise _missing("bond", type_id)
        k, r0 = coeffs.bond[type_id]
        top_bonds.append(TopBond(i=i, j=j, func=1, params=[r0, k]))

    top_angles: list[TopAngle] = []
    for type_id, i, j, k_idx in angle_rows:
        if type_id not in coeffs.angle:
            raise _missing("angle", type_id)
        k, theta0 = coeffs.angle[type_id]
        top_angles.append(TopAngle(i=i, j=j, k=k_idx, func=1, params=[theta0, k]))

    top_dihedrals: list[TopDihedral] = []
    for type_id, i, j, k_idx, l_idx in dihedral_rows:
        if type_id not in coeffs.dihedral:
            raise _missing("dihedral", type_id)
        top_dihedrals.append(
            TopDihedral(i=i, j=j, k=k_idx, atom_l=l_idx, func=3, params=coeffs.dihedral[type_id])
        )

    referenced_atom_types = {row[2] for row in atom_rows}
    for type_id in referenced_atom_types:
        if type_id not in coeffs.pair:
            raise _missing("pair", type_id)

    atom_types = {
        _at_type_name(type_id): AtomType(
            name=_at_type_name(type_id),
            mass=masses.get(type_id, 0.0),
            charge=0.0,
            ptype="A",
            sigma=coeffs.pair[type_id][1],
            epsilon=coeffs.pair[type_id][0],
        )
        for type_id in referenced_atom_types
    }

    mol_name = "AT_FRAGMENT"
    top = Topology(
        atom_types=atom_types,
        molecule_types={
            mol_name: MoleculeType(
                name=mol_name,
                nrexcl=0,
                atoms=top_atoms,
                bonds=top_bonds,
                angles=top_angles,
                dihedrals=top_dihedrals,
            )
        },
        molecules=[(mol_name, 1)],
    )

    return gro, top
