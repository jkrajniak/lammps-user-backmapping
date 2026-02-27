"""GROMACS .top/.itp topology file parser."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class AtomType:
    name: str
    mass: float
    charge: float
    ptype: str  # 'A' for atom, 'V' for virtual/CG
    sigma: float = 0.0  # nm
    epsilon: float = 0.0  # kJ/mol


@dataclass
class TopAtom:
    index: int
    type: str
    resid: int
    resname: str
    name: str
    charge_group: int
    charge: float
    mass: float


@dataclass
class TopBond:
    i: int
    j: int
    func: int
    params: list[float] = field(default_factory=list)


@dataclass
class TopAngle:
    i: int
    j: int
    k: int
    func: int
    params: list[float] = field(default_factory=list)


@dataclass
class TopDihedral:
    i: int
    j: int
    k: int
    atom_l: int  # 4th atom index (named to avoid E741 lint)
    func: int
    params: list[float] = field(default_factory=list)


@dataclass
class MoleculeType:
    name: str
    nrexcl: int
    atoms: list[TopAtom] = field(default_factory=list)
    bonds: list[TopBond] = field(default_factory=list)
    angles: list[TopAngle] = field(default_factory=list)
    dihedrals: list[TopDihedral] = field(default_factory=list)


@dataclass
class Topology:
    atom_types: dict[str, AtomType] = field(default_factory=dict)
    molecule_types: dict[str, MoleculeType] = field(default_factory=dict)
    molecules: list[tuple[str, int]] = field(default_factory=list)
    combination_rule: int = 2
    fudge_lj: float = 1.0
    fudge_qq: float = 1.0
    defaults_gen_pairs: str = "yes"


def parse_top(path: Path, include_dirs: list[Path] | None = None) -> Topology:
    """Parse a GROMACS .top or .itp file, resolving #include directives."""
    top = Topology()
    _parse_file(Path(path), top, include_dirs or [])
    return top


def _parse_file(path: Path, top: Topology, include_dirs: list[Path]) -> None:
    """Recursively parse a topology file."""
    if not path.exists():
        raise FileNotFoundError(f"Topology file not found: {path}")

    lines = _preprocess(path, include_dirs)
    section: str | None = None
    current_mol: MoleculeType | None = None

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith(";"):
            continue

        # Remove inline comments
        if ";" in stripped:
            stripped = stripped[: stripped.index(";")].strip()

        # Section headers
        if stripped.startswith("["):
            section = stripped.strip("[] ").lower()
            continue

        tokens = stripped.split()

        if section == "defaults":
            if len(tokens) >= 2:
                top.combination_rule = int(tokens[1])
            if len(tokens) >= 4:
                top.fudge_lj = float(tokens[3])
            if len(tokens) >= 5:
                top.fudge_qq = float(tokens[4])

        elif section == "atomtypes":
            _parse_atomtype(tokens, top)

        elif section == "moleculetype":
            if len(tokens) >= 2:
                current_mol = MoleculeType(name=tokens[0], nrexcl=int(tokens[1]))
                top.molecule_types[tokens[0]] = current_mol

        elif section == "atoms" and current_mol:
            _parse_atom(tokens, current_mol)

        elif section == "bonds" and current_mol:
            _parse_bond(tokens, current_mol)

        elif section == "angles" and current_mol:
            _parse_angle(tokens, current_mol)

        elif section == "dihedrals" and current_mol:
            _parse_dihedral(tokens, current_mol)

        elif section == "molecules":
            if len(tokens) >= 2:
                top.molecules.append((tokens[0], int(tokens[1])))


def _preprocess(path: Path, include_dirs: list[Path]) -> list[str]:
    """Read file and resolve #include directives, skip #ifdef blocks."""
    result: list[str] = []
    defines: dict[str, str] = {}
    skip_depth = 0

    for raw_line in path.read_text().splitlines():
        stripped = raw_line.strip()

        if stripped.startswith(("#ifdef", "#ifndef")):
            skip_depth += 1
            continue
        if stripped.startswith("#else"):
            continue
        if stripped.startswith("#endif"):
            if skip_depth > 0:
                skip_depth -= 1
            continue

        if skip_depth > 0:
            continue

        if stripped.startswith("#define"):
            parts = stripped.split(None, 2)
            if len(parts) >= 3:
                defines[parts[1]] = parts[2]
            continue

        if stripped.startswith("#include"):
            inc_file = stripped.split('"')[1] if '"' in stripped else stripped.split()[1]
            inc_path = _resolve_include(inc_file, path.parent, include_dirs)
            if inc_path:
                result.extend(_preprocess(inc_path, include_dirs))
            continue

        result.append(raw_line)

    return result


def _resolve_include(filename: str, base_dir: Path, include_dirs: list[Path]) -> Path | None:
    """Find an included file in base_dir or include_dirs."""
    candidate = base_dir / filename
    if candidate.exists():
        return candidate
    for d in include_dirs:
        candidate = d / filename
        if candidate.exists():
            return candidate
    return None


def _parse_atomtype(tokens: list[str], top: Topology) -> None:
    """Parse an atomtypes entry. Handles varying column counts."""
    if len(tokens) < 6:
        return

    # Common format: name  bond_type  mass  charge  ptype  sigma  epsilon
    # or:            name  at_num  mass  charge  ptype  sigma  epsilon
    name = tokens[0]

    # Detect format by finding ptype column (A, V, S, D)
    ptype_idx = -1
    for i in range(2, min(len(tokens), 6)):
        if tokens[i] in ("A", "V", "S", "D"):
            ptype_idx = i
            break

    if ptype_idx < 0:
        return

    mass = float(tokens[ptype_idx - 2])
    charge = float(tokens[ptype_idx - 1])
    ptype = tokens[ptype_idx]
    sigma = float(tokens[ptype_idx + 1]) if ptype_idx + 1 < len(tokens) else 0.0
    epsilon = float(tokens[ptype_idx + 2]) if ptype_idx + 2 < len(tokens) else 0.0

    top.atom_types[name] = AtomType(
        name=name,
        mass=mass,
        charge=charge,
        ptype=ptype,
        sigma=sigma,
        epsilon=epsilon,
    )


def _parse_atom(tokens: list[str], mol: MoleculeType) -> None:
    if len(tokens) < 7:
        return
    mol.atoms.append(
        TopAtom(
            index=int(tokens[0]),
            type=tokens[1],
            resid=int(tokens[2]),
            resname=tokens[3],
            name=tokens[4],
            charge_group=int(tokens[5]),
            charge=float(tokens[6]),
            mass=float(tokens[7]) if len(tokens) > 7 else 0.0,
        )
    )


def _parse_bond(tokens: list[str], mol: MoleculeType) -> None:
    if len(tokens) < 3:
        return
    params = [float(t) for t in tokens[3:]]
    mol.bonds.append(
        TopBond(
            i=int(tokens[0]),
            j=int(tokens[1]),
            func=int(tokens[2]),
            params=params,
        )
    )


def _parse_angle(tokens: list[str], mol: MoleculeType) -> None:
    if len(tokens) < 4:
        return
    params = [float(t) for t in tokens[4:]]
    mol.angles.append(
        TopAngle(
            i=int(tokens[0]),
            j=int(tokens[1]),
            k=int(tokens[2]),
            func=int(tokens[3]),
            params=params,
        )
    )


def _parse_dihedral(tokens: list[str], mol: MoleculeType) -> None:
    if len(tokens) < 5:
        return
    params = [float(t) for t in tokens[5:]]
    mol.dihedrals.append(
        TopDihedral(
            i=int(tokens[0]),
            j=int(tokens[1]),
            k=int(tokens[2]),
            atom_l=int(tokens[3]),
            func=int(tokens[4]),
            params=params,
        )
    )
