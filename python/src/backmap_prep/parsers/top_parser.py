"""GROMACS .top/.itp topology file parser."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class AtomType:
    name: str
    mass: float
    charge: float
    ptype: str  # 'A' for atom, 'V' for virtual/CG
    bond_type: str = ""  # GROMACS interaction class (CT, OS, HC, …)
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
    param_tokens: list[str] = field(default_factory=list)
    # func 9 with explicit parameters: (phase, k, n) of every consecutive line
    merged_terms: list[tuple[float, float, int]] = field(default_factory=list)


@dataclass
class TopPair:
    i: int
    j: int
    func: int
    params: list[float] = field(default_factory=list)


@dataclass
class DihedralTypeEntry:
    func: int
    params: list[float] = field(default_factory=list)
    # func 9: every (phase, k, n) line of the entry; consecutive lines with the
    # same atom types add terms, as in GROMACS.
    terms: list[list[float]] = field(default_factory=list)


@dataclass
class MoleculeType:
    name: str
    nrexcl: int
    atoms: list[TopAtom] = field(default_factory=list)
    bonds: list[TopBond] = field(default_factory=list)
    cross_bonds: list[TopBond] = field(default_factory=list)
    angles: list[TopAngle] = field(default_factory=list)
    cross_angles: list[TopAngle] = field(default_factory=list)
    dihedrals: list[TopDihedral] = field(default_factory=list)
    cross_dihedrals: list[TopDihedral] = field(default_factory=list)
    cross_pairs: list[TopPair] = field(default_factory=list)
    pairs: list[TopPair] = field(default_factory=list)


@dataclass
class Topology:
    atom_types: dict[str, AtomType] = field(default_factory=dict)
    dihedraltypes: dict[str, dict[str, dict[str, dict[str, DihedralTypeEntry]]]] = field(
        default_factory=dict
    )
    molecule_types: dict[str, MoleculeType] = field(default_factory=dict)
    molecules: list[tuple[str, int]] = field(default_factory=list)
    combination_rule: int = 2
    fudge_lj: float = 1.0
    fudge_qq: float = 1.0
    defaults_gen_pairs: str = "yes"
    has_defaults: bool = False  # True once a [ defaults ] section was read
    # [ bondtypes ] / [ angletypes ] (func 1): bond classes -> GROMACS params,
    # stored under both atom orders.
    bondtypes: dict[tuple[str, str], list[float]] = field(default_factory=dict)
    angletypes: dict[tuple[str, str, str], list[float]] = field(default_factory=dict)
    # [ angletypes ] func 5 (Urey-Bradley): theta0, k_theta, r13, k_UB
    angletypes_ub: dict[tuple[str, str, str], list[float]] = field(default_factory=dict)
    # [ nonbond_params ] func 1: (type_i, type_j) -> (c6/sigma, c12/epsilon) as written
    nonbond_params: dict[tuple[str, str], tuple[float, float]] = field(default_factory=dict)
    _last_dihedraltype: tuple[str, str, str, str, int] | None = None


def parse_top(
    path: Path,
    include_dirs: list[Path] | None = None,
    forcefield_dirs: list[Path] | None = None,
) -> Topology:
    """Parse a GROMACS .top or .itp file, resolving #include directives."""
    top = Topology()
    ff_dirs = forcefield_dirs or []
    _parse_file(Path(path), top, include_dirs or [], ff_dirs)
    if ff_dirs and not top.dihedraltypes:
        _merge_forcefield_types(top, ff_dirs)
    return top


def _merge_forcefield_types(top: Topology, forcefield_dirs: list[Path]) -> None:
    """Load OPLS dihedraltypes/atomtypes from bundled forcefield when TOP has no include."""
    for ff_dir in forcefield_dirs:
        ff_itp = ff_dir / "oplsaa.ff" / "forcefield.itp"
        if not ff_itp.is_file():
            continue
        ff_top = Topology()
        _parse_file(ff_itp, ff_top, [ff_dir / "oplsaa.ff"], forcefield_dirs)
        for name, entry in ff_top.atom_types.items():
            top.atom_types.setdefault(name, entry)
        for key, params in ff_top.bondtypes.items():
            top.bondtypes.setdefault(key, params)
        for akey, aparams in ff_top.angletypes.items():
            top.angletypes.setdefault(akey, aparams)
        for ukey, uparams in ff_top.angletypes_ub.items():
            top.angletypes_ub.setdefault(ukey, uparams)
        for i, j_map in ff_top.dihedraltypes.items():
            for j, k_map in j_map.items():
                for k, l_map in k_map.items():
                    for l, entry in l_map.items():
                        _store_dihedraltype(top, i, j, k, l, entry)
        return


def _parse_file(
    path: Path,
    top: Topology,
    include_dirs: list[Path],
    forcefield_dirs: list[Path],
) -> None:
    """Recursively parse a topology file."""
    if not path.exists():
        raise FileNotFoundError(f"Topology file not found: {path}")

    lines = _preprocess(path, include_dirs, forcefield_dirs)
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
            top.has_defaults = True
            if len(tokens) >= 2:
                top.combination_rule = int(tokens[1])
            if len(tokens) >= 3:
                top.defaults_gen_pairs = tokens[2].lower()
            if len(tokens) >= 4:
                top.fudge_lj = float(tokens[3])
            if len(tokens) >= 5:
                top.fudge_qq = float(tokens[4])

        elif section == "atomtypes":
            _parse_atomtype(tokens, top)

        elif section == "dihedraltypes":
            _parse_dihedraltype(tokens, top)

        elif section == "nonbond_params":
            if len(tokens) >= 5 and tokens[2] == "1":
                values = (float(tokens[3]), float(tokens[4]))
                top.nonbond_params[(tokens[0], tokens[1])] = values
                top.nonbond_params[(tokens[1], tokens[0])] = values

        elif section == "bondtypes":
            if len(tokens) >= 5 and tokens[2] == "1":
                params = [float(t) for t in tokens[3:5]]
                top.bondtypes[(tokens[0], tokens[1])] = params
                top.bondtypes[(tokens[1], tokens[0])] = params

        elif section == "angletypes":
            if len(tokens) >= 6 and tokens[3] == "1":
                params = [float(t) for t in tokens[4:6]]
                top.angletypes[(tokens[0], tokens[1], tokens[2])] = params
                top.angletypes[(tokens[2], tokens[1], tokens[0])] = params
            elif len(tokens) >= 8 and tokens[3] == "5":
                params = [float(t) for t in tokens[4:8]]
                top.angletypes_ub[(tokens[0], tokens[1], tokens[2])] = params
                top.angletypes_ub[(tokens[2], tokens[1], tokens[0])] = params

        elif section == "moleculetype":
            if len(tokens) >= 2:
                current_mol = MoleculeType(name=tokens[0], nrexcl=int(tokens[1]))
                top.molecule_types[tokens[0]] = current_mol

        elif section == "atoms" and current_mol:
            _parse_atom(tokens, current_mol)

        elif section == "bonds" and current_mol:
            _parse_bond(tokens, current_mol)

        elif section == "cross_bonds" and current_mol:
            _parse_cross_bond(tokens, current_mol)

        elif section == "angles" and current_mol:
            _parse_angle(tokens, current_mol)

        elif section == "cross_angles" and current_mol:
            _parse_cross_angle(tokens, current_mol)

        elif section == "dihedrals" and current_mol:
            _parse_dihedral(tokens, current_mol)

        elif section == "cross_dihedrals" and current_mol:
            _parse_cross_dihedral(tokens, current_mol)

        elif section == "cross_pairs" and current_mol:
            _parse_cross_pair(tokens, current_mol)

        elif section == "pairs" and current_mol:
            _parse_cross_pair(tokens, current_mol, target="pairs")

        elif section == "molecules":
            if len(tokens) >= 2:
                top.molecules.append((tokens[0], int(tokens[1])))


def _preprocess(
    path: Path,
    include_dirs: list[Path],
    forcefield_dirs: list[Path],
    defines: dict[str, str] | None = None,
) -> list[str]:
    """Read file and resolve #include directives, evaluating #ifdef/#ifndef/#else
    blocks against `defines` (shared across the whole recursive parse, so a
    #define earlier in this file or an included one is visible everywhere
    after it -- matching the C-preprocessor-like semantics GROMACS .itp files
    rely on).

    Each stack entry is (branch_condition_met, else_branch_active): #else
    flips which branch is active without losing the original condition (a
    naive skip-depth counter previously skipped BOTH branches of every
    #ifdef/#else/#endif unconditionally -- e.g. oplsaa.ff/ffnonbonded.itp's
    #ifdef HEAVY_H / #else / #endif around opls_116/opls_117 (water O/H)
    caused their LJ parameters to silently resolve to 0.0/0.0 regardless of
    HEAVY_H, see experiments/20260724_pet-pairwise-diagnostic-cgcg-collapse.md).
    """
    result: list[str] = []
    if defines is None:
        defines = {}
    # stack of (condition_met_for_this_branch, is_negated_ifndef)
    cond_stack: list[bool] = []

    def active() -> bool:
        return all(cond_stack)

    for raw_line in path.read_text().splitlines():
        stripped = raw_line.strip()

        if stripped.startswith("#ifdef"):
            macro = stripped.split(None, 1)[1].strip() if len(stripped.split()) > 1 else ""
            cond_stack.append(macro in defines)
            continue
        if stripped.startswith("#ifndef"):
            macro = stripped.split(None, 1)[1].strip() if len(stripped.split()) > 1 else ""
            cond_stack.append(macro not in defines)
            continue
        if stripped.startswith("#else"):
            if cond_stack:
                cond_stack[-1] = not cond_stack[-1]
            continue
        if stripped.startswith("#endif"):
            if cond_stack:
                cond_stack.pop()
            continue

        if not active():
            continue

        if stripped.startswith("#define"):
            parts = stripped.split(None, 2)
            if len(parts) >= 2:
                defines[parts[1]] = parts[2] if len(parts) >= 3 else ""
            continue

        if stripped.startswith("#include"):
            inc_file = stripped.split('"')[1] if '"' in stripped else stripped.split()[1]
            inc_path = _resolve_include(inc_file, path.parent, include_dirs, forcefield_dirs)
            if inc_path:
                result.extend(_preprocess(inc_path, include_dirs, forcefield_dirs, defines))
            continue

        result.append(raw_line)

    return result


def _resolve_include(
    filename: str,
    base_dir: Path,
    include_dirs: list[Path],
    forcefield_dirs: list[Path],
) -> Path | None:
    """Find an included file in base_dir, include_dirs, or forcefield_dirs."""
    absolute = Path(filename)
    if absolute.is_absolute() and absolute.exists():
        return absolute

    candidates: list[Path] = [base_dir / filename]
    candidates.extend(d / filename for d in include_dirs)
    candidates.extend(d / filename for d in forcefield_dirs)

    # Remap standard GROMACS install paths when forcefield is bundled locally.
    if absolute.is_absolute():
        for marker in ("/top/", "/share/gromacs/top/"):
            if marker in filename:
                tail = filename.split(marker, 1)[1]
                candidates.extend(d / tail for d in forcefield_dirs)
                candidates.extend(d / Path(tail).name for d in forcefield_dirs)

    for env_name in ("GMXDATA", "GROMACS_DATA"):
        env_root = os.environ.get(env_name)
        if not env_root:
            continue
        candidates.append(Path(env_root) / "top" / filename)
        candidates.append(Path(env_root) / "share" / "gromacs" / "top" / filename)

    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _store_dihedraltype(
    top: Topology,
    type_i: str,
    type_j: str,
    type_k: str,
    type_l: str,
    entry: DihedralTypeEntry,
) -> None:
    for quad in ((type_i, type_j, type_k, type_l), (type_l, type_k, type_j, type_i)):
        i, j, k, l = quad
        top.dihedraltypes.setdefault(i, {}).setdefault(j, {}).setdefault(k, {})[l] = entry


def _parse_dihedraltype(tokens: list[str], top: Topology) -> None:
    if len(tokens) < 5:
        return
    try:
        params = [float(t) for t in tokens[5:]]
    except ValueError:
        return
    func = int(tokens[4])
    key = (tokens[0], tokens[1], tokens[2], tokens[3], func)
    if func == 9 and top._last_dihedraltype == key:
        existing = top.dihedraltypes[key[0]][key[1]][key[2]][key[3]]
        existing.terms.append(params[:3])
        return
    entry = DihedralTypeEntry(func=func, params=params)
    if func == 9:
        entry.terms = [params[:3]]
    _store_dihedraltype(top, tokens[0], tokens[1], tokens[2], tokens[3], entry)
    top._last_dihedraltype = key


def _lookup_dihedraltype(
    dihedraltypes: dict[str, dict[str, dict[str, dict[str, DihedralTypeEntry]]]],
    type_i: str,
    type_j: str,
    type_k: str,
    type_l: str,
) -> DihedralTypeEntry | None:
    """Lookup dihedraltypes with GROMACS wildcard and reverse-order matching."""
    wildcard = "X"
    quads = (
        (type_i, type_j, type_k, type_l),
        (wildcard, type_j, type_k, type_l),
        (type_i, type_j, type_k, wildcard),
        (wildcard, type_j, type_k, wildcard),
        (type_l, type_k, type_j, type_i),
        (wildcard, type_k, type_j, type_i),
        (type_l, type_k, type_j, wildcard),
        (wildcard, type_k, type_j, wildcard),
    )
    for i, j, k, l in quads:
        entry = dihedraltypes.get(i, {}).get(j, {}).get(k, {}).get(l)
        if entry is not None:
            return entry
    return None


def lookup_dihedraltype_most_specific(
    dihedraltypes: dict[str, dict[str, dict[str, dict[str, DihedralTypeEntry]]]],
    quad: tuple[str, str, str, str],
    funcs: set[int],
) -> DihedralTypeEntry | None:
    """Best [ dihedraltypes ] match for ``quad`` among entries of ``funcs``.

    As in GROMACS, the entry with the fewest ``X`` wildcards wins (any position,
    either atom order); used for func 2/4/9, where wildcards also sit in the
    middle (e.g. improper ``OBL X X CL``). func 3 keeps ``_lookup_dihedraltype``.
    """
    best: DihedralTypeEntry | None = None
    best_wild = 5
    for mask in range(16):
        wild = bin(mask).count("1")
        if wild >= best_wild:
            continue
        for order in (quad, quad[::-1]):
            key = tuple("X" if mask >> n & 1 else order[n] for n in range(4))
            entry = dihedraltypes.get(key[0], {}).get(key[1], {}).get(key[2], {}).get(key[3])
            if entry is not None and entry.func in funcs:
                best, best_wild = entry, wild
                break
    return best


def _dihedral_type_names(
    dih: TopDihedral, atom_types: dict[int, str], atom_type_defs: dict[str, AtomType] | None
) -> list[tuple[str, str, str, str]]:
    names = tuple(atom_types.get(i) for i in (dih.i, dih.j, dih.k, dih.atom_l))
    if not all(names):
        raise ValueError(
            f"Cannot resolve dihedral types for atoms {dih.i}-{dih.j}-{dih.k}-{dih.atom_l}"
        )
    variants = [names]
    if atom_type_defs:
        mapped = tuple(_interaction_type(n, atom_type_defs) for n in names)
        if mapped != names:
            variants.append(mapped)
    return variants  # type: ignore[return-value]


def resolve_periodic_terms(
    dih: TopDihedral,
    atom_types: dict[int, str],
    dihedraltypes: dict[str, dict[str, dict[str, dict[str, DihedralTypeEntry]]]],
    atom_type_defs: dict[str, AtomType] | None = None,
) -> list[tuple[float, float, int]]:
    """(phase deg, k kJ/mol, n) terms of a func 4/9 dihedral, explicit or from types."""
    if dih.params:
        if len(dih.params) < 3:
            raise ValueError(f"func {dih.func} dihedral needs phase, k, n: {dih}")
        return [(dih.params[0], dih.params[1], int(dih.params[2]))]
    for quad in _dihedral_type_names(dih, atom_types, atom_type_defs):
        entry = lookup_dihedraltype_most_specific(dihedraltypes, quad, {dih.func})
        if entry is not None:
            rows = entry.terms or [entry.params[:3]]
            return [(r[0], r[1], int(r[2])) for r in rows]
    raise ValueError(
        f"Missing func-{dih.func} dihedraltypes entry for atoms "
        f"{dih.i}-{dih.j}-{dih.k}-{dih.atom_l}"
    )


def resolve_improper_harmonic(
    dih: TopDihedral,
    atom_types: dict[int, str],
    dihedraltypes: dict[str, dict[str, dict[str, dict[str, DihedralTypeEntry]]]],
    atom_type_defs: dict[str, AtomType] | None = None,
) -> tuple[float, float]:
    """(xi0 deg, k kJ/mol/rad^2) of a func 2 improper, explicit or from types."""
    if dih.params:
        if len(dih.params) < 2:
            raise ValueError(f"func 2 improper needs xi0, k: {dih}")
        return dih.params[0], dih.params[1]
    for quad in _dihedral_type_names(dih, atom_types, atom_type_defs):
        entry = lookup_dihedraltype_most_specific(dihedraltypes, quad, {2})
        if entry is not None:
            return entry.params[0], entry.params[1]
    raise ValueError(
        f"Missing func-2 dihedraltypes entry for improper {dih.i}-{dih.j}-{dih.k}-{dih.atom_l}"
    )


def resolve_dihedral_params(
    dih: TopDihedral,
    atom_types: dict[int, str],
    dihedraltypes: dict[str, dict[str, dict[str, dict[str, DihedralTypeEntry]]]],
    atom_type_defs: dict[str, AtomType] | None = None,
) -> tuple[int, list[float]]:
    """Return (func, params) for a dihedral, resolving from dihedraltypes if needed."""
    if dih.params:
        return dih.func, dih.params

    type_i = atom_types.get(dih.i)
    type_j = atom_types.get(dih.j)
    type_k = atom_types.get(dih.k)
    type_l = atom_types.get(dih.atom_l)
    if not all([type_i, type_j, type_k, type_l]):
        raise ValueError(
            f"Cannot resolve dihedral types for atoms {dih.i}-{dih.j}-{dih.k}-{dih.atom_l}"
        )

    quad_variants: list[tuple[str, str, str, str]] = [(type_i, type_j, type_k, type_l)]
    if atom_type_defs:
        mapped = tuple(_interaction_type(name, atom_type_defs) for name in quad_variants[0])
        if mapped != quad_variants[0]:
            quad_variants.append(mapped)

    for quad in quad_variants:
        entry = _lookup_dihedraltype(dihedraltypes, *quad)
        if entry is not None:
            return entry.func, [float(p) for p in entry.params]

    raise ValueError(
        "Missing dihedraltypes entry for "
        f"{type_i}-{type_j}-{type_k}-{type_l} (atoms {dih.i}-{dih.j}-{dih.k}-{dih.atom_l})"
    )


def _interaction_type(type_name: str, atom_type_defs: dict[str, AtomType]) -> str:
    """Map an atom type name to its GROMACS interaction class when available."""
    entry = atom_type_defs.get(type_name)
    if entry is not None and entry.bond_type:
        return entry.bond_type
    return type_name


# OPLS improper defines from oplsaa.ff/ffbonded.itp (func 1: phi0 deg, k kJ/mol, n)
OPLS_IMPROPER_DEFINES: dict[str, tuple[float, float, int]] = {
    "improper_O_C_X_Y": (180.0, 43.93200, 2),
    "improper_X_NO_ON_NO": (180.0, 43.93200, 2),
    "improper_N2_X_N2_N2": (180.0, 43.93200, 2),
    "improper_Z_N_X_Y": (180.0, 4.18400, 2),
    "improper_Z_CM_X_Y": (180.0, 62.76000, 2),
    "improper_Z_CA_X_Y": (180.0, 4.60240, 2),
}


def resolve_opls_improper_params(dih: TopDihedral) -> tuple[float, float, int] | None:
    """Resolve func-1 OPLS improper dihedral parameters from inline or define tokens."""
    if len(dih.params) >= 3:
        return dih.params[0], dih.params[1], int(dih.params[2])
    if dih.param_tokens:
        define = OPLS_IMPROPER_DEFINES.get(dih.param_tokens[0])
        if define is not None:
            return define
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
    bond_type = tokens[1] if ptype_idx >= 5 else ""

    top.atom_types[name] = AtomType(
        name=name,
        mass=mass,
        charge=charge,
        ptype=ptype,
        bond_type=bond_type,
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
    try:
        params = [float(t) for t in tokens[3:]]
    except ValueError:
        return
    mol.bonds.append(
        TopBond(
            i=int(tokens[0]),
            j=int(tokens[1]),
            func=int(tokens[2]),
            params=params,
        )
    )


def _parse_cross_bond(tokens: list[str], mol: MoleculeType) -> None:
    if len(tokens) < 2:
        return
    if len(tokens) == 2:
        mol.cross_bonds.append(
            TopBond(
                i=int(tokens[0]),
                j=int(tokens[1]),
                func=0,
                params=[],
            )
        )
        return
    try:
        params = [float(t) for t in tokens[3:]]
    except ValueError:
        return
    mol.cross_bonds.append(
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
    try:
        params = [float(t) for t in tokens[4:]]
    except ValueError:
        return
    mol.angles.append(
        TopAngle(
            i=int(tokens[0]),
            j=int(tokens[1]),
            k=int(tokens[2]),
            func=int(tokens[3]),
            params=params,
        )
    )


def _parse_cross_angle(tokens: list[str], mol: MoleculeType) -> None:
    if len(tokens) < 3:
        return
    if len(tokens) == 3:
        mol.cross_angles.append(
            TopAngle(
                i=int(tokens[0]),
                j=int(tokens[1]),
                k=int(tokens[2]),
                func=0,
                params=[],
            )
        )
        return
    try:
        params = [float(t) for t in tokens[4:]]
    except ValueError:
        return
    mol.cross_angles.append(
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
    param_tokens = tokens[5:]
    try:
        params = [float(t) for t in param_tokens]
    except ValueError:
        params = []
    mol.dihedrals.append(
        TopDihedral(
            i=int(tokens[0]),
            j=int(tokens[1]),
            k=int(tokens[2]),
            atom_l=int(tokens[3]),
            func=int(tokens[4]),
            params=params,
            param_tokens=param_tokens,
        )
    )


def _parse_cross_dihedral(tokens: list[str], mol: MoleculeType) -> None:
    if len(tokens) < 4:
        return
    if len(tokens) == 4:
        mol.cross_dihedrals.append(
            TopDihedral(
                i=int(tokens[0]),
                j=int(tokens[1]),
                k=int(tokens[2]),
                atom_l=int(tokens[3]),
                func=0,
                params=[],
            )
        )
        return
    param_tokens = tokens[5:] if len(tokens) > 5 else []
    try:
        params = [float(t) for t in param_tokens]
    except ValueError:
        params = []
    mol.cross_dihedrals.append(
        TopDihedral(
            i=int(tokens[0]),
            j=int(tokens[1]),
            k=int(tokens[2]),
            atom_l=int(tokens[3]),
            func=int(tokens[4]),
            params=params,
            param_tokens=param_tokens,
        )
    )


def _parse_cross_pair(tokens: list[str], mol: MoleculeType, target: str = "cross_pairs") -> None:
    """Parse a ``[ cross_pairs ]`` or (``target="pairs"``) ``[ pairs ]`` line.

    A line with only the two atom indices uses function 1, as in GROMACS.
    """
    if len(tokens) < 2:
        return
    if len(tokens) == 2:
        tokens = [*tokens, "1"]
    try:
        params = [float(t) for t in tokens[3:]]
    except ValueError:
        params = []
    getattr(mol, target).append(
        TopPair(
            i=int(tokens[0]),
            j=int(tokens[1]),
            func=int(tokens[2]),
            params=params,
        )
    )


def resolve_pair_lj_params(
    top: Topology,
    atom_i: TopAtom,
    atom_j: TopAtom,
    func: int,
    params: list[float],
) -> tuple[float, float]:
    """Resolve GROMACS pair LJ parameters to LAMMPS real units."""
    from .. import units

    if func != 1:
        raise ValueError(f"Unsupported cross-pair func {func} for {atom_i.name}-{atom_j.name}")
    if len(params) >= 2:
        sigma_nm, eps_kj = params[0], params[1]
        if top.combination_rule == 1:
            sigma_nm, eps_kj = units.c6c12_to_sigma_epsilon(sigma_nm, eps_kj)
        return units.distance(sigma_nm), units.energy(eps_kj)
    type_i = top.atom_types.get(atom_i.type)
    type_j = top.atom_types.get(atom_j.type)
    if type_i is None or type_j is None:
        raise ValueError(f"Unknown atom types for pair {atom_i.type}-{atom_j.type}")
    sigma, epsilon = units.lj_pair_params(
        type_i.epsilon,
        type_j.epsilon,
        type_i.sigma,
        type_j.sigma,
        combination_rule=top.combination_rule,
        fudge_lj=top.fudge_lj,
    )
    return sigma, epsilon
