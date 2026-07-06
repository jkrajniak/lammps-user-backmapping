"""Build LAMMPS `System` from a network-engine hybrid GRO/TOP pair."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from backmap_prep import units
from backmap_prep.builder import (
    AngleTypeInfo,
    AtomTypeInfo,
    BondTypeInfo,
    DihedralTypeInfo,
    LammpsAngle,
    LammpsAtom,
    LammpsBond,
    LammpsCrossPair,
    LammpsDihedral,
    PairTypeInfo,
    System,
)
from backmap_prep.network.pbc import prepare_network_coordinates, validate_bond_geometry
from backmap_prep.parsers import parse_gro, parse_top
from backmap_prep.parsers.top_parser import (
    resolve_dihedral_params,
    resolve_opls_improper_params,
    resolve_pair_lj_params,
)

if TYPE_CHECKING:
    from backmap_prep.parsers.gro_parser import GroAtom
    from backmap_prep.parsers.top_parser import (
        MoleculeType,
        TopAngle,
        TopAtom,
        TopBond,
        TopDihedral,
        Topology,
    )
    from backmap_prep.schema import CrossAngle, CrossBond, CrossDihedral, Settings


def _parse_float_tokens(params: str | None) -> list[float]:
    if not params:
        return []
    tokens: list[float] = []
    for token in params.split():
        try:
            tokens.append(float(token))
        except ValueError:
            continue
    return tokens


def _token_name(value: str) -> str:
    return value.split(":")[-1]


def _single_letter_cg(atom_type: str) -> bool:
    return len(atom_type) == 1 and atom_type.isalpha() and atom_type.isupper()


def _cg_types_from_settings(settings: Settings) -> set[str]:
    return {bead.type for mol in settings.molecules for bead in mol.beads if bead.type}


def _cross_bond_defaults(
    settings: Settings,
) -> dict[frozenset[str], tuple[str, list[float], str | None]]:
    defaults: dict[frozenset[str], tuple[str, list[float], str | None]] = {}
    for cross_bond in settings.cross_interactions.bonds:
        default_keyword = "cg" if cross_bond.cg_bonded else "at"
        default_params = _cross_bond_params(cross_bond)
        for pair in cross_bond.pairs:
            pair_key = frozenset((_token_name(pair[0]), _token_name(pair[1])))
            defaults[pair_key] = (default_keyword, default_params, cross_bond.table)
    return defaults


def _cross_angle_defaults(
    settings: Settings,
) -> dict[tuple[str, str, str], tuple[str, list[float]]]:
    defaults: dict[tuple[str, str, str], tuple[str, list[float]]] = {}
    for cross_angle in settings.cross_interactions.angles:
        keyword = "cg" if cross_angle.cg_bonded else "at"
        params = _cross_angle_params(cross_angle)
        for triple in cross_angle.triples:
            defaults[tuple(_token_name(token) for token in triple)] = (keyword, params)
    return defaults


def _cross_bond_params(cross_bond: CrossBond) -> list[float]:
    values = _parse_float_tokens(cross_bond.params)
    if len(values) >= 3:
        return [units.spring_bond(values[2]), units.distance(values[1])]
    return [0.0, 0.0]


def _cross_angle_params(cross_angle: CrossAngle) -> list[float]:
    values = _parse_float_tokens(cross_angle.params)
    if len(values) >= 3:
        return [units.spring_angle(values[2]), values[1]]
    return [0.0, 0.0]


def _topology_molecule(topology: Topology) -> MoleculeType:
    molecule = next(iter(topology.molecule_types.values()), None)
    if molecule is None:
        raise ValueError("Hybrid topology has no molecule definitions")
    return molecule


def _atom_type_info(
    topology: Topology,
    molecule: MoleculeType,
    cg_type_names: set[str],
) -> tuple[list[AtomTypeInfo], dict[str, int]]:
    type_ids: dict[str, int] = {}
    type_samples: dict[str, TopAtom] = {}
    for atom in molecule.atoms:
        if atom.type not in type_ids:
            type_ids[atom.type] = len(type_ids) + 1
            type_samples[atom.type] = atom

    atom_types: list[AtomTypeInfo] = []
    for atom_type_name, type_id in type_ids.items():
        atom_type = topology.atom_types.get(atom_type_name)
        sample = type_samples[atom_type_name]
        is_cg = atom_type_name in cg_type_names or _single_letter_cg(atom_type_name)
        sigma = units.sigma(atom_type.sigma) if atom_type else 0.0
        epsilon = units.epsilon(atom_type.epsilon) if atom_type else 0.0
        atom_types.append(
            AtomTypeInfo(
                type_id=type_id,
                name=atom_type_name,
                mass=atom_type.mass if atom_type else sample.mass,
                is_cg=is_cg,
                sigma=sigma,
                epsilon=epsilon,
            )
        )
    return atom_types, type_ids


def _backmap_mol_ids(
    molecule: MoleculeType,
    cg_type_names: set[str],
) -> dict[int, int]:
    """Assign LAMMPS mol_id so each CG bead and its AT atoms share one molecule.

    Network hybrid topologies list atoms as CG bead followed by its mapped AT
    atoms (see bakery ``prepare_hybrid`` output). ``fix backmap`` requires
    ``n_at % n_cg == 0`` per mol_id, which holds when each group has exactly
    one CG bead.
    """
    mol_ids: dict[int, int] = {}
    current_mol_id = 0
    for atom in molecule.atoms:
        is_cg = atom.type in cg_type_names or _single_letter_cg(atom.type)
        if is_cg:
            current_mol_id += 1
        if current_mol_id == 0:
            raise ValueError(
                f"Hybrid topology atom {atom.index} ({atom.name}) precedes the first CG bead"
            )
        mol_ids[atom.index] = current_mol_id
    return mol_ids


def _plain_cg_style(style: str) -> str:
    """Map backmap/* styles to plain LAMMPS styles for CG-only dynamics."""
    if style == "backmap/table":
        return "table"
    if style == "backmap/harmonic":
        return "harmonic"
    if style == "backmap/ryckaert":
        return "ryckaert"
    if style == "backmap/charmm":
        return "charmm"
    return style


def _cg_atoms(
    molecule: MoleculeType,
    gro_atoms: list[GroAtom],
    type_ids: dict[str, int],
) -> list[LammpsAtom]:
    """Build CG-only atoms (single network molecule, mol_id=1)."""
    if len(gro_atoms) < len(molecule.atoms):
        raise ValueError(
            f"CG GRO has {len(gro_atoms)} atoms but topology expects {len(molecule.atoms)}"
        )

    atoms: list[LammpsAtom] = []
    for top_atom, gro_atom in zip(molecule.atoms, gro_atoms, strict=False):
        atoms.append(
            LammpsAtom(
                atom_id=top_atom.index,
                mol_id=1,
                type_id=type_ids[top_atom.type],
                charge=top_atom.charge,
                x=units.distance(gro_atom.x),
                y=units.distance(gro_atom.y),
                z=units.distance(gro_atom.z),
                type_name=top_atom.type,
                is_cg=True,
            )
        )
    return atoms


def _hybrid_atoms(
    molecule: MoleculeType,
    gro_atoms: list[GroAtom],
    type_ids: dict[str, int],
    cg_type_names: set[str],
) -> list[LammpsAtom]:
    if len(gro_atoms) < len(molecule.atoms):
        raise ValueError(
            f"Hybrid GRO has {len(gro_atoms)} atoms but topology expects {len(molecule.atoms)}"
        )

    mol_ids = _backmap_mol_ids(molecule, cg_type_names)
    atoms: list[LammpsAtom] = []
    for top_atom, gro_atom in zip(molecule.atoms, gro_atoms, strict=False):
        is_cg = top_atom.type in cg_type_names or _single_letter_cg(top_atom.type)
        atoms.append(
            LammpsAtom(
                atom_id=top_atom.index,
                mol_id=mol_ids[top_atom.index],
                type_id=type_ids[top_atom.type],
                charge=top_atom.charge,
                x=units.distance(gro_atom.x),
                y=units.distance(gro_atom.y),
                z=units.distance(gro_atom.z),
                type_name=top_atom.type,
                is_cg=is_cg,
            )
        )
    return atoms


def _add_bond_type(
    system: System,
    bond_type_map: dict[tuple[str, str, float, float, str], int],
    *,
    style: str,
    keyword: str,
    params: list[float],
    table_file: str | None = None,
) -> int:
    table_token = table_file or ""
    key = (style, keyword, round(params[0], 8), round(params[1], 8), table_token)
    type_id = bond_type_map.get(key)
    if type_id is not None:
        return type_id

    type_id = len(system.bond_types) + 1
    system.bond_types.append(
        BondTypeInfo(
            type_id=type_id,
            style=style,
            keyword=keyword,
            params=params,
            table_file=table_file,
            table_keyword="ENTRY" if table_file else None,
        )
    )
    bond_type_map[key] = type_id
    return type_id


def _add_angle_type(
    system: System,
    angle_type_map: dict[tuple[str, str, str, str, float, float], int],
    *,
    style: str = "backmap/harmonic",
    keyword: str,
    params: list[float],
    table_file: str | None = None,
) -> int:
    p0 = round(params[0], 8) if params else 0.0
    p1 = round(params[1], 8) if len(params) > 1 else 0.0
    key = (style, keyword, table_file or "", p0, p1)
    type_id = angle_type_map.get(key)
    if type_id is not None:
        return type_id
    type_id = len(system.angle_types) + 1
    system.angle_types.append(
        AngleTypeInfo(
            type_id=type_id,
            style=style,
            keyword=keyword,
            params=params,
            table_file=table_file,
            table_keyword="ENTRY" if table_file else None,
        )
    )
    angle_type_map[key] = type_id
    return type_id


def _cross_dihedral_defaults(
    settings: Settings,
) -> dict[tuple[str, str, str, str], tuple[str, list[float], str | None]]:
    defaults: dict[tuple[str, str, str, str], tuple[str, list[float], str | None]] = {}
    for cross_dihedral in settings.cross_interactions.dihedrals:
        keyword = "cg" if cross_dihedral.cg_bonded else "at"
        params, table_file = _cross_dihedral_params(cross_dihedral)
        for quad in cross_dihedral.quadruples:
            defaults[tuple(_token_name(token) for token in quad)] = (
                keyword,
                params,
                table_file,
            )
    return defaults


def _cross_dihedral_params(cross_dihedral: CrossDihedral) -> tuple[list[float], str | None]:
    if cross_dihedral.table:
        return [], cross_dihedral.table
    values = _parse_float_tokens(cross_dihedral.params)
    if len(values) >= 6:
        return units.gromacs_rb_to_lammps(values[:6]), None
    return [0.0] * 6, None


def _pair_14_terms(
    system: System,
    molecule: MoleculeType,
    topology: Topology,
) -> None:
    atom_by_index = {atom.index: atom for atom in molecule.atoms}
    seen: set[tuple[int, int]] = set()
    for pair in molecule.cross_pairs:
        atom_i = atom_by_index.get(pair.i)
        atom_j = atom_by_index.get(pair.j)
        if atom_i is None or atom_j is None:
            continue
        sigma, epsilon = resolve_pair_lj_params(topology, atom_i, atom_j, pair.func, pair.params)
        i_id, j_id = pair.i, pair.j
        if i_id > j_id:
            i_id, j_id = j_id, i_id
        key = (i_id, j_id)
        if key in seen:
            continue
        seen.add(key)
        system.cross_pairs.append(
            LammpsCrossPair(i=i_id, j=j_id, sigma=sigma, epsilon=epsilon, keyword="at")
        )
    if system.cross_pairs:
        system.has_cross_pairs = True


def _forcefield_dirs(settings: Settings, base_dir: Path) -> list[Path]:
    dirs: list[Path] = []
    if settings.prep.forcefield_dir:
        ff_dir = Path(settings.prep.forcefield_dir)
        if not ff_dir.is_absolute():
            ff_dir = (base_dir / ff_dir).resolve()
        dirs.append(ff_dir)
    return dirs


def _is_intra_bead_dihedral(
    atoms_by_id: dict[int, LammpsAtom],
    i: int,
    j: int,
    k: int,
    l: int,
) -> bool:
    mol_ids = {atoms_by_id[idx].mol_id for idx in (i, j, k, l) if idx in atoms_by_id}
    return len(mol_ids) == 1


def _add_dihedral_type(
    system: System,
    dihedral_type_map: dict[tuple[str, str, str, tuple[float, ...], str], int],
    *,
    style: str,
    keyword: str,
    params: list[float],
    table_file: str | None = None,
) -> int:
    coeff_key = tuple(round(value, 8) for value in params[:6])
    if style == "charmm" or style == "backmap/charmm" or style == "harmonic":
        coeff_key = tuple(round(value, 8) for value in params[:3])
    key = (style, keyword, table_file or "", coeff_key)
    type_id = dihedral_type_map.get(key)
    if type_id is not None:
        return type_id
    type_id = len(system.dihedral_types) + 1
    system.dihedral_types.append(
        DihedralTypeInfo(
            type_id=type_id,
            style=style,
            keyword=keyword,
            params=params,
            table_file=table_file,
            table_keyword="ENTRY" if table_file else None,
        )
    )
    dihedral_type_map[key] = type_id
    return type_id


def _register_dihedral_table(
    system: System,
    xvg_name: str,
    search_dirs: list[Path],
) -> str | None:
    if _find_xvg(xvg_name, search_dirs) is None:
        return None
    table_out = Path(xvg_name).stem + ".table"
    if (xvg_name, table_out) not in system.dihedral_table_files:
        system.dihedral_table_files.append((xvg_name, table_out))
    return table_out


def _cg_dihedral_table_name(dih: TopDihedral, search_dirs: list[Path]) -> str | None:
    if dih.func != 8 or not dih.params:
        return None
    index = int(dih.params[0])
    candidate = f"table_d{index}.xvg"
    if _find_xvg(candidate, search_dirs):
        return candidate
    return None


def _dihedral_terms(
    system: System,
    molecule: MoleculeType,
    topology: Topology,
    dihedral_defaults: dict[tuple[str, str, str, str], tuple[str, list[float], str | None]],
    cg_type_names: set[str],
    search_dirs: list[Path],
    *,
    plain_cg: bool = False,
) -> list[LammpsDihedral]:
    atom_by_index = {atom.index: atom for atom in molecule.atoms}
    atoms_by_id = {atom.atom_id: atom for atom in system.atoms}
    atom_type_by_index = {atom.index: atom.type for atom in molecule.atoms}
    dihedrals: list[LammpsDihedral] = []
    dihedral_type_map: dict[tuple[str, str, str, tuple[float, ...], str], int] = {}
    all_dihedrals = [*molecule.dihedrals, *molecule.cross_dihedrals]

    for dihedral_index, dih in enumerate(all_dihedrals, start=1):
        atom_i = atom_by_index[dih.i]
        atom_j = atom_by_index[dih.j]
        atom_k = atom_by_index[dih.k]
        atom_l = atom_by_index[dih.atom_l]
        quad = (
            _token_name(atom_i.name),
            _token_name(atom_j.name),
            _token_name(atom_k.name),
            _token_name(atom_l.name),
        )
        is_cg_quad = all(
            atom.type in cg_type_names or _single_letter_cg(atom.type)
            for atom in (atom_i, atom_j, atom_k, atom_l)
        )
        keyword = "cg" if is_cg_quad else "at"
        style = "backmap/ryckaert"
        table_file: str | None = None
        params: list[float] = [0.0] * 6

        func = dih.func
        raw_params = dih.params
        if func == 3 and not raw_params:
            func, raw_params = resolve_dihedral_params(
                dih, atom_type_by_index, topology.dihedraltypes, topology.atom_types
            )

        if func == 8:
            style = "backmap/table"
            params = []
            keyword = "cg"
            xvg_name = _cg_dihedral_table_name(dih, search_dirs)
            if not xvg_name:
                tablenr = int(raw_params[0]) if raw_params else 0
                raise ValueError(f"Missing dihedral table table_d{tablenr}.xvg for func-8 dihedral")
            table_file = _register_dihedral_table(system, xvg_name, search_dirs)
        elif func == 3:
            params = units.gromacs_rb_to_lammps(raw_params)
            if _is_intra_bead_dihedral(atoms_by_id, dih.i, dih.j, dih.k, dih.atom_l):
                style = "ryckaert"
                keyword = ""
            else:
                style = "backmap/ryckaert"
                keyword = "at"
        elif func == 1:
            improper = resolve_opls_improper_params(dih)
            if improper is None:
                raise ValueError(
                    f"Unsupported func-1 dihedral for atoms "
                    f"{dih.i}-{dih.j}-{dih.k}-{dih.atom_l} ({quad})"
                )
            phi0, k_kj, multiplicity = improper
            params = units.gromacs_harmonic_dihedral(k_kj, multiplicity, phi0)
            if _is_intra_bead_dihedral(atoms_by_id, dih.i, dih.j, dih.k, dih.atom_l):
                style = "harmonic"
                keyword = ""
            else:
                style = "harmonic"
                keyword = ""
        elif func == 0:
            default = dihedral_defaults.get(quad)
            if default:
                keyword, params, table_file_src = default
                if table_file_src:
                    style = "backmap/table"
                    table_file = _register_dihedral_table(system, table_file_src, search_dirs)
                    params = []
                else:
                    style = "backmap/ryckaert"
            else:
                params = [0.0] * 6
        else:
            raise ValueError(
                f"Unsupported dihedral func {func} for atoms "
                f"{dih.i}-{dih.j}-{dih.k}-{dih.atom_l} ({quad})"
            )

        export_style = _plain_cg_style(style) if plain_cg else style
        export_keyword = "" if plain_cg else keyword
        type_id = _add_dihedral_type(
            system,
            dihedral_type_map,
            style=export_style,
            keyword=export_keyword,
            params=params,
            table_file=table_file,
        )
        dihedrals.append(
            LammpsDihedral(
                dihedral_id=dihedral_index,
                type_id=type_id,
                i=dih.i,
                j=dih.j,
                k=dih.k,
                l=dih.atom_l,
            )
        )
    return dihedrals


def _bond_terms(
    system: System,
    molecule: MoleculeType,
    bond_defaults: dict[frozenset[str], tuple[str, list[float], str | None]],
    cg_type_names: set[str],
    search_dirs: list[Path],
    *,
    plain_cg: bool = False,
) -> list[LammpsBond]:
    atom_by_index = {atom.index: atom for atom in molecule.atoms}
    bonds: list[LammpsBond] = []
    bond_type_map: dict[tuple[str, str, float, float, str], int] = {}
    all_bonds = [*molecule.bonds, *molecule.cross_bonds]
    for bond_index, bond in enumerate(all_bonds, start=1):
        atom_i = atom_by_index[bond.i]
        atom_j = atom_by_index[bond.j]
        is_cg_bond = (atom_i.type in cg_type_names or _single_letter_cg(atom_i.type)) and (
            atom_j.type in cg_type_names or _single_letter_cg(atom_j.type)
        )
        keyword = "cg" if is_cg_bond else "at"

        table_file: str | None = None
        if bond.func == 8:
            style = "backmap/table"
            params = [0.0, 0.0]
            xvg_name = _cg_bond_table_name(bond, search_dirs)
            if xvg_name:
                table_out = _register_bond_table(system, xvg_name, search_dirs)
                if table_out:
                    table_file = table_out
        elif bond.func == 0:
            style = "backmap/harmonic"
            default = bond_defaults.get(
                frozenset((_token_name(atom_i.name), _token_name(atom_j.name)))
            )
            if default:
                keyword, params, table_file_src = default
                if table_file_src:
                    style = "backmap/table"
                    table_file = _register_bond_table(system, table_file_src, search_dirs)
            else:
                params = [0.0, 0.0]
        elif len(bond.params) >= 2:
            style = "backmap/harmonic"
            params = [units.spring_bond(bond.params[1]), units.distance(bond.params[0])]
        else:
            style = "backmap/harmonic"
            default = bond_defaults.get(
                frozenset((_token_name(atom_i.name), _token_name(atom_j.name)))
            )
            if default:
                keyword, params, table_file_src = default
                if table_file_src:
                    style = "backmap/table"
                    table_file = _register_bond_table(system, table_file_src, search_dirs)
            else:
                params = [0.0, 0.0]

        export_style = _plain_cg_style(style) if plain_cg else style
        export_keyword = "" if plain_cg else keyword
        type_id = _add_bond_type(
            system,
            bond_type_map,
            style=export_style,
            keyword=export_keyword,
            params=params,
            table_file=table_file,
        )

        bonds.append(
            LammpsBond(
                bond_id=bond_index,
                type_id=type_id,
                i=bond.i,
                j=bond.j,
            )
        )
    return bonds


def _angle_terms(
    system: System,
    molecule: MoleculeType,
    angle_defaults: dict[tuple[str, str, str], tuple[str, list[float]]],
    cg_type_names: set[str],
    search_dirs: list[Path],
    *,
    plain_cg: bool = False,
) -> list[LammpsAngle]:
    atom_by_index = {atom.index: atom for atom in molecule.atoms}
    angles: list[LammpsAngle] = []
    angle_type_map: dict[tuple[str, str, str, str, float, float], int] = {}
    all_angles = [*molecule.angles, *molecule.cross_angles]
    for angle_index, angle in enumerate(all_angles, start=1):
        atom_i = atom_by_index[angle.i]
        atom_j = atom_by_index[angle.j]
        atom_k = atom_by_index[angle.k]
        is_cg_angle = all(
            atom.type in cg_type_names or _single_letter_cg(atom.type)
            for atom in (atom_i, atom_j, atom_k)
        )
        keyword = "cg" if is_cg_angle else "at"
        triple = (
            _token_name(atom_i.name),
            _token_name(atom_j.name),
            _token_name(atom_k.name),
        )
        style = "backmap/harmonic"
        table_file: str | None = None
        if angle.func == 8:
            style = "backmap/table"
            params: list[float] = []
            keyword = "cg"
            xvg_name = _cg_angle_table_name(angle, search_dirs)
            if not xvg_name:
                tablenr = int(angle.params[0]) if angle.params else 0
                raise ValueError(f"Missing angle table table_a{tablenr}.xvg for func-8 angle")
            table_file = _register_angle_table(system, xvg_name, search_dirs)
        elif angle.func == 0:
            default = angle_defaults.get(triple)
            params = default[1] if default else [0.0, 0.0]
            if default:
                keyword = default[0]
        elif len(angle.params) >= 2:
            params = [units.spring_angle(angle.params[1]), angle.params[0]]
        else:
            default = angle_defaults.get(triple)
            if default:
                keyword, params = default
            else:
                params = [0.0, 0.0]
        export_style = _plain_cg_style(style) if plain_cg else style
        export_keyword = "" if plain_cg else keyword
        type_id = _add_angle_type(
            system,
            angle_type_map,
            style=export_style,
            keyword=export_keyword,
            params=params,
            table_file=table_file,
        )
        angles.append(
            LammpsAngle(
                angle_id=angle_index,
                type_id=type_id,
                i=angle.i,
                j=angle.j,
                k=angle.k,
            )
        )
    return angles


def _find_xvg(name: str, search_dirs: list[Path]) -> str | None:
    for directory in search_dirs:
        if (directory / name).is_file():
            return name
    return None


def _register_bond_table(
    system: System,
    xvg_name: str,
    search_dirs: list[Path],
) -> str | None:
    if _find_xvg(xvg_name, search_dirs) is None:
        return None
    table_out = Path(xvg_name).stem + ".table"
    if (xvg_name, table_out) not in system.table_files:
        system.table_files.append((xvg_name, table_out))
    return table_out


def _cg_bond_table_name(bond: TopBond, search_dirs: list[Path]) -> str | None:
    if bond.func != 8 or not bond.params:
        return None
    index = int(bond.params[0])
    for candidate in (f"table_b{index}.xvg", f"table_a{index}.xvg"):
        if _find_xvg(candidate, search_dirs):
            return candidate
    return None


def _cg_angle_table_name(angle: TopAngle, search_dirs: list[Path]) -> str | None:
    if angle.func != 8 or not angle.params:
        return None
    index = int(angle.params[0])
    candidate = f"table_a{index}.xvg"
    if _find_xvg(candidate, search_dirs):
        return candidate
    return None


def _register_angle_table(
    system: System,
    xvg_name: str,
    search_dirs: list[Path],
) -> str | None:
    if _find_xvg(xvg_name, search_dirs) is None:
        return None
    table_out = Path(xvg_name).stem + ".table"
    if (xvg_name, table_out) not in system.angle_table_files:
        system.angle_table_files.append((xvg_name, table_out))
    return table_out


def _resolve_pair_tables(
    system: System,
    settings: Settings,
    search_dirs: list[Path],
) -> None:
    table_groups = set(settings.simulation.table_groups)
    if not table_groups:
        return
    for pair_type in system.pair_types:
        if pair_type.kind != "cg":
            continue
        type_i = system.atom_types[pair_type.itype - 1]
        type_j = system.atom_types[pair_type.jtype - 1]
        if type_i.name not in table_groups or type_j.name not in table_groups:
            continue
        name_a, name_b = sorted([type_i.name, type_j.name])
        for xvg_name in (f"table_{name_a}_{name_b}.xvg", f"table_{name_b}_{name_a}.xvg"):
            if _find_xvg(xvg_name, search_dirs) is None:
                continue
            table_out = Path(xvg_name).stem + ".table"
            pair_type.table_file = table_out
            pair_type.table_keyword = "ENTRY"
            if (xvg_name, table_out) not in system.pair_table_files:
                system.pair_table_files.append((xvg_name, table_out))
            break


def _pair_terms(atom_types: list[AtomTypeInfo]) -> list[PairTypeInfo]:
    pair_types: list[PairTypeInfo] = []
    for i, atom_type_i in enumerate(atom_types, start=1):
        for j, atom_type_j in enumerate(atom_types[i - 1 :], start=i):
            if atom_type_i.is_cg and atom_type_j.is_cg:
                pair_types.append(PairTypeInfo(itype=i, jtype=j, kind="cg"))
            elif not atom_type_i.is_cg and not atom_type_j.is_cg:
                sigma = (
                    0.5 * (atom_type_i.sigma + atom_type_j.sigma)
                    if atom_type_i.sigma > 0 and atom_type_j.sigma > 0
                    else 0.0
                )
                epsilon = (
                    (atom_type_i.epsilon * atom_type_j.epsilon) ** 0.5
                    if atom_type_i.epsilon > 0 and atom_type_j.epsilon > 0
                    else 0.0
                )
                pair_types.append(
                    PairTypeInfo(
                        itype=i,
                        jtype=j,
                        kind="atomistic",
                        sigma=sigma,
                        epsilon=epsilon,
                    )
                )
            else:
                pair_types.append(PairTypeInfo(itype=i, jtype=j, kind="none"))
    return pair_types


def build_system_from_cg(
    settings: Settings,
    settings_path: Path,
    table_search_dirs: list[Path] | None = None,
) -> System:
    """Build a pure CG LAMMPS system from ``cg_system`` GRO/TOP (network prep)."""
    if settings.cg_system is None:
        raise ValueError("cg_system is required for CG-only network build")

    from backmap_prep.schema import resolve_data_dir

    base_dir = resolve_data_dir(settings_path, settings)
    gro_path = (base_dir / settings.cg_system.coordinates).resolve()
    top_path = (base_dir / settings.cg_system.topology).resolve()

    gro_file = parse_gro(gro_path)
    top_file = parse_top(top_path, include_dirs=[base_dir])
    molecule = _topology_molecule(top_file)
    cg_type_names = {atom.type for atom in molecule.atoms}

    system = System()
    system.box = (
        units.distance(gro_file.box[0]),
        units.distance(gro_file.box[1]),
        units.distance(gro_file.box[2]),
    )

    atom_types, type_ids = _atom_type_info(top_file, molecule, cg_type_names)
    for atom_type in atom_types:
        atom_type.is_cg = True
    system.atom_types = atom_types
    system.atoms = _cg_atoms(molecule, gro_file.atoms, type_ids)
    system.cg_type_id = atom_types[0].type_id if atom_types else 0

    search_dirs = [base_dir, *(table_search_dirs or [])]
    system.bonds = _bond_terms(
        system,
        molecule,
        {},
        cg_type_names,
        search_dirs,
        plain_cg=True,
    )
    system.angles = _angle_terms(
        system,
        molecule,
        {},
        cg_type_names,
        search_dirs,
        plain_cg=True,
    )
    system.dihedrals = _dihedral_terms(
        system,
        molecule,
        top_file,
        {},
        cg_type_names,
        search_dirs,
        plain_cg=True,
    )
    system.pair_types = _pair_terms(system.atom_types)
    _resolve_pair_tables(system, settings, search_dirs)
    prepare_network_coordinates(system)
    cg_cut = units.distance(settings.simulation.cg_cutoff)
    lj_cut = units.distance(settings.simulation.lj_cutoff)
    validate_bond_geometry(system, max(lj_cut, cg_cut))
    return system


def build_system_from_hybrid(
    settings: Settings,
    base_dir: Path,
    gro_path: Path,
    top_path: Path,
    table_search_dirs: list[Path] | None = None,
    forcefield_dirs: list[Path] | None = None,
) -> System:
    """Build a LAMMPS `System` from network-engine hybrid GRO/TOP outputs."""
    ff_dirs = (
        forcefield_dirs if forcefield_dirs is not None else _forcefield_dirs(settings, base_dir)
    )
    gro_file = parse_gro(
        (base_dir / gro_path).resolve() if not gro_path.is_absolute() else gro_path
    )
    top_file = parse_top(
        (base_dir / top_path).resolve() if not top_path.is_absolute() else top_path,
        include_dirs=[base_dir],
        forcefield_dirs=ff_dirs,
    )
    molecule = _topology_molecule(top_file)
    cg_type_names = _cg_types_from_settings(settings)

    system = System()
    system.box = (
        units.distance(gro_file.box[0]),
        units.distance(gro_file.box[1]),
        units.distance(gro_file.box[2]),
    )

    atom_types, type_ids = _atom_type_info(top_file, molecule, cg_type_names)
    system.atom_types = atom_types
    system.atoms = _hybrid_atoms(molecule, gro_file.atoms, type_ids, cg_type_names)
    system.cg_type_id = next((atom_type.type_id for atom_type in atom_types if atom_type.is_cg), 0)

    bond_defaults = _cross_bond_defaults(settings)
    angle_defaults = _cross_angle_defaults(settings)
    dihedral_defaults = _cross_dihedral_defaults(settings)
    search_dirs = [base_dir, *(table_search_dirs or [])]
    system.bonds = _bond_terms(system, molecule, bond_defaults, cg_type_names, search_dirs)
    system.angles = _angle_terms(system, molecule, angle_defaults, cg_type_names, search_dirs)
    system.dihedrals = _dihedral_terms(
        system, molecule, top_file, dihedral_defaults, cg_type_names, search_dirs
    )
    _pair_14_terms(system, molecule, top_file)
    system.pair_types = _pair_terms(system.atom_types)
    _resolve_pair_tables(system, settings, search_dirs)
    system.has_cross_bonds = any(bond_type.keyword == "at" for bond_type in system.bond_types)
    system.has_cross_angles = any(angle_type.keyword == "at" for angle_type in system.angle_types)
    system.has_cross_dihedrals = any(
        dihedral_type.keyword == "at" for dihedral_type in system.dihedral_types
    )
    prepare_network_coordinates(system)
    lj_cut = units.distance(settings.simulation.lj_cutoff)
    cg_cut = units.distance(settings.simulation.cg_cutoff)
    validate_bond_geometry(system, max(lj_cut, cg_cut))
    return system
