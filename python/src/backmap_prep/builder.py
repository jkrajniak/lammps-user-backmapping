"""Build the hybrid CG+AT system from settings and source files."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

from . import units
from .parsers import parse_gro, parse_top

if TYPE_CHECKING:
    from .schema import Settings


@dataclass
class LammpsAtom:
    atom_id: int
    mol_id: int
    type_id: int
    charge: float
    x: float  # Angstrom
    y: float  # Angstrom
    z: float  # Angstrom
    type_name: str = ""
    is_cg: bool = False


@dataclass
class LammpsBond:
    bond_id: int
    type_id: int
    i: int
    j: int


@dataclass
class LammpsAngle:
    angle_id: int
    type_id: int
    i: int
    j: int
    k: int


@dataclass
class BondTypeInfo:
    type_id: int
    style: str  # "harmonic", "backmap/harmonic", "backmap/table"
    keyword: str  # "at" or "cg" for backmap/* styles, "" for static
    params: list[float]
    table_file: str | None = None
    table_keyword: str | None = None


@dataclass
class AngleTypeInfo:
    type_id: int
    style: str
    keyword: str
    params: list[float]


@dataclass
class AtomTypeInfo:
    type_id: int
    name: str
    mass: float
    is_cg: bool
    sigma: float = 0.0  # Angstrom
    epsilon: float = 0.0  # kcal/mol


@dataclass
class PairTypeInfo:
    """Which pair kind for each type pair: 'atomistic', 'cg', 'none'."""

    itype: int
    jtype: int
    kind: str
    sigma: float = 0.0
    epsilon: float = 0.0
    table_file: str | None = None
    table_keyword: str | None = None


@dataclass
class System:
    """Complete LAMMPS system representation."""

    atoms: list[LammpsAtom] = field(default_factory=list)
    bonds: list[LammpsBond] = field(default_factory=list)
    angles: list[LammpsAngle] = field(default_factory=list)
    atom_types: list[AtomTypeInfo] = field(default_factory=list)
    bond_types: list[BondTypeInfo] = field(default_factory=list)
    angle_types: list[AngleTypeInfo] = field(default_factory=list)
    pair_types: list[PairTypeInfo] = field(default_factory=list)
    box: tuple[float, float, float] = (0.0, 0.0, 0.0)  # Angstrom
    cg_type_id: int = 0
    has_cross_bonds: bool = False
    has_cross_angles: bool = False

    # Table files to convert
    table_files: list[tuple[str, str]] = field(default_factory=list)  # (src, dst) pairs


def build_system(settings: Settings, base_dir: Path) -> System:
    """Build the complete hybrid system from settings and source files."""
    sys = System()

    mol_def = settings.molecules[0]
    ident = mol_def.ident or mol_def.name

    assert isinstance(mol_def.source.coordinates, str)
    assert isinstance(mol_def.source.topology, str)
    at_gro = parse_gro(base_dir / mol_def.source.coordinates)
    at_top = parse_top(base_dir / mol_def.source.topology, include_dirs=[base_dir])

    cg_gro = parse_gro(base_dir / settings.cg_system.coordinates)
    cg_top = parse_top(base_dir / settings.cg_system.topology, include_dirs=[base_dir])

    sys.box = (
        units.distance(cg_gro.box[0]),
        units.distance(cg_gro.box[1]),
        units.distance(cg_gro.box[2]),
    )

    # Determine molecule counts from CG topology
    n_molecules = 0
    cg_mol_name = None
    for name, count in cg_top.molecules:
        if name in cg_top.molecule_types:
            cg_mol_name = name
            n_molecules = count
            break

    if n_molecules == 0 or cg_mol_name is None:
        raise ValueError("No molecules found in CG topology")

    # Build atom type map
    type_map: dict[str, int] = {}
    type_id_counter = 0

    # CG atom types first
    cg_mol = cg_top.molecule_types.get(cg_mol_name)
    if not cg_mol:
        raise ValueError(f"Molecule type '{cg_mol_name}' not found in CG topology")

    for atom in cg_mol.atoms:
        if atom.type not in type_map:
            type_id_counter += 1
            type_map[atom.type] = type_id_counter
            mass = atom.mass
            if mass == 0.0 and atom.type in cg_top.atom_types:
                mass = cg_top.atom_types[atom.type].mass

            sys.atom_types.append(
                AtomTypeInfo(
                    type_id=type_id_counter,
                    name=atom.type,
                    mass=mass,
                    is_cg=True,
                )
            )

    # AT atom types
    at_mol_name = ident
    at_mol = None
    for name, mol in at_top.molecule_types.items():
        if name == at_mol_name or name == mol_def.name:
            at_mol = mol
            at_mol_name = name
            break
    if at_mol is None:
        at_mol = next(iter(at_top.molecule_types.values()), None)
    if at_mol is None:
        raise ValueError("No molecule type found in AT topology")

    for atom in at_mol.atoms:
        if atom.type not in type_map:
            type_id_counter += 1
            type_map[atom.type] = type_id_counter
            mass = atom.mass
            if mass == 0.0 and atom.type in at_top.atom_types:
                mass = at_top.atom_types[atom.type].mass

            sig = 0.0
            eps = 0.0
            if atom.type in at_top.atom_types:
                sig = units.sigma(at_top.atom_types[atom.type].sigma)
                eps = units.epsilon(at_top.atom_types[atom.type].epsilon)

            sys.atom_types.append(
                AtomTypeInfo(
                    type_id=type_id_counter,
                    name=atom.type,
                    mass=mass,
                    is_cg=False,
                    sigma=sig,
                    epsilon=eps,
                )
            )

    # Track which CG type id is used for fix backmap
    for ati in sys.atom_types:
        if ati.is_cg:
            sys.cg_type_id = ati.type_id
            break

    # Build bead → atom-name mapping from YAML
    bead_atoms: dict[str, list[str]] = {}
    for bead in mol_def.beads:
        atom_names = []
        for ref in bead.atoms:
            parts = ref.split(":")
            atom_names.append(parts[-1])
        bead_atoms[bead.name] = atom_names

    # Build atom name → TopAtom index mapping for the AT molecule
    at_name_to_idx: dict[str, int] = {}
    for atom in at_mol.atoms:
        at_name_to_idx[atom.name] = atom.index

    # Build intra-CG bond set (bonds between atoms within the same bead)
    atoms_in_bead: dict[str, set[str]] = {}
    for bead_name, atom_names in bead_atoms.items():
        atoms_in_bead[bead_name] = set(atom_names)

    def is_intra_cg_bond(name_i: str, name_j: str) -> bool:
        for _bead_name, names in atoms_in_bead.items():
            if name_i in names and name_j in names:
                return True
        return False

    # Build bond types: static intra-CG bonds from AT topology
    bond_type_counter = 0

    # Collect unique AT bond parameter sets for intra-CG
    intra_bond_params: dict[tuple[str, float, float], int] = {}
    for bond in at_mol.bonds:
        ai = at_mol.atoms[bond.i - 1]
        aj = at_mol.atoms[bond.j - 1]
        if not is_intra_cg_bond(ai.name, aj.name):
            continue
        if bond.func == 1 and len(bond.params) >= 2:
            k_lammps = units.spring_bond(bond.params[1])
            r0_lammps = units.distance(bond.params[0])
            key = ("harmonic", round(k_lammps, 6), round(r0_lammps, 6))
            if key not in intra_bond_params:
                bond_type_counter += 1
                intra_bond_params[key] = bond_type_counter
                sys.bond_types.append(
                    BondTypeInfo(
                        type_id=bond_type_counter,
                        style="harmonic",
                        keyword="",
                        params=[k_lammps, r0_lammps],
                    )
                )

    # Cross-CG bond types from settings
    cross_bond_type_map: dict[str, int] = {}
    for cb in settings.cross_interactions.bonds:
        bond_type_counter += 1
        sys.has_cross_bonds = True
        params_tokens = cb.params.split()
        func_type = int(params_tokens[0])

        if cb.cg_bonded:
            if cb.table:
                style = "backmap/table"
                kw = "cg"
                table_out = Path(cb.table).stem + ".table"
                sys.bond_types.append(
                    BondTypeInfo(
                        type_id=bond_type_counter,
                        style=style,
                        keyword=kw,
                        params=[],
                        table_file=table_out,
                        table_keyword="ENTRY",
                    )
                )
                sys.table_files.append((cb.table, table_out))
            else:
                style = "backmap/harmonic"
                kw = "cg"
                k_lammps = (
                    units.spring_bond(float(params_tokens[2])) if len(params_tokens) > 2 else 0.0
                )
                r0_lammps = (
                    units.distance(float(params_tokens[1])) if len(params_tokens) > 1 else 0.0
                )
                sys.bond_types.append(
                    BondTypeInfo(
                        type_id=bond_type_counter,
                        style=style,
                        keyword=kw,
                        params=[k_lammps, r0_lammps],
                    )
                )
        else:
            if func_type == 1 and len(params_tokens) >= 3:
                style = "backmap/harmonic"
                kw = "at"
                r0_lammps = units.distance(float(params_tokens[1]))
                k_lammps = units.spring_bond(float(params_tokens[2]))
                sys.bond_types.append(
                    BondTypeInfo(
                        type_id=bond_type_counter,
                        style=style,
                        keyword=kw,
                        params=[k_lammps, r0_lammps],
                    )
                )
            elif func_type == 8 and cb.table:
                style = "backmap/table"
                kw = "at"
                table_out = Path(cb.table).stem + ".table"
                sys.bond_types.append(
                    BondTypeInfo(
                        type_id=bond_type_counter,
                        style=style,
                        keyword=kw,
                        params=[],
                        table_file=table_out,
                        table_keyword="ENTRY",
                    )
                )
                sys.table_files.append((cb.table, table_out))
            else:
                style = "backmap/harmonic"
                kw = "at"
                sys.bond_types.append(
                    BondTypeInfo(
                        type_id=bond_type_counter,
                        style=style,
                        keyword=kw,
                        params=[0.0, 0.0],
                    )
                )

        # Store the cross bond type index for pair lookup
        for pair in cb.pairs:
            key_str = f"{pair[0]}_{pair[1]}"
            cross_bond_type_map[key_str] = bond_type_counter

    # Cross-CG angle types from settings
    angle_type_counter = 0
    intra_angle_params: dict[tuple[str, float, float], int] = {}
    for angle in at_mol.angles:
        ai = at_mol.atoms[angle.i - 1]
        aj = at_mol.atoms[angle.j - 1]
        ak = at_mol.atoms[angle.k - 1]
        all_in_same = False
        for _bead_name, names in atoms_in_bead.items():
            if ai.name in names and aj.name in names and ak.name in names:
                all_in_same = True
                break
        if not all_in_same:
            continue
        if angle.func == 1 and len(angle.params) >= 2:
            k_lammps = units.spring_angle(angle.params[1])
            theta0 = angle.params[0]
            key = ("harmonic", round(k_lammps, 6), round(theta0, 4))
            if key not in intra_angle_params:
                angle_type_counter += 1
                intra_angle_params[key] = angle_type_counter
                sys.angle_types.append(
                    AngleTypeInfo(
                        type_id=angle_type_counter,
                        style="harmonic",
                        keyword="",
                        params=[k_lammps, theta0],
                    )
                )

    for ca in settings.cross_interactions.angles:
        angle_type_counter += 1
        sys.has_cross_angles = True
        params_tokens = ca.params.split()

        kw = "cg" if ca.cg_bonded else "at"

        if len(params_tokens) >= 3:
            theta0 = float(params_tokens[1])
            k_lammps = units.spring_angle(float(params_tokens[2]))
        else:
            theta0 = 0.0
            k_lammps = 0.0

        sys.angle_types.append(
            AngleTypeInfo(
                type_id=angle_type_counter,
                style="backmap/harmonic",
                keyword=kw,
                params=[k_lammps, theta0],
            )
        )

    # Build pair type classification
    n_types = type_id_counter
    for i in range(1, n_types + 1):
        for j in range(i, n_types + 1):
            ti = sys.atom_types[i - 1]
            tj = sys.atom_types[j - 1]
            if ti.is_cg and tj.is_cg:
                sys.pair_types.append(
                    PairTypeInfo(
                        itype=i,
                        jtype=j,
                        kind="cg",
                    )
                )
            elif not ti.is_cg and not tj.is_cg:
                sig_ij = 0.5 * (ti.sigma + tj.sigma) if ti.sigma > 0 and tj.sigma > 0 else 0.0
                eps_ij = (
                    (ti.epsilon * tj.epsilon) ** 0.5 if ti.epsilon > 0 and tj.epsilon > 0 else 0.0
                )
                sys.pair_types.append(
                    PairTypeInfo(
                        itype=i,
                        jtype=j,
                        kind="atomistic",
                        sigma=sig_ij,
                        epsilon=eps_ij,
                    )
                )
            else:
                sys.pair_types.append(
                    PairTypeInfo(
                        itype=i,
                        jtype=j,
                        kind="none",
                    )
                )

    # Build atoms and bonds for each molecule instance
    atom_id = 0
    bond_id = 0
    angle_id = 0
    cg_atom_count = len(cg_mol.atoms)

    for mol_idx in range(n_molecules):
        mol_id = mol_idx + 1
        mol_atom_offset = atom_id

        # CG atom positions from CG .gro
        cg_start = mol_idx * cg_atom_count
        cg_local_to_global: dict[int, int] = {}

        for ci, cg_atom in enumerate(cg_mol.atoms):
            atom_id += 1
            gro_idx = cg_start + ci
            if gro_idx < len(cg_gro.atoms):
                ga = cg_gro.atoms[gro_idx]
                x = units.distance(ga.x)
                y = units.distance(ga.y)
                z = units.distance(ga.z)
            else:
                x = y = z = 0.0

            sys.atoms.append(
                LammpsAtom(
                    atom_id=atom_id,
                    mol_id=mol_id,
                    type_id=type_map[cg_atom.type],
                    charge=cg_atom.charge,
                    x=x,
                    y=y,
                    z=z,
                    type_name=cg_atom.type,
                    is_cg=True,
                )
            )
            cg_local_to_global[cg_atom.index] = atom_id

        # AT atoms: place relative to their CG bead's position
        at_local_to_global: dict[int, int] = {}
        for bead_idx, bead in enumerate(mol_def.beads):
            cg_atom_gro_idx = cg_start + bead_idx
            if cg_atom_gro_idx < len(cg_gro.atoms):
                cx = units.distance(cg_gro.atoms[cg_atom_gro_idx].x)
                cy = units.distance(cg_gro.atoms[cg_atom_gro_idx].y)
                cz = units.distance(cg_gro.atoms[cg_atom_gro_idx].z)
            else:
                cx = cy = cz = 0.0

            for atom_ref in bead.atoms:
                parts = atom_ref.split(":")
                atom_name = parts[-1]
                at_idx = at_name_to_idx.get(atom_name)
                if at_idx is None:
                    continue
                at_atom = at_mol.atoms[at_idx - 1]

                # Position from AT .gro (single molecule template)
                at_gro_atom = None
                for ga in at_gro.atoms:
                    if ga.name == atom_name:
                        at_gro_atom = ga
                        break

                if at_gro_atom:
                    dx = units.distance(at_gro_atom.x) - units.distance(at_gro.atoms[0].x)
                    dy = units.distance(at_gro_atom.y) - units.distance(at_gro.atoms[0].y)
                    dz = units.distance(at_gro_atom.z) - units.distance(at_gro.atoms[0].z)
                    x = cx + dx
                    y = cy + dy
                    z = cz + dz
                else:
                    x, y, z = cx, cy, cz

                atom_id += 1
                sys.atoms.append(
                    LammpsAtom(
                        atom_id=atom_id,
                        mol_id=mol_id,
                        type_id=type_map[at_atom.type],
                        charge=at_atom.charge,
                        x=x,
                        y=y,
                        z=z,
                        type_name=at_atom.type,
                        is_cg=False,
                    )
                )
                at_local_to_global[at_idx] = atom_id

        # Intra-CG bonds from AT topology
        for bond in at_mol.bonds:
            ai = at_mol.atoms[bond.i - 1]
            aj = at_mol.atoms[bond.j - 1]
            if not is_intra_cg_bond(ai.name, aj.name):
                continue
            if bond.func == 1 and len(bond.params) >= 2:
                k_lammps = units.spring_bond(bond.params[1])
                r0_lammps = units.distance(bond.params[0])
                key = ("harmonic", round(k_lammps, 6), round(r0_lammps, 6))
                bt_id = intra_bond_params.get(key)
                if bt_id is None:
                    continue
                gi = at_local_to_global.get(bond.i)
                gj = at_local_to_global.get(bond.j)
                if gi and gj:
                    bond_id += 1
                    sys.bonds.append(
                        LammpsBond(
                            bond_id=bond_id,
                            type_id=bt_id,
                            i=gi,
                            j=gj,
                        )
                    )

        # Cross-CG bonds from settings
        for cb_idx, cb in enumerate(settings.cross_interactions.bonds):
            bt_offset = len(intra_bond_params)
            bt_id = bt_offset + cb_idx + 1
            for pair in cb.pairs:
                name_i = pair[0].split(":")[-1]
                name_j = pair[1].split(":")[-1]

                if cb.cg_bonded:
                    gi = None
                    gj = None
                    for _ci, cg_atom in enumerate(cg_mol.atoms):
                        if cg_atom.name == name_i:
                            gi = cg_local_to_global.get(cg_atom.index)
                        if cg_atom.name == name_j:
                            gj = cg_local_to_global.get(cg_atom.index)
                    # Also try bead names
                    if gi is None:
                        for bi, bead in enumerate(mol_def.beads):
                            if bead.name == name_i:
                                gi = mol_atom_offset + bi + 1
                            if bead.name == name_j:
                                gj = mol_atom_offset + bi + 1
                else:
                    at_idx_i = at_name_to_idx.get(name_i)
                    at_idx_j = at_name_to_idx.get(name_j)
                    gi = at_local_to_global.get(at_idx_i) if at_idx_i is not None else None
                    gj = at_local_to_global.get(at_idx_j) if at_idx_j is not None else None

                if gi and gj:
                    bond_id += 1
                    sys.bonds.append(
                        LammpsBond(
                            bond_id=bond_id,
                            type_id=bt_id,
                            i=gi,
                            j=gj,
                        )
                    )

        # Intra-CG angles from AT topology
        for angle in at_mol.angles:
            ai = at_mol.atoms[angle.i - 1]
            aj = at_mol.atoms[angle.j - 1]
            ak = at_mol.atoms[angle.k - 1]
            all_in_same = False
            for _bead_name, names in atoms_in_bead.items():
                if ai.name in names and aj.name in names and ak.name in names:
                    all_in_same = True
                    break
            if not all_in_same:
                continue
            if angle.func == 1 and len(angle.params) >= 2:
                k_lammps = units.spring_angle(angle.params[1])
                theta0 = angle.params[0]
                key = ("harmonic", round(k_lammps, 6), round(theta0, 4))
                at_id = intra_angle_params.get(key)
                if at_id is None:
                    continue
                gi = at_local_to_global.get(angle.i)
                gj = at_local_to_global.get(angle.j)
                gk = at_local_to_global.get(angle.k)
                if gi and gj and gk:
                    angle_id += 1
                    sys.angles.append(
                        LammpsAngle(
                            angle_id=angle_id,
                            type_id=at_id,
                            i=gi,
                            j=gj,
                            k=gk,
                        )
                    )

        # Cross-CG angles from settings
        for ca_idx, ca in enumerate(settings.cross_interactions.angles):
            at_offset = len(intra_angle_params)
            at_id = at_offset + ca_idx + 1
            for triple in ca.triples:
                triple_names = [t.split(":")[-1] for t in triple]
                idx_i = at_name_to_idx.get(triple_names[0])
                idx_j = at_name_to_idx.get(triple_names[1])
                idx_k = at_name_to_idx.get(triple_names[2])
                if idx_i is None or idx_j is None or idx_k is None:
                    continue
                gi = at_local_to_global.get(idx_i)
                gj = at_local_to_global.get(idx_j)
                gk = at_local_to_global.get(idx_k)
                if gi and gj and gk:
                    angle_id += 1
                    sys.angles.append(
                        LammpsAngle(
                            angle_id=angle_id,
                            type_id=at_id,
                            i=gi,
                            j=gj,
                            k=gk,
                        )
                    )

    return sys
