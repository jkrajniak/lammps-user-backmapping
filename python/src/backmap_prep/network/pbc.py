"""Periodic boundary helpers for network hybrid LAMMPS export."""

from __future__ import annotations

import math
from collections import defaultdict, deque
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from backmap_prep.builder import LammpsAtom, LammpsBond, System


def _place_near_reference(
    ref: LammpsAtom, mobile: LammpsAtom, box: tuple[float, float, float]
) -> None:
    """Move *mobile* to the periodic image nearest *ref* (in-place)."""
    bx, by, bz = box
    for attr, length in (("x", bx), ("y", by), ("z", bz)):
        if length <= 0:
            continue
        ref_val = getattr(ref, attr)
        mobile_val = getattr(mobile, attr)
        delta = mobile_val - ref_val
        if delta > length / 2:
            setattr(mobile, attr, mobile_val - length)
        elif delta < -length / 2:
            setattr(mobile, attr, mobile_val + length)


def _bond_adjacency(bonds: list[LammpsBond]) -> dict[int, list[int]]:
    adjacency: dict[int, list[int]] = defaultdict(list)
    for bond in bonds:
        adjacency[bond.i].append(bond.j)
        adjacency[bond.j].append(bond.i)
    return adjacency


def _atoms_by_mol_id(atoms: list[LammpsAtom]) -> dict[int, list[LammpsAtom]]:
    grouped: dict[int, list[LammpsAtom]] = defaultdict(list)
    for atom in atoms:
        grouped[atom.mol_id].append(atom)
    return grouped


def _mol_id_by_atom_id(atoms: list[LammpsAtom]) -> dict[int, int]:
    return {atom.atom_id: atom.mol_id for atom in atoms}


def _intra_mol_bonds(
    bonds: list[LammpsBond],
    mol_id_by_atom: dict[int, int],
) -> list[LammpsBond]:
    return [
        bond
        for bond in bonds
        if mol_id_by_atom.get(bond.i) == mol_id_by_atom.get(bond.j)
        and mol_id_by_atom.get(bond.i) is not None
    ]


def _unwrap_component(
    root_id: int,
    component: set[int],
    by_id: dict[int, LammpsAtom],
    adjacency: dict[int, list[int]],
    box: tuple[float, float, float],
) -> None:
    component_queue: deque[int] = deque([root_id])
    component_seen = {root_id}
    while component_queue:
        current_id = component_queue.popleft()
        current = by_id[current_id]
        for neighbor_id in adjacency[current_id]:
            if neighbor_id in component and neighbor_id not in component_seen:
                component_seen.add(neighbor_id)
                _place_near_reference(current, by_id[neighbor_id], box)
                component_queue.append(neighbor_id)


def unwrap_bond_graph(
    atoms: list[LammpsAtom], bonds: list[LammpsBond], box: tuple[float, float, float]
) -> None:
    """Unwrap coordinates so bonded neighbors are close in Cartesian space."""
    if not atoms or not bonds:
        return
    bx, by, bz = box
    if min(bx, by, bz) <= 0:
        return

    by_id = {atom.atom_id: atom for atom in atoms}
    adjacency = _bond_adjacency(bonds)
    visited: set[int] = set()

    for root_id in sorted(by_id):
        if root_id in visited:
            continue
        component: set[int] = set()
        queue: deque[int] = deque([root_id])
        while queue:
            current_id = queue.popleft()
            if current_id in component:
                continue
            component.add(current_id)
            visited.add(current_id)
            for neighbor_id in adjacency[current_id]:
                if neighbor_id not in component:
                    queue.append(neighbor_id)

        component_root = min(component)
        _unwrap_component(component_root, component, by_id, adjacency, box)

    for _ in range(16):
        for bond in bonds:
            _place_near_reference(by_id[bond.i], by_id[bond.j], box)


def _unwrap_molecule(
    mol_atoms: list[LammpsAtom],
    mol_bonds: list[LammpsBond],
    box: tuple[float, float, float],
) -> None:
    if not mol_atoms or not mol_bonds:
        return
    by_id = {atom.atom_id: atom for atom in mol_atoms}
    adjacency = _bond_adjacency(mol_bonds)
    cg_roots = sorted(atom.atom_id for atom in mol_atoms if atom.is_cg)
    root_id = cg_roots[0] if cg_roots else min(by_id)
    component = set(by_id)
    _unwrap_component(root_id, component, by_id, adjacency, box)
    for _ in range(8):
        for bond in mol_bonds:
            _place_near_reference(by_id[bond.i], by_id[bond.j], box)


def _translate_molecule(
    atoms: list[LammpsAtom],
    mol_id: int,
    shift_x: float,
    shift_y: float,
    shift_z: float,
) -> None:
    for atom in atoms:
        if atom.mol_id != mol_id:
            continue
        atom.x += shift_x
        atom.y += shift_y
        atom.z += shift_z


def _shift_molecule_near_reference(
    atoms: list[LammpsAtom],
    mol_id: int,
    ref: LammpsAtom,
    mobile: LammpsAtom,
    box: tuple[float, float, float],
) -> None:
    """Translate every atom in *mol_id* so *mobile* sits near *ref*."""
    before = (mobile.x, mobile.y, mobile.z)
    _place_near_reference(ref, mobile, box)
    shift_x = mobile.x - before[0]
    shift_y = mobile.y - before[1]
    shift_z = mobile.z - before[2]
    if shift_x == 0.0 and shift_y == 0.0 and shift_z == 0.0:
        return
    mobile.x, mobile.y, mobile.z = before
    _translate_molecule(atoms, mol_id, shift_x, shift_y, shift_z)


def _inter_mol_bonds(
    bonds: list[LammpsBond],
    mol_id_by_atom: dict[int, int],
) -> list[LammpsBond]:
    return [
        bond
        for bond in bonds
        if mol_id_by_atom.get(bond.i) != mol_id_by_atom.get(bond.j)
        and mol_id_by_atom.get(bond.i) is not None
    ]


def _place_mols_by_spanning_tree(
    atoms: list[LammpsAtom],
    bonds: list[LammpsBond],
    box: tuple[float, float, float],
) -> None:
    """Place molecule fragments using a spanning tree over inter-mol bonds."""
    by_id = {atom.atom_id: atom for atom in atoms}
    mol_id_by_atom = _mol_id_by_atom_id(atoms)
    inter_bonds = _inter_mol_bonds(bonds, mol_id_by_atom)
    if not inter_bonds:
        return

    mol_adjacency: dict[int, list[tuple[int, LammpsBond]]] = defaultdict(list)
    for bond in inter_bonds:
        mol_i = mol_id_by_atom[bond.i]
        mol_j = mol_id_by_atom[bond.j]
        mol_adjacency[mol_i].append((mol_j, bond))
        mol_adjacency[mol_j].append((mol_i, bond))

    root_mol = min(mol_adjacency)
    placed = {root_mol}
    queue: deque[int] = deque([root_mol])
    while queue:
        current_mol = queue.popleft()
        for neighbor_mol, bond in mol_adjacency[current_mol]:
            if neighbor_mol in placed:
                continue
            placed.add(neighbor_mol)
            queue.append(neighbor_mol)
            atom_i = by_id[bond.i]
            atom_j = by_id[bond.j]
            if mol_id_by_atom[bond.i] == neighbor_mol:
                _shift_molecule_near_reference(atoms, neighbor_mol, atom_j, atom_i, box)
            else:
                _shift_molecule_near_reference(atoms, neighbor_mol, atom_i, atom_j, box)


def _min_image_vector_to(
    ref: tuple[float, float, float],
    mobile: tuple[float, float, float],
    box: tuple[float, float, float],
) -> tuple[float, float, float]:
    bx, by, bz = box
    dx = mobile[0] - ref[0]
    dy = mobile[1] - ref[1]
    dz = mobile[2] - ref[2]
    if bx > 0:
        if dx > bx / 2:
            dx -= bx
        elif dx < -bx / 2:
            dx += bx
    if by > 0:
        if dy > by / 2:
            dy -= by
        elif dy < -by / 2:
            dy += by
    if bz > 0:
        if dz > bz / 2:
            dz -= bz
        elif dz < -bz / 2:
            dz += bz
    return (dx, dy, dz)


def _assign_image_flags_by_bond_tree(
    atoms: list[LammpsAtom],
    bonds: list[LammpsBond],
    box: tuple[float, float, float],
) -> None:
    """Assign LAMMPS image flags by unwrapping along a bond spanning tree."""
    if not atoms:
        return
    bx, by, bz = box
    by_id = {atom.atom_id: atom for atom in atoms}
    unwrapped: dict[int, tuple[float, float, float]] = {
        atom.atom_id: (atom.x, atom.y, atom.z) for atom in atoms
    }
    adjacency = _bond_adjacency(bonds)
    root_id = min(by_id)
    queue: deque[int] = deque([root_id])
    visited = {root_id}

    while queue:
        current_id = queue.popleft()
        current_pos = unwrapped[current_id]
        for neighbor_id in adjacency[current_id]:
            if neighbor_id in visited:
                continue
            visited.add(neighbor_id)
            queue.append(neighbor_id)
            delta = _min_image_vector_to(current_pos, unwrapped[neighbor_id], box)
            unwrapped[neighbor_id] = (
                current_pos[0] + delta[0],
                current_pos[1] + delta[1],
                current_pos[2] + delta[2],
            )

    for atom in atoms:
        ux, uy, uz = unwrapped[atom.atom_id]
        atom.ix = int(math.floor(ux / bx)) if bx > 0 else 0
        atom.iy = int(math.floor(uy / by)) if by > 0 else 0
        atom.iz = int(math.floor(uz / bz)) if bz > 0 else 0
        atom.x = ux - atom.ix * bx if bx > 0 else ux
        atom.y = uy - atom.iy * by if by > 0 else uy
        atom.z = uz - atom.iz * bz if bz > 0 else uz


def _fold_molecules_with_shared_images(
    atoms: list[LammpsAtom],
    box: tuple[float, float, float],
) -> None:
    """Fold each molecule rigidly into the primary cell using one image triple per mol."""
    bx, by, bz = box
    for mol_atoms in _atoms_by_mol_id(atoms).values():
        if not mol_atoms:
            continue
        cg_atoms = [atom for atom in mol_atoms if atom.is_cg]
        ref = cg_atoms[0] if cg_atoms else min(mol_atoms, key=lambda atom: atom.atom_id)
        ix = int(math.floor(ref.x / bx)) if bx > 0 else 0
        iy = int(math.floor(ref.y / by)) if by > 0 else 0
        iz = int(math.floor(ref.z / bz)) if bz > 0 else 0
        for atom in mol_atoms:
            atom.ix = ix
            atom.iy = iy
            atom.iz = iz
            atom.x = atom.x - ix * bx if bx > 0 else atom.x
            atom.y = atom.y - iy * by if by > 0 else atom.y
            atom.z = atom.z - iz * bz if bz > 0 else atom.z


def _shorten_inter_mol_bonds(
    atoms: list[LammpsAtom],
    bonds: list[LammpsBond],
    box: tuple[float, float, float],
    *,
    max_passes: int = 32,
    target_length: float | None = None,
) -> None:
    """Translate whole molecules so cross-mol bonds are short in Cartesian space."""
    by_id = {atom.atom_id: atom for atom in atoms}
    mol_id_by_atom = _mol_id_by_atom_id(atoms)
    threshold = target_length if target_length is not None else min(box) / 2

    for _ in range(max_passes):
        long_bonds: list[tuple[float, LammpsBond]] = []
        for bond in bonds:
            atom_i = by_id.get(bond.i)
            atom_j = by_id.get(bond.j)
            if atom_i is None or atom_j is None:
                continue
            if mol_id_by_atom.get(bond.i) == mol_id_by_atom.get(bond.j):
                continue
            length = math.dist((atom_i.x, atom_i.y, atom_i.z), (atom_j.x, atom_j.y, atom_j.z))
            if length > threshold:
                long_bonds.append((length, bond))

        if not long_bonds:
            return

        long_bonds.sort(reverse=True)
        for _, bond in long_bonds:
            atom_i = by_id[bond.i]
            atom_j = by_id[bond.j]
            mol_i = mol_id_by_atom[bond.i]
            mol_j = mol_id_by_atom[bond.j]
            if mol_j >= mol_i:
                _shift_molecule_near_reference(atoms, mol_j, atom_i, atom_j, box)
            else:
                _shift_molecule_near_reference(atoms, mol_i, atom_j, atom_i, box)


def fold_atoms_with_images(atoms: list[LammpsAtom], box: tuple[float, float, float]) -> None:
    """Fold coordinates into the primary cell and store LAMMPS image flags."""
    bx, by, bz = box
    for atom in atoms:
        x, y, z = atom.x, atom.y, atom.z
        ix = int(math.floor(x / bx)) if bx > 0 else 0
        iy = int(math.floor(y / by)) if by > 0 else 0
        iz = int(math.floor(z / bz)) if bz > 0 else 0
        atom.ix = ix
        atom.iy = iy
        atom.iz = iz
        atom.x = x - ix * bx if bx > 0 else x
        atom.y = y - iy * by if by > 0 else y
        atom.z = z - iz * bz if bz > 0 else z


def _unwrapped_position(
    atom: LammpsAtom, box: tuple[float, float, float]
) -> tuple[float, float, float]:
    bx, by, bz = box
    return (
        atom.x + atom.ix * bx if bx > 0 else atom.x,
        atom.y + atom.iy * by if by > 0 else atom.y,
        atom.z + atom.iz * bz if bz > 0 else atom.z,
    )


def _min_image_displacement(
    pos_i: tuple[float, float, float],
    pos_j: tuple[float, float, float],
    box: tuple[float, float, float],
) -> float:
    bx, by, bz = box
    dx = pos_j[0] - pos_i[0]
    dy = pos_j[1] - pos_i[1]
    dz = pos_j[2] - pos_i[2]
    if bx > 0:
        if dx > bx / 2:
            dx -= bx
        elif dx < -bx / 2:
            dx += bx
    if by > 0:
        if dy > by / 2:
            dy -= by
        elif dy < -by / 2:
            dy += by
    if bz > 0:
        if dz > bz / 2:
            dz -= bz
        elif dz < -bz / 2:
            dz += bz
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _min_image_distance(
    atom_i: LammpsAtom,
    atom_j: LammpsAtom,
    box: tuple[float, float, float],
) -> float:
    return _min_image_displacement(
        _unwrapped_position(atom_i, box),
        _unwrapped_position(atom_j, box),
        box,
    )


def max_bond_length(
    atoms: list[LammpsAtom], bonds: list[LammpsBond], box: tuple[float, float, float]
) -> float:
    """Longest bonded minimum-image distance (Angstrom)."""
    if not bonds:
        return 0.0
    by_id = {atom.atom_id: atom for atom in atoms}
    longest = 0.0
    for bond in bonds:
        atom_i = by_id.get(bond.i)
        atom_j = by_id.get(bond.j)
        if atom_i is None or atom_j is None:
            continue
        length = _min_image_distance(atom_i, atom_j, box)
        longest = max(longest, length)
    return longest


def max_euclidean_bond_length(atoms: list[LammpsAtom], bonds: list[LammpsBond]) -> float:
    """Longest direct bonded distance in primary-cell coordinates (Angstrom)."""
    if not bonds:
        return 0.0
    by_id = {atom.atom_id: atom for atom in atoms}
    longest = 0.0
    for bond in bonds:
        atom_i = by_id.get(bond.i)
        atom_j = by_id.get(bond.j)
        if atom_i is None or atom_j is None:
            continue
        length = math.dist((atom_i.x, atom_i.y, atom_i.z), (atom_j.x, atom_j.y, atom_j.z))
        longest = max(longest, length)
    return longest


def _bond_geometry_limit(interaction_cutoff_ang: float) -> float:
    return max(20.0, 2.0 * interaction_cutoff_ang)


def assign_network_image_flags(system: System) -> None:
    """Fold unwrapped hybrid coordinates and assign bond-tree image flags."""
    _assign_image_flags_by_bond_tree(system.atoms, system.bonds, system.box)


def prepare_network_coordinates(system: System) -> None:
    """Prepare hybrid coordinates for LAMMPS bonded communication."""
    mol_id_by_atom = _mol_id_by_atom_id(system.atoms)
    intra_bonds = _intra_mol_bonds(system.bonds, mol_id_by_atom)
    atoms_by_mol = _atoms_by_mol_id(system.atoms)

    for mol_id in sorted(atoms_by_mol):
        mol_atoms = atoms_by_mol[mol_id]
        mol_atom_ids = {atom.atom_id for atom in mol_atoms}
        mol_bonds = [
            bond for bond in intra_bonds if bond.i in mol_atom_ids and bond.j in mol_atom_ids
        ]
        _unwrap_molecule(mol_atoms, mol_bonds, system.box)

    _place_mols_by_spanning_tree(system.atoms, system.bonds, system.box)
    _shorten_inter_mol_bonds(
        system.atoms,
        system.bonds,
        system.box,
        max_passes=128,
        target_length=15.0,
    )
    _assign_image_flags_by_bond_tree(system.atoms, system.bonds, system.box)


def validate_bond_geometry(system: System, interaction_cutoff_ang: float) -> None:
    """Fail fast when bonded pairs are unphysical after PBC preparation."""
    limit = _bond_geometry_limit(interaction_cutoff_ang)
    max_bond = max_bond_length(system.atoms, system.bonds, system.box)
    if max_bond > limit:
        raise ValueError(
            f"Bonded pair spans {max_bond:.1f} Å (min-image) after PBC prep (limit {limit:.1f} Å); "
            "check hybrid coordinates or bond topology"
        )
