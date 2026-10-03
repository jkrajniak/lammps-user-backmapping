"""Tests for network periodic boundary preparation."""

from __future__ import annotations

import math

import pytest

from backmap_prep.builder import LammpsAngle, LammpsAtom, LammpsBond, System
from backmap_prep.network.pbc import (
    fold_atoms_with_images,
    max_bond_length,
    max_interaction_extent,
    prepare_network_coordinates,
)


def test_fold_preserves_minimum_image_bond() -> None:
    box = (50.0, 50.0, 50.0)
    atoms = [
        LammpsAtom(1, 1, 1, 0.0, 1.0, 2.0, 3.0, "A", True),
        LammpsAtom(2, 2, 1, 0.0, 49.0, 2.0, 3.0, "B", True),
    ]
    bonds = [LammpsBond(1, 1, 1, 2)]
    fold_atoms_with_images(atoms, box)
    assert max_bond_length(atoms, bonds, box) < 5.0
    assert all(0.0 <= atom.x < box[0] for atom in atoms)


def test_prepare_network_coordinates_keeps_bonds_short_in_cycle() -> None:
    box = (40.0, 40.0, 40.0)
    atoms = [
        LammpsAtom(1, 1, 1, 0.0, 1.0, 1.0, 1.0, "A", True),
        LammpsAtom(2, 1, 1, 0.0, 39.0, 1.0, 1.0, "B", True),
        LammpsAtom(3, 1, 1, 0.0, 39.0, 39.0, 1.0, "C", True),
        LammpsAtom(4, 1, 1, 0.0, 1.0, 39.0, 1.0, "D", True),
    ]
    bonds = [
        LammpsBond(1, 1, 1, 2),
        LammpsBond(2, 1, 2, 3),
        LammpsBond(3, 1, 3, 4),
        LammpsBond(4, 1, 4, 1),
    ]
    system = System(atoms=atoms, bonds=bonds, box=box, has_cross_bonds=True)
    prepare_network_coordinates(system)
    assert max_bond_length(system.atoms, system.bonds, system.box) < 5.0


def test_prepare_network_coordinates_shortens_cross_mol_bond() -> None:
    box = (50.0, 50.0, 50.0)
    atoms = [
        LammpsAtom(1, 1, 1, 0.0, 5.0, 5.0, 5.0, "A", True),
        LammpsAtom(2, 1, 1, 0.0, 1.0, 5.0, 5.0, "B", False),
        LammpsAtom(3, 2, 1, 0.0, 49.0, 5.0, 5.0, "C", True),
        LammpsAtom(4, 2, 1, 0.0, 48.0, 5.0, 5.0, "D", False),
    ]
    bonds = [
        LammpsBond(1, 1, 1, 2),
        LammpsBond(2, 1, 3, 4),
        LammpsBond(3, 1, 2, 3),
    ]
    system = System(atoms=atoms, bonds=bonds, box=box, has_cross_bonds=True)
    prepare_network_coordinates(system)
    assert max_bond_length(system.atoms, system.bonds, system.box) < 5.0


def test_fold_sets_image_flags_for_out_of_box_coords() -> None:
    atoms = [LammpsAtom(1, 1, 1, 0.0, 55.0, -2.0, 3.0, "A", True)]
    fold_atoms_with_images(atoms, (50.0, 50.0, 50.0))
    assert atoms[0].ix == 1
    assert atoms[0].iy == -1
    assert atoms[0].x == pytest.approx(5.0)
    assert atoms[0].y == pytest.approx(48.0)


def _max_flag_unwrapped_bond(system: System) -> float:
    bx, by, bz = system.box
    by_id = {atom.atom_id: atom for atom in system.atoms}

    def unwrapped(atom: LammpsAtom) -> tuple[float, float, float]:
        return (atom.x + atom.ix * bx, atom.y + atom.iy * by, atom.z + atom.iz * bz)

    return max(
        math.dist(unwrapped(by_id[bond.i]), unwrapped(by_id[bond.j])) for bond in system.bonds
    )


def test_prepare_network_coordinates_flags_consistent_in_every_component() -> None:
    """Image flags must follow the bonds in every connected component.

    A hybrid melt molecule has at least two components (the CG chain and the
    AT chain are not bonded to each other), and a melt has many molecules.
    Walking the bond graph from one root only left the other components with
    per-atom flags, so bonded atoms a bond length apart in the cell got
    flags one box vector apart.
    """
    box = (20.0, 20.0, 20.0)
    n_at = 24  # 12 two-atom beads; the AT chain crosses the x boundary
    atoms = []
    bonds = []
    for mol in (1, 2):
        base = 100 * (mol - 1)
        y = 5.0 + 8.0 * (mol - 1)
        atoms.append(LammpsAtom(base + 1, mol, 1, 0.0, 10.0, y, 5.0, "A", True))
        atoms.append(LammpsAtom(base + 2, mol, 1, 0.0, 14.0, y, 5.0, "B", True))
        for k in range(n_at):
            x = (15.0 + 1.2 * k) % box[0]
            atoms.append(LammpsAtom(base + 3 + k, mol, 2, 0.0, x, y + 1.0, 5.0, "C", False))
        bonds.append(LammpsBond(len(bonds) + 1, 1, base + 1, base + 2))
        # backmap-prep order: all intra-bead bonds first, then the cross-bead ones
        for k in range(0, n_at, 2):
            bonds.append(LammpsBond(len(bonds) + 1, 2, base + 3 + k, base + 4 + k))
        for k in range(1, n_at - 1, 2):
            bonds.append(LammpsBond(len(bonds) + 1, 3, base + 3 + k, base + 4 + k))
    system = System(atoms=atoms, bonds=bonds, box=box, has_cross_bonds=True)

    prepare_network_coordinates(system)

    assert max_bond_length(system.atoms, system.bonds, system.box) < 4.5
    assert _max_flag_unwrapped_bond(system) < 4.5
    assert all(0.0 <= atom.x < box[0] for atom in system.atoms)


def test_interaction_extent_uses_minimum_image_not_folded_coordinates() -> None:
    """A bond across the boundary has a small extent, not about a box length.

    The communication cutoff used to include the distance between folded file
    coordinates, which inflated it to ~85-100 A for any system with a bond
    crossing the box boundary.
    """
    box = (40.0, 40.0, 40.0)
    atoms = [
        LammpsAtom(1, 1, 1, 0.0, 39.5, 5.0, 5.0, "A", True),
        LammpsAtom(2, 1, 2, 0.0, 0.5, 5.0, 5.0, "C", False),
        LammpsAtom(3, 1, 2, 0.0, 1.5, 5.0, 5.0, "C", False),
    ]
    bonds = [LammpsBond(1, 1, 2, 3)]
    angles = [LammpsAngle(1, 1, 1, 2, 3)]
    system = System(atoms=atoms, bonds=bonds, angles=angles, box=box)

    extent = max_interaction_extent(system)

    assert extent == pytest.approx(2.0)
