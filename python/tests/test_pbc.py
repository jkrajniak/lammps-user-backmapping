"""Tests for network periodic boundary preparation."""

from __future__ import annotations

import pytest

from backmap_prep.builder import LammpsAtom, LammpsBond, System
from backmap_prep.network.pbc import (
    fold_atoms_with_images,
    max_bond_length,
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
