"""Unit tests for hybrid dihedral intra/cross-bead classification."""

from __future__ import annotations

from backmap_prep.builder import LammpsAtom
from backmap_prep.network.lammps_builder import _is_intra_bead_dihedral


def _atom(atom_id: int, mol_id: int) -> LammpsAtom:
    return LammpsAtom(
        atom_id=atom_id,
        mol_id=mol_id,
        type_id=1,
        charge=0.0,
        x=0.0,
        y=0.0,
        z=0.0,
        type_name="X",
        is_cg=False,
    )


def test_intra_bead_requires_all_four_atoms() -> None:
    atoms = {1: _atom(1, 10), 2: _atom(2, 10), 3: _atom(3, 10)}
    assert not _is_intra_bead_dihedral(atoms, 1, 2, 3, 4)


def test_intra_bead_same_mol_id() -> None:
    atoms = {i: _atom(i, 7) for i in (1, 2, 3, 4)}
    assert _is_intra_bead_dihedral(atoms, 1, 2, 3, 4)


def test_cross_bead_different_mol_ids() -> None:
    atoms = {
        1: _atom(1, 1),
        2: _atom(2, 1),
        3: _atom(3, 2),
        4: _atom(4, 2),
    }
    assert not _is_intra_bead_dihedral(atoms, 1, 2, 3, 4)
