"""1-4 pairs come from the bond graph, once each (GROMACS convention).

bakery's network [ cross_pairs ] lists are incomplete and contain duplicates,
so the listed lines only supply explicit parameters.
"""

from __future__ import annotations

from backmap_prep.network.lammps_builder import _topological_14_pairs
from backmap_prep.parsers.top_parser import MoleculeType, TopBond


def _molecule(bonds: list[tuple[int, int]], cross: list[tuple[int, int]] = ()) -> MoleculeType:
    mol = MoleculeType(name="MOL", nrexcl=3)
    mol.bonds = [TopBond(i, j, 1) for i, j in bonds]
    mol.cross_bonds = [TopBond(i, j, 1) for i, j in cross]
    return mol


def test_chain_has_one_pair_per_dihedral() -> None:
    mol = _molecule([(1, 2), (2, 3), (3, 4), (4, 5)])
    assert _topological_14_pairs(mol, {1, 2, 3, 4, 5}) == {(1, 4), (2, 5)}


def test_pairs_closer_through_a_ring_are_not_1_4() -> None:
    # Five-membered ring: every pair is at most two bonds apart.
    ring = _molecule([(1, 2), (2, 3), (3, 4), (4, 5), (5, 1)])
    assert _topological_14_pairs(ring, {1, 2, 3, 4, 5}) == set()


def test_cross_bonds_count_and_cg_atoms_are_skipped() -> None:
    # 1-2 intra, 2-3 crosslink, 3-4 intra; atom 9 is a CG bead bonded to 1.
    mol = _molecule([(1, 2), (3, 4), (9, 1)], cross=[(2, 3)])
    assert _topological_14_pairs(mol, {1, 2, 3, 4}) == {(1, 4)}
