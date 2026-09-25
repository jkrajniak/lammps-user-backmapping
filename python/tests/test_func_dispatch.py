"""Bonds and angles are converted by GROMACS function type, never by parameter count.

Any two-parameter angle used to become LAMMPS harmonic, so MARTINI G96 angles
(func 2) and Urey-Bradley angles (func 5) were converted wrongly without error.
"""

from __future__ import annotations

import pytest

from backmap_prep import units
from backmap_prep.builder import System
from backmap_prep.network.lammps_builder import _angle_terms, _bond_terms
from backmap_prep.parsers.top_parser import MoleculeType, TopAngle, TopAtom, TopBond


def _molecule() -> MoleculeType:
    mol = MoleculeType(name="MOL", nrexcl=3)
    mol.atoms = [TopAtom(i, "opls_135", 1, "MOL", f"C{i}", i, 0.0, 12.0) for i in (1, 2, 3)]
    return mol


def test_harmonic_bond_and_angle_convert() -> None:
    mol = _molecule()
    mol.bonds = [TopBond(1, 2, 1, [0.153, 224262.4])]
    mol.angles = [TopAngle(1, 2, 3, 1, [112.7, 488.273])]
    system = System()
    _bond_terms(system, mol, {}, set(), [])
    _angle_terms(system, mol, {}, set(), [])
    assert system.bond_types[0].params == pytest.approx(
        [units.spring_bond(224262.4), units.distance(0.153)]
    )
    assert system.angle_types[0].params == pytest.approx([units.spring_angle(488.273), 112.7])


@pytest.mark.parametrize("func", [2, 5, 6])
def test_other_angle_functions_are_an_error(func: int) -> None:
    mol = _molecule()
    mol.angles = [TopAngle(1, 2, 3, func, [108.0, 21.5, 0.25, 1000.0])]
    with pytest.raises(ValueError, match=rf"Unsupported angle func {func}"):
        _angle_terms(System(), mol, {}, set(), [])


@pytest.mark.parametrize("func", [2, 3, 6])
def test_other_bond_functions_are_an_error(func: int) -> None:
    mol = _molecule()
    mol.bonds = [TopBond(1, 2, func, [0.47, 1250.0, 2.0])]
    with pytest.raises(ValueError, match=rf"Unsupported bond func {func}"):
        _bond_terms(System(), mol, {}, set(), [])
