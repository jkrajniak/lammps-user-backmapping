"""GROMACS func-1 dihedrals, V = k[1 + cos(n phi - phi_s)]: the phase must survive the mapping."""

from __future__ import annotations

import pytest

from backmap_prep import units
from backmap_prep.builder import System
from backmap_prep.network.lammps_builder import _dihedral_terms
from backmap_prep.parsers.top_parser import MoleculeType, TopAtom, TopDihedral, Topology


def _style_params(phase: float) -> tuple[str, list[float]]:
    mol = MoleculeType(name="MOL", nrexcl=3)
    mol.atoms = [TopAtom(i, "C", 1, "MOL", f"C{i}", i, 0.0, 12.0) for i in range(1, 5)]
    mol.dihedrals = [TopDihedral(1, 2, 3, 4, 1, [phase, 4.6, 3])]
    system = System()
    _dihedral_terms(system, mol, Topology(), {}, set(), [])
    (t,) = system.dihedral_types
    return t.style, t.params


@pytest.mark.parametrize(("phase", "sign"), [(0.0, 1.0), (180.0, -1.0), (-180.0, -1.0)])
def test_phase_0_or_180_stays_harmonic(phase: float, sign: float) -> None:
    style, params = _style_params(phase)
    assert style == "backmap/harmonic"
    assert params == pytest.approx([units.energy(4.6), sign, 3.0])


@pytest.mark.parametrize("phase", [60.0, -90.0, 120.0])
def test_other_phases_become_fourier_with_the_phase(phase: float) -> None:
    style, params = _style_params(phase)
    assert style == "backmap/fourier"
    assert params == pytest.approx([1.0, units.energy(4.6), 3.0, phase])


def test_harmonic_conversion_rejects_a_phase_it_cannot_hold() -> None:
    with pytest.raises(ValueError, match="no phase 60"):
        units.gromacs_harmonic_dihedral(4.6, 3, 60.0)
