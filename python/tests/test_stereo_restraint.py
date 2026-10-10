"""Stereocenter detection and the backmapping-time restraint (``simulation.stereo_restraint``)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

from backmap_prep.builder import AtomTypeInfo, LammpsAtom, LammpsBond, System
from backmap_prep.schema import SimulationParams
from backmap_prep.stereo import IDEAL_IMPROPER_DEG, find_stereocenters, write_stereo_files

if TYPE_CHECKING:
    from pathlib import Path

# tetrahedron directions
TET = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float) / np.sqrt(3)


def _molecule(neighbor_types: list[int], mirror: bool = False) -> System:
    """A center (type 1) with four neighbors of the given types; type 10 is hydrogen."""
    directions = TET.copy()
    if mirror:
        directions[:, 0] *= -1
    system = System(box=(50.0, 50.0, 50.0))
    system.atom_types = [
        AtomTypeInfo(1, "C", 12.0, False),
        AtomTypeInfo(10, "H", 1.008, False),
        AtomTypeInfo(11, "F", 19.0, False),
        AtomTypeInfo(12, "Cl", 35.5, False),
        AtomTypeInfo(13, "Br", 79.9, False),
    ]
    system.atoms.append(LammpsAtom(1, 1, 1, 0.0, 25.0, 25.0, 25.0))
    for n, (t, d) in enumerate(zip(neighbor_types, directions, strict=True), 2):
        system.atoms.append(LammpsAtom(n, 1, t, 0.0, *(25.0 + 1.5 * d)))
        system.bonds.append(LammpsBond(n - 1, 1, 1, n))
    return system


def test_chiral_center_gets_ideal_improper_with_the_sign_of_its_handedness() -> None:
    right = find_stereocenters(_molecule([10, 11, 12, 13]))
    left = find_stereocenters(_molecule([10, 11, 12, 13], mirror=True))
    assert len(right) == len(left) == 1
    assert right[0][3] == left[0][3] == 1
    assert abs(right[0][4]) == pytest.approx(IDEAL_IMPROPER_DEG)
    assert np.sign(right[0][4]) == -np.sign(left[0][4])


def test_hydrogen_is_left_out() -> None:
    a, b, c, _, _ = find_stereocenters(_molecule([10, 11, 12, 13]))[0]
    assert 2 not in (a, b, c)


def test_center_with_equivalent_neighbors_is_not_a_stereocenter() -> None:
    assert find_stereocenters(_molecule([10, 10, 12, 13])) == []
    assert find_stereocenters(_molecule([11, 11, 11, 13])) == []


def test_environment_distinguishes_neighbors_of_the_same_type() -> None:
    """Two CH2 neighbors of the same type are inequivalent when their substituents differ."""
    system = _molecule([10, 1, 1, 11])
    # extend the two same-type neighbors (atoms 3 and 4) by different atoms
    system.atoms.append(LammpsAtom(6, 1, 12, 0.0, 25.0, 25.0, 30.0))
    system.atoms.append(LammpsAtom(7, 1, 13, 0.0, 25.0, 30.0, 25.0))
    system.bonds.append(LammpsBond(5, 1, 3, 6))
    system.bonds.append(LammpsBond(6, 1, 4, 7))
    centers = [c for c in find_stereocenters(system) if c[3] == 1]
    assert len(centers) == 1


def test_files_hold_one_restraint_per_center(tmp_path: Path) -> None:
    n = write_stereo_files(tmp_path, _molecule([10, 11, 12, 13]), 83.68)
    assert n == 1
    fix = (tmp_path / "stereo.lmp").read_text().splitlines()[-1].split()
    assert fix[:4] == ["fix", "stereo1", "all", "restrain"]
    assert float(fix[-3]) == pytest.approx(83.68 / 4.184, rel=1e-3)
    assert "unfix stereo1" in (tmp_path / "stereo_release.lmp").read_text()


def test_restraint_is_off_by_default_and_must_be_positive() -> None:
    assert SimulationParams().stereo_restraint is None
    with pytest.raises(ValueError, match="stereo_restraint"):
        SimulationParams(stereo_restraint=-1.0)
