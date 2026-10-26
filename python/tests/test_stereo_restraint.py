"""Stereocenter detection and the backmapping-time restraint (``simulation.stereo_restraint``)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

from backmap_prep.builder import AtomTypeInfo, LammpsAtom, LammpsBond, System
from backmap_prep.schema import SimulationParams
from backmap_prep.stereo import (
    IDEAL_IMPROPER_DEG,
    _dihedral,
    find_stereocenters,
    restrain_target,
    write_stereo_files,
)

if TYPE_CHECKING:
    from pathlib import Path

# tetrahedron directions
TET = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float) / np.sqrt(3)


def _molecule(
    neighbor_types: list[int], mirror: bool = False, first: np.ndarray | None = None
) -> System:
    """A center (type 1) with four neighbors of the given types; type 10 is hydrogen.

    ``first`` replaces the direction of the first neighbor.
    """
    directions = TET.copy()
    if mirror:
        directions[:, 0] *= -1
    if first is not None:
        directions[0] = first
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
    assert len(right) == len(left) == 2
    for r, m in zip(right, left, strict=True):
        assert r[3] == m[3] == 1
        assert abs(r[4]) == pytest.approx(IDEAL_IMPROPER_DEG)
        assert np.sign(r[4]) == -np.sign(m[4])


def test_dihedral_has_the_iupac_sign() -> None:
    """Looking along 2 -> 3, the bond to 1 turns clockwise onto the bond to 4: +90 degrees."""
    points = [np.array(p, dtype=float) for p in [(1, 0, 0), (0, 0, 0), (0, 0, 1), (0, 1, 1)]]
    assert _dihedral(*points) == pytest.approx(90.0)


@pytest.mark.parametrize("mirror", [False, True])
def test_target_is_the_improper_of_the_current_geometry(mirror: bool) -> None:
    system = _molecule([10, 11, 12, 13], mirror=mirror)
    pos = {a.atom_id: np.array([a.x, a.y, a.z]) for a in system.atoms}
    for a, b, c, center, phi0 in find_stereocenters(system):
        assert phi0 == pytest.approx(_dihedral(pos[a], pos[b], pos[c], pos[center]), abs=1e-3)


@pytest.mark.parametrize(
    ("phi", "target"),
    [
        (-IDEAL_IMPROPER_DEG, 180.0 - IDEAL_IMPROPER_DEG),
        (IDEAL_IMPROPER_DEG, IDEAL_IMPROPER_DEG - 180.0),
    ],
)
def test_restrain_target_puts_the_minimum_of_fix_restrain_at_phi(phi: float, target: float) -> None:
    assert restrain_target(phi) == pytest.approx(target)
    # fix restrain dihedral: E = K[1 + cos(phi - phi0)], zero at phi = phi0 + 180
    assert 1.0 + np.cos(np.radians(phi - restrain_target(phi))) == pytest.approx(0.0, abs=1e-12)


def test_heavy_quartet_leaves_hydrogen_out_and_the_second_holds_all_four_neighbors() -> None:
    heavy, second = find_stereocenters(_molecule([10, 11, 12, 13]))
    assert set(heavy[:3]) == {3, 4, 5}
    assert set(heavy[:3]) | set(second[:3]) == {2, 3, 4, 5}
    assert second[2] == 2


def test_hydrogen_crossing_to_the_heavy_neighbors_flips_only_the_second_quartet() -> None:
    """The failure of a heavy-atom restraint: the hydrogen, not the center, crosses over."""
    kept = _molecule([10, 11, 12, 13])
    crossed = _molecule([10, 11, 12, 13], first=-0.6 * TET[0])
    pos = {a.atom_id: np.array([a.x, a.y, a.z]) for a in crossed.atoms}
    for n, (a, b, c, center, phi0) in enumerate(find_stereocenters(kept)):
        phi = _dihedral(pos[a], pos[b], pos[c], pos[center])
        assert (np.sign(phi) == np.sign(phi0)) == (n == 0)


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
    centers = {c[3] for c in find_stereocenters(system)}
    assert 1 in centers


def test_files_hold_two_restraints_per_center(tmp_path: Path) -> None:
    n = write_stereo_files(tmp_path, _molecule([10, 11, 12, 13]), 83.68)
    assert n == 1
    fixes = [line.split() for line in (tmp_path / "stereo.lmp").read_text().splitlines()[1:]]
    quartets = find_stereocenters(_molecule([10, 11, 12, 13]))
    assert len(fixes) == len(quartets) == 2
    for k, (fix, quartet) in enumerate(zip(fixes, quartets, strict=True), 1):
        assert fix[:5] == ["fix", f"stereo{k}", "all", "restrain", "dihedral"]
        assert [int(x) for x in fix[5:9]] == list(quartet[:4])
        assert float(fix[-3]) == pytest.approx(83.68 / 4.184, rel=1e-3)
        assert float(fix[-1]) == pytest.approx(restrain_target(quartet[4]), abs=1e-4)
    release = (tmp_path / "stereo_release.lmp").read_text().splitlines()
    assert release[1:] == ["unfix stereo1", "unfix stereo2"]


def test_restraint_is_off_by_default_and_must_be_positive() -> None:
    assert SimulationParams().stereo_restraint is None
    with pytest.raises(ValueError, match="stereo_restraint"):
        SimulationParams(stereo_restraint=-1.0)
