"""Backmapping-time LJ core on attracting 1-4 pairs without LJ (``simulation.core_guard``)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from backmap_prep.builder import AtomTypeInfo, LammpsAtom, LammpsCrossPair, System
from backmap_prep.core_guard import unguarded_pairs, write_core_guard_files
from backmap_prep.schema import SimulationParams

if TYPE_CHECKING:
    from pathlib import Path


def _cooh() -> System:
    """The atoms of a COOH group and its ring carbon, with OPLS-AA charges of the PET example.

    1 O (OH), 2 C, 3 O (=O), 4 H (acid, no LJ), 5 ring C.
    """
    system = System(box=(30.0, 30.0, 30.0))
    system.atom_types = [AtomTypeInfo(1, "O", 15.9994, False), AtomTypeInfo(2, "H", 1.008, False)]
    charges = {1: -0.56, 2: 0.56, 3: -0.48, 4: 0.43, 5: -0.13}
    for i, q in charges.items():
        system.atoms.append(LammpsAtom(i, 1, 2 if i == 4 else 1, q, float(i), 0.0, 0.0))
    system.has_cross_pairs = True
    system.cross_pairs = [
        LammpsCrossPair(3, 4, 0.0, 0.0, qq_scale=0.5),  # carbonyl O / acid H: guarded
        LammpsCrossPair(4, 5, 0.0, 0.0, qq_scale=0.5),  # acid H / negative ring C: guarded
        LammpsCrossPair(1, 5, 3.0, 0.1, qq_scale=0.5),  # has LJ
        LammpsCrossPair(2, 4, 0.0, 0.0, qq_scale=0.5),  # repulsive Coulomb
        LammpsCrossPair(3, 1, 0.0, 0.0, qq_scale=0.0),  # no Coulomb
    ]
    return system


def test_only_attracting_pairs_without_lj_are_guarded() -> None:
    assert [(p.i, p.j) for p in unguarded_pairs(_cooh())] == [(3, 4), (4, 5)]


def test_zero_sigma_or_zero_epsilon_counts_as_no_lj() -> None:
    system = _cooh()
    system.cross_pairs = [
        LammpsCrossPair(3, 4, 2.0, 0.0, qq_scale=0.5),
        LammpsCrossPair(4, 5, 0.0, 0.2, qq_scale=0.5),
    ]
    assert len(unguarded_pairs(system)) == 2


def test_files_and_units(tmp_path: Path) -> None:
    assert write_core_guard_files(tmp_path, _cooh(), 0.15, 0.2092) == 2
    lines = (tmp_path / "core_guard.dat").read_text().splitlines()
    assert lines[0] == "2"
    # LJ only (four columns): 0.15 nm = 1.5 A, 0.2092 kJ/mol = 0.05 kcal/mol
    i, j, sigma, epsilon = lines[1].split()
    assert (int(i), int(j)) == (3, 4)
    assert float(sigma) == pytest.approx(1.5)
    assert float(epsilon) == pytest.approx(0.05)
    fix = (tmp_path / "core_guard.lmp").read_text().splitlines()[1].split()
    assert fix == ["fix", "core_guard", "all", "backmap/pairs", "at", "file", "core_guard.dat"]
    assert "unfix core_guard" in (tmp_path / "core_guard_release.lmp").read_text()


def test_no_files_without_unguarded_pairs(tmp_path: Path) -> None:
    system = _cooh()
    system.cross_pairs = system.cross_pairs[2:]
    assert write_core_guard_files(tmp_path, system, 0.15, 0.2092) == 0
    assert not list(tmp_path.iterdir())


def test_guard_is_on_by_default_and_must_be_positive() -> None:
    params = SimulationParams()
    assert params.core_guard == pytest.approx(0.15)
    assert params.core_guard_epsilon == pytest.approx(0.2092)
    assert SimulationParams(core_guard=None).core_guard is None
    with pytest.raises(ValueError, match="core_guard"):
        SimulationParams(core_guard=-0.15)
    with pytest.raises(ValueError, match="core_guard_epsilon"):
        SimulationParams(core_guard_epsilon=0.0)
