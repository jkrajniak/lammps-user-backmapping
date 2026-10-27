"""Ring detection and the backmapping-time ring guard (``simulation.ring_guard``)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

from backmap_prep.builder import AtomTypeInfo, LammpsAtom, LammpsBond, System
from backmap_prep.ring_guard import find_rings, write_ring_guard_files
from backmap_prep.schema import SimulationParams

if TYPE_CHECKING:
    from pathlib import Path


def _system(n_atoms: int, bonds: list[tuple[int, int]], cg: set[int] | None = None) -> System:
    """Atoms 1..n_atoms on a line (positions do not matter for detection), the given bonds."""
    cg = cg or set()
    system = System(box=(50.0, 50.0, 50.0))
    system.atom_types = [AtomTypeInfo(1, "C", 12.011, False), AtomTypeInfo(2, "A", 72.0, True)]
    for i in range(1, n_atoms + 1):
        system.atoms.append(
            LammpsAtom(i, 1, 2 if i in cg else 1, 0.0, float(i), 0.0, 0.0, is_cg=i in cg)
        )
    for n, (i, j) in enumerate(bonds, 1):
        system.bonds.append(LammpsBond(n, 1, i, j))
    return system


def _cycle(atoms: list[int]) -> list[tuple[int, int]]:
    return [(a, atoms[(k + 1) % len(atoms)]) for k, a in enumerate(atoms)]


def test_benzene_ring_in_ring_order_from_the_lowest_id() -> None:
    # ring 1-2-3-4-5-6 with substituents 7 on 1 and 8 on 4
    system = _system(8, [*_cycle([3, 1, 2, 6, 5, 4]), (1, 7), (4, 8)])
    rings = find_rings(system)
    assert len(rings) == 1
    ring = rings[0]
    assert ring[0] == 1
    assert set(ring) == {1, 2, 3, 4, 5, 6}
    # consecutive atoms are bonded
    bonded = {frozenset(b) for b in _cycle([3, 1, 2, 6, 5, 4])}
    assert all(frozenset((ring[k], ring[(k + 1) % 6])) in bonded for k in range(6))


@pytest.mark.parametrize("size", [5, 6, 7])
def test_rings_of_five_to_seven_atoms_are_guarded(size: int) -> None:
    assert [set(r) for r in find_rings(_system(size, _cycle(list(range(1, size + 1)))))] == [
        set(range(1, size + 1))
    ]


@pytest.mark.parametrize("size", [3, 4, 8])
def test_small_and_large_rings_are_not_guarded(size: int) -> None:
    assert find_rings(_system(size, _cycle(list(range(1, size + 1))))) == []


def test_fused_rings_give_one_ring_each_and_not_the_envelope() -> None:
    # naphthalene skeleton: rings 1-2-3-4-5-6 and 5-6-7-8-9-10 share the bond 5-6
    bonds = [*_cycle([1, 2, 3, 4, 5, 6]), (5, 10), (10, 9), (9, 8), (8, 7), (7, 6)]
    rings = find_rings(_system(10, bonds))
    assert sorted(set(r) for r in rings) == sorted([{1, 2, 3, 4, 5, 6}, {5, 6, 7, 8, 9, 10}])


def test_chain_has_no_rings() -> None:
    assert find_rings(_system(6, [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)])) == []


def test_cg_atoms_are_ignored() -> None:
    assert find_rings(_system(6, _cycle([1, 2, 3, 4, 5, 6]), cg={6})) == []


def test_files_and_units(tmp_path: Path) -> None:
    system = _system(12, [*_cycle([1, 2, 3, 4, 5, 6]), *_cycle([7, 8, 9, 10, 11, 12])])
    assert write_ring_guard_files(tmp_path, system, 0.2, 41840.0, 0.12) == 2
    rows = [
        line.split("|")[0].split()
        for line in (tmp_path / "rings.txt").read_text().splitlines()
        if not line.startswith("#")
    ]
    assert [set(map(int, r)) for r in rows] == [set(range(1, 7)), set(range(7, 13))]
    fix = (tmp_path / "ring_guard.lmp").read_text().splitlines()[1].split()
    assert fix[:4] == ["fix", "ring_guard", "at_atoms", "backmap/ringguard"]
    # 41840 kJ/mol/nm^2 = 100 kcal/mol/A^2, 0.2 nm = 2 A
    assert float(fix[4]) == pytest.approx(100.0)
    assert float(fix[5]) == pytest.approx(2.0)
    assert float(fix[6]) == pytest.approx(1.2)
    assert fix[7] == "rings.txt"
    assert "unfix ring_guard" in (tmp_path / "ring_guard_release.lmp").read_text()


def test_no_files_without_rings(tmp_path: Path) -> None:
    assert write_ring_guard_files(tmp_path, _system(3, [(1, 2), (2, 3)]), 0.2, 41840.0) == 0
    assert not list(tmp_path.iterdir())


def test_guard_is_off_by_default_and_must_be_positive() -> None:
    params = SimulationParams()
    assert params.ring_guard is None
    assert params.ring_guard_k == pytest.approx(418400.0)
    assert params.ring_guard_bonded == pytest.approx(0.12)
    with pytest.raises(ValueError, match="ring_guard"):
        SimulationParams(ring_guard=-0.2)
    with pytest.raises(ValueError, match="ring_guard_k"):
        SimulationParams(ring_guard=0.2, ring_guard_k=0.0)
    with pytest.raises(ValueError, match="ring_guard_bonded"):
        SimulationParams(ring_guard=0.2, ring_guard_bonded=0.25)


def test_positions_do_not_change_detection() -> None:
    system = _system(6, _cycle([1, 2, 3, 4, 5, 6]))
    rng = np.random.default_rng(1)
    for atom in system.atoms:
        atom.x, atom.y, atom.z = rng.uniform(0.0, 50.0, 3)
    assert len(find_rings(system)) == 1


def test_atoms_bonded_to_a_ring_are_listed_after_the_bar(tmp_path: Path) -> None:
    # ring 1-6, hydrogen 7 on 1, CH2 8 on 4 with N 9 on 8: 7 and 8 are excluded, N 9 is not
    system = _system(9, [*_cycle([1, 2, 3, 4, 5, 6]), (1, 7), (4, 8), (8, 9)])
    assert write_ring_guard_files(tmp_path, system, 0.2, 41840.0) == 1
    line = (tmp_path / "rings.txt").read_text().splitlines()[1]
    ring, excluded = line.split("|")
    assert set(map(int, ring.split())) == set(range(1, 7))
    assert sorted(map(int, excluded.split())) == [7, 8]
