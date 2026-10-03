"""GROMACS Urey-Bradley angles, multi-term periodic dihedrals and harmonic impropers.

Slipids / CHARMM36 lipids use angle func 5, dihedral func 9 (several lines per
atom-type quadruple, wildcards) and improper func 2 (wildcards in the middle).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from backmap_prep import units
from backmap_prep.builder import System
from backmap_prep.network.lammps_builder import _angle_terms, _dihedral_terms, _improper_terms
from backmap_prep.parsers.top_parser import (
    AtomType,
    MoleculeType,
    TopAngle,
    TopAtom,
    TopDihedral,
    Topology,
    _parse_file,
    lookup_dihedraltype_most_specific,
)
from backmap_prep.writers import write_lammps_data

if TYPE_CHECKING:
    from pathlib import Path

FF = """
[ angletypes ]
CTL2 CTL2 CTL2   5   113.60   488.273   0.25610   9338.72
[ dihedraltypes ]
CTL2 CTL2 CTL2 CTL2   9   0.0   0.54392   2
CTL2 CTL2 CTL2 CTL2   9   180.0 0.30962   3
CTL2 CTL2 CTL2 CTL2   9   60.0  0.10     6
X    CTL2 CTL2 X      9   0.0   0.8368    3
OBL  X    X    CTL2   2   0.00  836.8
"""


def _topology(tmp_path: Path) -> Topology:
    itp = tmp_path / "ff.itp"
    itp.write_text(FF)
    top = Topology()
    _parse_file(itp, top, [tmp_path], [])
    top.atom_types = {
        "CTL2": AtomType("CTL2", 12.011, 0.0, "A"),
        "OBL": AtomType("OBL", 15.999, 0.0, "A"),
    }
    return top


def _molecule(types: list[str]) -> MoleculeType:
    mol = MoleculeType(name="MOL", nrexcl=3)
    mol.atoms = [
        TopAtom(i + 1, t, 1, "MOL", f"A{i + 1}", i + 1, 0.0, 12.0) for i, t in enumerate(types)
    ]
    return mol


def test_consecutive_func9_lines_are_one_entry(tmp_path: Path) -> None:
    top = _topology(tmp_path)
    entry = lookup_dihedraltype_most_specific(top.dihedraltypes, ("CTL2",) * 4, {9})
    assert entry is not None
    assert entry.terms == [[0.0, 0.54392, 2.0], [180.0, 0.30962, 3.0], [60.0, 0.10, 6.0]]


def test_most_specific_match_wins_and_middle_wildcards_resolve(tmp_path: Path) -> None:
    top = _topology(tmp_path)
    exact = lookup_dihedraltype_most_specific(top.dihedraltypes, ("CTL2",) * 4, {9})
    wild = lookup_dihedraltype_most_specific(top.dihedraltypes, ("OBL", "CTL2", "CTL2", "OBL"), {9})
    improper = lookup_dihedraltype_most_specific(
        top.dihedraltypes, ("OBL", "CTL2", "OBL", "CTL2"), {2}
    )
    assert exact is not None
    assert len(exact.terms) == 3
    assert wild is not None
    assert wild.terms == [[0.0, 0.8368, 3.0]]
    assert improper is not None
    assert improper.params == [0.0, 836.8]


def test_urey_bradley_angle(tmp_path: Path) -> None:
    top = _topology(tmp_path)
    mol = _molecule(["CTL2"] * 3)
    mol.angles = [TopAngle(1, 2, 3, 5)]
    system = System()
    _angle_terms(system, mol, {}, set(), [], topology=top)
    (t,) = system.angle_types
    assert t.style == "backmap/charmm"
    assert t.params == pytest.approx(
        [units.spring_angle(488.273), 113.60, units.spring_bond(9338.72), units.distance(0.25610)]
    )


def test_func9_dihedral_becomes_fourier(tmp_path: Path) -> None:
    top = _topology(tmp_path)
    mol = _molecule(["CTL2"] * 4)
    mol.dihedrals = [TopDihedral(1, 2, 3, 4, 9)]
    system = System()
    _dihedral_terms(system, mol, top, {}, set(), [])
    (t,) = system.dihedral_types
    assert t.style == "backmap/fourier"
    assert t.params == pytest.approx(
        [
            3,
            units.energy(0.54392),
            2,
            0.0,
            units.energy(0.30962),
            3,
            180.0,
            units.energy(0.10),
            6,
            60.0,
        ]
    )


def test_explicit_func9_lines_merge(tmp_path: Path) -> None:
    top = _topology(tmp_path)
    mol = _molecule(["CTL2"] * 4)
    mol.dihedrals = [
        TopDihedral(1, 2, 3, 4, 9, [0.0, 1.0, 1]),
        TopDihedral(1, 2, 3, 4, 9, [180.0, 2.0, 2]),
    ]
    system = System()
    dihedrals = _dihedral_terms(system, mol, top, {}, set(), [])
    assert len(dihedrals) == 1
    assert system.dihedral_types[0].params[0] == 2


def test_improper_written_to_data_and_forcefield(tmp_path: Path) -> None:
    top = _topology(tmp_path)
    mol = _molecule(["OBL", "CTL2", "OBL", "CTL2"])
    mol.dihedrals = [TopDihedral(1, 2, 3, 4, 2)]
    system = System()
    assert _dihedral_terms(system, mol, top, {}, set(), []) == []
    system.impropers = _improper_terms(system, mol, top, set())
    (t,) = system.improper_types
    assert t.params == pytest.approx([units.spring_angle(836.8), 0.0])
    data = tmp_path / "x.data"
    write_lammps_data(system, data)
    text = data.read_text()
    assert "1 impropers" in text
    assert "1 improper types" in text
    assert "Impropers\n\n1 1 1 2 3 4" in text
