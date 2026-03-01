"""Tests for backmap_prep.builder — system construction from settings and source files."""

from __future__ import annotations

from textwrap import dedent
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

import pytest

from backmap_prep.builder import (
    AtomTypeInfo,
    BondTypeInfo,
    LammpsAngle,
    LammpsAtom,
    LammpsBond,
    PairTypeInfo,
    System,
    build_system,
)
from backmap_prep.schema import Settings


def _write_cg_files(base: Path, n_mols: int = 1) -> None:
    """Write minimal CG .gro and .top files for n_mols identical molecules."""
    atoms = []
    for m in range(n_mols):
        resid = m + 1
        atoms.append(f"{resid:5d}MOL   B1{m * 2 + 1:5d}   0.500   0.500   0.500")
        atoms.append(f"{resid:5d}MOL   B2{m * 2 + 2:5d}   0.600   0.600   0.600")

    gro = f"CG system\n{len(atoms):5d}\n" + "\n".join(atoms) + "\n   5.00000   5.00000   5.00000\n"
    (base / "cg.gro").write_text(gro)

    top = dedent("""\
        [ defaults ]
        1  2  yes  0.5  0.8333

        [ atomtypes ]
        CGA  72.0  0.0  V  0.47  3.5
        CGB  72.0  0.0  V  0.47  3.5

        [ moleculetype ]
        MOL  3

        [ atoms ]
        1  CGA  1  MOL  B1  1  0.0  72.0
        2  CGB  1  MOL  B2  1  0.0  72.0

        [ bonds ]
        1  2  1  0.350  5000.0

        [ molecules ]
    """)
    top += f"MOL {n_mols}\n"
    (base / "cg.top").write_text(top)


def _write_at_files(base: Path) -> None:
    """Write minimal AT .gro and .top files for a single-molecule template."""
    gro = dedent("""\
        AT template
            4
            1MOL  C1    1   0.500   0.500   0.500
            1MOL  C2    2   0.510   0.500   0.500
            1MOL  C3    3   0.600   0.600   0.600
            1MOL  C4    4   0.610   0.600   0.600
           5.00000   5.00000   5.00000
    """)
    (base / "at.gro").write_text(gro)

    top = dedent("""\
        [ defaults ]
        1  2  yes  0.5  0.8333

        [ atomtypes ]
        CH2  14.0  0.0  A  0.395  0.382
        CH3  15.0  0.0  A  0.395  0.815

        [ moleculetype ]
        TestMol  3

        [ atoms ]
        1  CH2  1  MOL  C1  1  0.0  14.0
        2  CH3  1  MOL  C2  1  0.0  15.0
        3  CH2  1  MOL  C3  1  0.0  14.0
        4  CH3  1  MOL  C4  1  0.0  15.0

        [ bonds ]
        1  2  1  0.154  200000.0
        3  4  1  0.154  200000.0

        [ angles ]
        1  2  3  1  111.0  500.0
    """)
    (base / "at.top").write_text(top)


def _make_settings(
    base: Path,
    *,
    cross_bonds: list | None = None,
    cross_angles: list | None = None,
) -> Settings:
    """Build a Settings object pointing at files in base."""
    d: dict = {
        "molecules": [
            {
                "name": "TestMol",
                "source": {"coordinates": "at.gro", "topology": "at.top"},
                "beads": [
                    {"name": "B1", "type": "CGA", "atoms": ["C1", "C2"]},
                    {"name": "B2", "type": "CGB", "atoms": ["C3", "C4"]},
                ],
            }
        ],
        "cg_system": {"coordinates": "cg.gro", "topology": "cg.top"},
    }
    if cross_bonds or cross_angles:
        d["cross_interactions"] = {}
        if cross_bonds:
            d["cross_interactions"]["bonds"] = cross_bonds
        if cross_angles:
            d["cross_interactions"]["angles"] = cross_angles
    return Settings(**d)


class TestDataclasses:
    def test_lammps_atom(self) -> None:
        a = LammpsAtom(atom_id=1, mol_id=1, type_id=2, charge=0.5, x=1.0, y=2.0, z=3.0)
        assert a.is_cg is False
        assert a.type_name == ""

    def test_lammps_bond(self) -> None:
        b = LammpsBond(bond_id=1, type_id=1, i=1, j=2)
        assert b.bond_id == 1

    def test_lammps_angle(self) -> None:
        a = LammpsAngle(angle_id=1, type_id=1, i=1, j=2, k=3)
        assert a.k == 3

    def test_system_defaults(self) -> None:
        s = System()
        assert s.atoms == []
        assert s.bonds == []
        assert s.box == (0.0, 0.0, 0.0)
        assert s.cg_type_id == 0
        assert not s.has_cross_bonds

    def test_bond_type_info(self) -> None:
        bt = BondTypeInfo(type_id=1, style="harmonic", keyword="", params=[100.0, 1.54])
        assert bt.table_file is None

    def test_atom_type_info(self) -> None:
        at = AtomTypeInfo(type_id=1, name="CH2", mass=14.0, is_cg=False, sigma=3.95, epsilon=0.382)
        assert at.sigma == 3.95

    def test_pair_type_info(self) -> None:
        pt = PairTypeInfo(itype=1, jtype=2, kind="none")
        assert pt.table_file is None


class TestBuildSystem:
    def test_basic_build(self, tmp_path: Path) -> None:
        _write_cg_files(tmp_path, n_mols=1)
        _write_at_files(tmp_path)
        settings = _make_settings(tmp_path)
        system = build_system(settings, tmp_path)

        assert len(system.atoms) > 0
        assert len(system.atom_types) > 0
        assert system.box[0] == pytest.approx(50.0)

    def test_cg_types_first(self, tmp_path: Path) -> None:
        _write_cg_files(tmp_path)
        _write_at_files(tmp_path)
        settings = _make_settings(tmp_path)
        system = build_system(settings, tmp_path)

        cg_types = [t for t in system.atom_types if t.is_cg]
        at_types = [t for t in system.atom_types if not t.is_cg]
        assert all(c.type_id < a.type_id for c in cg_types for a in at_types)

    def test_cg_type_id_set(self, tmp_path: Path) -> None:
        _write_cg_files(tmp_path)
        _write_at_files(tmp_path)
        settings = _make_settings(tmp_path)
        system = build_system(settings, tmp_path)

        assert system.cg_type_id == 1

    def test_atom_count_per_molecule(self, tmp_path: Path) -> None:
        _write_cg_files(tmp_path, n_mols=1)
        _write_at_files(tmp_path)
        settings = _make_settings(tmp_path)
        system = build_system(settings, tmp_path)

        cg_atoms = [a for a in system.atoms if a.is_cg]
        at_atoms = [a for a in system.atoms if not a.is_cg]
        assert len(cg_atoms) == 2
        assert len(at_atoms) == 4

    def test_multiple_molecules(self, tmp_path: Path) -> None:
        _write_cg_files(tmp_path, n_mols=3)
        _write_at_files(tmp_path)
        settings = _make_settings(tmp_path)
        system = build_system(settings, tmp_path)

        mol_ids = {a.mol_id for a in system.atoms}
        assert mol_ids == {1, 2, 3}

    def test_intra_bonds_created(self, tmp_path: Path) -> None:
        _write_cg_files(tmp_path)
        _write_at_files(tmp_path)
        settings = _make_settings(tmp_path)
        system = build_system(settings, tmp_path)

        assert len(system.bonds) > 0
        harmonic_types = [bt for bt in system.bond_types if bt.style == "harmonic"]
        assert len(harmonic_types) > 0

    def test_cross_bonds(self, tmp_path: Path) -> None:
        _write_cg_files(tmp_path)
        _write_at_files(tmp_path)
        settings = _make_settings(
            tmp_path,
            cross_bonds=[
                {
                    "params": "1 0.154 200000.0",
                    "pairs": [["C2", "C3"]],
                }
            ],
        )
        system = build_system(settings, tmp_path)
        assert system.has_cross_bonds

        backmap_bond_types = [bt for bt in system.bond_types if bt.style == "backmap/harmonic"]
        assert len(backmap_bond_types) > 0

    def test_pair_type_classification(self, tmp_path: Path) -> None:
        _write_cg_files(tmp_path)
        _write_at_files(tmp_path)
        settings = _make_settings(tmp_path)
        system = build_system(settings, tmp_path)

        kinds = {pt.kind for pt in system.pair_types}
        assert "cg" in kinds
        assert "atomistic" in kinds
        assert "none" in kinds

    def test_missing_cg_molecule_raises(self, tmp_path: Path) -> None:
        gro = "Empty\n    0\n   5.00000   5.00000   5.00000\n"
        (tmp_path / "cg.gro").write_text(gro)
        top = dedent("""\
            [ moleculetype ]
            OTHER  3
            [ atoms ]
            1  X  1  OTHER  A1  1  0.0  10.0
            [ molecules ]
        """)
        (tmp_path / "cg.top").write_text(top)
        _write_at_files(tmp_path)
        settings = _make_settings(tmp_path)
        with pytest.raises(ValueError, match="No molecules found"):
            build_system(settings, tmp_path)
