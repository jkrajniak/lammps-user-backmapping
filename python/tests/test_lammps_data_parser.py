"""Tests for backmap_prep.parsers.lammps_data_parser."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

import pytest

from backmap_prep.parsers.lammps_data_parser import parse_cg_system, parse_lammps_data

SAMPLE_DATA = """\
LAMMPS data file — test

3 atoms
1 bonds
0 angles
0 dihedrals
0 impropers

2 atom types
1 bond types
0 angle types
0 dihedral types
0 improper types

0.0 50.000000 xlo xhi
0.0 50.000000 ylo yhi
0.0 50.000000 zlo zhi

Masses

1 29.062000 # A (CG)
2 14.027000 # CH2

Atoms # full

1 1 1 0.000000 10.500000 20.300000 30.100000 0 0 0
2 1 2 0.000000 11.000000 21.000000 31.000000 1 0 -1
3 1 2 0.000000 12.000000 22.000000 32.000000

Bonds

1 1 2 3
"""


class TestParseLammpsData:
    def test_parses_atoms(self, tmp_path: Path) -> None:
        p = tmp_path / "test.data"
        p.write_text(SAMPLE_DATA)
        result = parse_lammps_data(p)
        assert len(result.atoms) == 3
        assert result.atoms[0].x == 10.5
        assert result.atoms[0].mol_id == 1
        assert result.atoms[0].type_id == 1

    def test_parses_image_flags(self, tmp_path: Path) -> None:
        p = tmp_path / "test.data"
        p.write_text(SAMPLE_DATA)
        result = parse_lammps_data(p)
        assert result.atoms[1].ix == 1
        assert result.atoms[1].iy == 0
        assert result.atoms[1].iz == -1
        assert result.atoms[2].ix == 0  # default when not present

    def test_parses_bonds(self, tmp_path: Path) -> None:
        p = tmp_path / "test.data"
        p.write_text(SAMPLE_DATA)
        result = parse_lammps_data(p)
        assert len(result.bonds) == 1
        assert result.bonds[0].i == 2
        assert result.bonds[0].j == 3

    def test_parses_box(self, tmp_path: Path) -> None:
        p = tmp_path / "test.data"
        p.write_text(SAMPLE_DATA)
        result = parse_lammps_data(p)
        assert result.box[0] == 50.0
        assert result.box[1] == 50.0
        assert result.box[2] == 50.0

    def test_parses_type_counts(self, tmp_path: Path) -> None:
        p = tmp_path / "test.data"
        p.write_text(SAMPLE_DATA)
        result = parse_lammps_data(p)
        assert result.n_atom_types == 2
        assert result.n_bond_types == 1


CG_SYSTEM_DATA = """\
LAMMPS data file -- pure CG for equilibration

6 atoms
0 bonds
0 angles
0 dihedrals
0 impropers

2 atom types
0 bond types
0 angle types

0.0 57.537600 xlo xhi
0.0 57.537600 ylo yhi
0.0 57.537600 zlo zhi

Masses

1 29.062000
2 28.054000

Atoms # full

1 1 1 0.000000 1.0 2.0 3.0
2 1 2 0.500000 4.0 5.0 6.0
3 1 1 0.000000 7.0 8.0 9.0
4 2 1 0.000000 10.0 11.0 12.0
5 2 2 0.500000 13.0 14.0 15.0
6 2 1 0.000000 16.0 17.0 18.0
"""

CG_SYSTEM_DATA_NO_MASSES = CG_SYSTEM_DATA.replace("Masses\n\n1 29.062000\n2 28.054000\n\n", "")

CG_SYSTEM_DATA_NO_ATOMS = """\
LAMMPS data file -- no atoms section

0.0 10.0 xlo xhi
0.0 10.0 ylo yhi
0.0 10.0 zlo zhi

Masses

1 29.062000
"""

CG_SYSTEM_DATA_UNEVEN = """\
LAMMPS data file -- 5 atoms, not a multiple of the 3-atom first block

0.0 10.0 xlo xhi
0.0 10.0 ylo yhi
0.0 10.0 zlo zhi

Masses

1 29.062000

Atoms # full

1 1 1 0.0 1.0 1.0 1.0
2 1 1 0.0 2.0 2.0 2.0
3 1 1 0.0 3.0 3.0 3.0
4 2 1 0.0 4.0 4.0 4.0
5 2 1 0.0 5.0 5.0 5.0
"""


class TestParseCgSystem:
    def test_parses_box(self, tmp_path: Path) -> None:
        p = tmp_path / "cg_system.data"
        p.write_text(CG_SYSTEM_DATA)
        gro, _top = parse_cg_system(p)
        assert gro.box == (57.5376, 57.5376, 57.5376)

    def test_parses_full_atom_list_positionally(self, tmp_path: Path) -> None:
        p = tmp_path / "cg_system.data"
        p.write_text(CG_SYSTEM_DATA)
        gro, _top = parse_cg_system(p)
        assert len(gro.atoms) == 6
        assert (gro.atoms[0].x, gro.atoms[0].y, gro.atoms[0].z) == (1.0, 2.0, 3.0)
        assert (gro.atoms[5].x, gro.atoms[5].y, gro.atoms[5].z) == (16.0, 17.0, 18.0)

    def test_derives_molecule_template_and_count(self, tmp_path: Path) -> None:
        p = tmp_path / "cg_system.data"
        p.write_text(CG_SYSTEM_DATA)
        _gro, top = parse_cg_system(p)
        assert top.molecules == [("CG", 2)]
        mol = top.molecule_types["CG"]
        assert [a.type for a in mol.atoms] == ["1", "2", "1"]
        assert mol.atoms[1].charge == pytest.approx(0.5)

    def test_type_mass_from_masses_section(self, tmp_path: Path) -> None:
        p = tmp_path / "cg_system.data"
        p.write_text(CG_SYSTEM_DATA)
        _gro, top = parse_cg_system(p)
        mol = top.molecule_types["CG"]
        assert mol.atoms[0].mass == pytest.approx(29.062)
        assert mol.atoms[1].mass == pytest.approx(28.054)
        assert top.atom_types["1"].mass == pytest.approx(29.062)

    def test_missing_masses_section_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "cg_system.data"
        p.write_text(CG_SYSTEM_DATA_NO_MASSES)
        with pytest.raises(ValueError, match="Masses"):
            parse_cg_system(p)

    def test_missing_atoms_section_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "cg_system.data"
        p.write_text(CG_SYSTEM_DATA_NO_ATOMS)
        with pytest.raises(ValueError, match="Atoms"):
            parse_cg_system(p)

    def test_uneven_molecule_grouping_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "cg_system.data"
        p.write_text(CG_SYSTEM_DATA_UNEVEN)
        with pytest.raises(ValueError, match="not a multiple"):
            parse_cg_system(p)
