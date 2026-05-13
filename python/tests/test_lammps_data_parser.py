"""Tests for backmap_prep.parsers.lammps_data_parser."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

from backmap_prep.parsers.lammps_data_parser import parse_lammps_data

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
