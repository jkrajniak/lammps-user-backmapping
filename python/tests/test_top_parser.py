"""Tests for backmap_prep.parsers.top_parser — GROMACS topology parsing."""

from __future__ import annotations

from textwrap import dedent
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

import pytest

from backmap_prep.parsers.top_parser import (
    Topology,
    parse_top,
)


class TestParseTop:
    def test_parse_basic(self, top_file: Path) -> None:
        top = parse_top(top_file)
        assert isinstance(top, Topology)
        assert "MOL" in top.molecule_types
        assert len(top.molecules) == 1

    def test_defaults_section(self, top_file: Path) -> None:
        top = parse_top(top_file)
        assert top.combination_rule == 2
        assert top.fudge_lj == pytest.approx(0.5)
        assert top.fudge_qq == pytest.approx(0.8333)

    def test_atomtypes(self, top_file: Path) -> None:
        top = parse_top(top_file)
        assert "CG1" in top.atom_types
        assert "AT1" in top.atom_types
        cg = top.atom_types["CG1"]
        assert cg.mass == pytest.approx(72.0)
        assert cg.ptype == "A"
        assert cg.sigma == pytest.approx(0.47)
        assert cg.epsilon == pytest.approx(3.5)

    def test_molecule_atoms(self, top_file: Path) -> None:
        top = parse_top(top_file)
        mol = top.molecule_types["MOL"]
        assert len(mol.atoms) == 3
        assert mol.atoms[0].type == "CG1"
        assert mol.atoms[0].name == "B1"
        assert mol.atoms[1].mass == pytest.approx(12.0)

    def test_bonds(self, top_file: Path) -> None:
        top = parse_top(top_file)
        mol = top.molecule_types["MOL"]
        assert len(mol.bonds) == 1
        b = mol.bonds[0]
        assert b.i == 2
        assert b.j == 3
        assert b.func == 1
        assert b.params == pytest.approx([0.154, 200000.0])

    def test_angles(self, top_file: Path) -> None:
        top = parse_top(top_file)
        mol = top.molecule_types["MOL"]
        assert len(mol.angles) == 1
        a = mol.angles[0]
        assert a.i == 1
        assert a.j == 2
        assert a.k == 3
        assert a.params == pytest.approx([120.0, 500.0])

    def test_molecules_count(self, top_file: Path) -> None:
        top = parse_top(top_file)
        assert top.molecules[0] == ("MOL", 1)

    def test_nrexcl(self, top_file: Path) -> None:
        top = parse_top(top_file)
        assert top.molecule_types["MOL"].nrexcl == 3


class TestParseTopInclude:
    def test_include_directive(self, tmp_path: Path) -> None:
        itp = tmp_path / "frag.itp"
        itp.write_text(
            dedent("""\
            [ atomtypes ]
            X1  10.0  0.0  A  0.30  0.5
        """)
        )
        main = tmp_path / "main.top"
        main.write_text(
            dedent("""\
            #include "frag.itp"

            [ moleculetype ]
            M1  3

            [ atoms ]
            1  X1  1  M1  A1  1  0.0  10.0

            [ molecules ]
            M1 2
        """)
        )
        top = parse_top(main)
        assert "X1" in top.atom_types
        assert top.molecules[0] == ("M1", 2)

    def test_include_dirs(self, tmp_path: Path) -> None:
        inc_dir = tmp_path / "includes"
        inc_dir.mkdir()
        itp = inc_dir / "types.itp"
        itp.write_text(
            dedent("""\
            [ atomtypes ]
            Y1  14.0  0.0  A  0.35  0.6
        """)
        )
        main = tmp_path / "main.top"
        main.write_text(
            dedent("""\
            #include "types.itp"
            [ moleculetype ]
            M2  2
            [ atoms ]
            1  Y1  1  M2  A1  1  0.0  14.0
            [ molecules ]
            M2 1
        """)
        )
        top = parse_top(main, include_dirs=[inc_dir])
        assert "Y1" in top.atom_types


class TestPreprocessor:
    def test_ifdef_skipped(self, tmp_path: Path) -> None:
        content = dedent("""\
            [ atomtypes ]
            R1  16.0  0.0  A  0.30  0.5
            #ifdef POSRES
            SHOULD_NOT_APPEAR  99.0  0.0  A  0.0  0.0
            #endif
            [ moleculetype ]
            MX 3
            [ atoms ]
            1  R1  1  MX  A1  1  0.0  16.0
            [ molecules ]
            MX 1
        """)
        p = tmp_path / "test.top"
        p.write_text(content)
        top = parse_top(p)
        assert "SHOULD_NOT_APPEAR" not in top.atom_types
        assert "R1" in top.atom_types

    def test_inline_comments_stripped(self, tmp_path: Path) -> None:
        content = dedent("""\
            [ atomtypes ]
            Z1  20.0  0.0  A  0.40  1.0 ; this is a comment

            [ moleculetype ]
            MC  3

            [ atoms ]
            1  Z1  1  MC  A1  1  0.0  20.0

            [ molecules ]
            MC  1
        """)
        p = tmp_path / "test.top"
        p.write_text(content)
        top = parse_top(p)
        assert "Z1" in top.atom_types
        assert top.atom_types["Z1"].epsilon == pytest.approx(1.0)


class TestDihedrals:
    def test_dihedrals_parsed(self, tmp_path: Path) -> None:
        content = dedent("""\
            [ moleculetype ]
            DH 3
            [ atoms ]
            1  C  1  DH  A1  1  0.0  12.0
            2  C  1  DH  A2  1  0.0  12.0
            3  C  1  DH  A3  1  0.0  12.0
            4  C  1  DH  A4  1  0.0  12.0
            [ dihedrals ]
            1  2  3  4  1  180.0  10.0  2.0
            [ molecules ]
            DH 1
        """)
        p = tmp_path / "test.top"
        p.write_text(content)
        top = parse_top(p)
        mol = top.molecule_types["DH"]
        assert len(mol.dihedrals) == 1
        d = mol.dihedrals[0]
        assert d.i == 1
        assert d.j == 2
        assert d.k == 3
        assert d.atom_l == 4
        assert d.func == 1
        assert d.params == pytest.approx([180.0, 10.0, 2.0])


class TestMissingFile:
    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="Topology file not found"):
            parse_top(tmp_path / "nonexistent.top")
