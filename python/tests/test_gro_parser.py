"""Tests for backmap_prep.parsers.gro_parser — GROMACS .gro file parsing."""

from __future__ import annotations

from textwrap import dedent
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

import pytest

from backmap_prep.parsers.gro_parser import GroAtom, GroFile, parse_gro


class TestParseGro:
    def test_parse_basic(self, gro_file: Path) -> None:
        result = parse_gro(gro_file)
        assert isinstance(result, GroFile)
        assert result.title == "Test system"
        assert len(result.atoms) == 3

    def test_atom_fields(self, gro_file: Path) -> None:
        result = parse_gro(gro_file)
        a1 = result.atoms[0]
        assert a1.resid == 1
        assert a1.resname == "RES"
        assert a1.name == "C1"
        assert a1.index == 1
        assert a1.x == pytest.approx(0.100)
        assert a1.y == pytest.approx(0.200)
        assert a1.z == pytest.approx(0.300)

    def test_velocities_present(self, gro_file: Path) -> None:
        result = parse_gro(gro_file)
        a1 = result.atoms[0]
        assert a1.vx == pytest.approx(0.01)
        assert a1.vy == pytest.approx(0.02)
        assert a1.vz == pytest.approx(0.03)

    def test_velocities_absent(self, gro_file: Path) -> None:
        result = parse_gro(gro_file)
        a2 = result.atoms[1]
        assert a2.vx == pytest.approx(0.0)
        assert a2.vy == pytest.approx(0.0)
        assert a2.vz == pytest.approx(0.0)

    def test_box_dimensions(self, gro_file: Path) -> None:
        result = parse_gro(gro_file)
        assert result.box == pytest.approx((5.0, 5.0, 5.0))

    def test_invalid_file_too_short(self, tmp_path: Path) -> None:
        p = tmp_path / "bad.gro"
        p.write_text("Title\n")
        with pytest.raises(ValueError, match=r"Invalid .gro file"):
            parse_gro(p)

    def test_single_atom(self, tmp_path: Path) -> None:
        content = dedent("""\
            Single atom
                1
                1MOL  A1    1   1.234   5.678   9.012
               3.00000   4.00000   5.00000
        """)
        p = tmp_path / "single.gro"
        p.write_text(content)
        result = parse_gro(p)
        assert len(result.atoms) == 1
        assert result.atoms[0].name == "A1"
        assert result.box == pytest.approx((3.0, 4.0, 5.0))


class TestGroAtomDataclass:
    def test_defaults(self) -> None:
        a = GroAtom(resid=1, resname="RES", name="C1", index=1, x=0.0, y=0.0, z=0.0)
        assert a.vx == 0.0
        assert a.vy == 0.0
        assert a.vz == 0.0
