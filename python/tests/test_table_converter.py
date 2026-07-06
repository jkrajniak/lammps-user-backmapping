"""Tests for backmap_prep.table_converter — XVG to LAMMPS table conversion."""

from __future__ import annotations

from textwrap import dedent
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

import pytest

from backmap_prep import units
from backmap_prep.table_converter import (
    _convert_angle_xvg,
    _convert_dihedral_xvg,
    _convert_xvg,
    convert_tables,
)


class TestConvertXvg:
    def test_basic_conversion(self, xvg_file: Path, tmp_path: Path) -> None:
        dst = tmp_path / "output.table"
        _convert_xvg(xvg_file, dst)
        assert dst.exists()

    def test_output_format(self, xvg_file: Path, tmp_path: Path) -> None:
        dst = tmp_path / "output.table"
        _convert_xvg(xvg_file, dst)
        content = dst.read_text()
        assert "ENTRY" in content
        assert "N 3" in content

    def test_unit_conversion_applied(self, xvg_file: Path, tmp_path: Path) -> None:
        dst = tmp_path / "output.table"
        _convert_xvg(xvg_file, dst)
        lines = dst.read_text().splitlines()

        data_lines = [line for line in lines if line and not line.startswith(("#", "E", "N"))]
        first_data = data_lines[0].split()
        assert first_data[0] == "1"
        assert float(first_data[1]) == pytest.approx(units.distance(0.10))
        assert float(first_data[2]) == pytest.approx(units.energy(10.0))
        assert float(first_data[3]) == pytest.approx(units.force(-200.0))

    def test_comments_skipped(self, tmp_path: Path) -> None:
        content = dedent("""\
            # comment
            @ grace metadata
            @TYPE xy
            0.10  1.0  -10.0
        """)
        src = tmp_path / "data.xvg"
        src.write_text(content)
        dst = tmp_path / "data.table"
        _convert_xvg(src, dst)
        assert "N 1" in dst.read_text()

    def test_empty_xvg_raises(self, tmp_path: Path) -> None:
        src = tmp_path / "empty.xvg"
        src.write_text("# only comments\n@ metadata\n")
        dst = tmp_path / "empty.table"
        with pytest.raises(ValueError, match="No data found"):
            _convert_xvg(src, dst)


class TestConvertAngleXvg:
    def test_angle_degrees_preserved(self, tmp_path: Path) -> None:
        src = tmp_path / "table_a1.xvg"
        src.write_text("0.0  10.0  -200.0\n90.0  5.0  -100.0\n")
        dst = tmp_path / "table_a1.table"
        _convert_angle_xvg(src, dst)
        lines = [ln for ln in dst.read_text().splitlines() if ln and ln[0].isdigit()]
        first = lines[0].split()
        assert float(first[1]) == pytest.approx(0.0)
        assert float(first[2]) == pytest.approx(units.energy(10.0))
        assert float(first[3]) == pytest.approx(units.angular_force(-200.0))

    def test_angle_table_via_convert_tables(self, tmp_path: Path) -> None:
        from backmap_prep.builder import System
        from backmap_prep.schema import Settings

        xvg = tmp_path / "table_a1.xvg"
        xvg.write_text("0.0  1.0  -10.0\n180.0  0.5  -5.0\n")

        settings = Settings(
            molecules=[
                {
                    "name": "M",
                    "source": {"coordinates": "a.gro", "topology": "a.top"},
                    "beads": [{"name": "B", "type": "C", "atoms": ["X"]}],
                }
            ],
            cg_system={"coordinates": "c.gro", "topology": "c.top"},
        )
        system = System()
        system.angle_table_files = [("table_a1.xvg", "table_a1.table")]

        result = convert_tables(system, settings, tmp_path)
        assert len(result) == 1
        content = result[0].read_text()
        assert "180.00000000" in content


class TestConvertDihedralXvg:
    def test_dihedral_degrees_preserved(self, tmp_path: Path) -> None:
        src = tmp_path / "table_d1.xvg"
        src.write_text("-180.0  10.0  -200.0\n0.0  5.0  -100.0\n")
        dst = tmp_path / "table_d1.table"
        _convert_dihedral_xvg(src, dst)
        lines = [ln for ln in dst.read_text().splitlines() if ln and ln[0].isdigit()]
        first = lines[0].split()
        assert float(first[1]) == pytest.approx(-180.0)
        assert float(first[2]) == pytest.approx(units.energy(10.0))
        assert float(first[3]) == pytest.approx(units.angular_force(-200.0))


class TestConvertTables:
    def test_no_table_files(self, tmp_path: Path) -> None:
        from backmap_prep.builder import System
        from backmap_prep.schema import Settings

        settings_dict = {
            "molecules": [
                {
                    "name": "M",
                    "source": {"coordinates": "a.gro", "topology": "a.top"},
                    "beads": [{"name": "B", "type": "C", "atoms": ["X"]}],
                }
            ],
            "cg_system": {"coordinates": "c.gro", "topology": "c.top"},
        }
        settings = Settings(**settings_dict)
        system = System()
        result = convert_tables(system, settings, tmp_path)
        assert result == []

    def test_xvg_table_converted(self, tmp_path: Path) -> None:
        from backmap_prep.builder import System
        from backmap_prep.schema import Settings

        xvg = tmp_path / "bond.xvg"
        xvg.write_text("0.10  1.0  -10.0\n0.20  2.0  -20.0\n")

        settings_dict = {
            "molecules": [
                {
                    "name": "M",
                    "source": {"coordinates": "a.gro", "topology": "a.top"},
                    "beads": [{"name": "B", "type": "C", "atoms": ["X"]}],
                }
            ],
            "cg_system": {"coordinates": "c.gro", "topology": "c.top"},
        }
        settings = Settings(**settings_dict)
        system = System()
        system.table_files = [("bond.xvg", "bond.table")]

        result = convert_tables(system, settings, tmp_path)
        assert len(result) == 1
        assert result[0].name == "bond.table"

    def test_missing_source_skipped(self, tmp_path: Path) -> None:
        from backmap_prep.builder import System
        from backmap_prep.schema import Settings

        settings_dict = {
            "molecules": [
                {
                    "name": "M",
                    "source": {"coordinates": "a.gro", "topology": "a.top"},
                    "beads": [{"name": "B", "type": "C", "atoms": ["X"]}],
                }
            ],
            "cg_system": {"coordinates": "c.gro", "topology": "c.top"},
        }
        settings = Settings(**settings_dict)
        system = System()
        system.table_files = [("missing.xvg", "missing.table")]

        result = convert_tables(system, settings, tmp_path)
        assert result == []

    def test_table_format_copied(self, tmp_path: Path) -> None:
        from backmap_prep.builder import System
        from backmap_prep.schema import Settings

        src = tmp_path / "existing.table"
        src.write_text("ENTRY\nN 2\n\n1 1.0 0.5 -0.1\n2 2.0 0.3 -0.05\n")

        settings_dict = {
            "molecules": [
                {
                    "name": "M",
                    "source": {"coordinates": "a.gro", "topology": "a.top"},
                    "beads": [{"name": "B", "type": "C", "atoms": ["X"]}],
                }
            ],
            "cg_system": {"coordinates": "c.gro", "topology": "c.top"},
        }
        settings = Settings(**settings_dict)
        system = System()
        system.table_files = [("existing.table", "copy.table")]

        result = convert_tables(system, settings, tmp_path)
        assert len(result) == 1
        assert (tmp_path / "copy.table").read_text() == src.read_text()
