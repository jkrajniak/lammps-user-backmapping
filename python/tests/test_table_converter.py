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
    extend_pair_table_to_zero,
    parse_lammps_pair_table,
    write_lammps_pair_table,
)


class TestExtendPairTableToZero:
    def test_already_at_floor_unchanged(self) -> None:
        r = [1.0e-4, 0.1, 0.2]
        e = [1e6, 100.0, 10.0]
        f = [1e7, 1000.0, 50.0]
        r2, e2, f2 = extend_pair_table_to_zero(r, e, f)
        assert r2 == r
        assert e2 == e
        assert f2 == f

    def test_bumps_illegal_zero_lower_bound(self) -> None:
        r = [0.0, 0.1, 0.2]
        e = [1e6, 100.0, 10.0]
        f = [1e7, 1000.0, 50.0]
        r2, e2, f2 = extend_pair_table_to_zero(r, e, f)
        assert r2[0] == pytest.approx(1.0e-4)
        assert e2 == e
        assert f2 == f

    def test_extends_below_rmin_with_steep_wall(self) -> None:
        r = [0.20, 0.40, 0.60]
        e = [1000.0, 100.0, 10.0]
        f = [5000.0, 400.0, 20.0]
        r2, e2, f2 = extend_pair_table_to_zero(r, e, f)
        assert r2[0] == pytest.approx(1.0e-4)
        assert r2[0] > 0.0
        assert r2[-3:] == r
        assert e2[0] > e[0]
        assert f2[0] > 0.0
        # Continuity at original r_min
        idx = r2.index(0.20)
        assert e2[idx] == pytest.approx(1000.0)
        assert f2[idx] == pytest.approx(5000.0)

    def test_roundtrip_parse_write(self, tmp_path: Path) -> None:
        path = tmp_path / "t.table"
        write_lammps_pair_table(path, [1.0e-4, 0.1], [1e6, 10.0], [1e7, 5.0], source="test")
        r, e, f = parse_lammps_pair_table(path)
        assert r[0] == pytest.approx(1.0e-4)
        assert e[1] == pytest.approx(10.0)


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
        _convert_xvg(xvg_file, dst, extend_core=True)
        lines = dst.read_text().splitlines()

        data_lines = [line for line in lines if line and not line.startswith(("#", "E", "N"))]
        # Pair tables are extended to r=0; find the original first XVG point (0.10 nm → 1 Å)
        match = None
        for line in data_lines:
            parts = line.split()
            if abs(float(parts[1]) - units.distance(0.10)) < 1e-9:
                match = parts
                break
        assert match is not None
        assert float(match[2]) == pytest.approx(units.energy(10.0))
        assert float(match[3]) == pytest.approx(units.force(-200.0))
        assert float(data_lines[0].split()[1]) == pytest.approx(1.0e-4)

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
