"""Tests for backmap_prep.table_converter — XVG to LAMMPS table conversion."""

from __future__ import annotations

import math
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

    def test_wall_energy_capped_for_large_rmin_to_rfloor_gap(self) -> None:
        # Regression test: table_A_A.table-shaped input (r_min=0.02, far from
        # r_floor=1e-4) previously produced an uncapped (r_min/r)**12 wall of
        # ~1.4e31 kcal/mol at the floor -- large enough to poison LAMMPS's
        # pair_style table spline and destabilize a real simulation within a
        # few steps (found running examples/pe-lammps/ on a real LAMMPS
        # build). The wall must now saturate at a bounded ceiling instead.
        r = [0.02, 0.04, 0.06]
        e = [3334.63867905, 3213.94571597, 3097.24728647]
        f = [6129.10296633, 5934.78481442, 5729.64978872]
        r2, e2, f2 = extend_pair_table_to_zero(r, e, f)
        assert r2[0] == pytest.approx(1.0e-4)
        assert e2[0] < 1.0e8, f"wall energy {e2[0]:.3e} not bounded (was ~1.4e31 before the fix)"
        assert e2[0] > e[0]  # still a real repulsive wall, just not exploding
        assert f2[0] > 0.0
        # Continuity at original r_min (the cap only changes the region below r_cap)
        idx = r2.index(0.02)
        assert e2[idx] == pytest.approx(3334.63867905)
        assert f2[idx] == pytest.approx(6129.10296633)

    def test_extension_points_are_log_spaced_not_one_far_knot(self) -> None:
        # Regression test: a single extension point at r_floor=1e-4, jumping
        # straight to r_min=0.02 (a 200x gap versus the table's own uniform
        # 0.02 spacing), was found to corrupt LAMMPS's cubic-spline fit for
        # pair_style table near the low-r region even with the energy cap in
        # place (confirmed by reconstructing the spline directly against the
        # melamine `table_A_A.table` evidence -- see
        # examples/pe-lammps/README.md). Every consecutive pair of r values,
        # including the seam into the table's real data, must now stay within
        # a bounded ratio of each other.
        r = [0.02, 0.04, 0.06]
        e = [3334.63867905, 3213.94571597, 3097.24728647]
        f = [6129.10296633, 5934.78481442, 5729.64978872]
        r2, _e2, _f2 = extend_pair_table_to_zero(r, e, f)
        assert r2[0] == pytest.approx(1.0e-4)
        assert len(r2) > len(r) + 1, "expected more than one extension point for a 200x gap"
        # Only the extension zone (up to and including the seam into r_min=0.02)
        # is bounded by the ratio cap -- the original table's own point spacing
        # beyond r_min (e.g. 0.02 -> 0.04, a 2x ratio near the origin) is
        # untouched, pre-existing IBI data, not something this fix controls.
        extension_zone = [r for r in r2 if r <= 0.02 + 1e-12]
        ratios = [extension_zone[i + 1] / extension_zone[i] for i in range(len(extension_zone) - 1)]
        assert max(ratios) <= 1.2 + 1e-9, f"largest consecutive ratio {max(ratios):.3f} exceeds cap"

    def test_roundtrip_parse_write(self, tmp_path: Path) -> None:
        path = tmp_path / "t.table"
        write_lammps_pair_table(path, [1.0e-4, 0.1], [1e6, 10.0], [1e7, 5.0], source="test")
        r, e, _f = parse_lammps_pair_table(path)
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

    def test_degenerate_force_column_falls_back_to_numerical_gradient(self, tmp_path: Path) -> None:
        """Regression test for the PET/Dacron 2017 ESPResSo++ export bug:
        7-column tables where the force column is all zero but the energy
        column clearly varies (e.g. a decaying repulsive wall)."""
        src = tmp_path / "broken.xvg"
        # r(nm), f, -f', g, -g', h(V, kJ/mol), -h'(F, always 0 -- the bug)
        rows = []
        for i in range(20):
            r_nm = 0.02 + i * 0.02
            v_kj = 1.0e5 * math.exp(-r_nm / 0.05)
            rows.append(f"{r_nm:.6f} 0 0 0 0 {v_kj:.6f} 0")
        src.write_text("\n".join(rows) + "\n")

        dst = tmp_path / "broken.table"
        _convert_xvg(src, dst)
        lines = [
            line
            for line in dst.read_text().splitlines()
            if line and not line.startswith(("#", "E", "N"))
        ]
        forces = [float(line.split()[3]) for line in lines]
        # The wall is repulsive and decaying: F = -dV/dr must be positive
        # (pushes atoms apart) and largest near the steep short-range part.
        assert any(abs(f) > 1.0e-6 for f in forces)
        assert forces[0] > 0.0
        assert forces[0] > forces[-1]

    def test_healthy_force_column_is_not_overwritten(self, xvg_file: Path, tmp_path: Path) -> None:
        """A table with a real, populated force column must be passed
        through unchanged (no accidental re-derivation)."""
        dst = tmp_path / "output.table"
        _convert_xvg(xvg_file, dst)
        lines = dst.read_text().splitlines()
        data_lines = [line for line in lines if line and not line.startswith(("#", "E", "N"))]
        first_data = data_lines[0].split()
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
        assert float(first[3]) == pytest.approx(units.energy(-200.0))

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
        assert float(first[3]) == pytest.approx(units.energy(-200.0))


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

    def test_pair_table_format_is_pure_passthrough(self, tmp_path: Path) -> None:
        # Regression test: pair .table sources used to be re-parsed and run
        # through extend_pair_table_to_zero() even when already native
        # format, which corrupted LAMMPS's pair_style table spline enough to
        # crash a real simulation (see extend_pair_table_to_zero's docs and
        # openspec history). A .table source is already native LAMMPS
        # format -- copy it unchanged, exactly like bond/angle/dihedral
        # .table sources already do, not re-extend it.
        from backmap_prep.builder import System
        from backmap_prep.schema import Settings

        src = tmp_path / "table_A_A.table"
        src.write_text("ENTRY\nN 3\n\n1 0.02 1000.0 5000.0\n2 0.04 100.0 400.0\n3 0.06 10.0 20.0\n")

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
        system.pair_table_files = [("table_A_A.table", "table_1_1.table")]

        result = convert_tables(system, settings, tmp_path)
        assert len(result) == 1
        assert (tmp_path / "table_1_1.table").read_text() == src.read_text()
