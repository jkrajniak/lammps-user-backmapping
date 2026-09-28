"""Tests for backmap_prep.parsers.lammps_script_parser."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

import math

import pytest

from backmap_prep.parsers.lammps_script_parser import opls_to_ryckaert, parse_at_fragment_script

BOUNDED_SCRIPT = """\
units real
atom_style full

pair_style lj/cut 14.0
pair_coeff 1 1 0.207 3.748
pair_coeff 1 2 0.18 3.8
pair_coeff 2 2 0.118 3.905

bond_style harmonic
bond_coeff 1 800.0 1.53

angle_style harmonic
angle_coeff 1 126.6 111.0

dihedral_style ryckaert
dihedral_coeff 1 9.27 12.15 -13.12 -3.06 26.24 -13.49

neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes
special_bonds lj 0.0 0.0 0.0 coul 0.0 0.0 0.0
fix integrate all nvt temp 298.0 298.0 100.0
thermo 1000
run 10000
write_data out.data
"""


class TestParseAtFragmentScript:
    def test_parses_bond_coefficients(self, tmp_path: Path) -> None:
        p = tmp_path / "in.frag"
        p.write_text(BOUNDED_SCRIPT)
        coeffs = parse_at_fragment_script(p)
        assert coeffs.bond == {1: (800.0, 1.53)}

    def test_parses_angle_coefficients(self, tmp_path: Path) -> None:
        p = tmp_path / "in.frag"
        p.write_text(BOUNDED_SCRIPT)
        coeffs = parse_at_fragment_script(p)
        assert coeffs.angle == {1: (126.6, 111.0)}

    def test_parses_dihedral_coefficients(self, tmp_path: Path) -> None:
        p = tmp_path / "in.frag"
        p.write_text(BOUNDED_SCRIPT)
        coeffs = parse_at_fragment_script(p)
        assert coeffs.dihedral == {1: [9.27, 12.15, -13.12, -3.06, 26.24, -13.49]}

    def test_diagonal_pair_coefficients_only(self, tmp_path: Path) -> None:
        p = tmp_path / "in.frag"
        p.write_text(BOUNDED_SCRIPT)
        coeffs = parse_at_fragment_script(p)
        assert coeffs.pair == {1: (0.207, 3.748), 2: (0.118, 3.905)}

    def test_ignores_unrelated_commands(self, tmp_path: Path) -> None:
        p = tmp_path / "in.frag"
        p.write_text(BOUNDED_SCRIPT)
        # Should not raise despite fix/run/thermo/special_bonds/neighbor present.
        parse_at_fragment_script(p)

    def test_missing_units_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "in.frag"
        p.write_text("bond_style harmonic\nbond_coeff 1 1.0 1.0\n")
        with pytest.raises(ValueError, match="units real"):
            parse_at_fragment_script(p)

    def test_wrong_units_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "in.frag"
        p.write_text("units metal\nbond_style harmonic\nbond_coeff 1 1.0 1.0\n")
        with pytest.raises(ValueError, match="units real"):
            parse_at_fragment_script(p)

    def test_unsupported_bond_style_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "in.frag"
        p.write_text("units real\nbond_style morse\n")
        with pytest.raises(ValueError, match="unsupported bond_style"):
            parse_at_fragment_script(p)

    def test_unsupported_angle_style_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "in.frag"
        p.write_text("units real\nangle_style cosine\n")
        with pytest.raises(ValueError, match="unsupported angle_style"):
            parse_at_fragment_script(p)

    def test_unsupported_dihedral_style_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "in.frag"
        p.write_text("units real\ndihedral_style quadratic\n")
        with pytest.raises(ValueError, match="unsupported dihedral_style"):
            parse_at_fragment_script(p)

    def test_improper_style_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "in.frag"
        p.write_text("units real\nimproper_style harmonic\n")
        with pytest.raises(ValueError, match="unsupported style directive"):
            parse_at_fragment_script(p)

    def test_pair_style_any_value_accepted(self, tmp_path: Path) -> None:
        p = tmp_path / "in.frag"
        p.write_text("units real\npair_style lj/cut/coul/long 10.0\npair_coeff 1 1 0.1 3.0\n")
        coeffs = parse_at_fragment_script(p)
        assert coeffs.pair == {1: (0.1, 3.0)}


def _opls_energy(k: tuple[float, float, float, float], phi: float) -> float:
    k1, k2, k3, k4 = k
    return 0.5 * (
        k1 * (1 + math.cos(phi))
        + k2 * (1 - math.cos(2 * phi))
        + k3 * (1 + math.cos(3 * phi))
        + k4 * (1 - math.cos(4 * phi))
    )


class TestOplsDihedral:
    @pytest.mark.parametrize(
        "k", [(1.3, -0.05, 0.2, 0.1), (0.0, 0.0, 0.3, 0.0), (1.74, -0.157, 0.279, 0.0)]
    )
    def test_ryckaert_equals_opls_energy(self, k: tuple[float, float, float, float]) -> None:
        c = opls_to_ryckaert(*k)
        for step in range(73):
            phi = math.radians(-180.0 + 5.0 * step)
            rb = sum(cn * math.cos(phi) ** n for n, cn in enumerate(c))
            assert rb == pytest.approx(_opls_energy(k, phi), abs=1e-12)

    def test_trans_is_zero(self) -> None:
        c = opls_to_ryckaert(1.3, -0.05, 0.2, 0.1)
        assert sum(cn * (-1.0) ** n for n, cn in enumerate(c)) == pytest.approx(0.0, abs=1e-12)

    def test_parses_opls_fragment(self, tmp_path: Path) -> None:
        p = tmp_path / "in.frag"
        p.write_text("units real\ndihedral_style opls\ndihedral_coeff 2 1.3 -0.05 0.2 0.1\n")
        coeffs = parse_at_fragment_script(p)
        assert coeffs.dihedral == {2: opls_to_ryckaert(1.3, -0.05, 0.2, 0.1)}

    def test_malformed_opls_line_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "in.frag"
        p.write_text("units real\ndihedral_style opls\ndihedral_coeff 1 1.3 -0.05 0.2\n")
        with pytest.raises(ValueError, match="malformed dihedral_coeff"):
            parse_at_fragment_script(p)
