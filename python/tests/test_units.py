"""Tests for backmap_prep.units — GROMACS → LAMMPS real unit conversions."""

from __future__ import annotations

import pytest

from backmap_prep import units


class TestDistanceConversion:
    def test_one_nm_to_angstrom(self) -> None:
        assert units.distance(1.0) == pytest.approx(10.0)

    def test_zero(self) -> None:
        assert units.distance(0.0) == 0.0

    def test_fractional(self) -> None:
        assert units.distance(0.154) == pytest.approx(1.54)


class TestEnergyConversion:
    def test_one_kj_to_kcal(self) -> None:
        assert units.energy(1.0) == pytest.approx(0.239006)

    def test_zero(self) -> None:
        assert units.energy(0.0) == 0.0

    def test_negative(self) -> None:
        assert units.energy(-4.184) == pytest.approx(-4.184 * 0.239006)


class TestForceConversion:
    def test_force(self) -> None:
        expected = 0.239006 / 10.0
        assert units.force(1.0) == pytest.approx(expected)

    def test_zero(self) -> None:
        assert units.force(0.0) == 0.0


class TestTimeConversion:
    def test_one_ps_to_fs(self) -> None:
        assert units.time(1.0) == pytest.approx(1000.0)

    def test_fractional(self) -> None:
        assert units.time(0.001) == pytest.approx(1.0)


class TestSpringBondConversion:
    def test_spring_bond(self) -> None:
        expected = 0.239006 / 100.0
        assert units.spring_bond(1.0) == pytest.approx(expected)


class TestSpringAngleConversion:
    def test_spring_angle(self) -> None:
        assert units.spring_angle(1.0) == pytest.approx(0.239006)


class TestSigmaEpsilon:
    def test_sigma_is_distance(self) -> None:
        assert units.sigma(0.34) == pytest.approx(3.4)

    def test_epsilon_is_energy(self) -> None:
        assert units.epsilon(1.0) == pytest.approx(0.239006)


@pytest.mark.parametrize(
    ("func", "input_val", "expected"),
    [
        (units.distance, 2.5, 25.0),
        (units.energy, 10.0, 2.39006),
        (units.time, 0.5, 500.0),
        (units.sigma, 0.5, 5.0),
        (units.epsilon, 0.5, 0.119503),
    ],
)
def test_parametrized_conversions(func: object, input_val: float, expected: float) -> None:
    assert func(input_val) == pytest.approx(expected)
