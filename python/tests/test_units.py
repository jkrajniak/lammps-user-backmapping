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
        # GROMACS E=(k/2)x^2 -> LAMMPS E=Kx^2: unit conversion, then halved.
        expected = (0.239006 / 100.0) / 2.0
        assert units.spring_bond(1.0) == pytest.approx(expected)


class TestSpringAngleConversion:
    def test_spring_angle(self) -> None:
        assert units.spring_angle(1.0) == pytest.approx(0.239006 / 2.0)


class TestGromacsRbConversion:
    def test_sign_flip_and_energy(self) -> None:
        gromacs = [10.0, -5.0, 2.0, 0.0, 0.0, 0.0]
        lammps = units.gromacs_rb_to_lammps(gromacs)
        assert lammps[0] == pytest.approx(units.energy(10.0))
        assert lammps[1] == pytest.approx(-units.energy(-5.0))
        assert lammps[2] == pytest.approx(units.energy(2.0))


class TestLjPairParams:
    def test_opls_combination_rule2_with_fudge(self) -> None:
        sigma, epsilon = units.lj_pair_params(
            0.5,
            0.8,
            0.34,
            0.40,
            combination_rule=2,
            fudge_lj=0.5,
        )
        assert sigma == pytest.approx(0.5 * (units.sigma(0.34) + units.sigma(0.40)))
        assert epsilon == pytest.approx(0.5 * (units.epsilon(0.5) * units.epsilon(0.8)) ** 0.5)

    def test_geometric_combination_rule(self) -> None:
        sigma, epsilon = units.lj_pair_params(
            1.0,
            1.0,
            0.3,
            0.5,
            combination_rule=1,
            fudge_lj=1.0,
        )
        assert sigma == pytest.approx(units.sigma((0.3 * 0.5) ** 0.5))
        assert epsilon == pytest.approx(units.epsilon(1.0))


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
