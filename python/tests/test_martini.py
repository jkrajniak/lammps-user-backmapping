"""MARTINI CG tables (pair LJ + reaction field, G96 angles) and split exclusions."""

from __future__ import annotations

import math

import pytest

from backmap_prep.martini import (
    ONE_4PI_EPS0,
    cg_type_charges,
    g96_angle_xvg,
    lj_sigma_epsilon,
    pair_xvg,
)
from backmap_prep.parsers.top_parser import AtomType, MoleculeType, TopAtom, Topology
from backmap_prep.schema import CGNonbonded, SimulationParams
from backmap_prep.writers import _cg_special, _special_bonds_split

NB = CGNonbonded(cutoff=1.1, epsilon_r=15.0, epsilon_rf=0.0)


def _rows(text: str) -> list[list[float]]:
    return [[float(v) for v in ln.split()] for ln in text.splitlines() if not ln.startswith("#")]


def test_pair_table_matches_gromacs_formulas_and_vanishes_at_cutoff() -> None:
    sigma, eps, qi, qj = 0.47, 3.5, 1.0, -1.0
    rows = _rows(pair_xvg(sigma, eps, qi, qj, NB))
    rc = NB.cutoff
    k_rf, c_rf = 1.0 / (2 * rc**3), 1.5 / rc

    def v(r: float) -> float:
        s6 = (sigma / r) ** 6
        lj = 4 * eps * (s6 * s6 - s6) - 4 * eps * ((sigma / rc) ** 12 - (sigma / rc) ** 6)
        return lj + ONE_4PI_EPS0 * qi * qj / 15.0 * (1 / r + k_rf * r * r - c_rf)

    r, e, _ = next(row for row in rows if abs(row[0] - 0.5) < 1e-9)
    assert e == pytest.approx(v(r), rel=1e-9)
    assert rows[-1][0] == pytest.approx(rc)
    assert rows[-1][1] == pytest.approx(0.0, abs=1e-8)
    # force column = -dV/dr (fine-step derivative of the analytic potential)
    h = 1e-6
    assert rows[0][0] == pytest.approx(0.236)  # table starts at sigma / 2
    for r_probe in (0.3, 0.42, 0.55, 0.8, 1.05):
        r = next(row[0] for row in rows if row[0] >= r_probe)
        k = next(i for i, row in enumerate(rows) if row[0] == r)
        assert rows[k][2] == pytest.approx(-(v(r + h) - v(r - h)) / (2 * h), rel=1e-7)


def test_g96_angle_table_energy_and_force_per_degree() -> None:
    rows = _rows(g96_angle_xvg(120.0, 35.0))
    assert rows[0][0] == 0.0
    assert rows[-1][0] == 180.0
    th, e, _ = rows[200]
    c0 = math.cos(math.radians(120.0))
    assert e == pytest.approx(0.5 * 35.0 * (math.cos(math.radians(th)) - c0) ** 2)

    def v(theta_deg: float) -> float:
        return 0.5 * 35.0 * (math.cos(math.radians(theta_deg)) - c0) ** 2

    h = 1e-6  # degrees
    for k in range(20, 340, 31):
        th = rows[k][0]
        assert rows[k][2] == pytest.approx(-(v(th + h) - v(th - h)) / (2 * h), rel=1e-6, abs=1e-12)


def test_nonbond_params_take_precedence_and_rule_1_converts() -> None:
    top = Topology()
    top.combination_rule = 2
    top.atom_types = {"A": AtomType("A", 72, 0, "A", sigma=0.4, epsilon=1.0)}
    top.nonbond_params = {("A", "A"): (0.47, 4.5)}
    assert lj_sigma_epsilon(top, "A", "A") == (0.47, 4.5)
    top.combination_rule = 1
    sigma, eps = 0.47, 4.5
    c6, c12 = 4 * eps * sigma**6, 4 * eps * sigma**12
    top.nonbond_params = {("A", "A"): (c6, c12)}
    assert lj_sigma_epsilon(top, "A", "A") == pytest.approx((sigma, eps))


def test_one_charge_per_cg_type() -> None:
    top = Topology()
    mol = MoleculeType(name="M", nrexcl=1)
    mol.atoms = [
        TopAtom(1, "Q1", 1, "M", "A", 1, 1.0, 72),
        TopAtom(2, "Q1", 1, "M", "B", 2, 0.5, 72),
    ]
    top.molecule_types = {"M": mol}
    with pytest.raises(ValueError, match="one charge per type"):
        cg_type_charges(top)


def test_split_exclusions_martini() -> None:
    assert _special_bonds_split(3, 1) == (
        "special_bonds lj 0.0 1.0e-100 1.0e-100 coul 0.0 1.0e-100 1.0e-100\n\n"
    )
    sim = SimulationParams(exclusion_nrexcl=3, exclusion_nrexcl_cg=1)
    assert _cg_special(sim) == " cg_special 0.0 1.0 1.0"
    assert _cg_special(SimulationParams(exclusion_nrexcl=3)) == ""
