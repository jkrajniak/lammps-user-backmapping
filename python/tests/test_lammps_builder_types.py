"""Atom types of the network LAMMPS builder: CG beads, masses, and atomistic pair parameters."""

from __future__ import annotations

from types import SimpleNamespace
from typing import Any

import pytest

from backmap_prep import units
from backmap_prep.builder import AtomTypeInfo
from backmap_prep.network.lammps_builder import (
    _atom_type_info,
    _backmap_mol_ids,
    _cg_type_names,
    _pair_terms,
)
from backmap_prep.parsers.top_parser import AtomType, MoleculeType, TopAtom, Topology


def _settings(*bead_types: str) -> Any:
    beads = [SimpleNamespace(type=t) for t in bead_types]
    return SimpleNamespace(molecules=[SimpleNamespace(beads=beads)])


def _molecule(types: list[str], masses: list[float] | None = None) -> MoleculeType:
    mol = MoleculeType(name="M", nrexcl=3)
    for n, t in enumerate(types, 1):
        mass = masses[n - 1] if masses else 0.0
        mol.atoms.append(TopAtom(n, t, 1, "M", f"{t}{n}", n, 0.0, mass))
    return mol


# A GROMOS-like hybrid: CG bead "CGA", then atomistic C, O, H
GROMOS = ["CGA", "C", "O", "H", "CGA", "C", "O", "H"]


def test_single_letter_atomistic_types_are_not_beads_when_the_settings_name_the_beads() -> None:
    mol = _molecule(GROMOS)
    names = _cg_type_names(_settings("CGA"), mol)
    assert names == {"CGA"}
    assert _backmap_mol_ids(mol, names) == {1: 1, 2: 1, 3: 1, 4: 1, 5: 2, 6: 2, 7: 2, 8: 2}
    infos, _ = _atom_type_info(Topology(), mol, names)
    assert {i.name for i in infos if i.is_cg} == {"CGA"}


def test_single_letter_fallback_warns_when_the_settings_name_no_beads() -> None:
    with pytest.warns(UserWarning, match="single-letter"):
        assert _cg_type_names(_settings(), _molecule(["A", "c1", "c2", "B", "c3"])) == {"A", "B"}


def _masses(mol: MoleculeType, top: Topology) -> dict[str, float]:
    infos, _ = _atom_type_info(top, mol, {"CGA"})
    return {i.name: i.mass for i in infos}


def test_atoms_mass_takes_precedence_over_the_atomtypes_mass() -> None:
    top = Topology(
        atom_types={"C": AtomType("C", 12.011, 0.0, "A"), "H": AtomType("H", 1.008, 0.0, "A")}
    )
    # hydrogen mass repartitioning written in [ atoms ]; C keeps its type mass
    mol = _molecule(["CGA", "C", "H"], masses=[72.0, 0.0, 3.024])
    assert _masses(mol, top) == pytest.approx({"CGA": 72.0, "C": 12.011, "H": 3.024})


def test_masses_kept_only_in_atoms_are_not_zero() -> None:
    mol = _molecule(["CGA", "C"], masses=[72.0, 12.011])
    assert _masses(mol, Topology()) == pytest.approx({"CGA": 72.0, "C": 12.011})


def test_atoms_of_one_type_with_different_masses_are_an_error() -> None:
    mol = _molecule(["CGA", "H", "H"], masses=[72.0, 1.008, 3.024])
    with pytest.raises(ValueError, match="type H have different masses"):
        _masses(mol, Topology())


def test_cg_beads_keep_the_mass_of_the_cg_force_field() -> None:
    """The [ atoms ] mass of a bead is the sum of its atoms and differs between reacted beads."""
    top = Topology(atom_types={"CGA": AtomType("CGA", 72.0, 0.0, "A")})
    mol = _molecule(["CGA", "C", "CGA", "C"], masses=[87.166, 12.011, 86.158, 12.011])
    assert _masses(mol, top) == pytest.approx({"CGA": 72.0, "C": 12.011})


def _at(type_id: int, name: str, sigma_nm: float, eps_kj: float) -> AtomTypeInfo:
    return AtomTypeInfo(
        type_id, name, 12.0, False, sigma=units.sigma(sigma_nm), epsilon=units.epsilon(eps_kj)
    )


@pytest.mark.parametrize("rule", [2, 3])
def test_nonbond_params_override_the_mixed_atomistic_pair(rule: int) -> None:
    types = [_at(1, "OA", 0.30, 0.8), _at(2, "NT", 0.33, 0.6)]
    nbfix = {("OA", "NT"): (0.29, 1.5), ("NT", "OA"): (0.29, 1.5)}
    mixed = {(p.itype, p.jtype): p for p in _pair_terms(types, rule)}
    fixed = {(p.itype, p.jtype): p for p in _pair_terms(types, rule, nbfix)}
    assert fixed[1, 2].sigma == pytest.approx(units.sigma(0.29))
    assert fixed[1, 2].epsilon == pytest.approx(units.epsilon(1.5))
    assert mixed[1, 2].sigma != pytest.approx(fixed[1, 2].sigma)
    for key in ((1, 1), (2, 2)):
        assert (fixed[key].sigma, fixed[key].epsilon) == (mixed[key].sigma, mixed[key].epsilon)


def test_nonbond_params_in_c6_c12_for_combination_rule_1() -> None:
    types = [_at(1, "OA", 0.30, 0.8), _at(2, "NT", 0.33, 0.6)]
    c6, c12 = 2.0e-3, 1.0e-6
    (p,) = [p for p in _pair_terms(types, 1, {("OA", "NT"): (c6, c12)}) if p.itype != p.jtype]
    assert p.sigma == pytest.approx(units.sigma((c12 / c6) ** (1 / 6)))
    assert p.epsilon == pytest.approx(units.epsilon(c6 * c6 / (4 * c12)))
