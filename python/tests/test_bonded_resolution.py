"""Bonds and angles given without parameters resolve from [ bondtypes ] / [ angletypes ].

bakery writes PET crosslink bonds as ``i j`` in [ cross_bonds ]; GROMACS takes
their parameters from the force field's bondtypes. They used to become K = 0,
so the crosslinks held nothing together.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import pytest

from backmap_prep import units
from backmap_prep.network.lammps_builder import _angletype_params, _bondtype_params
from backmap_prep.parsers.top_parser import (
    AtomType,
    MoleculeType,
    TopAtom,
    Topology,
    _parse_cross_pair,
)

REPO = Path(__file__).resolve().parents[2]
PET_DATA = REPO.parent / "paper-reverse-mapping-polymer-networks/preparation/dacron/backmapping"


def _atom(index: int, type_name: str) -> TopAtom:
    return TopAtom(index, type_name, 1, "MOL", f"X{index}", index, 0.0, 12.0)


def _topology() -> Topology:
    top = Topology()
    top.atom_types = {
        "opls_1": AtomType("opls_1", 12.0, 0.0, "A", bond_type="CT"),
        "opls_2": AtomType("opls_2", 16.0, 0.0, "A", bond_type="OS"),
    }
    top.bondtypes = {("CT", "OS"): [0.141, 267776.0], ("OS", "CT"): [0.141, 267776.0]}
    top.angletypes = {("CT", "OS", "CT"): [109.5, 418.4]}
    return top


def test_bond_resolves_from_bondtypes_by_class() -> None:
    k, r0 = _bondtype_params(_topology(), _atom(1, "opls_1"), _atom(2, "opls_2"))
    assert k == pytest.approx(units.spring_bond(267776.0))
    assert k == pytest.approx(320.0, rel=1e-4)
    assert r0 == pytest.approx(1.41)


def test_angle_resolves_from_angletypes_by_class() -> None:
    atoms = (_atom(1, "opls_1"), _atom(2, "opls_2"), _atom(3, "opls_1"))
    k, theta = _angletype_params(_topology(), atoms)
    assert k == pytest.approx(units.spring_angle(418.4))
    assert theta == pytest.approx(109.5)


def test_unresolvable_bond_is_an_error_not_zero() -> None:
    with pytest.raises(ValueError, match=r"no \[ bondtypes \] entry"):
        _bondtype_params(_topology(), _atom(1, "opls_1"), _atom(2, "opls_1"))


def test_two_field_pair_line_uses_function_1() -> None:
    mol = MoleculeType(name="MOL", nrexcl=3)
    _parse_cross_pair(["2", "9565"], mol)
    assert len(mol.cross_pairs) == 1
    assert mol.cross_pairs[0].func == 1


@pytest.mark.integration
@pytest.mark.skipif(not PET_DATA.is_dir(), reason="PET paper-data bundle not present")
def test_pet_has_no_zero_at_bonded_types(tmp_path: Path) -> None:
    import yaml

    from backmap_prep.cli import main

    example = REPO / "examples/pet/large"
    raw = yaml.safe_load((example / "settings.v2.yaml").read_text())
    for key in ("data_dir", "tables_dir", "forcefield_dir"):
        if raw["prep"].get(key):
            raw["prep"][key] = str((example / raw["prep"][key]).resolve())
    for key in ("angles_file", "dihedrals_file"):
        if raw["cross_interactions"].get(key):
            raw["cross_interactions"][key] = str(example / raw["cross_interactions"][key])
    settings = tmp_path / "settings.v2.yaml"
    settings.write_text(yaml.safe_dump(raw))
    old = os.getcwd()
    try:
        os.chdir(tmp_path)
        assert main(["build", str(settings)]) == 0
    finally:
        os.chdir(old)
    ff = (tmp_path / "pet.ff.lmp").read_text()
    zero = [ln for ln in ff.splitlines() if re.match(r"^(bond|angle)_coeff \d+ \S+ at 0 0$", ln)]
    assert zero == []
