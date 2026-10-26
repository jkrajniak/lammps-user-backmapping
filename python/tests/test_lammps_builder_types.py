"""Atom types of the network LAMMPS builder: which types are CG beads, and their masses."""

from __future__ import annotations

from types import SimpleNamespace
from typing import Any

import pytest

from backmap_prep.network.lammps_builder import (
    _atom_type_info,
    _backmap_mol_ids,
    _cg_type_names,
)
from backmap_prep.parsers.top_parser import MoleculeType, TopAtom, Topology


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
