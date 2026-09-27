"""GROMOS-style 1-4 pairs: [ pairtypes ] give the parameters of listed [ pairs ].

With gen-pairs no, a listed pair without its own parameters takes the
pairtypes entry of its atom types (C6/C12 under combination rule 1), not a
fudge-scaled normal LJ.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from backmap_prep import units
from backmap_prep.parsers.top_parser import parse_top, resolve_pair_lj_params

if TYPE_CHECKING:
    from pathlib import Path

TOP = """[ defaults ]
1 1 no 1.0 1.0
[ atomtypes ]
CH2 6 14.027 0.0 A 0.0074684164 3.3965584e-05
CH3 6 15.035 0.0 A 0.0096138025 2.6646244e-05
[ pairtypes ]
CH2 CH2 1 0.0047238129 4.7419261e-06
CH3 CH2 1 0.0056894694 5.347702e-06
[ moleculetype ]
BUT 3
[ atoms ]
1 CH3 1 BUT C1 1 0.0 15.035
2 CH2 1 BUT C2 2 0.0 14.027
3 CH2 1 BUT C3 3 0.0 14.027
4 CH2 1 BUT C4 4 0.0 14.027
[ pairs ]
1 4 1
"""


def test_pairtype_parameters_are_used(tmp_path: Path) -> None:
    path = tmp_path / "t.top"
    path.write_text(TOP)
    top = parse_top(path)
    assert top.defaults_gen_pairs == "no"
    mol = top.molecule_types["BUT"]
    a1, a4 = mol.atoms[0], mol.atoms[3]
    sigma, epsilon = resolve_pair_lj_params(top, a1, a4, 1, [])
    s_nm, e_kj = units.c6c12_to_sigma_epsilon(0.0056894694, 5.347702e-06)
    assert sigma == pytest.approx(units.distance(s_nm))
    assert epsilon == pytest.approx(units.energy(e_kj))
