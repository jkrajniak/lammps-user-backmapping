"""The AT-only force field covers the bonded extent of its frame.

LAMMPS finds a bonded partner among the ghosts within the communication
cutoff (pair cutoff + skin by default). In a stretched, unrelaxed frame a
bond longer than that made it pair the far periodic image of the partner, so
E_bond was wrong (found by the GROMACS energy check of pe_10 on the build
frame); at-system now measures the extent and sets comm_modify cutoff.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from backmap_prep.at_system import bonded_extent

if TYPE_CHECKING:
    from pathlib import Path

DATA = """test

3 atoms
1 atom types
2 bonds
1 bond types
1 angles
1 angle types

0.0 40.0 xlo xhi
0.0 40.0 ylo yhi
0.0 40.0 zlo zhi

Masses

1 12.0

Atoms # full

1 1 1 0.0 39.5 10.0 10.0 0 0 0
2 1 1 0.0 1.0 10.0 10.0 0 0 0
3 1 1 0.0 {x3} 10.0 10.0 0 0 0

Bonds

1 1 1 2
2 1 2 3

Angles

1 1 1 2 3
"""


@pytest.mark.parametrize(("x3", "want"), [(2.5, 3.0), (19.0, 19.5)])
def test_bonded_extent_min_image(tmp_path: Path, x3: float, want: float) -> None:
    path = tmp_path / "at.data"
    path.write_text(DATA.format(x3=x3))
    # bond 1-2 crosses the boundary (1.5 A by minimum image), the angle spans 1..3
    assert bonded_extent(path) == pytest.approx(want)
