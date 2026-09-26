"""pair_style backmap honours special-bond factors; cg_special sets them for CG pairs.

The pair style used to call its sub-styles with factor 1 for every listed pair,
which was only right because every generated input excludes 1-2..1-4 fully.
A MARTINI CG model excludes only bonded neighbours (nrexcl = 1) while the AT
force field excludes up to 1-4, so CG pairs need their own special weights.

Skipped unless ``BACKMAP_LMP`` points to a LAMMPS binary with the backmap package.
"""

from __future__ import annotations

import itertools
import math
import os
import re
import subprocess
from pathlib import Path

import pytest

LMP_ENV = "BACKMAP_LMP"
EPS, SIG, CUT = 0.3, 3.2, 12.0
# One molecule: CG beads 1-4 (type 1), AT atoms 5-8 (type 2), one AT atom per
# bead, both chains bonded linearly. fix backmap puts each bead on its atom.
XYZ = [(10.0, 10.0, 10.0), (13.5, 11.0, 10.0), (16.5, 10.0, 11.0), (19.0, 12.0, 10.5)]


def _data() -> str:
    atoms = [f"{k + 1} 1 1 0.0 {x} {y} {z}" for k, (x, y, z) in enumerate(XYZ)]
    atoms += [f"{k + 5} 1 2 0.0 {x} {y} {z}" for k, (x, y, z) in enumerate(XYZ)]
    bonds = ["1 1 1 2", "2 1 2 3", "3 1 3 4", "4 1 5 6", "5 1 6 7", "6 1 7 8"]
    return "\n".join(
        [
            "pair special test",
            "",
            "8 atoms",
            "6 bonds",
            "2 atom types",
            "1 bond types",
            "",
            "0.0 40.0 xlo xhi",
            "0.0 40.0 ylo yhi",
            "0.0 40.0 zlo zhi",
            "",
            "Masses",
            "",
            "1 12.0",
            "2 12.0",
            "",
            "Atoms # full",
            "",
            *atoms,
            "",
            "Bonds",
            "",
            *bonds,
            "",
        ]
    )


def _lj(r: float) -> float:
    return 4.0 * EPS * ((SIG / r) ** 12 - (SIG / r) ** 6)


def _expected(weights: tuple[float, float, float]) -> float:
    """Chain energy with special weights (1-2, 1-3, 1-4) by topological distance."""
    total = 0.0
    for a, b in itertools.combinations(range(4), 2):
        w = weights[b - a - 1]
        total += w * _lj(math.dist(XYZ[a], XYZ[b]))
    return total


def _run(tmp_path: Path, special: str, cg_special: str, lam: float) -> float:
    lmp = os.environ.get(LMP_ENV)
    if not lmp or not Path(lmp).is_file():
        pytest.skip(f"set {LMP_ENV} to a LAMMPS binary built with the backmap package")
    work = tmp_path / f"run_{abs(hash((special, cg_special, lam)))}"
    work.mkdir()
    (work / "test.data").write_text(_data())
    (work / "in.test").write_text(
        f"""units real
atom_style full
boundary p p p
read_data test.data
pair_style backmap {CUT} lj/cut {CUT} {CUT} lj/cut {CUT} {cg_special}
pair_coeff 1 1 cg {EPS} {SIG}
pair_coeff 1 2 none
pair_coeff 2 2 atomistic {EPS} {SIG}
bond_style zero
bond_coeff 1
special_bonds {special}
fix bm all backmap cg_type 1 alpha 0.0001 lambda0 {lam}
thermo_style custom step evdwl
run 0
print "RESULT evdwl=$(evdwl:%.12e)"
"""
    )
    proc = subprocess.run(
        [lmp, "-in", "in.test", "-log", "log.test", "-screen", "none"],
        cwd=work,
        capture_output=True,
        text=True,
        check=False,
    )
    log = (work / "log.test").read_text() if (work / "log.test").exists() else ""
    assert proc.returncode == 0, (log + proc.stderr)[-2000:]
    match = re.search(r"^RESULT evdwl=(\S+)", log, re.MULTILINE)
    assert match, log[-2000:]
    return float(match.group(1))


TINY = "lj 0.0 1.0e-100 1.0e-100 coul 0.0 1.0e-100 1.0e-100"


@pytest.mark.integration
@pytest.mark.parametrize(
    ("special", "cg_special", "lam", "weights"),
    [
        # MARTINI-style: CG keeps 1-3 and 1-4, AT excludes them.
        (TINY, "cg_special 0.0 1.0 1.0", 0.0, (0.0, 1.0, 1.0)),
        (TINY, "cg_special 0.0 1.0 1.0", 1.0, (0.0, 0.0, 0.0)),
        # Plain special factors now reach the sub-styles (1-4 at one half).
        ("lj 0.0 0.0 0.5 coul 0.0 0.0 0.5", "", 0.0, (0.0, 0.0, 0.5)),
        ("lj 0.0 0.0 0.5 coul 0.0 0.0 0.5", "", 1.0, (0.0, 0.0, 0.5)),
        # Full exclusion, the setting of every existing example.
        ("lj 0.0 0.0 0.0 coul 0.0 0.0 0.0", "", 0.0, (0.0, 0.0, 0.0)),
    ],
)
def test_special_weights(
    tmp_path: Path,
    special: str,
    cg_special: str,
    lam: float,
    weights: tuple[float, float, float],
) -> None:
    got = _run(tmp_path, special, cg_special, lam)
    want = _expected(weights)
    assert got == pytest.approx(want, rel=1e-10, abs=1e-12)
