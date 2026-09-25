"""LAMMPS regression tests for fix backmap/pairs (1-4 LJ and Coulomb).

Skipped unless ``BACKMAP_LMP`` points to a LAMMPS binary built with the
backmap package.
"""

from __future__ import annotations

import math
import os
import re
import subprocess
from pathlib import Path

import pytest

LMP_ENV = "BACKMAP_LMP"
QQRD2E_REAL = 332.06371  # LAMMPS units real
SIGMA, EPSILON, QQ_SCALE = 3.0, 0.2, 0.5
Q2, Q5 = 0.4, -0.4

# One CG bead (type 1) and a four-atom AT chain (type 2), atoms 2-5; atoms 2
# and 5 are 1-4 partners 3.6 A apart along x.
DATA = """\
fix backmap/pairs test

5 atoms
3 bonds
2 atom types
1 bond types

0.0 40.0 xlo xhi
0.0 40.0 ylo yhi
0.0 40.0 zlo zhi

Masses

1 56.0
2 14.0

Atoms # full

1 1 1 0.0 21.8 20.0 20.0
2 1 2 {q2} 20.0 20.0 20.0
3 1 2 0.0 21.2 20.0 20.0
4 1 2 0.0 22.4 20.0 20.0
5 1 2 {q5} 23.6 20.0 20.0

Bonds

1 1 2 3
2 1 3 4
3 1 4 5
"""

INPUT = """\
units real
atom_style full
boundary p p p
read_data test.data
pair_style backmap 12.0 lj/cut/coul/cut 12.0 12.0 12.0 lj/cut 12.0
pair_coeff 1 1 none
pair_coeff 1 2 none
pair_coeff 2 2 atomistic 0.0 1.0
bond_style harmonic
bond_coeff 1 0.0 1.2
special_bonds lj 0.0 0.0 0.0 coul 0.0 0.0 0.0
fix bm all backmap cg_type 1 alpha 0.0001 lambda0 1.0
fix pairs all backmap/pairs at file pairs.dat cut 12.0
thermo_style custom step pe f_pairs f_pairs[1] f_pairs[2]
run 0
variable f2x equal fx[2]
variable f5x equal fx[5]
print "RESULT pe=$(pe:%.12e) total=$(f_pairs:%.12e) lj=$(f_pairs[1]:%.12e) coul=$(f_pairs[2]:%.12e) f2x=$(v_f2x:%.12e) f5x=$(v_f5x:%.12e)"
"""


def _lmp() -> Path:
    path = os.environ.get(LMP_ENV)
    if not path or not Path(path).is_file():
        pytest.skip(f"set {LMP_ENV} to a LAMMPS binary built with the backmap package")
    return Path(path)


def _run(tmp_path: Path, pairs_line: str, shift: float = 0.0) -> dict[str, float]:
    data = DATA.format(q2=Q2, q5=Q5)
    if shift:
        # Translate the molecule along x and wrap it into the 40 A box.
        lines = []
        section = None
        for line in data.splitlines():
            if line.strip() in ("Atoms # full", "Bonds", "Masses"):
                section = line.strip()
            parts = line.split()
            if section == "Atoms # full" and len(parts) == 7:
                parts[4] = f"{(float(parts[4]) + shift) % 40.0:.4f}"
                line = " ".join(parts)
            lines.append(line)
        data = "\n".join(lines) + "\n"
    (tmp_path / "test.data").write_text(data)
    (tmp_path / "pairs.dat").write_text(f"1\n{pairs_line}\n")
    (tmp_path / "in.test").write_text(INPUT)
    proc = subprocess.run(
        [str(_lmp()), "-in", "in.test", "-log", "log.test", "-screen", "none"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    log = (tmp_path / "log.test").read_text()
    assert proc.returncode == 0, log[-2000:]
    match = re.search(r"^RESULT (.*)$", log, re.MULTILINE)
    assert match, log[-2000:]
    return {k: float(v) for k, v in (kv.split("=") for kv in match.group(1).split())}


def _expected(r: float, qq_scale: float) -> tuple[float, float, float]:
    """(E_LJ, E_Coul, force on atom 5 along +x) for the 2-5 pair at distance r."""
    sr6 = (SIGMA / r) ** 6
    e_lj = 4.0 * EPSILON * (sr6 * sr6 - sr6)
    e_coul = QQRD2E_REAL * qq_scale * Q2 * Q5 / r
    de_dr = -24.0 * EPSILON * (2.0 * sr6 * sr6 - sr6) / r - e_coul / r
    return e_lj, e_coul, -de_dr


@pytest.mark.integration
def test_lj_and_coulomb_14_energy_and_force(tmp_path: Path) -> None:
    res = _run(tmp_path, f"2 5 {SIGMA} {EPSILON} {QQ_SCALE}")
    e_lj, e_coul, f5 = _expected(3.6, QQ_SCALE)

    assert res["lj"] == pytest.approx(e_lj, rel=1e-10)
    assert res["coul"] == pytest.approx(e_coul, rel=1e-10)
    assert res["total"] == pytest.approx(e_lj + e_coul, rel=1e-10)
    # Counted in the potential energy from the first force evaluation (setup).
    assert res["pe"] == pytest.approx(e_lj + e_coul, rel=1e-10)
    assert res["f5x"] == pytest.approx(f5, rel=1e-10)
    assert res["f2x"] == pytest.approx(-f5, rel=1e-10)
    assert math.isfinite(res["pe"])


@pytest.mark.integration
def test_four_column_pairs_file_is_lj_only(tmp_path: Path) -> None:
    res = _run(tmp_path, f"2 5 {SIGMA} {EPSILON}")
    e_lj, _, _ = _expected(3.6, 0.0)

    assert res["lj"] == pytest.approx(e_lj, rel=1e-10)
    assert res["coul"] == 0.0


@pytest.mark.integration
def test_pair_across_periodic_boundary(tmp_path: Path) -> None:
    """The 1-4 pair spans the x boundary; both atoms are owned by the one rank.

    The closest image of the partner is then a periodic ghost. Its force must
    still land on the owned atom and the energy must be counted once.
    """
    res = _run(tmp_path, f"2 5 {SIGMA} {EPSILON} {QQ_SCALE}", shift=18.9)
    e_lj, e_coul, f5 = _expected(3.6, QQ_SCALE)

    assert res["lj"] == pytest.approx(e_lj, rel=1e-10)
    assert res["coul"] == pytest.approx(e_coul, rel=1e-10)
    assert res["f5x"] == pytest.approx(f5, rel=1e-10)
    assert res["f2x"] == pytest.approx(-f5, rel=1e-10)
