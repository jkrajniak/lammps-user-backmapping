"""fix backmap/pairs without fix backmap: explicit 1-4 pairs in a plain AT system.

The AT-only continuation of a backmapped frame lists its 1-4 pairs explicitly
when they are not the fudge-scaled normal LJ (GROMOS [ pairtypes ]); there is no
fix backmap there, so every pair has weight 1.

Skipped unless ``BACKMAP_LMP`` points to a LAMMPS binary with the backmap package.
"""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

import pytest

LMP_ENV = "BACKMAP_LMP"
SIGMA, EPS, QQ = 3.4, 0.25, 0.5
R = 4.1
Q1, Q2 = 0.3, -0.2
QQRD2E = 332.06371  # real units


@pytest.mark.integration
def test_pairs_without_fix_backmap(tmp_path: Path) -> None:
    lmp = os.environ.get(LMP_ENV)
    if not lmp or not Path(lmp).is_file():
        pytest.skip(f"set {LMP_ENV} to a LAMMPS binary built with the backmap package")
    (tmp_path / "d.data").write_text(
        f"""plain pairs

2 atoms
1 atom types

0.0 30.0 xlo xhi
0.0 30.0 ylo yhi
0.0 30.0 zlo zhi

Masses

1 12.0

Atoms # full

1 1 1 {Q1} 10.0 10.0 10.0
2 1 1 {Q2} {10.0 + R} 10.0 10.0
"""
    )
    (tmp_path / "pairs.dat").write_text(f"1\n1 2 {SIGMA} {EPS} {QQ}\n")
    (tmp_path / "in.t").write_text(
        """units real
atom_style full
boundary p p p
read_data d.data
pair_style zero 10.0
pair_coeff * *
fix pairs14 all backmap/pairs at file pairs.dat
fix_modify pairs14 energy yes
thermo_style custom step pe
run 0
print "RESULT $(pe:%.12g)"
"""
    )
    proc = subprocess.run(
        [lmp, "-in", "in.t", "-log", "log.t", "-screen", "none"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    log = (tmp_path / "log.t").read_text() if (tmp_path / "log.t").exists() else ""
    assert proc.returncode == 0, (log + proc.stderr)[-2000:]
    got = float(re.search(r"^RESULT (\S+)", log, re.MULTILINE).group(1))
    lj = 4 * EPS * ((SIGMA / R) ** 12 - (SIGMA / R) ** 6)
    coul = QQRD2E * QQ * Q1 * Q2 / R
    assert got == pytest.approx(lj + coul, rel=1e-6)
