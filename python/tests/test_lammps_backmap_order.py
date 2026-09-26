"""fix backmap gives the same trajectory whether defined before or after the integrator.

The COM update ran in initial_integrate(), so with fix backmap defined before
the integration fix the beads followed the AT positions of the previous step.
It now runs in post_integrate(), after every fix has moved its atoms.

Skipped unless ``BACKMAP_LMP`` points to a LAMMPS binary with the backmap package.
"""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

import pytest

DATA = """order test

6 atoms
2 atom types

0.0 40.0 xlo xhi
0.0 40.0 ylo yhi
0.0 40.0 zlo zhi

Masses

1 24.0
2 12.0

Atoms # full

1 1 1 0.0 10.0 10.0 10.0
2 1 1 0.0 15.0 10.0 10.0
3 1 2 0.0 10.0 10.6 10.0
4 1 2 0.0 10.0 9.4 10.0
5 1 2 0.0 15.0 10.0 10.6
6 1 2 0.0 15.0 10.0 9.4

Velocities

1 0 0 0
2 0 0 0
3 0.002 -0.001 0.003
4 -0.002 0.001 0.0
5 0.0 0.003 -0.002
6 0.001 0.0 0.002
"""

FIX_BM = "fix bm all backmap cg_type 1 alpha 0.01 lambda0 0.0\nfix_modify bm active yes"
FIX_NVE = "fix integrate all nve"


def _positions(tmp_path: Path, first: str, second: str) -> list[float]:
    lmp = os.environ.get("BACKMAP_LMP")
    if not lmp or not Path(lmp).is_file():
        pytest.skip("set BACKMAP_LMP to a LAMMPS binary built with the backmap package")
    (tmp_path / "test.data").write_text(DATA)
    xs = " ".join(f"$(x[{i}]:%.14e) $(y[{i}]:%.14e) $(z[{i}]:%.14e)" for i in range(1, 7))
    (tmp_path / "in.test").write_text(
        f"""units real
atom_style full
boundary p p p
read_data test.data
pair_style backmap 12.0 lj/cut 12.0 12.0 lj/cut 12.0
pair_coeff 1 1 cg 0.8 4.5
pair_coeff 1 2 none
pair_coeff 2 2 atomistic 0.1 3.0
{first}
{second}
timestep 1.0
run 20
print "RESULT {xs}"
"""
    )
    proc = subprocess.run(
        [lmp, "-in", "in.test", "-log", "log.test", "-screen", "none"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    log = (tmp_path / "log.test").read_text() if (tmp_path / "log.test").exists() else ""
    assert proc.returncode == 0, (log + proc.stderr)[-2000:]
    m = re.search(r"^RESULT (.*)$", log, re.MULTILINE)
    assert m, log[-2000:]
    return [float(v) for v in m.group(1).split()]


@pytest.mark.integration
def test_trajectory_independent_of_fix_order(tmp_path: Path) -> None:
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    before = _positions(tmp_path / "a", FIX_BM, FIX_NVE)
    after = _positions(tmp_path / "b", FIX_NVE, FIX_BM)
    assert before == pytest.approx(after, rel=0, abs=1e-12)
    # the atoms did move (the test is not trivially equal)
    assert abs(before[6] - 10.0) > 1e-3
