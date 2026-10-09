"""fix backmap distributes CG forces at setup; `peratom full` reports them.

setup() did not call post_force(), so after the force evaluation that opens
every run (and run 0) the CG forces stayed on the beads and the AT atoms got
none: the first velocity half-kick of each run missed them.

Skipped unless ``BACKMAP_LMP`` points to a LAMMPS binary with the backmap package.
"""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

import pytest

EPS, SIG = 0.8, 4.5
R = 5.0  # bead-bead distance along x
# Beads 1, 2 (type 1); atoms 3, 4 -> bead 1 and 5, 6 -> bead 2 (type 2, equal
# masses, placed symmetrically about their bead so each bead is at its COM).
XYZ = {
    1: (10.0, 10.0, 10.0),
    2: (10.0 + R, 10.0, 10.0),
    3: (10.0, 10.6, 10.0),
    4: (10.0, 9.4, 10.0),
    5: (10.0 + R, 10.0, 10.6),
    6: (10.0 + R, 10.0, 9.4),
}


def _run(tmp_path: Path) -> dict[str, float]:
    lmp = os.environ.get("BACKMAP_LMP")
    if not lmp or not Path(lmp).is_file():
        pytest.skip("set BACKMAP_LMP to a LAMMPS binary built with the backmap package")
    atoms = "\n".join(f"{i} 1 {1 if i <= 2 else 2} 0.0 {x} {y} {z}" for i, (x, y, z) in XYZ.items())
    (tmp_path / "test.data").write_text(
        f"""setup test

6 atoms
2 atom types

0.0 40.0 xlo xhi
0.0 40.0 ylo yhi
0.0 40.0 zlo zhi

Masses

1 24.0
2 12.0

Atoms # full

{atoms}
"""
    )
    fx = " ".join(f"fx{i}=$(fx[{i}]:%.12e)" for i in XYZ)
    col = " ".join(f"c5_{i}=$(v_c5[{i}]:%.12e) c2_{i}=$(v_c2[{i}]:%.12e)" for i in XYZ)
    (tmp_path / "in.test").write_text(
        f"""units real
atom_style full
boundary p p p
read_data test.data
pair_style backmap 12.0 lj/cut 12.0 12.0 lj/cut 12.0
pair_coeff 1 1 cg {EPS} {SIG}
pair_coeff 1 2 none
pair_coeff 2 2 atomistic 0.0 1.0
fix bm all backmap cg_type 1 alpha 0.0001 lambda0 0.0 peratom full
variable c5 atom f_bm[5]
variable c2 atom f_bm[2]
run 0
print "RESULT {fx} {col}"
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
    return {k: float(v) for k, v in (kv.split("=") for kv in m.group(1).split())}


@pytest.mark.integration
def test_cg_force_distributed_at_setup(tmp_path: Path) -> None:
    sr6 = (SIG / R) ** 6
    f_on_2 = 24.0 * EPS * (2.0 * sr6 * sr6 - sr6) / R  # +x on bead 2 (repulsive)
    res = _run(tmp_path)
    for bead in (1, 2):
        assert res[f"fx{bead}"] == pytest.approx(0.0, abs=1e-12)
    for atom, sign in ((3, -1), (4, -1), (5, 1), (6, 1)):
        assert res[f"fx{atom}"] == pytest.approx(sign * 0.5 * f_on_2, rel=1e-10)
        # per-atom array: redistributed force (col 5) and bead position x (col 2)
        assert res[f"c5_{atom}"] == pytest.approx(sign * 0.5 * f_on_2, rel=1e-10)
        assert res[f"c2_{atom}"] == pytest.approx(10.0 if atom < 5 else 10.0 + R)
    # beads report their CG force before redistribution
    assert res["c5_2"] == pytest.approx(f_on_2, rel=1e-10)
    assert res["c5_1"] == pytest.approx(-f_on_2, rel=1e-10)
