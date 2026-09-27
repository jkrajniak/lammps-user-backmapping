"""Intra-bead AT terms keep full weight wherever the bead sits in the box.

fix backmap mapped only one copy of each AT atom to its bead, so a bonded
term whose partner in the bond list was another periodic image or a ghost on
another rank was classified as inter-bead and weighted by lambda instead of 1:
an intra-bead bond across the periodic boundary was switched off at lambda = 0,
and the result depended on the domain decomposition. The minimizer also never
built the bead map (no MIN_PRE_FORCE), so its first evaluation weighted every
intra-bead term as inter-bead.

One molecule: CG bead (type 1) and two AT atoms (type 2) joined by one
intra-bead backmap/harmonic bond; the molecule is placed inside the box,
across the periodic boundary and across the processor boundary of a 2-rank
run. The bond energy must be K (r - r0)^2 in every case.

Skipped unless ``BACKMAP_LMP`` points to a LAMMPS binary with the backmap package.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

LMP_ENV = "BACKMAP_LMP"
L = 20.0
K, R0 = 100.0, 1.5
X_BEAD, X_A, X_B = 19.9, 19.2, 20.6  # bond length 1.4 A, crosses x = L when unshifted
E_FULL = K * ((X_B - X_A) - R0) ** 2

DATA = f"""same-bead image test

3 atoms
2 atom types
1 bonds
1 bond types

0.0 {L} xlo xhi
0.0 {L} ylo yhi
0.0 {L} zlo zhi

Masses

1 28.0
2 14.0

Atoms # full

1 1 1 0.0 {X_BEAD} 10.0 10.0
2 1 2 0.0 {X_A} 10.0 10.0
3 1 2 0.0 {X_B - L:.4f} 10.0 10.0

Bonds

1 1 2 3
"""

# shift along x: 0 = across the periodic boundary, -9.5 = inside the box and
# across the x = 10 processor boundary of a "processors 2 1 1" run, -5 = inside
# one subdomain.
PLACEMENTS = {"periodic": 0.0, "processor": -9.5, "inside": -5.0}


def _lmp() -> str:
    lmp = os.environ.get(LMP_ENV)
    if not lmp or not Path(lmp).is_file():
        pytest.skip(f"set {LMP_ENV} to a LAMMPS binary built with the backmap package")
    return lmp


def _run(tmp: Path, shift: float, lam: float, np_: int, command: str) -> float:
    lmp = _lmp()
    tmp.mkdir(parents=True, exist_ok=True)
    (tmp / "d.data").write_text(DATA)
    procs = "processors 2 1 1" if np_ == 2 else ""
    (tmp / "in.t").write_text(
        f"""units real
atom_style full
boundary p p p
{procs}
read_data d.data
displace_atoms all move {shift} 0 0 units box
pair_style zero 5.0
pair_coeff * *
bond_style backmap/harmonic
bond_coeff 1 at {K} {R0}
special_bonds lj 0 0 0
comm_modify cutoff 8.0
fix bm all backmap cg_type 1 alpha 0.0001 lambda0 {lam}
fix hold all setforce 0.0 0.0 0.0
thermo_style custom step ebond
thermo_modify format float %.12g
{command}
"""
    )
    cmd = [lmp, "-in", "in.t", "-log", "log.t", "-screen", "none"]
    if np_ > 1:
        mpirun = shutil.which("mpirun")
        if mpirun is None:
            pytest.skip("mpirun not available")
        cmd = [mpirun, "-np", str(np_), *cmd]
    proc = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True, check=False)
    log = (tmp / "log.t").read_text() if (tmp / "log.t").exists() else ""
    assert proc.returncode == 0, (log + proc.stderr)[-2000:]
    # first thermo line after the header: step 0 of the run or minimization
    match = re.search(r"^\s+Step\s+E_bond\s*\n\s*0\s+(\S+)", log, re.MULTILINE)
    assert match, log[-2000:]
    return float(match.group(1))


@pytest.mark.integration
@pytest.mark.parametrize("np_", [1, 2])
@pytest.mark.parametrize("lam", [0.0, 0.5])
@pytest.mark.parametrize("where", sorted(PLACEMENTS))
def test_intra_bead_bond_full_weight(tmp_path: Path, where: str, lam: float, np_: int) -> None:
    got = _run(tmp_path, PLACEMENTS[where], lam, np_, "run 0")
    assert got == pytest.approx(E_FULL, rel=1e-10)


@pytest.mark.integration
def test_minimize_builds_bead_map(tmp_path: Path) -> None:
    """The minimizer's first energy evaluation sees the intra-bead bond at full weight."""
    got = _run(tmp_path, PLACEMENTS["inside"], 0.0, 1, "minimize 1.0e-30 1.0e-30 1 1")
    assert got == pytest.approx(E_FULL, rel=1e-10)
