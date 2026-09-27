"""A bead's CG force reaches its AT atoms on every domain decomposition.

fix backmap forward-communicates each bead's force and AT mass sum to its
ghost copies. A ghost that is forwarded again (2-D/3-D brick or tiled
decompositions, where a diagonal neighbour receives it in two hops) was
packed from its own partial force and zero mass sum instead of the values it
had received, so AT atoms whose bead reached them that way got no CG force
(found by the static MPI parity check, review item B1: dodecane, 4 ranks, RCB).

One molecule: a bead (type 1) near the corner of the lower-left subdomain of
a 2 x 2 x 1 grid and its two AT atoms (type 2, masses 12 and 16) in the
diagonally opposite subdomain. A constant force on the bead (fix addforce,
defined before fix backmap) must be split over the AT atoms by mass.

Skipped unless ``BACKMAP_LMP`` points to a LAMMPS binary with the backmap package.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

LMP_ENV = "BACKMAP_LMP"
F_BEAD = (3.0, -2.0, 1.5)
M_AT = {2: 12.0, 3: 16.0}

DATA = """cg force multihop test

3 atoms
3 atom types
1 bonds
1 bond types

0.0 20.0 xlo xhi
0.0 20.0 ylo yhi
0.0 20.0 zlo zhi

Masses

1 28.0
2 12.0
3 16.0

Atoms # full

1 1 1 0.0 9.9 9.9 10.0
2 1 2 0.0 10.3 10.4 10.0
3 1 3 0.0 10.6 10.9 10.0

Bonds

1 1 2 3
"""

DECOMPOSITIONS = {
    "serial": (1, ""),
    "brick_2x2": (4, "processors 2 2 1"),
    "tiled_rcb": (4, "comm_style tiled"),
}


def _run(tmp: Path, np_: int, decomp: str) -> dict[int, np.ndarray]:
    lmp = os.environ.get(LMP_ENV)
    if not lmp or not Path(lmp).is_file():
        pytest.skip(f"set {LMP_ENV} to a LAMMPS binary built with the backmap package")
    (tmp / "d.data").write_text(DATA)
    balance = "balance 1.0 rcb" if "tiled" in decomp else ""
    (tmp / "in.t").write_text(
        f"""units real
atom_style full
boundary p p p
{decomp}
read_data d.data
{balance}
pair_style zero 5.0
pair_coeff * *
bond_style zero
bond_coeff 1
comm_modify cutoff 6.0
group bead type 1
fix push bead addforce {F_BEAD[0]} {F_BEAD[1]} {F_BEAD[2]}
fix bm all backmap cg_type 1 alpha 0.0001 lambda0 0.0
dump d all custom 1 f.dump id type fx fy fz
dump_modify d sort id format float %.15g
run 0
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
    lines = (tmp / "f.dump").read_text().splitlines()
    head = next(i for i, ln in enumerate(lines) if ln.startswith("ITEM: ATOMS"))
    rows = np.loadtxt(lines[head + 1 :], ndmin=2)
    return {int(r[1]): r[2:5] for r in rows}


@pytest.mark.integration
@pytest.mark.parametrize("decomposition", sorted(DECOMPOSITIONS))
def test_cg_force_reaches_at_atoms(tmp_path: Path, decomposition: str) -> None:
    np_, decomp = DECOMPOSITIONS[decomposition]
    forces = _run(tmp_path, np_, decomp)
    m_bead = sum(M_AT.values())
    for atom_type, mass in M_AT.items():
        want = np.array(F_BEAD) * mass / m_bead
        np.testing.assert_allclose(forces[atom_type], want, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(forces[1], 0.0, atol=1e-12)
