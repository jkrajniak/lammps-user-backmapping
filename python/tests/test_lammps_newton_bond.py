"""Bonded backmap styles give the same forces with ``newton_bond`` off, across MPI ranks.

With ``newton on off`` a bond or angle whose atoms are on two ranks is computed on both; each
rank must add force only to its owned atoms. Six styles added it to the ghost copies as well,
and the reverse communication of ``newton_pair`` then doubled the cross-rank bonded force
(code review F11).

A four-atom chain straddles the boundary of a 2 x 1 x 1 decomposition; the forces of each style
pair with two ranks and ``newton on off`` must equal the serial ones.

Skipped unless ``BACKMAP_LMP`` points to a LAMMPS binary with the backmap package, and for the
parallel runs ``mpirun``.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

LMP_ENV = "BACKMAP_LMP"
XYZ = [(13.8, 10.0, 10.0), (14.8, 10.8, 10.2), (15.6, 10.1, 10.5), (16.7, 10.6, 11.0)]
CASES = {
    "harmonic": ("backmap/harmonic", "1 at 300.0 1.2", "backmap/harmonic", "1 at 60.0 100.0"),
    "gromos_cosine_squared": (
        "backmap/gromos",
        "1 at 42.7 1.13",
        "backmap/cosine/squared",
        "1 at 63.3 111.0",
    ),
    "table": (
        "backmap/table linear 2001",
        "1 at bond.table B",
        "backmap/table linear 1801",
        "1 at angle.table A",
    ),
}
RUNS = {"serial": (1, "on"), "np2_newton_on": (2, "on"), "np2_newton_bond_off": (2, "on off")}


def _tables(d: Path) -> None:
    with (d / "bond.table").open("w") as fh:
        fh.write("B\nN 2001\n\n")
        for i in range(2001):
            r = 0.5 + 3.0 * i / 2000
            fh.write(f"{i + 1} {r:.8f} {50 * (r - 1.2) ** 2:.12g} {-100 * (r - 1.2):.12g}\n")
    with (d / "angle.table").open("w") as fh:
        fh.write("A\nN 1801\n\n")
        for i in range(1801):
            t = 180.0 * i / 1800
            fh.write(f"{i + 1} {t:.8f} {0.002 * (t - 110) ** 2:.12g} {-0.004 * (t - 110):.12g}\n")


def _data() -> str:
    lines = ["newton_bond test", "", "4 atoms", "2 atom types", "3 bonds", "1 bond types"]
    lines += ["2 angles", "1 angle types", "", "0.0 30.0 xlo xhi", "0.0 20.0 ylo yhi"]
    lines += ["0.0 20.0 zlo zhi", "", "Masses", "", "1 12.0", "2 12.0", "", "Atoms # full", ""]
    lines += [f"{k} {k} 1 0.0 {x} {y} {z}" for k, (x, y, z) in enumerate(XYZ, 1)]
    lines += ["", "Bonds", "", "1 1 1 2", "2 1 2 3", "3 1 3 4"]
    lines += ["", "Angles", "", "1 1 1 2 3", "2 1 2 3 4"]
    return "\n".join(lines) + "\n"


def _forces(tmp: Path, case: str, run: str) -> tuple[float, np.ndarray]:
    lmp = os.environ.get(LMP_ENV)
    if not lmp or not Path(lmp).is_file():
        pytest.skip(f"set {LMP_ENV} to a LAMMPS binary built with the backmap package")
    nprocs, newton = RUNS[run]
    bond, bond_coeff, angle, angle_coeff = CASES[case]
    d = tmp / f"{case}_{run}"
    d.mkdir()
    _tables(d)
    (d / "d.data").write_text(_data())
    # the atoms belong to no bead (CG type 2 has no atoms): lambda 1, weight 1
    (d / "in.t").write_text(
        f"""units real
atom_style full
boundary p p p
newton {newton}
processors {nprocs} 1 1
atom_modify map array
read_data d.data
pair_style zero 8.0
pair_coeff * *
bond_style {bond}
bond_coeff {bond_coeff}
angle_style {angle}
angle_coeff {angle_coeff}
special_bonds lj 0.0 0.0 0.0 coul 0.0 0.0 0.0
fix bm all backmap cg_type 2 alpha 0.0001 lambda0 1.0
dump d all custom 1 f.dump id fx fy fz
dump_modify d sort id format float %.15g
thermo_style custom step pe
run 0
print "RESULT $(pe:%.15g)"
"""
    )
    cmd = [lmp, "-in", "in.t", "-log", "log.t", "-screen", "none"]
    if nprocs > 1:
        mpirun = shutil.which("mpirun")
        if mpirun is None:
            pytest.skip("mpirun not available")
        cmd = [mpirun, "-np", str(nprocs), *cmd]
    proc = subprocess.run(cmd, cwd=d, capture_output=True, text=True, check=False)
    log = (d / "log.t").read_text() if (d / "log.t").exists() else ""
    assert proc.returncode == 0, (log + proc.stderr)[-2000:]
    energy = float(next(ln.split()[1] for ln in log.splitlines() if ln.startswith("RESULT ")))
    lines = (d / "f.dump").read_text().splitlines()
    head = next(i for i, ln in enumerate(lines) if ln.startswith("ITEM: ATOMS"))
    return energy, np.loadtxt(lines[head + 1 :], ndmin=2)[:, 1:4]


@pytest.mark.integration
@pytest.mark.parametrize("case", sorted(CASES))
def test_two_ranks_with_newton_bond_off_reproduce_the_serial_forces(
    tmp_path: Path, case: str
) -> None:
    e_ref, f_ref = _forces(tmp_path, case, "serial")
    assert np.abs(f_ref).max() > 1.0
    for run in ("np2_newton_on", "np2_newton_bond_off"):
        energy, forces = _forces(tmp_path, case, run)
        assert energy == pytest.approx(e_ref, rel=1e-10), run
        np.testing.assert_allclose(
            forces, f_ref, rtol=1e-9, atol=1e-9 * np.abs(f_ref).max(), err_msg=run
        )
