"""LAMMPS regression tests for energies reported by the backmap styles.

These run a LAMMPS binary built with the backmap package. They are skipped
unless ``BACKMAP_LMP`` points to that binary.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

LMP_ENV = "BACKMAP_LMP"
_REPO = Path(__file__).resolve().parents[2]
_DODECANE = _REPO / "examples" / "dodecane" / "large"

_AT_PAIR_COEFFS = """\
pair_coeff 3 3 {tag}0.207266 3.748000
pair_coeff 3 4 {tag}0.156387 3.826500
pair_coeff 4 4 {tag}0.117997 3.905000
"""

_HYBRID = """\
units real
atom_style full
boundary p p p
read_data dodecane.data
pair_style backmap 14.00 lj/cut/coul/cut 14.00 9.00 14.00 table linear 1000
pair_coeff 1 1 cg table_A_A.table ENTRY
pair_coeff 1 2 cg table_A_B.table ENTRY
pair_coeff 1 3 none
pair_coeff 1 4 none
pair_coeff 2 2 cg table_B_B.table ENTRY
pair_coeff 2 3 none
pair_coeff 2 4 none
{at_pairs}\
bond_style hybrid harmonic backmap/harmonic backmap/table linear 1000
bond_coeff 1 harmonic 800.000883 1.530000
bond_coeff 2 backmap/table cg table_b1.table ENTRY
bond_coeff 3 backmap/table cg table_b1.table ENTRY
bond_coeff 4 backmap/harmonic at 800.000883 1.530000
angle_style backmap/harmonic
angle_coeff 1 at 126.673180 111.0000
special_bonds lj 0.0 0.0 0.0 coul 0.0 0.0 0.0
fix bm all backmap cg_type 1 2 alpha 0.0001 lambda0 1.0
compute pa all pe/atom pair
compute pa_sum all reduce sum c_pa
thermo_style custom step evdwl ecoul c_pa_sum
run 0
print "RESULT evdwl=$(evdwl:%.10e) ecoul=$(ecoul:%.10e) peratom=$(c_pa_sum:%.10e)"
"""

_AT_ONLY = """\
units real
atom_style full
boundary p p p
read_data dodecane.data
group cg type 1 2
delete_atoms group cg bond yes mol no
pair_style lj/cut 14.00
pair_coeff * * 0.0 1.0
{at_pairs}\
bond_style harmonic
bond_coeff * 800.000883 1.530000
angle_style harmonic
angle_coeff 1 126.673180 111.0000
special_bonds lj 0.0 0.0 0.0 coul 0.0 0.0 0.0
thermo_style custom step evdwl
run 0
print "RESULT evdwl=$(evdwl:%.10e)"
"""


def _lmp() -> Path:
    path = os.environ.get(LMP_ENV)
    if not path or not Path(path).is_file():
        pytest.skip(f"set {LMP_ENV} to a LAMMPS binary built with the backmap package")
    return Path(path)


def _run(lmp: Path, workdir: Path, script: str) -> dict[str, float]:
    (workdir / "in.test").write_text(script)
    proc = subprocess.run(
        [str(lmp), "-in", "in.test", "-log", "log.test", "-screen", "none"],
        cwd=workdir,
        capture_output=True,
        text=True,
        check=False,
    )
    log = (workdir / "log.test").read_text()
    assert proc.returncode == 0, log[-2000:]
    match = re.search(r"^RESULT (.*)$", log, re.MULTILINE)
    assert match, log[-2000:]
    return {k: float(v) for k, v in (kv.split("=") for kv in match.group(1).split())}


@pytest.fixture
def dodecane_dir(tmp_path: Path) -> Path:
    shutil.copy(_DODECANE / "dodecane.data", tmp_path)
    for table in _DODECANE.glob("table_*.table"):
        shutil.copy(table, tmp_path)
    return tmp_path


@pytest.mark.integration
def test_pair_backmap_energy_matches_plain_lj_at_lambda_one(dodecane_dir: Path) -> None:
    """At lambda = 1 the hybrid pair energy must equal the plain AT pair energy.

    Regression for a double tally in PairBackmap::compute(), which reported
    exactly twice the pair energy while the forces were correct.
    """
    lmp = _lmp()
    hybrid = _run(
        lmp,
        dodecane_dir,
        _HYBRID.format(at_pairs=_AT_PAIR_COEFFS.format(tag="atomistic ")),
    )
    plain = _run(lmp, dodecane_dir, _AT_ONLY.format(at_pairs=_AT_PAIR_COEFFS.format(tag="")))

    assert hybrid["evdwl"] == pytest.approx(plain["evdwl"], rel=1e-10)
    assert hybrid["peratom"] == pytest.approx(hybrid["evdwl"] + hybrid["ecoul"], rel=1e-10)
