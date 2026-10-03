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

_HYBRID = """\
units real
atom_style full
boundary p p p
read_data dodecane.data
include dodecane.ff.lmp
{fix_bm}
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
group cg type {cg_types}
delete_atoms group cg bond yes mol no
pair_style lj/cut {cutoff}
pair_coeff * * 0.0 1.0
{at_pairs}\
bond_style zero
bond_coeff *
angle_style zero
angle_coeff *
dihedral_style zero
dihedral_coeff *
special_bonds lj 0.0 0.0 0.0 coul 0.0 0.0 0.0
thermo_style custom step evdwl
run 0
print "RESULT evdwl=$(evdwl:%.10e)"
"""


def _example_inputs(workdir: Path) -> tuple[str, str, str, str]:
    """fix backmap line at lambda = 1, CG types, AT cutoff and plain AT pair coeffs.

    Taken from the example's generated files, so the test follows the example
    when it is regenerated (type numbering, coefficients).
    """
    ff = (workdir / "dodecane.ff.lmp").read_text()
    bm = (workdir / "dodecane.backmap.lmp").read_text()
    fix_line = next(ln for ln in bm.splitlines() if ln.startswith("fix bm all backmap"))
    fix_bm = re.sub(r"lambda0 \S+", "lambda0 1.0", fix_line)
    cg_types = re.search(r"cg_type ((?:\d+ ?)+)", fix_line).group(1).strip()
    cutoff = re.search(r"^pair_style backmap \S+ lj/cut/coul/cut (\S+)", ff, re.MULTILINE).group(1)
    at_pairs = "".join(
        f"pair_coeff {i} {j} {eps} {sig}\n"
        for i, j, eps, sig in re.findall(
            r"^pair_coeff (\d+) (\d+) atomistic (\S+) (\S+)\s*$", ff, re.MULTILINE
        )
    )
    assert at_pairs, "no atomistic pair_coeff lines in dodecane.ff.lmp"
    return fix_bm, cg_types, cutoff, at_pairs


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
    for name in ("dodecane.data", "dodecane.ff.lmp", "dodecane.backmap.lmp"):
        shutil.copy(_DODECANE / name, tmp_path)
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
    fix_bm, cg_types, cutoff, at_pairs = _example_inputs(dodecane_dir)
    hybrid = _run(lmp, dodecane_dir, _HYBRID.format(fix_bm=fix_bm))
    plain = _run(
        lmp,
        dodecane_dir,
        _AT_ONLY.format(cg_types=cg_types, cutoff=cutoff, at_pairs=at_pairs),
    )

    assert hybrid["evdwl"] == pytest.approx(plain["evdwl"], rel=1e-10)
    assert hybrid["peratom"] == pytest.approx(hybrid["evdwl"] + hybrid["ecoul"], rel=1e-10)
