"""pair_style backmap runs its sub-styles' init_style().

It used to skip it (to keep the sub-styles' neighbor requests out), so state
that only init_style() sets was never set: lj/cut's rRESPA cutoff pointer stayed
uninitialized and init_one() dereferenced it, which segfaulted in about one
run in five, depending on the memory layout of the binary.

The charge check of lj/cut/coul/cut lives in init_style(), so it shows
deterministically whether init_style() ran. The repeated lj/cut run covers the
crash itself.

Skipped unless ``BACKMAP_LMP`` points to a LAMMPS binary with the backmap package.
"""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

import pytest

LMP_ENV = "BACKMAP_LMP"

DATA = """substyle init test

4 atoms
2 atom types

0.0 30.0 xlo xhi
0.0 30.0 ylo yhi
0.0 30.0 zlo zhi

Masses

1 12.0
2 12.0

Atoms # {style}

{atoms}
"""

XYZ = [(10.0, 10.0, 10.0), (14.0, 10.0, 10.0)]


def _atoms(style: str) -> str:
    rows = []
    for k, (x, y, z) in enumerate(XYZ):
        for off, typ in ((0, 1), (2, 2)):
            tag = k + 1 + off
            q = " 0.0" if style == "full" else ""
            rows.append(f"{tag} {k + 1} {typ}{q} {x} {y} {z}")
    return "\n".join(sorted(rows, key=lambda r: int(r.split()[0])))


def _run(tmp_path: Path, style: str, at_style: str) -> subprocess.CompletedProcess:
    lmp = os.environ.get(LMP_ENV)
    if not lmp or not Path(lmp).is_file():
        pytest.skip(f"set {LMP_ENV} to a LAMMPS binary built with the backmap package")
    (tmp_path / "test.data").write_text(DATA.format(style=style, atoms=_atoms(style)))
    (tmp_path / "in.test").write_text(
        f"""units real
atom_style {style}
boundary p p p
read_data test.data
pair_style backmap 10.0 {at_style} 10.0 lj/cut 10.0
pair_coeff 1 1 cg 0.3 3.2
pair_coeff 1 2 none
pair_coeff 2 2 atomistic 0.3 3.2
fix bm all backmap cg_type 1 alpha 0.0001 lambda0 0.0
thermo_style custom step evdwl
run 0
print "RESULT evdwl=$(evdwl:%.12e)"
"""
    )
    return subprocess.run(
        [lmp, "-in", "in.test", "-log", "log.test", "-screen", "none"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )


@pytest.mark.integration
def test_substyle_init_style_runs(tmp_path: Path) -> None:
    """lj/cut/coul/cut's own init_style() error reaches the user."""
    proc = _run(tmp_path, "bond", "lj/cut/coul/cut")
    log = (tmp_path / "log.test").read_text() if (tmp_path / "log.test").exists() else ""
    assert proc.returncode != 0, log[-2000:]
    assert "requires atom attribute q" in log + proc.stdout + proc.stderr, log[-2000:]


@pytest.mark.integration
def test_lj_cut_substyle_init_is_stable(tmp_path: Path) -> None:
    """lj/cut sub-style: no crash from an uninitialized rRESPA pointer (20 runs)."""
    for k in range(20):
        work = tmp_path / f"run_{k}"
        work.mkdir()
        proc = _run(work, "full", "lj/cut")
        log = (work / "log.test").read_text() if (work / "log.test").exists() else ""
        assert proc.returncode == 0, f"run {k}: rc {proc.returncode}\n{log[-2000:]}"
        assert re.search(r"^RESULT evdwl=", log, re.MULTILINE), log[-2000:]


@pytest.mark.integration
@pytest.mark.parametrize(
    ("at_style", "message"),
    [
        ("lj/cut/coul/long 10.0", "needs a long-range solver"),
        ("lj/cut/tip4p/long 1 2 1 1 0.1 10.0", "does not implement single()"),
    ],
)
def test_unsupported_substyle_is_rejected(tmp_path: Path, at_style: str, message: str) -> None:
    """Sub-styles without single() or needing kspace are refused with a clear message."""
    lmp = os.environ.get(LMP_ENV)
    if not lmp or not Path(lmp).is_file():
        pytest.skip(f"set {LMP_ENV} to a LAMMPS binary built with the backmap package")
    (tmp_path / "test.data").write_text(DATA.format(style="full", atoms=_atoms("full")))
    (tmp_path / "in.test").write_text(
        f"""units real
atom_style full
boundary p p p
read_data test.data
pair_style backmap 10.0 {at_style} 10.0 lj/cut 10.0
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
    assert proc.returncode != 0
    assert message in log + proc.stdout + proc.stderr, log[-2000:]
