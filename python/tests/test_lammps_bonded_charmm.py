"""LAMMPS regression tests for backmap/charmm, backmap/fourier and improper backmap/harmonic.

Each style must equal the stock LAMMPS style (angle charmm, dihedral fourier,
improper harmonic) times the backmap weight: 1 within a bead, lambda for an
`at` term across beads, 1 - lambda for a `cg` term.

Skipped unless ``BACKMAP_LMP`` points to a LAMMPS binary built with the
backmap, MOLECULE and EXTRA-MOLECULE packages.
"""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

import pytest

LMP_ENV = "BACKMAP_LMP"

# Four AT atoms in a non-planar chain (A).
AT_XYZ = [(10.0, 10.0, 10.0), (11.5, 10.0, 10.0), (12.0, 11.4, 10.0), (13.4, 11.6, 10.9)]

ANGLE = ("K theta0 K_ub r_ub", "35.0 109.5 12.0 2.4")
DIHEDRAL = "3 0.8 1 0.0 0.35 2 180.0 0.2 6 60.0"
IMPROPER = "20.0 5.0"


def _data(n_cg: int) -> str:
    """One molecule: n_cg CG beads (type 1) then 4 AT atoms (type 2)."""
    lines = [
        "backmap charmm-type styles test",
        "",
        f"{n_cg + 4} atoms",
        "1 angles",
        "1 dihedrals",
        "1 impropers",
        "2 atom types",
        "1 angle types",
        "1 dihedral types",
        "1 improper types",
        "",
        "0.0 30.0 xlo xhi",
        "0.0 30.0 ylo yhi",
        "0.0 30.0 zlo zhi",
        "",
        "Masses",
        "",
        "1 24.0",
        "2 12.0",
        "",
        "Atoms # full",
        "",
    ]
    for i in range(n_cg):
        lines.append(f"{i + 1} 1 1 0.0 11.5 10.5 10.2")
    ids = [n_cg + 1 + k for k in range(4)]
    for tag, (x, y, z) in zip(ids, AT_XYZ, strict=True):
        lines.append(f"{tag} 1 2 0.0 {x} {y} {z}")
    a, b, c, d = ids
    lines += [
        "",
        "Angles",
        "",
        f"1 1 {a} {b} {c}",
        "",
        "Dihedrals",
        "",
        f"1 1 {a} {b} {c} {d}",
        "",
        "Impropers",
        "",
        f"1 1 {a} {b} {c} {d}",
        "",
    ]
    return "\n".join(lines)


def _input(styles: str, fix: str, n_cg: int) -> str:
    forces = " ".join(f"f{k}{c}=$(f{c}[{n_cg + 1 + k}]:%.12e)" for k in range(4) for c in "xyz")
    return f"""units real
atom_style full
boundary p p p
read_data test.data
pair_style zero 5.0
pair_coeff * *
{styles}
{fix}
thermo_style custom step eangle edihed eimp
run 0
print "RESULT eangle=$(eangle:%.12e) edihed=$(edihed:%.12e) eimp=$(eimp:%.12e) {forces}"
"""


def _lmp() -> Path:
    path = os.environ.get(LMP_ENV)
    if not path or not Path(path).is_file():
        pytest.skip(f"set {LMP_ENV} to a LAMMPS binary built with the backmap package")
    return Path(path)


def _run(tmp_path: Path, n_cg: int, styles: str, fix: str) -> dict[str, float]:
    (tmp_path / "test.data").write_text(_data(n_cg))
    (tmp_path / "in.test").write_text(_input(styles, fix, n_cg))
    proc = subprocess.run(
        [str(_lmp()), "-in", "in.test", "-log", "log.test", "-screen", "none"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    log = (tmp_path / "log.test").read_text() if (tmp_path / "log.test").exists() else ""
    assert proc.returncode == 0, (log + proc.stdout + proc.stderr)[-3000:]
    match = re.search(r"^RESULT (.*)$", log, re.MULTILINE)
    assert match, log[-2000:]
    return {k: float(v) for k, v in (kv.split("=") for kv in match.group(1).split())}


def _stock(tmp_path: Path) -> dict[str, float]:
    styles = f"""angle_style charmm
angle_coeff 1 {ANGLE[1]}
dihedral_style fourier
dihedral_coeff 1 {DIHEDRAL}
improper_style harmonic
improper_coeff 1 {IMPROPER}"""
    return _run(tmp_path / "stock", 1, styles, "")


def _backmap(tmp_path: Path, n_cg: int, tag: str, lam: float) -> dict[str, float]:
    styles = f"""angle_style backmap/charmm
angle_coeff 1 {tag} {ANGLE[1]}
dihedral_style backmap/fourier
dihedral_coeff 1 {tag} {DIHEDRAL}
improper_style backmap/harmonic
improper_coeff 1 {tag} {IMPROPER}"""
    fix = f"fix bm all backmap cg_type 1 alpha 0.0001 lambda0 {lam}"
    return _run(tmp_path / f"bm_{n_cg}_{tag}_{lam}", n_cg, styles, fix)


@pytest.mark.integration
@pytest.mark.parametrize(
    ("n_cg", "tag", "lam", "weight"),
    [
        (1, "at", 0.3, 1.0),  # all four atoms in one bead: never weighted
        (2, "at", 1.0, 1.0),  # two beads of two atoms: across beads
        (2, "at", 0.5, 0.5),
        (2, "cg", 0.25, 0.75),
    ],
)
def test_equals_stock_style_times_weight(
    tmp_path: Path, n_cg: int, tag: str, lam: float, weight: float
) -> None:
    for sub in ("stock", f"bm_{n_cg}_{tag}_{lam}"):
        (tmp_path / sub).mkdir()
    ref = _stock(tmp_path)
    res = _backmap(tmp_path, n_cg, tag, lam)
    for term in ("eangle", "edihed", "eimp"):
        assert ref[term] != 0.0
        assert res[term] == pytest.approx(weight * ref[term], rel=1e-10), term
    for key in (f"f{k}{c}" for k in range(4) for c in "xyz"):
        assert res[key] == pytest.approx(weight * ref[key], rel=1e-9, abs=1e-12), key
