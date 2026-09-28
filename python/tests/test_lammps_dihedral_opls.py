"""LAMMPS regression tests for dihedral_style backmap/opls.

The style must equal the stock LAMMPS dihedral_style opls times the backmap
weight: 1 within a bead, lambda for an `at` term across beads, 1 - lambda for
a `cg` term. A data-file round trip (write_data, read back) must keep K1..K4.

Skipped unless ``BACKMAP_LMP`` points to a LAMMPS binary built with the
backmap and MOLECULE packages.
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

# OPLS-AA-like K1..K4 (kcal/mol), all four terms nonzero and of both signs.
OPLS = "1.3 -0.05 0.2 0.1"


def _data(n_cg: int) -> str:
    """One molecule: n_cg CG beads (type 1) then 4 AT atoms (type 2)."""
    lines = [
        "backmap opls test",
        "",
        f"{n_cg + 4} atoms",
        "1 dihedrals",
        "2 atom types",
        "1 dihedral types",
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
        "Dihedrals",
        "",
        f"1 1 {a} {b} {c} {d}",
        "",
    ]
    return "\n".join(lines)


def _input(styles: str, fix: str, n_cg: int, extra: str = "") -> str:
    forces = " ".join(f"f{k}{c}=$(f{c}[{n_cg + 1 + k}]:%.12e)" for k in range(4) for c in "xyz")
    return f"""units real
atom_style full
boundary p p p
read_data test.data
pair_style zero 5.0
pair_coeff * *
{styles}
{fix}
thermo_style custom step edihed
run 0
print "RESULT edihed=$(edihed:%.12e) {forces}"
{extra}
"""


def _lmp() -> Path:
    path = os.environ.get(LMP_ENV)
    if not path or not Path(path).is_file():
        pytest.skip(f"set {LMP_ENV} to a LAMMPS binary built with the backmap package")
    return Path(path)


def _run(tmp_path: Path, n_cg: int, styles: str, fix: str, extra: str = "") -> dict[str, float]:
    (tmp_path / "test.data").write_text(_data(n_cg))
    (tmp_path / "in.test").write_text(_input(styles, fix, n_cg, extra))
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
    styles = f"dihedral_style opls\ndihedral_coeff 1 {OPLS}"
    return _run(tmp_path / "stock", 1, styles, "")


def _fix(lam: float) -> str:
    return f"fix bm all backmap cg_type 1 alpha 0.0001 lambda0 {lam}"


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
def test_equals_stock_opls_times_weight(
    tmp_path: Path, n_cg: int, tag: str, lam: float, weight: float
) -> None:
    sub = f"bm_{n_cg}_{tag}_{lam}"
    for d in ("stock", sub):
        (tmp_path / d).mkdir()
    ref = _stock(tmp_path)
    styles = f"dihedral_style backmap/opls\ndihedral_coeff 1 {tag} {OPLS}"
    res = _run(tmp_path / sub, n_cg, styles, _fix(lam))
    assert ref["edihed"] != 0.0
    assert res["edihed"] == pytest.approx(weight * ref["edihed"], rel=1e-10)
    for key in (f"f{k}{c}" for k in range(4) for c in "xyz"):
        assert res[key] == pytest.approx(weight * ref[key], rel=1e-9, abs=1e-12), key


@pytest.mark.integration
def test_write_data_keeps_opls_coefficients(tmp_path: Path) -> None:
    (tmp_path / "w").mkdir()
    styles = f"dihedral_style backmap/opls\ndihedral_coeff 1 at {OPLS}"
    _run(tmp_path / "w", 1, styles, _fix(1.0), extra="write_data out.data")
    text = (tmp_path / "w" / "out.data").read_text()
    after_header = text.split("Dihedral Coeffs", 1)[1].split("\n", 1)[1]
    first = next(ln for ln in after_header.splitlines() if ln.strip())
    assert first.split() == ["1", "at", *OPLS.split()]
