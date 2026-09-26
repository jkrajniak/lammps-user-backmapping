"""LAMMPS regression tests: bonded backmap coefficients survive write_data/read_data.

A style that sets ``writedata = 1`` without implementing ``write_data()``
leaves an empty ``* Coeffs`` section that ``read_data`` rejects. Each test
writes a data file, reads it back without re-specifying the bonded
coefficients, and requires the same energies.

Skipped unless ``BACKMAP_LMP`` points to a LAMMPS binary built with the
backmap package.
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
_PE = _REPO / "examples" / "pe" / "large"

_PAIR = """\
pair_style backmap 14.00 lj/cut/coul/cut 14.00 9.00 14.00 table linear 1000
pair_coeff 1 1 cg table_A_A.table ENTRY
pair_coeff 1 2 cg table_A_B.table ENTRY
pair_coeff 1 3 none
pair_coeff 1 4 none
pair_coeff 2 2 cg table_B_B.table ENTRY
pair_coeff 2 3 none
pair_coeff 2 4 none
pair_coeff 3 3 atomistic 0.207266 3.748000
pair_coeff 3 4 atomistic 0.156387 3.826500
pair_coeff 4 4 atomistic 0.117997 3.905000
"""

_BONDED_COEFFS = """\
bond_coeff 1 at 316.682950 1.530000
bond_coeff 2 cg 50.0 4.0
bond_coeff 3 cg 50.0 4.0
bond_coeff 4 at 316.682950 1.530000
angle_coeff 1 cg 10.0 150.0
angle_coeff 2 cg 10.0 150.0
angle_coeff 3 at 62.141560 111.0000
{dihedral_coeffs}\
"""

_DIHEDRALS = {
    "backmap/ryckaert": (
        "dihedral_coeff 1 cg 0.1 -0.2 0.3 0.4 -0.05 -0.6\n"
        "dihedral_coeff 2 at 0.717018 -2.217713 2.905285 3.135783 -0.731287 -6.271589\n"
    ),
    "backmap/harmonic": "dihedral_coeff 1 cg 1.5 1 3\ndihedral_coeff 2 at 0.7 -1 2\n",
}

_HYBRID = """\
units real
atom_style full
boundary p p p
bond_style backmap/harmonic
angle_style backmap/harmonic
dihedral_style {dihedral_style}
read_data {data}
{pair}\
{coeffs}\
special_bonds lj 0.0 0.0 0.0 coul 0.0 0.0 0.0
fix bm all backmap cg_type 1 2 alpha 0.0001 lambda0 0.5
thermo_style custom step pe ebond eangle edihed evdwl
run 0
print "RESULT pe=$(pe:%.12e) ebond=$(ebond:%.12e) eangle=$(eangle:%.12e) edihed=$(edihed:%.12e)"
{tail}\
"""

_AT_ONLY = """\
units real
atom_style full
boundary p p p
pair_style lj/cut 14.00
bond_style harmonic
angle_style harmonic
dihedral_style ryckaert
read_data {data}
{setup}\
special_bonds lj 0.0 0.0 0.0 coul 0.0 0.0 0.0
thermo_style custom step pe edihed
run 0
print "RESULT pe=$(pe:%.12e) edihed=$(edihed:%.12e)"
{tail}\
"""

_AT_ONLY_SETUP = """\
group cg type 1 2
delete_atoms group cg bond yes mol no
pair_coeff * * 0.0 1.0
pair_coeff 3 3 0.207266 3.748000
pair_coeff 3 4 0.156387 3.826500
pair_coeff 4 4 0.117997 3.905000
bond_coeff * 316.682950 1.530000
angle_coeff * 62.141560 111.0000
dihedral_coeff * 0.717018 -2.217713 2.905285 3.135783 -0.731287 -6.271589
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


def _empty_coeff_sections(data: Path) -> list[str]:
    lines = [line.strip() for line in data.read_text().splitlines()]
    empty = []
    for idx, line in enumerate(lines):
        if not line.split("#")[0].strip().endswith("Coeffs"):
            continue
        following = next((nxt for nxt in lines[idx + 1 :] if nxt), "")
        if not following[:1].isdigit():
            empty.append(line)
    return empty


@pytest.fixture
def pe_dir(tmp_path: Path) -> Path:
    shutil.copy(_PE / "pe.data", tmp_path)
    for table in ("table_A_A", "table_A_B", "table_B_B"):
        shutil.copy(_PE / f"{table}.table", tmp_path)
    return tmp_path


@pytest.mark.integration
@pytest.mark.parametrize("dihedral_style", sorted(_DIHEDRALS))
def test_backmap_bonded_coeffs_round_trip(pe_dir: Path, dihedral_style: str) -> None:
    lmp = _lmp()
    coeffs = _BONDED_COEFFS.format(dihedral_coeffs=_DIHEDRALS[dihedral_style])
    first = _run(
        lmp,
        pe_dir,
        _HYBRID.format(
            dihedral_style=dihedral_style,
            data="pe.data",
            pair=_PAIR,
            coeffs=coeffs,
            tail="write_data written.data\n",
        ),
    )
    written = pe_dir / "written.data"
    assert _empty_coeff_sections(written) == []

    second = _run(
        lmp,
        pe_dir,
        _HYBRID.format(
            dihedral_style=dihedral_style,
            data="written.data",
            pair=_PAIR,
            coeffs="",
            tail="",
        ),
    )
    for key in ("pe", "ebond", "eangle", "edihed"):
        assert second[key] == pytest.approx(first[key], rel=1e-12, abs=1e-9), key


@pytest.mark.integration
def test_ryckaert_coeffs_round_trip(pe_dir: Path) -> None:
    lmp = _lmp()
    first = _run(
        lmp,
        pe_dir,
        _AT_ONLY.format(data="pe.data", setup=_AT_ONLY_SETUP, tail="write_data written.data\n"),
    )
    written = pe_dir / "written.data"
    assert _empty_coeff_sections(written) == []

    second = _run(
        lmp,
        pe_dir,
        _AT_ONLY.format(data="written.data", setup="", tail=""),
    )
    assert first["edihed"] != 0.0
    assert second["edihed"] == pytest.approx(first["edihed"], rel=1e-12)
    # delete_atoms renumbers atom IDs, so the reread sums the (overlapping,
    # ~1e10 kcal/mol) start-frame pair energy in a different order.
    assert second["pe"] == pytest.approx(first["pe"], rel=1e-9)
