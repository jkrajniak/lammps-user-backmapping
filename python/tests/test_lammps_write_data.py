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

# Coefficients by role; the type numbers come from the example's generated
# force-field file, so the test follows the example when it is regenerated.
_BOND = {"at": "at 316.682950 1.530000", "cg": "cg 50.0 4.0"}
_ANGLE = {"at": "at 62.141560 111.0000", "cg": "cg 10.0 150.0"}
_DIHEDRALS = {
    "backmap/ryckaert": {
        "cg": "cg 0.1 -0.2 0.3 0.4 -0.05 -0.6",
        "at": "at 0.717018 -2.217713 2.905285 3.135783 -0.731287 -6.271589",
    },
    "backmap/harmonic": {"cg": "cg 1.5 1 3", "at": "at 0.7 -1 2"},
}


def _roles(ff: str, kind: str) -> dict[int, str]:
    """Type -> "at" or "cg" for one bonded kind, from ``<kind>_coeff`` lines."""
    roles = {}
    for t, role in re.findall(rf"^{kind}_coeff (\d+) \S+ (at|cg)\b", ff, re.MULTILINE):
        roles[int(t)] = role
    return roles


def _example(workdir: Path) -> dict[str, str]:
    ff = (workdir / "pe.ff.lmp").read_text()
    bm = (workdir / "pe.backmap.lmp").read_text()
    fix_line = next(ln for ln in bm.splitlines() if ln.startswith("fix bm all backmap"))
    pair = "".join(
        ln + "\n" for ln in ff.splitlines() if ln.startswith(("pair_style", "pair_coeff"))
    )
    at_pairs = "".join(
        f"pair_coeff {i} {j} {eps} {sig}\n"
        for i, j, eps, sig in re.findall(
            r"^pair_coeff (\d+) (\d+) atomistic (\S+) (\S+)\s*$", ff, re.MULTILINE
        )
    )
    return {
        "ff": ff,
        "pair": pair,
        "at_pairs": at_pairs,
        "cg_types": re.search(r"cg_type ((?:\d+ ?)+)", fix_line).group(1).strip(),
        "fix_bm": re.sub(r"lambda0 \S+", "lambda0 0.5", fix_line),
    }


def _bonded_coeffs(ff: str, dihedral_style: str) -> str:
    lines = []
    for kind, table in (
        ("bond", _BOND),
        ("angle", _ANGLE),
        ("dihedral", _DIHEDRALS[dihedral_style]),
    ):
        roles = _roles(ff, kind)
        assert roles, f"no {kind}_coeff lines in pe.ff.lmp"
        lines += [f"{kind}_coeff {t} {table[role]}" for t, role in sorted(roles.items())]
    return "\n".join(lines) + "\n"


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
{fix_bm}
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

# Stock lj/cut and harmonic write their Coeffs with %g and only the i-i pair
# coefficients (cross terms are re-mixed on read), so they are given again on
# the re-read; only the Dihedral Coeffs are taken from the written file.
_AT_ONLY_COEFFS = """\
pair_coeff * * 0.0 1.0
{at_pairs}\
bond_coeff * 316.682950 1.530000
angle_coeff * 62.141560 111.0000
"""

_AT_ONLY_SETUP = """\
group cg type {cg_types}
delete_atoms group cg bond yes mol no
{coeffs}\
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
    for name in ("pe.data", "pe.ff.lmp", "pe.backmap.lmp"):
        shutil.copy(_PE / name, tmp_path)
    for table in ("table_A_A", "table_A_B", "table_B_B"):
        shutil.copy(_PE / f"{table}.table", tmp_path)
    return tmp_path


@pytest.mark.integration
@pytest.mark.parametrize("dihedral_style", sorted(_DIHEDRALS))
def test_backmap_bonded_coeffs_round_trip(pe_dir: Path, dihedral_style: str) -> None:
    lmp = _lmp()
    ex = _example(pe_dir)
    coeffs = _bonded_coeffs(ex["ff"], dihedral_style)
    first = _run(
        lmp,
        pe_dir,
        _HYBRID.format(
            dihedral_style=dihedral_style,
            data="pe.data",
            pair=ex["pair"],
            fix_bm=ex["fix_bm"],
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
            pair=ex["pair"],
            fix_bm=ex["fix_bm"],
            coeffs="",
            tail="",
        ),
    )
    for key in ("pe", "ebond", "eangle", "edihed"):
        assert second[key] == pytest.approx(first[key], rel=1e-12, abs=1e-9), key


@pytest.mark.integration
def test_ryckaert_coeffs_round_trip(pe_dir: Path) -> None:
    lmp = _lmp()
    ex = _example(pe_dir)
    coeffs = _AT_ONLY_COEFFS.format(at_pairs=ex["at_pairs"])
    setup = _AT_ONLY_SETUP.format(cg_types=ex["cg_types"], coeffs=coeffs)
    first = _run(
        lmp,
        pe_dir,
        _AT_ONLY.format(data="pe.data", setup=setup, tail="write_data written.data\n"),
    )
    written = pe_dir / "written.data"
    assert _empty_coeff_sections(written) == []

    second = _run(
        lmp,
        pe_dir,
        _AT_ONLY.format(data="written.data", setup=coeffs, tail=""),
    )
    assert first["edihed"] != 0.0
    assert second["edihed"] == pytest.approx(first["edihed"], rel=1e-12)
    # delete_atoms renumbers atom IDs, so the reread sums the (overlapping,
    # ~1e10 kcal/mol) start-frame pair energy in a different order.
    assert second["pe"] == pytest.approx(first["pe"], rel=1e-9)
