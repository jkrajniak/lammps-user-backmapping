"""AT-only force field and data files derived from a hybrid system."""

from __future__ import annotations

import os
import re
import shutil
from pathlib import Path

import pytest

from backmap_prep.cli import main
from backmap_prep.schema import load_settings

EXAMPLE = Path(__file__).resolve().parents[2] / "examples" / "pe"


@pytest.fixture
def pe_dir(tmp_path: Path) -> Path:
    dst = tmp_path / "pe"
    shutil.copytree(EXAMPLE, dst)
    return dst


def _run(workdir: Path, *args: str) -> None:
    old = os.getcwd()
    try:
        os.chdir(workdir)
        assert main([*args]) == 0
    finally:
        os.chdir(old)


def _counts(path: Path) -> dict[str, int]:
    text = path.read_text()
    return {
        key: int(m.group(1))
        for key in ("atoms", "bonds", "angles", "dihedrals", "atom types")
        if (m := re.search(rf"^(\d+) {key}$", text, re.MULTILINE))
    }


def test_at_forcefield_matches_hybrid_coefficients(pe_dir: Path) -> None:
    settings = pe_dir / "settings.yaml"
    _run(pe_dir, "build", str(settings))
    at_ff = (pe_dir / "pe.at.ff.lmp").read_text()
    hybrid_ff = (pe_dir / "pe.ff.lmp").read_text()

    assert "pair_style lj/cut/coul/cut" in at_ff
    assert "backmap/" not in at_ff  # no backmap styles
    # Same AT coefficients as the hybrid, now under plain styles.
    for pattern in (r"316\.6826004 1\.53", r"62\.1414914 111", r"2\.217710325 -2\.905282027"):
        assert re.search(pattern, hybrid_ff), pattern
        assert re.search(pattern, at_ff), pattern


def test_at_system_from_frame_and_reference_share_numbering(pe_dir: Path) -> None:
    settings = pe_dir / "settings.yaml"
    _run(pe_dir, "build", str(settings))
    _run(pe_dir, "at-system", str(settings), "--from", str(pe_dir / "pe.data"))
    _run(pe_dir, "at-system", str(settings), "--reference", "10", "--seed", "7")

    n_mol = 10
    frame = _counts(pe_dir / "pe_at.data")
    reference = _counts(pe_dir / "pe_at_ref.data")
    expected = {"atoms": 100 * n_mol, "bonds": 99 * n_mol, "angles": 98 * n_mol}
    expected |= {"dihedrals": 97 * n_mol, "atom types": 2}
    assert frame == expected
    assert reference == expected
    assert load_settings(settings).output.prefix == "pe"


def test_at_system_needs_one_mode(pe_dir: Path) -> None:
    settings = pe_dir / "settings.yaml"
    old = os.getcwd()
    try:
        os.chdir(pe_dir)
        assert main(["at-system", str(settings)]) == 1
    finally:
        os.chdir(old)
