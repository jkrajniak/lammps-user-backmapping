"""Integration tests — verify backmap-prep generates correct output for each example."""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import ClassVar

import pytest

from backmap_prep.cli import main

EXAMPLES_DIR = Path(__file__).resolve().parents[2] / "examples"

EXAMPLES = ["dodecane", "dodecane-lammps-cg", "pe", "pe4", "pe_10", "pe_aa", "melamine"]


@pytest.fixture(params=EXAMPLES)
def example_workdir(request, tmp_path: Path) -> Path:
    name = request.param
    src = EXAMPLES_DIR / name
    if not src.exists():
        pytest.skip(f"Example {name} not found at {src}")
    dst = tmp_path / name
    shutil.copytree(src, dst)
    ff_src = EXAMPLES_DIR / "forcefield"
    if ff_src.is_dir():
        shutil.copytree(ff_src, tmp_path / "forcefield", dirs_exist_ok=True)
    return dst


class TestExampleGeneration:
    def test_backmap_prep_succeeds(self, example_workdir: Path) -> None:
        settings = example_workdir / "settings.yaml"
        assert settings.exists(), f"Missing settings.yaml in {example_workdir}"
        old_cwd = os.getcwd()
        try:
            os.chdir(example_workdir)
            result = main([str(settings)])
        finally:
            os.chdir(old_cwd)
        assert result == 0, f"backmap-prep failed for {example_workdir.name}"

    def test_data_file_generated(self, example_workdir: Path) -> None:
        settings = example_workdir / "settings.yaml"
        old_cwd = os.getcwd()
        try:
            os.chdir(example_workdir)
            main([str(settings)])
        finally:
            os.chdir(old_cwd)
        data_files = list(example_workdir.glob("*.data"))
        assert len(data_files) >= 1, "Expected at least 1 .data file"

    def test_input_script_generated(self, example_workdir: Path) -> None:
        settings = example_workdir / "settings.yaml"
        old_cwd = os.getcwd()
        try:
            os.chdir(example_workdir)
            main([str(settings)])
        finally:
            os.chdir(old_cwd)
        in_files = list(example_workdir.glob("in.*"))
        in_files = [f for f in in_files if not f.name.endswith(".yaml")]
        assert len(in_files) >= 1, f"Expected at least 1 input script, found {len(in_files)}"

    def test_data_file_has_atoms_and_bonds(self, example_workdir: Path) -> None:
        settings = example_workdir / "settings.yaml"
        old_cwd = os.getcwd()
        try:
            os.chdir(example_workdir)
            main([str(settings)])
        finally:
            os.chdir(old_cwd)
        data_files = list(example_workdir.glob("*.data"))
        content = data_files[0].read_text()
        assert "atoms" in content
        assert "bonds" in content
        assert "Masses" in content
        assert "Atoms" in content


class TestDeterministicOutput:
    @pytest.mark.parametrize("example_name", EXAMPLES)
    def test_identical_runs(self, example_name: str, tmp_path: Path) -> None:
        src = EXAMPLES_DIR / example_name
        if not src.exists():
            pytest.skip(f"Example {example_name} not found")

        ff_src = EXAMPLES_DIR / "forcefield"
        if ff_src.is_dir():
            shutil.copytree(ff_src, tmp_path / "forcefield", dirs_exist_ok=True)
        run1 = tmp_path / "run1"
        run2 = tmp_path / "run2"
        shutil.copytree(src, run1)
        shutil.copytree(src, run2)

        for workdir in (run1, run2):
            old_cwd = os.getcwd()
            try:
                os.chdir(workdir)
                main([str(workdir / "settings.yaml")])
            finally:
                os.chdir(old_cwd)

        for f1 in run1.glob("*.data"):
            f2 = run2 / f1.name
            assert f2.exists()
            assert f1.read_text() == f2.read_text(), (
                f"Non-deterministic output for {example_name}/{f1.name}"
            )

        for f1 in run1.glob("in.*"):
            if f1.name.endswith(".yaml"):
                continue
            f2 = run2 / f1.name
            assert f2.exists()
            assert f1.read_text() == f2.read_text(), (
                f"Non-deterministic output for {example_name}/{f1.name}"
            )


class TestLammpsNativeCgParity:
    """dodecane-lammps-cg must build to the same hybrid system as dodecane.

    dodecane-lammps-cg/dodecane_cg.data was derived (see its README) from a
    hybrid build of dodecane/'s current cg_conf.gro/topol_cg.top, so a fresh
    build of both examples should differ only in the cosmetic type-name
    convention (A/B vs the LAMMPS numeric type IDs "1"/"2") documented there.
    """

    _NAME_SUBS: ClassVar[dict[str, str]] = {
        "# A (CG)": "# 1 (CG)",
        "# B (CG)": "# 2 (CG)",
        "table_A_A": "table_1_1",
        "table_A_B": "table_1_2",
        "table_B_B": "table_2_2",
    }

    def _normalize(self, text: str) -> str:
        for old, new in self._NAME_SUBS.items():
            text = text.replace(old, new)
        return text

    def _build(self, name: str, tmp_path: Path) -> Path:
        src = EXAMPLES_DIR / name
        dst = tmp_path / name
        shutil.copytree(src, dst)
        old_cwd = os.getcwd()
        try:
            os.chdir(dst)
            main([str(dst / "settings.yaml")])
        finally:
            os.chdir(old_cwd)
        return dst

    def test_hybrid_data_matches_after_type_name_normalization(self, tmp_path: Path) -> None:
        gromacs_dir = self._build("dodecane", tmp_path)
        lammps_dir = self._build("dodecane-lammps-cg", tmp_path)

        gromacs_data = self._normalize((gromacs_dir / "dodecane.data").read_text())
        lammps_data = (lammps_dir / "dodecane.data").read_text()
        assert gromacs_data == lammps_data

    def test_input_script_matches_after_type_name_normalization(self, tmp_path: Path) -> None:
        gromacs_dir = self._build("dodecane", tmp_path)
        lammps_dir = self._build("dodecane-lammps-cg", tmp_path)

        gromacs_in = self._normalize((gromacs_dir / "in.dodecane").read_text())
        lammps_in = (lammps_dir / "in.dodecane").read_text()
        assert gromacs_in == lammps_in
