"""Integration tests — verify backmap-prep generates correct output for each example."""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest

from backmap_prep.cli import main

EXAMPLES_DIR = Path(__file__).resolve().parents[2] / "examples"

EXAMPLES = ["dodecane", "pe", "pe4", "pe_10", "pe_aa", "melamine"]


@pytest.fixture(params=EXAMPLES)
def example_workdir(request, tmp_path: Path) -> Path:
    name = request.param
    src = EXAMPLES_DIR / name
    if not src.exists():
        pytest.skip(f"Example {name} not found at {src}")
    dst = tmp_path / name
    shutil.copytree(src, dst)
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
