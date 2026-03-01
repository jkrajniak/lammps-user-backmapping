"""Tests for backmap_prep.cli — command-line interface."""

from __future__ import annotations

from textwrap import dedent
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

from backmap_prep.cli import main


def _write_full_example(base: Path) -> Path:
    """Write all files needed for a minimal end-to-end CLI run."""
    cg_gro = dedent("""\
        CG system
            1
            1MOL  B1    1   0.500   0.500   0.500
           5.00000   5.00000   5.00000
    """)
    (base / "cg.gro").write_text(cg_gro)

    cg_top = dedent("""\
        [ atomtypes ]
        CG1  72.0  0.0  V  0.47  3.5

        [ moleculetype ]
        MOL  3

        [ atoms ]
        1  CG1  1  MOL  B1  1  0.0  72.0

        [ molecules ]
        MOL 1
    """)
    (base / "cg.top").write_text(cg_top)

    at_gro = dedent("""\
        AT template
            2
            1MOL  C1    1   0.500   0.500   0.500
            1MOL  C2    2   0.510   0.500   0.500
           5.00000   5.00000   5.00000
    """)
    (base / "at.gro").write_text(at_gro)

    at_top = dedent("""\
        [ atomtypes ]
        CH2  14.0  0.0  A  0.395  0.382

        [ moleculetype ]
        TestMol  3

        [ atoms ]
        1  CH2  1  MOL  C1  1  0.0  14.0
        2  CH2  1  MOL  C2  1  0.0  14.0

        [ bonds ]
        1  2  1  0.154  200000.0
    """)
    (base / "at.top").write_text(at_top)

    import yaml

    settings = {
        "molecules": [
            {
                "name": "TestMol",
                "source": {"coordinates": "at.gro", "topology": "at.top"},
                "beads": [{"name": "B1", "type": "CG1", "atoms": ["C1", "C2"]}],
            }
        ],
        "cg_system": {"coordinates": "cg.gro", "topology": "cg.top"},
        "output": {"prefix": "test_out"},
    }
    settings_path = base / "settings.yaml"
    settings_path.write_text(yaml.dump(settings))
    return settings_path


class TestCLI:
    def test_missing_file_returns_1(self, tmp_path: Path) -> None:
        result = main([str(tmp_path / "nonexistent.yaml")])
        assert result == 1

    def test_successful_run(self, tmp_path: Path) -> None:
        settings_path = _write_full_example(tmp_path)
        result = main([str(settings_path)])
        assert result == 0
        assert (tmp_path / "test_out.data").exists()
        assert (tmp_path / "in.test_out").exists()

    def test_output_prefix_override(self, tmp_path: Path) -> None:
        settings_path = _write_full_example(tmp_path)
        result = main([str(settings_path), "--output-prefix", "custom"])
        assert result == 0
        assert (tmp_path / "custom.data").exists()
        assert (tmp_path / "in.custom").exists()

    def test_data_file_content(self, tmp_path: Path) -> None:
        settings_path = _write_full_example(tmp_path)
        main([str(settings_path)])
        content = (tmp_path / "test_out.data").read_text()
        assert "atoms" in content
        assert "Masses" in content
        assert "Atoms # full" in content

    def test_input_file_content(self, tmp_path: Path) -> None:
        settings_path = _write_full_example(tmp_path)
        main([str(settings_path)])
        content = (tmp_path / "in.test_out").read_text()
        assert "units real" in content
        assert "fix bm all backmap" in content
