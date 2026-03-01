"""Shared test fixtures for backmap-prep tests."""

from __future__ import annotations

from textwrap import dedent
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

import pytest


@pytest.fixture
def gro_content() -> str:
    """Minimal 3-atom GRO file content (fixed-width GRO format)."""
    lines = [
        "Test system",
        "    3",
        "    1RES     C1    1   0.100   0.200   0.300  0.0100  0.0200  0.0300",
        "    1RES     C2    2   0.400   0.500   0.600",
        "    1RES     C3    3   0.700   0.800   0.900  0.0400  0.0500  0.0600",
        "   5.00000   5.00000   5.00000",
    ]
    return "\n".join(lines) + "\n"


@pytest.fixture
def gro_file(tmp_path: Path, gro_content: str) -> Path:
    """Write a minimal GRO file and return its path."""
    p = tmp_path / "test.gro"
    p.write_text(gro_content)
    return p


@pytest.fixture
def top_content() -> str:
    """Minimal topology file content with one molecule type."""
    return dedent("""\
        ; minimal topology
        [ defaults ]
        1  2  yes  0.5  0.8333

        [ atomtypes ]
        CG1  72.0  0.0  A  0.47  3.5
        AT1  12.0  0.0  A  0.34  0.4

        [ moleculetype ]
        MOL  3

        [ atoms ]
        1  CG1  1  MOL  B1  1  0.0  72.0
        2  AT1  1  MOL  C1  1  0.0  12.0
        3  AT1  1  MOL  C2  1  0.0  12.0

        [ bonds ]
        2  3  1  0.154  200000.0

        [ angles ]
        1  2  3  1  120.0  500.0

        [ molecules ]
        MOL 1
    """)


@pytest.fixture
def top_file(tmp_path: Path, top_content: str) -> Path:
    """Write a minimal TOP file and return its path."""
    p = tmp_path / "test.top"
    p.write_text(top_content)
    return p


@pytest.fixture
def xvg_content() -> str:
    """Minimal GROMACS .xvg file content."""
    return dedent("""\
        # Comment line
        @ title "Bond potential"
        @ xaxis label "r (nm)"
        0.10  10.0  -200.0
        0.15  5.0   -100.0
        0.20  1.0   -20.0
    """)


@pytest.fixture
def xvg_file(tmp_path: Path, xvg_content: str) -> Path:
    """Write a minimal XVG file and return its path."""
    p = tmp_path / "table.xvg"
    p.write_text(xvg_content)
    return p


@pytest.fixture
def minimal_settings_dict() -> dict:
    """Minimal valid settings dictionary."""
    return {
        "molecules": [
            {
                "name": "TestMol",
                "source": {
                    "coordinates": "at.gro",
                    "topology": "at.top",
                },
                "beads": [
                    {
                        "name": "B1",
                        "type": "CG1",
                        "atoms": ["C1", "C2"],
                    }
                ],
            }
        ],
        "cg_system": {
            "coordinates": "cg.gro",
            "topology": "cg.top",
        },
    }


@pytest.fixture
def settings_yaml(tmp_path: Path, minimal_settings_dict: dict) -> Path:
    """Write a minimal YAML settings file and return its path."""
    import yaml

    p = tmp_path / "settings.yaml"
    p.write_text(yaml.dump(minimal_settings_dict))
    return p
