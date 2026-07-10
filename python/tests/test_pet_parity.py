"""Tier A parity test for PET (Dacron) network hybrid build.

Verifies that backmap-prep's v2 YAML path reproduces the bakery reference
hyb_conf.gro (byte-identical) and hyb_topol.top (matching section counts)
for the 3-species PET polyester network.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

from backmap_prep.network.api import build_hybrid_gromacs
from backmap_prep.network.compare_topology import compare_topology_files
from backmap_prep.schema import load_settings

PET_ENV = "BACKMAP_PETE_DIR"
# python/tests/test_pet_parity.py -> parents[3] = user-backmapping (main checkout)
_HUB = Path(__file__).resolve().parents[3]
_DEFAULT_PET = (
    _HUB / "paper-reverse-mapping-polymer-networks" / "preparation" / "dacron" / "backmapping"
)


def pete_data_dir() -> Path | None:
    """Return the PET/Dacron paper-data backmapping dir if available."""
    override = os.environ.get(PET_ENV)
    if override:
        path = Path(override)
        return path if path.is_dir() else None
    if _DEFAULT_PET.is_dir() and (_DEFAULT_PET / "settings.xml").is_file():
        return _DEFAULT_PET
    return None


def _section_lines(top_path: Path, section: str) -> list[str]:
    """Return the non-empty, non-comment lines of a GROMACS topology section."""
    lines = top_path.read_text().splitlines()
    out: list[str] = []
    on = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("["):
            on = stripped == f"[ {section} ]"
            continue
        if on and stripped and not stripped.startswith(";"):
            out.append(stripped)
    return out


@pytest.mark.integration
def test_pet_build_hybrid_v2_parity(tmp_path: Path) -> None:
    """v2 YAML reproduces bakery hyb_conf.gro + hyb_topol.top for PET/Dacron."""
    data_dir = pete_data_dir()
    if data_dir is None:
        pytest.skip(f"PET/Dacron fixtures not found (set {PET_ENV} or clone paper-data repo)")

    v2_yaml = (
        Path(__file__).resolve().parents[2] / "examples" / "pet" / "large" / "settings.v2.yaml"
    )
    if not v2_yaml.is_file():
        pytest.skip("examples/pet/large/settings.v2.yaml missing")

    settings = load_settings(v2_yaml)

    ref_gro = data_dir / "hyb_conf.gro"
    ref_top = data_dir / "hyb_topol.top"
    if not ref_gro.is_file() or not ref_top.is_file():
        pytest.skip("bakery hyb_conf.gro / hyb_topol.top reference missing")
    # Back up to tmp_path (outside data_dir) so the build's own _*.1_ backup
    # numbering does not trip on our backup file (files_io.py:int(split[-1])).
    backup_gro = tmp_path / "hyb_conf.gro.ref"
    backup_top = tmp_path / "hyb_topol.top.ref"
    shutil.copy2(ref_gro, backup_gro)
    shutil.copy2(ref_top, backup_top)

    try:
        result = build_hybrid_gromacs(settings, base_dir=data_dir, chain_rng_seed=42)
        assert result.n_atoms == 33757, result.n_atoms
        # Coordinates byte-identical to bakery reference.
        assert result.coordinates_path.read_bytes() == backup_gro.read_bytes()
        # Topology: section counts match and the [ bonds ] section is line-identical
        # (the physical bonds, including ester linkages, are the same). The
        # compare_topology_files cross_bonds *set* can differ by atom-ID
        # bookkeeping while the [ bonds ] section is identical, so we assert on
        # the real section identity + counts.
        parity = compare_topology_files(str(backup_top), str(result.topology_path))
        for section in ("atoms", "bonds", "angles", "dihedrals", "pairs"):
            assert parity[section], f"{section} section keys differ: {parity}"
        assert _section_lines(backup_top, "bonds") == _section_lines(
            result.topology_path, "bonds"
        ), "[ bonds ] section lines differ"
        assert len(_section_lines(backup_top, "angles")) == len(
            _section_lines(result.topology_path, "angles")
        )
        assert len(_section_lines(backup_top, "dihedrals")) == len(
            _section_lines(result.topology_path, "dihedrals")
        )
    finally:
        # Restore the archive to pristine from the tmp backup; remove build artifacts.
        shutil.copy2(backup_gro, ref_gro)
        shutil.copy2(backup_top, ref_top)
        for artifact in (
            "_hyb_conf.gro.1_",
            "_hyb_topol.top.1_",
            "at_hyb_topol.top",
            "cross_angles_hyb_topol.dat",
            "cross_bonds_hyb_topol.dat",
            "cross_dihedrals_hyb_topol.dat",
            "graph_before_cross_bonds.pck",
        ):
            (data_dir / artifact).unlink(missing_ok=True)
        # Restore tracked files the build regenerates so the archive stays pristine.
        repo = data_dir.parents[2]
        for name in ("missing_definitions.txt", "exclusion_hyb_topol.list"):
            subprocess.run(
                ["git", "checkout", "--", str((data_dir / name).relative_to(repo))],
                cwd=repo,
                check=False,
                capture_output=True,
            )
