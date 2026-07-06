"""Public API for Phase 3 network (GROMACS hybrid) backmap preparation."""

from __future__ import annotations

import logging
import os
import random
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from backmap_prep.network.lammps_builder import build_system_from_hybrid
from backmap_prep.schema import resolve_data_dir, resolve_forcefield_dir, resolve_tables_dir

from .bakery.structures import BackmapperSettings2
from .v2_loader import has_native_network_config, settings_to_xml_root

if TYPE_CHECKING:
    from backmap_prep.builder import System
    from backmap_prep.schema import Settings

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class HybridBuildResult:
    """Paths and counts from a successful hybrid GROMACS build."""

    coordinates_path: Path
    topology_path: Path
    n_atoms: int
    missing_definitions_path: Path | None


@dataclass(frozen=True)
class NetworkLammpsBuildResult:
    """LAMMPS-ready network build output plus intermediate hybrid files."""

    system: System
    coordinates_path: Path
    topology_path: Path
    n_atoms: int
    missing_definitions_path: Path | None


def build_hybrid_gromacs(
    settings_source: Path | Settings,
    *,
    base_dir: Path | None = None,
    allow_no_bonds: bool = False,
    chain_rng_seed: int | None = None,
) -> HybridBuildResult:
    """Run the bakery network engine to produce hybrid GRO and topology files.

    Accepts either a bakery ``settings.xml`` path or a native v2 ``Settings`` object.
    All relative paths are resolved from ``base_dir`` (defaults to the XML parent or
    the caller-provided data directory for v2 YAML).
    """
    if chain_rng_seed is not None:
        random.seed(chain_rng_seed)

    if isinstance(settings_source, Path):
        settings_xml = settings_source.resolve()
        work_dir = (base_dir or settings_xml.parent).resolve()
        xml_root = None
        label = settings_xml.name
    else:
        settings = settings_source
        if not has_native_network_config(settings):
            raise ValueError("Settings object is not a complete native network configuration")
        work_dir = (base_dir or Path.cwd()).resolve()
        xml_root = settings_to_xml_root(settings)
        label = "settings.v2.yaml"

    previous_cwd = Path.cwd()
    try:
        os.chdir(work_dir)
        logger.info("Building hybrid system from %s (cwd=%s)", label, work_dir)
        if xml_root is not None:
            backmapper = BackmapperSettings2(
                xml_root,
                allow_no_bonds=allow_no_bonds,
            )
        else:
            backmapper = BackmapperSettings2(
                str(settings_xml),
                allow_no_bonds=allow_no_bonds,
            )
        backmapper.prepare_hybrid()
    finally:
        os.chdir(previous_cwd)

    coord_name = backmapper.hybrid_configuration["file"].file_name
    topol_name = backmapper.hyb_topology.file_name
    coordinates_path = work_dir / coord_name
    topology_path = work_dir / topol_name

    missing_path = work_dir / "missing_definitions.txt"
    if not missing_path.is_file():
        missing_path = None

    n_atoms = len(backmapper.hybrid_configuration["file"].atoms)
    return HybridBuildResult(
        coordinates_path=coordinates_path,
        topology_path=topology_path,
        n_atoms=n_atoms,
        missing_definitions_path=missing_path,
    )


def build_network_lammps(settings: Settings, settings_path: Path) -> NetworkLammpsBuildResult:
    """Build hybrid GRO/TOP for network systems and map them to a LAMMPS `System`."""
    work_dir = resolve_data_dir(settings_path, settings)
    hybrid = build_hybrid_gromacs(
        settings,
        base_dir=work_dir,
        allow_no_bonds=settings.prep.allow_no_bonds,
        chain_rng_seed=settings.prep.chain_rng_seed,
    )
    forcefield_dir = resolve_forcefield_dir(settings_path, settings)
    forcefield_dirs = [forcefield_dir] if forcefield_dir is not None else []
    system = build_system_from_hybrid(
        settings=settings,
        base_dir=work_dir,
        gro_path=hybrid.coordinates_path,
        top_path=hybrid.topology_path,
        table_search_dirs=[
            d for d in [resolve_tables_dir(settings_path, settings)] if d is not None
        ],
        forcefield_dirs=forcefield_dirs,
    )
    return NetworkLammpsBuildResult(
        system=system,
        coordinates_path=hybrid.coordinates_path,
        topology_path=hybrid.topology_path,
        n_atoms=hybrid.n_atoms,
        missing_definitions_path=hybrid.missing_definitions_path,
    )
