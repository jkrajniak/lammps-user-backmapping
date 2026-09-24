"""Public API for Phase 3 network (GROMACS hybrid) backmap preparation."""

from __future__ import annotations

import logging
import os
import random
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from backmap_prep.network.lammps_builder import build_system_from_hybrid
from backmap_prep.parsers.top_parser import parse_top
from backmap_prep.schema import (
    resolve_bakery_xml,
    resolve_data_dir,
    resolve_forcefield_dir,
    resolve_tables_dir,
)

from .bakery.structures import BackmapperSettings2
from .lammps_sources import materialize_lammps_sources
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


def check_molecule_names(settings: Settings, work_dir: Path) -> None:
    """Fail early when a CG residue has no molecule definition of the same name.

    The hybrid builder matches AT fragments to CG residues by name; a mismatch
    otherwise surfaces deep inside bakery as "could not find correct fragment".
    """
    if settings.cg_system is None or settings.cg_system.topology is None:
        return
    top_path = work_dir / settings.cg_system.topology
    if not top_path.is_file():
        return
    cg_top = parse_top(top_path, include_dirs=[top_path.parent])
    residues = {
        atom.resname
        for name, _ in cg_top.molecules
        if name in cg_top.molecule_types
        for atom in cg_top.molecule_types[name].atoms
    }
    defined = {mol.name for mol in settings.molecules}
    missing = sorted(residues - defined)
    if missing:
        raise ValueError(
            f"CG residue(s) {missing} in {settings.cg_system.topology} have no molecule "
            f"definition: molecules[].name must equal the CG residue name "
            f"(defined: {sorted(defined)})"
        )


def _absolute(path: str, data_dir: Path) -> str:
    candidate = Path(path)
    if candidate.is_absolute() or not (data_dir / candidate).exists():
        return path
    return str((data_dir / candidate).resolve())


def absolutize_inputs(settings: Settings, data_dir: Path) -> Settings:
    """Copy of settings with every input file path made absolute against ``data_dir``.

    The hybrid build then runs in the output directory and only reads from
    ``data_dir`` (which may be a published data archive).
    """
    new = settings.model_copy(deep=True)
    for mol in new.molecules:
        for attr in ("coordinates", "topology"):
            value = getattr(mol.source, attr)
            if isinstance(value, str):
                setattr(mol.source, attr, _absolute(value, data_dir))
            elif value:
                for entry in value:
                    entry.file = _absolute(entry.file, data_dir)
    if new.cg_system is not None:
        if new.cg_system.coordinates:
            new.cg_system.coordinates = _absolute(new.cg_system.coordinates, data_dir)
        if new.cg_system.topology:
            new.cg_system.topology = _absolute(new.cg_system.topology, data_dir)
    new.hybrid.topology_includes = [
        _absolute(inc, data_dir) for inc in new.hybrid.topology_includes
    ]
    return new


def build_hybrid_gromacs(
    settings_source: Path | Settings,
    *,
    base_dir: Path | None = None,
    output_dir: Path | None = None,
    allow_no_bonds: bool = False,
    chain_rng_seed: int | None = None,
) -> HybridBuildResult:
    """Run the bakery network engine to produce hybrid GRO and topology files.

    Accepts either a bakery ``settings.xml`` path or a native v2 ``Settings`` object.
    All relative paths are resolved from ``base_dir`` (defaults to the XML parent or
    the caller-provided data directory for v2 YAML). For v2 settings, outputs go
    to ``output_dir`` (default ``base_dir``); ``base_dir`` is only read. A bakery
    XML runs in its own directory, which must be writable.
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
        data_dir = (base_dir or Path.cwd()).resolve()
        work_dir = (output_dir or data_dir).resolve()
        work_dir.mkdir(parents=True, exist_ok=True)
        xml_root = settings_to_xml_root(absolutize_inputs(settings, data_dir))
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
    output_dir = settings_path.parent.resolve()
    if settings.prep.bakery_xml:
        # Bakery settings.xml passthrough (e.g. melamine_network): the v2
        # Settings object intentionally carries no native molecules/cg_system/
        # hybrid block in this mode, so has_native_network_config() is False
        # and build_hybrid_gromacs() must be driven from the XML path instead
        # -- mirrors cli.py's _cmd_build_hybrid bakery_xml branch.
        xml_path = resolve_bakery_xml(settings_path, settings)
        hybrid = build_hybrid_gromacs(
            xml_path,
            base_dir=xml_path.parent,
            allow_no_bonds=settings.prep.allow_no_bonds,
            chain_rng_seed=settings.prep.chain_rng_seed,
        )
    else:
        settings = materialize_lammps_sources(settings, work_dir, output_dir)
        check_molecule_names(settings, work_dir)
        hybrid = build_hybrid_gromacs(
            settings,
            base_dir=work_dir,
            output_dir=output_dir,
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
            d for d in [output_dir, resolve_tables_dir(settings_path, settings)] if d is not None
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
