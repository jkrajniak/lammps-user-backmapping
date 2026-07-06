"""Rebuild network hybrid LAMMPS systems from equilibrated CG coordinates."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from backmap_prep import units
from backmap_prep.network.api import NetworkLammpsBuildResult, build_hybrid_gromacs
from backmap_prep.network.lammps_builder import build_system_from_hybrid
from backmap_prep.network.pbc import (
    assign_network_image_flags,
    max_bond_length,
    validate_bond_geometry,
)
from backmap_prep.parsers import parse_gro, parse_lammps_data
from backmap_prep.schema import resolve_data_dir, resolve_forcefield_dir, resolve_tables_dir

if TYPE_CHECKING:
    from backmap_prep.builder import System
    from backmap_prep.schema import Settings

BeadKey = tuple[int, str, str]


def _bead_key(resid: int, resname: str, name: str) -> BeadKey:
    return (resid, resname.strip(), name.strip())


def load_cg_position_map(
    cg_frame: Path,
    *,
    template_gro: Path | None = None,
) -> tuple[dict[BeadKey, tuple[float, float, float]], tuple[float, float, float]]:
    """Load unwrapped CG bead positions (Å) keyed by (resid, resname, bead name)."""
    path = cg_frame.resolve()
    suffix = path.suffix.lower()

    if suffix == ".gro":
        gro = parse_gro(path)
        position_map = {
            _bead_key(atom.resid, atom.resname, atom.name): (
                units.distance(atom.x),
                units.distance(atom.y),
                units.distance(atom.z),
            )
            for atom in gro.atoms
        }
        box = (
            units.distance(gro.box[0]),
            units.distance(gro.box[1]),
            units.distance(gro.box[2]),
        )
        return position_map, box

    if suffix == ".data":
        if template_gro is None:
            raise ValueError("template GRO is required when CG frame is LAMMPS .data")
        frame = parse_lammps_data(path)
        template = parse_gro(template_gro.resolve())
        template_by_id = {atom.index: atom for atom in template.atoms}
        bx, by, bz = frame.box
        position_map: dict[BeadKey, tuple[float, float, float]] = {}
        for atom in frame.atoms:
            meta = template_by_id.get(atom.atom_id)
            if meta is None:
                raise ValueError(f"Template GRO missing atom id {atom.atom_id}")
            position_map[_bead_key(meta.resid, meta.resname, meta.name)] = (
                atom.x + atom.ix * bx if bx > 0 else atom.x,
                atom.y + atom.iy * by if by > 0 else atom.y,
                atom.z + atom.iz * bz if bz > 0 else atom.z,
            )
        return position_map, frame.box

    raise ValueError(f"Unsupported CG frame format: {path} (use .data or .gro)")


def _unwrap_all_atoms(system: System) -> None:
    """Convert folded coords + image flags to unwrapped Cartesian positions."""
    bx, by, bz = system.box
    for atom in system.atoms:
        atom.x += atom.ix * bx if bx > 0 else 0.0
        atom.y += atom.iy * by if by > 0 else 0.0
        atom.z += atom.iz * bz if bz > 0 else 0.0
        atom.ix = 0
        atom.iy = 0
        atom.iz = 0


def reposition_hybrid_from_cg(
    system: System,
    cg_positions: dict[BeadKey, tuple[float, float, float]],
    hybrid_gro: Path,
) -> None:
    """Rigidly translate each backmap group so its CG bead matches equilibrated coords."""
    _unwrap_all_atoms(system)
    gro = parse_gro(hybrid_gro.resolve())
    gro_by_index = {atom.index: atom for atom in gro.atoms}

    cg_atoms = [atom for atom in system.atoms if atom.is_cg]
    if len(cg_atoms) != len(cg_positions):
        raise ValueError(
            f"CG frame has {len(cg_positions)} beads but hybrid expects {len(cg_atoms)}"
        )

    atoms_by_mol: dict[int, list] = {}
    for atom in system.atoms:
        atoms_by_mol.setdefault(atom.mol_id, []).append(atom)

    for cg_atom in cg_atoms:
        meta = gro_by_index.get(cg_atom.atom_id)
        if meta is None:
            raise ValueError(f"Hybrid GRO missing atom id {cg_atom.atom_id}")
        key = _bead_key(meta.resid, meta.resname, meta.name)
        if key not in cg_positions:
            raise ValueError(f"CG frame missing bead {key!r}")
        new_x, new_y, new_z = cg_positions[key]
        delta = (new_x - cg_atom.x, new_y - cg_atom.y, new_z - cg_atom.z)
        if delta == (0.0, 0.0, 0.0):
            continue
        for atom in atoms_by_mol[cg_atom.mol_id]:
            atom.x += delta[0]
            atom.y += delta[1]
            atom.z += delta[2]


def rebuild_network_lammps(
    settings: Settings,
    settings_path: Path,
    cg_frame: Path,
) -> NetworkLammpsBuildResult:
    """Build hybrid topology and place AT fragments using equilibrated CG coordinates."""
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
            directory
            for directory in [resolve_tables_dir(settings_path, settings)]
            if directory is not None
        ],
        forcefield_dirs=forcefield_dirs,
    )

    template_gro: Path | None = None
    if settings.cg_system is not None:
        template_gro = (work_dir / settings.cg_system.coordinates).resolve()

    cg_positions, box = load_cg_position_map(cg_frame, template_gro=template_gro)
    if any(value > 0 for value in box):
        system.box = box
    reposition_hybrid_from_cg(system, cg_positions, hybrid.coordinates_path)
    assign_network_image_flags(system)
    lj_cut = units.distance(settings.simulation.lj_cutoff)
    cg_cut = units.distance(settings.simulation.cg_cutoff)
    validate_bond_geometry(system, max(lj_cut, cg_cut))
    system.write_image_flags = True

    max_bond = max_bond_length(system.atoms, system.bonds, system.box)
    print(
        f"Rebuild PBC ok: max bonded distance {max_bond:.2f} Å, "
        f"{sum(1 for atom in system.atoms if atom.ix or atom.iy or atom.iz)} atoms with image flags"
    )

    return NetworkLammpsBuildResult(
        system=system,
        coordinates_path=hybrid.coordinates_path,
        topology_path=hybrid.topology_path,
        n_atoms=hybrid.n_atoms,
        missing_definitions_path=hybrid.missing_definitions_path,
    )
