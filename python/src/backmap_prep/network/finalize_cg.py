"""Post-process equilibrated CG LAMMPS frames into consistent PBC representation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from backmap_prep import units
from backmap_prep.builder import LammpsAtom
from backmap_prep.network.lammps_builder import build_system_from_cg
from backmap_prep.network.pbc import (
    max_bond_length,
    prepare_network_coordinates,
    validate_bond_geometry,
)
from backmap_prep.parsers import parse_gro, parse_lammps_data

if TYPE_CHECKING:
    from backmap_prep.builder import System
    from backmap_prep.schema import Settings


@dataclass(frozen=True)
class FinalizeCgResult:
    system: System
    max_bond_ang: float
    nonzero_image_flags: int


def _unwrapped_angstrom(
    atom: LammpsAtom, box: tuple[float, float, float]
) -> tuple[float, float, float]:
    bx, by, bz = box
    return (
        atom.x + atom.ix * bx if bx > 0 else atom.x,
        atom.y + atom.iy * by if by > 0 else atom.y,
        atom.z + atom.iz * bz if bz > 0 else atom.z,
    )


def finalize_cg_from_equil(
    settings: Settings,
    settings_path: Path,
    equil_path: Path,
    table_search_dirs: list[Path] | None = None,
) -> FinalizeCgResult:
    """Apply bond-graph PBC fix to an equilibrated CG LAMMPS data frame."""
    system = build_system_from_cg(settings, settings_path, table_search_dirs=table_search_dirs)
    frame = parse_lammps_data(equil_path)
    if len(frame.atoms) != len(system.atoms):
        raise ValueError(
            f"Equil frame has {len(frame.atoms)} atoms but CG topology expects {len(system.atoms)}"
        )

    if any(value > 0 for value in frame.box):
        system.box = frame.box

    bx, by, bz = system.box
    frame_by_id = {atom.atom_id: atom for atom in frame.atoms}
    for atom in system.atoms:
        frame_atom = frame_by_id.get(atom.atom_id)
        if frame_atom is None:
            raise ValueError(f"Equil frame missing atom id {atom.atom_id}")
        atom.x = frame_atom.x + frame_atom.ix * bx if bx > 0 else frame_atom.x
        atom.y = frame_atom.y + frame_atom.iy * by if by > 0 else frame_atom.y
        atom.z = frame_atom.z + frame_atom.iz * bz if bz > 0 else frame_atom.z
        atom.ix = 0
        atom.iy = 0
        atom.iz = 0

    prepare_network_coordinates(system)
    cg_cut = units.distance(settings.simulation.cg_cutoff)
    lj_cut = units.distance(settings.simulation.lj_cutoff)
    validate_bond_geometry(system, max(lj_cut, cg_cut))

    system.write_image_flags = True
    max_bond = max_bond_length(system.atoms, system.bonds, system.box)
    nonzero_flags = sum(1 for atom in system.atoms if atom.ix or atom.iy or atom.iz)
    return FinalizeCgResult(
        system=system,
        max_bond_ang=max_bond,
        nonzero_image_flags=nonzero_flags,
    )


def write_cg_gro(
    system: System,
    template_gro: Path,
    output_path: Path,
    *,
    unwrapped: bool = True,
) -> None:
    """Write CG coordinates as GROMACS .gro (nm), preserving residue/atom names."""
    template = parse_gro(template_gro)
    template_by_id = {atom.index: atom for atom in template.atoms}
    bx, by, bz = system.box

    lines = [template.title, str(len(system.atoms))]
    for atom in sorted(system.atoms, key=lambda item: item.atom_id):
        meta = template_by_id.get(atom.atom_id)
        if meta is None:
            raise ValueError(f"Template GRO missing atom id {atom.atom_id}")
        if unwrapped:
            x_ang, y_ang, z_ang = _unwrapped_angstrom(atom, system.box)
        else:
            x_ang, y_ang, z_ang = atom.x, atom.y, atom.z
        x_nm = x_ang / units.NM_TO_ANGSTROM
        y_nm = y_ang / units.NM_TO_ANGSTROM
        z_nm = z_ang / units.NM_TO_ANGSTROM
        lines.append(
            f"{meta.resid:5d}{meta.resname:>5}{meta.name:>5}{atom.atom_id:5d}"
            f"{x_nm:8.3f}{y_nm:8.3f}{z_nm:8.3f}"
        )

    box_nm = (
        bx / units.NM_TO_ANGSTROM,
        by / units.NM_TO_ANGSTROM,
        bz / units.NM_TO_ANGSTROM,
    )
    lines.append(f"{box_nm[0]:10.5f}{box_nm[1]:10.5f}{box_nm[2]:10.5f}")
    output_path.write_text("\n".join(lines) + "\n")
