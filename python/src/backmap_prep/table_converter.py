"""Convert GROMACS .xvg tabulated potentials to LAMMPS .table format."""

from __future__ import annotations

from typing import TYPE_CHECKING

from . import units

if TYPE_CHECKING:
    from pathlib import Path

    from .builder import System
    from .schema import Settings


def convert_tables(system: System, settings: Settings, out_dir: Path) -> list[Path]:
    """Convert all referenced .xvg tables to LAMMPS .table format."""
    converted: list[Path] = []

    for src_name, dst_name in system.table_files:
        src_path = out_dir / src_name
        dst_path = out_dir / dst_name

        if not src_path.exists():
            # Try relative to settings directory
            continue

        suffix = src_path.suffix.lower()
        if suffix == ".xvg":
            _convert_xvg(src_path, dst_path)
            converted.append(dst_path)
        elif suffix == ".table":
            # Already LAMMPS format, skip conversion
            if src_path != dst_path:
                import shutil

                shutil.copy2(src_path, dst_path)
            converted.append(dst_path)

    return converted


def _convert_xvg(src: Path, dst: Path, is_pair: bool = True) -> None:
    """Convert a GROMACS .xvg file to LAMMPS table format.

    GROMACS .xvg columns: r(nm), V(kJ/mol), F(kJ/(mol·nm))
    LAMMPS table columns: index, r(Å), energy(kcal/mol), force(kcal/(mol·Å))
    """
    r_vals: list[float] = []
    e_vals: list[float] = []
    f_vals: list[float] = []

    for line in src.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "@")):
            continue
        tokens = line.split()
        if len(tokens) < 3:
            continue
        r_nm = float(tokens[0])
        v_kj = float(tokens[1])
        f_kj = float(tokens[2])

        r_vals.append(units.distance(r_nm))
        e_vals.append(units.energy(v_kj))
        f_vals.append(units.force(f_kj))

    if not r_vals:
        raise ValueError(f"No data found in {src}")

    n = len(r_vals)
    keyword = "ENTRY"

    with open(dst, "w") as f:
        f.write(f"# Converted from {src.name} by backmap-prep\n")
        f.write("# GROMACS units (nm, kJ/mol) → LAMMPS real (Å, kcal/mol)\n\n")
        f.write(f"{keyword}\n")
        f.write(f"N {n}\n\n")

        for i in range(n):
            f.write(f"{i + 1} {r_vals[i]:.8f} {e_vals[i]:.8f} {f_vals[i]:.8f}\n")
