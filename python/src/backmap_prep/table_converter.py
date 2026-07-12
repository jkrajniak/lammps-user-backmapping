"""Convert GROMACS .xvg tabulated potentials to LAMMPS .table format."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

from . import units

if TYPE_CHECKING:
    from pathlib import Path

    from .builder import System
    from .schema import Settings


def convert_tables(
    system: System,
    settings: Settings,
    out_dir: Path,
    extra_dirs: list[Path] | None = None,
) -> list[Path]:
    """Convert all referenced .xvg tables to LAMMPS .table format."""
    converted: list[Path] = []
    search_dirs = [out_dir, *(extra_dirs or [])]

    all_tables = (
        [(src, dst, "bond", False) for src, dst in system.table_files]
        + [(src, dst, "angle", False) for src, dst in system.angle_table_files]
        + [(src, dst, "dihedral", False) for src, dst in system.dihedral_table_files]
        + [(src, dst, "pair", True) for src, dst in system.pair_table_files]
    )

    for src_name, dst_name, kind, skip_zero in all_tables:
        src_path: Path | None = None
        for directory in search_dirs:
            candidate = directory / src_name
            if candidate.is_file():
                src_path = candidate
                break
        if src_path is None:
            continue

        dst_path = out_dir / dst_name

        suffix = src_path.suffix.lower()
        if suffix == ".xvg":
            if kind == "angle":
                _convert_angle_xvg(src_path, dst_path)
            elif kind == "dihedral":
                _convert_dihedral_xvg(src_path, dst_path)
            else:
                _convert_xvg(src_path, dst_path, skip_zero=skip_zero)
            converted.append(dst_path)
        elif suffix == ".table":
            if src_path != dst_path:
                import shutil

                shutil.copy2(src_path, dst_path)
            converted.append(dst_path)

    return converted


def _numerical_gradient(x_vals: list[float], y_vals: list[float]) -> list[float]:
    """Central-difference dy/dx (forward/backward at the endpoints)."""
    n = len(x_vals)
    grad = [0.0] * n
    if n < 2:
        return grad
    grad[0] = (y_vals[1] - y_vals[0]) / (x_vals[1] - x_vals[0])
    grad[-1] = (y_vals[-1] - y_vals[-2]) / (x_vals[-1] - x_vals[-2])
    for i in range(1, n - 1):
        grad[i] = (y_vals[i + 1] - y_vals[i - 1]) / (x_vals[i + 1] - x_vals[i - 1])
    return grad


def _is_force_column_degenerate(f_vals: list[float], e_vals: list[float]) -> bool:
    """True when the force column is (near-)zero while energy clearly varies.

    Some legacy ESPResSo++ -> GROMACS .xvg exports (dated 2017-04, found in
    the PET/Dacron table set) never populated the derivative/force column,
    leaving it all zeros while the energy column is a real, varying
    potential. GROMACS's own mdrun re-splines forces from the energy column
    internally, masking the bug; a LAMMPS `pair_style table` reads the force
    column literally, so a degenerate column silently produces zero
    nonbonded force everywhere.
    """
    max_f = max((abs(v) for v in f_vals), default=0.0)
    max_e = max((abs(v) for v in e_vals), default=0.0)
    if max_e < 1.0e-8:
        return False  # nothing to differentiate against; leave as-is
    return max_f < 1.0e-6 * max_e


def _convert_xvg(src: Path, dst: Path, skip_zero: bool = False) -> None:
    """Convert a GROMACS .xvg file to LAMMPS table format.

    Handles two formats:
      3-column: r(nm), V(kJ/mol), F(kJ/(mol·nm))
      7-column (energygrp-table): r, f, -f', g, -g', h(V), -h'(F)

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

        if len(tokens) >= 7:
            v_kj = float(tokens[5])
            f_kj = float(tokens[6])
        else:
            v_kj = float(tokens[1])
            f_kj = float(tokens[2])

        r_ang = units.distance(r_nm)
        if skip_zero and r_ang <= 0.0:
            continue

        r_vals.append(r_ang)
        e_vals.append(units.energy(v_kj))
        f_vals.append(units.force(f_kj))

    if not r_vals:
        raise ValueError(f"No data found in {src}")

    if _is_force_column_degenerate(f_vals, e_vals):
        # F(r) = -dV/dr; energies/positions are already in LAMMPS units.
        f_vals = [-g for g in _numerical_gradient(r_vals, e_vals)]

    _write_table_file(dst, src.name, r_vals, e_vals, f_vals, x_axis="distance")


def _convert_angle_xvg(src: Path, dst: Path) -> None:
    """Convert GROMACS angle table .xvg to LAMMPS angle table format.

    Column 0: angle in degrees (0–180)
    Column 1: V(kJ/mol)
    Column 2: -dV/dθ in kJ/(mol·rad)
    """
    theta_vals: list[float] = []
    e_vals: list[float] = []
    f_vals: list[float] = []

    for line in src.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "@")):
            continue
        tokens = line.split()
        if len(tokens) < 3:
            continue
        theta_deg = float(tokens[0])
        v_kj = float(tokens[1])
        f_kj_per_rad = float(tokens[2])

        theta_vals.append(theta_deg)
        e_vals.append(units.energy(v_kj))
        f_vals.append(units.angular_force(f_kj_per_rad))

    if not theta_vals:
        raise ValueError(f"No data found in {src}")

    if _is_force_column_degenerate(f_vals, e_vals):
        theta_rad = [math.radians(t) for t in theta_vals]
        f_vals = [-g for g in _numerical_gradient(theta_rad, e_vals)]

    _write_table_file(dst, src.name, theta_vals, e_vals, f_vals, x_axis="angle")


def _convert_dihedral_xvg(src: Path, dst: Path) -> None:
    """Convert GROMACS dihedral table .xvg to LAMMPS dihedral table format.

    Column 0: dihedral angle in degrees (-180..180)
    Column 1: V(kJ/mol)
    Column 2: -dV/dφ in kJ/(mol·rad)
    """
    phi_vals: list[float] = []
    e_vals: list[float] = []
    f_vals: list[float] = []

    for line in src.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "@")):
            continue
        tokens = line.split()
        if len(tokens) < 3:
            continue
        phi_deg = float(tokens[0])
        v_kj = float(tokens[1])
        f_kj_per_rad = float(tokens[2])

        phi_vals.append(phi_deg)
        e_vals.append(units.energy(v_kj))
        f_vals.append(units.angular_force(f_kj_per_rad))

    if not phi_vals:
        raise ValueError(f"No data found in {src}")

    if _is_force_column_degenerate(f_vals, e_vals):
        phi_rad = [math.radians(p) for p in phi_vals]
        f_vals = [-g for g in _numerical_gradient(phi_rad, e_vals)]

    _write_table_file(dst, src.name, phi_vals, e_vals, f_vals, x_axis="dihedral")


def _write_table_file(
    dst: Path,
    src_name: str,
    x_vals: list[float],
    e_vals: list[float],
    f_vals: list[float],
    *,
    x_axis: str,
) -> None:
    n = len(x_vals)
    keyword = "ENTRY"
    unit_note = "degrees" if x_axis in ("angle", "dihedral") else "Å"

    with open(dst, "w") as f:
        f.write(f"# Converted from {src_name} by backmap-prep\n")
        f.write(f"# GROMACS units → LAMMPS real (x={unit_note}, kcal/mol)\n\n")
        f.write(f"{keyword}\n")
        f.write(f"N {n}\n\n")

        for i in range(n):
            f.write(f"{i + 1} {x_vals[i]:.8f} {e_vals[i]:.8f} {f_vals[i]:.8f}\n")
