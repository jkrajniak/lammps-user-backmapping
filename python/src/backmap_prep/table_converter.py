"""Convert GROMACS .xvg tabulated potentials to LAMMPS .table format."""

from __future__ import annotations

import logging
import math
from typing import TYPE_CHECKING

from . import units

if TYPE_CHECKING:
    from pathlib import Path

    from .builder import System
    from .schema import Settings

logger = logging.getLogger(__name__)

# Soft-core wall for r < table r_min (LAMMPS aborts below inner cutoff; E++ clamps).
# LAMMPS pair_table requires rlo > 0 ("Invalid pair table lower boundary" if rlo <= 0).
_PAIR_CORE_POWER = 12.0
_PAIR_R_FLOOR = 1.0e-4  # Å — smallest allowed table lower bound
# Ceiling on the extrapolated wall energy. An uncapped (r_min/r)^12 power law,
# extended all the way down to _PAIR_R_FLOOR, blows up for tables whose usable
# data doesn't already start close to the floor (e.g. r_min=0.02 Å gives
# (0.02/0.0001)**12 ~ 4e27, producing a ~1e31 kcal/mol table entry) -- a value
# LAMMPS's pair_style table spline does not handle gracefully near a poisoned
# entry, causing simulations to blow up if minimization or a lambda ramp ever
# samples that close-contact region. 1e6 kcal/mol reuses this module's own
# existing "sufficiently repulsive" convention (the e_ref fallback below) --
# still many orders of magnitude past any physically meaningful MD energy
# scale, but small enough to stay numerically well-behaved.
_PAIR_WALL_ENERGY_CAP = 1.0e6  # kcal/mol
# Maximum allowed r[i+1]/r[i] ratio for the extension points prepended below
# r_min. A table whose real data is uniformly spaced (e.g. delta=0.02 Å) but
# whose r_min sits far above r_floor (e.g. 0.02 Å vs. a 1e-4 Å floor, a 200x
# gap) used to get exactly one extension point at r_floor, stepped by the
# table's own uniform delta -- i.e. a single point 200x further from r_min
# than any other adjacent pair in the table. Even with the energy cap above,
# that one badly-conditioned knot is enough to corrupt LAMMPS's cubic-spline
# fit for pair_style table near the low-r region (confirmed by reconstructing
# the spline directly: the corruption is confined to roughly r_floor..0.2 Å,
# but a real run's energy minimizer resolving an initial bad contact can
# land exactly there -- examples/pe-lammps/README.md has the full trail).
# Filling the gap with several log-spaced points instead, each within this
# ratio of its neighbor, keeps the extension smooth enough for LAMMPS's
# spline fit regardless of how large r_min/r_floor is.
_PAIR_EXT_MAX_RATIO = 1.2


def extend_pair_table_to_zero(
    r_vals: list[float],
    e_vals: list[float],
    f_vals: list[float],
    *,
    power: float = _PAIR_CORE_POWER,
    r_floor: float = _PAIR_R_FLOOR,
) -> tuple[list[float], list[float], list[float]]:
    """Prepend a steep repulsive wall down to a tiny positive r (not literally 0).

    LAMMPS ``pair_style table`` errors if any pair distance is below the table
    inner cutoff, and also rejects ``rlo <= 0``. ESPResSo++ clamps to the first
    bin instead. Extending from ``r_floor`` (> 0) to the original ``r_min`` with a
    continuous wall matching E,F at ``r_min`` avoids the LAMMPS abort while
    keeping the original IBI table for ``r >= r_min``.
    """
    if not r_vals:
        raise ValueError("empty pair table")
    if r_floor <= 0.0:
        raise ValueError("r_floor must be > 0 (LAMMPS rejects rlo <= 0)")

    r_min = r_vals[0]
    if r_min <= r_floor + 1e-15:
        # Already covers the floor (or illegal <=0 — bump first point).
        if r_min > 0.0:
            return list(r_vals), list(e_vals), list(f_vals)
        fixed_r = [r_floor, *r_vals[1:]]
        return fixed_r, list(e_vals), list(f_vals)

    e0 = e_vals[0]
    # Prefer a positive (repulsive) core; if e0 is tiny/negative, use a floor.
    e_ref = e0 if e0 > 1.0 else 1.0e6

    # Log-spaced extension points from r_floor up to (but excluding) r_min, so
    # that no two consecutive r values -- including the seam between the last
    # extension point and r_min itself -- differ by more than
    # _PAIR_EXT_MAX_RATIO. See the constant's docstring for why a single
    # far-spaced knot is a problem even with the energy cap.
    ratio_needed = r_min / r_floor
    n_extra = max(1, math.ceil(math.log(ratio_needed) / math.log(_PAIR_EXT_MAX_RATIO)))
    if n_extra > 100_000:
        raise ValueError("pair table core extension produced too many points")
    extra_r = [r_floor * (ratio_needed ** (i / n_extra)) for i in range(n_extra)]

    # Below r_cap, switch from the power-law wall to a constant-force (linear
    # energy) wall, continuous in both E and F at r_cap, so the extrapolated
    # energy saturates at _PAIR_WALL_ENERGY_CAP instead of diverging as
    # r -> r_floor. r_cap is where the pure power law would cross the cap;
    # clamped to [r_floor, r_min] so it's a no-op when the power law never
    # reaches the cap before r_floor (r_min already close to r_floor) and
    # never extends the linear segment above the table's own original data.
    r_cap = r_min * (e_ref / _PAIR_WALL_ENERGY_CAP) ** (1.0 / power)
    r_cap = min(max(r_cap, r_floor), r_min)
    e_at_cap = e_ref * (r_min / r_cap) ** power
    f_at_cap = power * e_ref * (r_min**power) / (r_cap ** (power + 1.0))

    def wall_e(r: float) -> float:
        if r >= r_cap:
            return float(e_ref * (r_min / r) ** power)
        return float(e_at_cap + f_at_cap * (r_cap - r))

    def wall_f(r: float) -> float:
        # F = -dE/dr for E = e_ref * (r_min/r)^n  → positive (repulsive)
        if r >= r_cap:
            return float(power * e_ref * (r_min**power) / (r ** (power + 1.0)))
        return float(f_at_cap)

    new_r = extra_r + list(r_vals)
    new_e = [wall_e(r) for r in extra_r] + list(e_vals)
    new_f = [wall_f(r) for r in extra_r] + list(f_vals)
    return new_r, new_e, new_f


def parse_lammps_pair_table(path: Path) -> tuple[list[float], list[float], list[float]]:
    """Parse a LAMMPS pair ``.table`` file into r, E, F columns."""
    r_vals: list[float] = []
    e_vals: list[float] = []
    f_vals: list[float] = []
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith(("#", "ENTRY", "N ")):
            continue
        if not s[0].isdigit():
            continue
        p = s.split()
        if len(p) < 4:
            continue
        r_vals.append(float(p[1]))
        e_vals.append(float(p[2]))
        f_vals.append(float(p[3]))
    if not r_vals:
        raise ValueError(f"No data found in {path}")
    return r_vals, e_vals, f_vals


def write_lammps_pair_table(
    dst: Path,
    r_vals: list[float],
    e_vals: list[float],
    f_vals: list[float],
    *,
    source: str = "",
) -> None:
    """Write a LAMMPS pair ``.table`` file."""
    _write_table_file(dst, source or dst.name, r_vals, e_vals, f_vals, x_axis="distance")


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
                _convert_xvg(src_path, dst_path, skip_zero=skip_zero, extend_core=(kind == "pair"))
            converted.append(dst_path)
        elif suffix == ".table":
            # A .table file is already native LAMMPS format -- pass it through
            # unchanged for every kind, including pair. Pair tables used to be
            # re-parsed and run through extend_pair_table_to_zero() even when
            # already in .table format, on the theory that this made a missing
            # low-r wall safe to add. In practice, inserting a single
            # extension point far below the table's existing r spacing (e.g.
            # r_floor=1e-4 prepended to a table starting at r_min=0.02, a 200x
            # gap versus its uniform 0.02 internal spacing) corrupts LAMMPS's
            # cubic-spline fit for pair_style table badly enough to
            # destabilize a real simulation within a few steps, independent
            # of how large the extrapolated wall energy is (found running
            # examples/pe-lammps/ on a real LAMMPS build: capping the wall
            # energy at a sane ceiling did not fix the crash; using the
            # untouched source .table file did). If a .table file's r_min is
            # genuinely too large for the pair distances a run will sample,
            # LAMMPS's own "Pair distance < table inner cutoff" error is a
            # clearer, more actionable failure than a silently corrupted fit.
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


def _convert_xvg(
    src: Path,
    dst: Path,
    skip_zero: bool = False,
    *,
    extend_core: bool = False,
) -> None:
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

    if extend_core and r_vals[0] > _PAIR_R_FLOOR:
        logger.warning(
            "Pair table %s starts at r_min=%.4f Å. Extending down to r=%.1e Å "
            "with a steep repulsive wall (LAMMPS requires coverage and rlo>0; "
            "E++ clamps).",
            src.name,
            r_vals[0],
            _PAIR_R_FLOOR,
        )
        r_vals, e_vals, f_vals = extend_pair_table_to_zero(r_vals, e_vals, f_vals)

    _write_table_file(dst, src.name, r_vals, e_vals, f_vals, x_axis="distance")


def _convert_angle_xvg(src: Path, dst: Path) -> None:
    """Convert GROMACS angle table .xvg to LAMMPS angle table format.

    Column 0: angle in degrees (0-180)
    Column 1: V(kJ/mol)
    Column 2: -dV/dtheta in kJ/(mol*deg) -- already a per-degree derivative,
        matching this file's own degree-based column 0 (confirmed by
        checking that dV/dtheta computed directly from columns 0-1 matches
        column 2 only when theta is treated as degrees, not radians; see
        decisions/2026-07-23-fix-angular-dihedral-table-force-units.md).
        Needs only the plain energy-unit conversion, no extra angle-unit
        rescaling.
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
        f_kj_per_deg = float(tokens[2])

        theta_vals.append(theta_deg)
        e_vals.append(units.energy(v_kj))
        f_vals.append(units.energy(f_kj_per_deg))

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
    Column 2: -dV/dphi in kJ/(mol*deg) -- already a per-degree derivative,
        matching this file's own degree-based column 0. Needs only the
        plain energy-unit conversion, no extra angle-unit rescaling (see
        decisions/2026-07-23-fix-angular-dihedral-table-force-units.md).
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
        f_kj_per_deg = float(tokens[2])

        phi_vals.append(phi_deg)
        e_vals.append(units.energy(v_kj))
        f_vals.append(units.energy(f_kj_per_deg))

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
        if x_axis == "distance":
            f.write(
                "# GROMACS units → LAMMPS real (x=Å, kcal/mol); "
                "pair core extended to r_floor>0 if needed\n\n"
            )
        else:
            f.write(f"# GROMACS units → LAMMPS real (x={unit_note}, kcal/mol)\n\n")
        f.write(f"{keyword}\n")
        f.write(f"N {n}\n\n")

        for i in range(n):
            f.write(f"{i + 1} {x_vals[i]:.8f} {e_vals[i]:.8f} {f_vals[i]:.8f}\n")
