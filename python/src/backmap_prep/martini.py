"""MARTINI CG force field terms written as GROMACS-format tables.

The CG side of a backmapping run evaluates CG-CG pairs and CG angles through
tables. For a MARTINI model those tables are generated here from the CG
topology, reproducing what GROMACS computes with the standard MARTINI settings:

- LJ from ``[ nonbond_params ]`` (else ``[ atomtypes ]`` and the combination
  rule), cut at ``cutoff`` with the potential shifted to zero there
  (``vdw-modifier = Potential-shift``);
- reaction-field Coulomb (``coulombtype = reaction-field``) with ``epsilon_r``
  and ``epsilon_rf`` (0 = infinity), which is zero at the cut-off by
  construction;
- G96 angles (func 2), E = 1/2 k (cos theta - cos theta0)^2.

Tables are written as 3-column ``.xvg`` (x, V, -dV/dx) in GROMACS units and go
through the same converter as every other table.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

    from .parsers.top_parser import Topology
    from .schema import CGNonbonded

__all__ = ["cg_type_charges", "g96_angle_xvg", "lj_sigma_epsilon", "pair_xvg"]

# GROMACS electric conversion factor 1/(4 pi eps0), kJ mol^-1 nm e^-2
ONE_4PI_EPS0 = 138.935458


def lj_sigma_epsilon(top: Topology, type_i: str, type_j: str) -> tuple[float, float]:
    """(sigma nm, epsilon kJ/mol) of a CG type pair as GROMACS resolves it."""
    pair = top.nonbond_params.get((type_i, type_j))
    if pair is not None:
        a, b = pair
    else:
        ti, tj = top.atom_types.get(type_i), top.atom_types.get(type_j)
        if ti is None or tj is None:
            raise ValueError(f"no LJ parameters for CG types {type_i}-{type_j}")
        if top.combination_rule == 2:
            a, b = 0.5 * (ti.sigma + tj.sigma), math.sqrt(ti.epsilon * tj.epsilon)
        else:
            a, b = math.sqrt(ti.sigma * tj.sigma), math.sqrt(ti.epsilon * tj.epsilon)
    if top.combination_rule == 1:  # a = C6, b = C12
        if a <= 0.0 or b <= 0.0:
            return 0.0, 0.0
        sigma = (b / a) ** (1.0 / 6.0)
        return sigma, a * a / (4.0 * b)
    return a, b


def cg_type_charges(top: Topology) -> dict[str, float]:
    """One charge per CG atom type (MARTINI charges are a property of the bead type)."""
    charges: dict[str, float] = {}
    for mol in top.molecule_types.values():
        for atom in mol.atoms:
            known = charges.setdefault(atom.type, atom.charge)
            if abs(known - atom.charge) > 1e-9:
                raise ValueError(
                    f"CG type {atom.type} carries charges {known} and {atom.charge}; "
                    "per-type-pair tables need one charge per type"
                )
    return charges


def pair_xvg(sigma: float, epsilon: float, q_i: float, q_j: float, nb: CGNonbonded) -> str:
    """Pair table r, V(r), -dV/dr (nm, kJ/mol) up to the cut-off."""
    rc = nb.cutoff
    if nb.epsilon_rf == 0.0:  # infinity
        k_rf = 1.0 / (2.0 * rc**3)
    else:
        k_rf = (nb.epsilon_rf - nb.epsilon_r) / ((2.0 * nb.epsilon_rf + nb.epsilon_r) * rc**3)
    c_rf = 1.0 / rc + k_rf * rc * rc
    qq = ONE_4PI_EPS0 * q_i * q_j / nb.epsilon_r

    def lj(r: float) -> tuple[float, float]:
        if epsilon == 0.0 or sigma == 0.0:
            return 0.0, 0.0
        s6 = (sigma / r) ** 6
        return 4.0 * epsilon * (s6 * s6 - s6), 24.0 * epsilon * (2.0 * s6 * s6 - s6) / r

    shift = lj(rc)[0] if nb.vdw_modifier == "potential-shift" else 0.0
    n = round(rc / nb.spacing)
    lines = [
        f"# MARTINI pair table: sigma={sigma:.10g} nm epsilon={epsilon:.10g} kJ/mol "
        f"qi={q_i:g} qj={q_j:g} rc={rc:g} eps_r={nb.epsilon_r:g} eps_rf={nb.epsilon_rf:g} "
        f"vdw-modifier={nb.vdw_modifier}",
    ]
    for step in range(1, n + 1):
        r = step * nb.spacing
        v_lj, f_lj = lj(r)
        v = v_lj - shift + qq * (1.0 / r + k_rf * r * r - c_rf)
        f = f_lj + qq * (1.0 / (r * r) - 2.0 * k_rf * r)
        lines.append(f"{r:.6f} {v:.10e} {f:.10e}")
    return "\n".join(lines) + "\n"


def g96_angle_xvg(theta0: float, k: float) -> str:
    """G96 angle table: theta (deg), V (kJ/mol), -dV/dtheta (kJ/mol/deg), 0-180 in 0.5 deg."""
    c0 = math.cos(math.radians(theta0))
    lines = [f"# G96 angle: theta0={theta0:g} deg k={k:g} kJ/mol"]
    for step in range(361):
        theta = 0.5 * step
        c = math.cos(math.radians(theta))
        s = math.sin(math.radians(theta))
        v = 0.5 * k * (c - c0) ** 2
        # dV/dtheta = k (c - c0) (-sin theta) per radian; per degree times pi/180
        f = k * (c - c0) * s * math.pi / 180.0
        lines.append(f"{theta:.1f} {v:.10e} {f:.10e}")
    return "\n".join(lines) + "\n"


def write_if_changed(path: Path, text: str) -> None:
    if not path.is_file() or path.read_text() != text:
        path.write_text(text)
