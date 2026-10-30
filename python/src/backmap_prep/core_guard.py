"""Backmapping-time LJ core on the 1-4 pairs that attract each other with no LJ repulsion.

OPLS-AA gives polar hydrogens no LJ (opls_270, the acid H, has sigma = epsilon = 0), so a 1-4
pair of such a hydrogen and a negative atom has only the scaled Coulomb attraction. Nothing keeps
the two apart once the bonded terms that normally do give way. In the PET/Dacron example the
O=C-OH angle is soft (200 kJ/mol/rad^2 in the published topology), and during the hot early stages
of the backmapping the acid H falls onto its own carbonyl O (0.006 A, O-C-O 45 deg) in 91 to 100
of the 243 COOH groups. The capped forces hold the collapsed pairs through the ramp. When they come
apart, the released Coulomb energy heats the acid H and, after the cap release, can crash the run
(``research/experiments/20261027_ring-guard-prototype.md``).

``simulation.core_guard`` gives each such pair a small LJ core (a second ``fix backmap/pairs``,
LJ only) from the first minimization until the hybrid structure is written. At the normal 1-4
distance of the PET pair (0.22 nm) it adds -0.02 kcal/mol; at 0.1 nm it adds +24.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from . import units

if TYPE_CHECKING:
    from pathlib import Path

    from .builder import LammpsCrossPair, System


def unguarded_pairs(system: System) -> list[LammpsCrossPair]:
    """1-4 pairs with Coulomb attraction (scale > 0, opposite charges) and no LJ term."""
    charge = {a.atom_id: a.charge for a in system.atoms}
    return [
        p
        for p in system.cross_pairs
        if p.qq_scale > 0.0
        and (p.sigma == 0.0 or p.epsilon == 0.0)
        and charge[p.i] * charge[p.j] < 0.0
    ]


def write_core_guard_files(
    directory: Path, system: System, sigma_nm: float, epsilon_kj: float
) -> int:
    """Write ``core_guard.dat`` and the fix that applies it; return the number of pairs.

    ``core_guard.lmp`` defines ``fix backmap/pairs`` on the pairs of :func:`unguarded_pairs`
    with an LJ term only, sigma from nm and epsilon from kJ/mol. It keeps the lambda weighting of
    the 1-4 pairs (``at``). ``core_guard_release.lmp`` removes it. Nothing is written when no
    pair needs a core.
    """
    pairs = unguarded_pairs(system)
    if not pairs:
        return 0
    sigma = units.distance(sigma_nm)
    epsilon = units.energy(epsilon_kj)
    lines = [f"{len(pairs)}"]
    lines += [f"{p.i} {p.j} {sigma:.6g} {epsilon:.6g}" for p in pairs]
    (directory / "core_guard.dat").write_text("\n".join(lines) + "\n")
    (directory / "core_guard.lmp").write_text(
        "# core guard: LJ core on the 1-4 pairs with Coulomb attraction and no LJ"
        " (backmap-prep, simulation.core_guard)\n"
        "fix core_guard all backmap/pairs at file core_guard.dat\n"
    )
    (directory / "core_guard_release.lmp").write_text(
        "# release of the core guard\nunfix core_guard\n"
    )
    return len(pairs)
