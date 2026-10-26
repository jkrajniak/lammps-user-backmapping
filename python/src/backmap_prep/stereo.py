"""Stereocenters of the atomistic part and the backmapping-time restraint that keeps them.

A fragment placed on a CG bead carries the handedness of the stereocenters inside it, but a
center whose neighbors lie in two beads has its handedness set by the CG geometry, and the
atomistic force field has no term that protects it (an inversion needs about one bond-angle
distortion). ``simulation.stereo_restraint`` holds such centers at the handedness they have in
the hybrid structure during the backmapping stages; the restraint is released before
production.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from . import units

if TYPE_CHECKING:
    from pathlib import Path

    from .builder import System

MORGAN_ITERATIONS = 6
HYDROGEN_MASS = 1.5
# Improper dihedral of the neighbors around the center of an ideal tetrahedron:
# arctan(1/sqrt(2)).
IDEAL_IMPROPER_DEG = 35.2644
# fix restrain dihedral evaluates K[1 + cos(phi - phi0)] with the phi0 it is given, so its energy
# minimum lies at phi0 + 180 deg, while its documentation puts the minimum at phi0 (LAMMPS
# stable_22Jul2025, stable_30Sep2026 and develop of 2026-10-25). The written target is shifted by
# this offset; test_lammps_stereo_restraint.py fails if LAMMPS changes the convention.
RESTRAIN_MINIMUM_OFFSET_DEG = 180.0


def _dihedral(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    """Dihedral angle in degrees with the IUPAC sign, the sign LAMMPS uses."""
    b0, b1, b2 = p1 - p0, p2 - p1, p3 - p2
    n1, n2 = np.cross(b0, b1), np.cross(b1, b2)
    m = np.cross(b1 / np.linalg.norm(b1), n1)
    return float(np.degrees(np.arctan2(np.dot(m, n2), np.dot(n1, n2))))


def restrain_target(phi_deg: float) -> float:
    """The phi0 of ``fix restrain dihedral`` whose energy minimum lies at ``phi_deg``."""
    return (phi_deg - RESTRAIN_MINIMUM_OFFSET_DEG + 180.0) % 360.0 - 180.0


def _classes(types: dict[int, int], adjacency: dict[int, list[int]]) -> dict[int, int]:
    """Morgan-type equivalence classes: atoms with the same type and the same environment."""
    label = dict(types)
    for _ in range(MORGAN_ITERATIONS):
        keys = {a: (label[a], tuple(sorted(label[n] for n in adjacency[a]))) for a in label}
        ids: dict[tuple[int, tuple[int, ...]], int] = {}
        label = {a: ids.setdefault(k, len(ids)) for a, k in keys.items()}
    return label


def find_stereocenters(system: System) -> list[tuple[int, int, int, int, float]]:
    """Atoms with four bonded atomistic neighbors that are all inequivalent.

    Returns two quartets ``(a, b, x, center, phi0)`` per center: a and b are the two heavy
    neighbors of lowest class, x is the third heavy neighbor in the first quartet and the fourth
    neighbor (hydrogen when there is one) in the second, and phi0 is the target improper dihedral
    of the quartet in degrees. Its magnitude is the ideal tetrahedral value and its sign is the
    handedness in ``system`` (IUPAC sign). The first quartet alone keeps the center on its side of
    the heavy neighbors but lets the hydrogen cross to that side under a stiff restraint, which
    leaves a nearly planar center that inverts once the restraint is released.
    """
    mass = {t.type_id: t.mass for t in system.atom_types}
    pos = {a.atom_id: np.array([a.x, a.y, a.z]) for a in system.atoms if not a.is_cg}
    atype = {a.atom_id: a.type_id for a in system.atoms if not a.is_cg}
    adjacency: dict[int, list[int]] = {i: [] for i in pos}
    for bond in system.bonds:
        if bond.i in pos and bond.j in pos:
            adjacency[bond.i].append(bond.j)
            adjacency[bond.j].append(bond.i)
    cls = _classes(atype, adjacency)
    box = np.array(system.box, dtype=float)

    def image(vec: np.ndarray) -> np.ndarray:
        return vec - box * np.round(vec / box) if np.all(box > 0) else vec

    out: list[tuple[int, int, int, int, float]] = []
    for center in sorted(pos):
        nbrs = adjacency[center]
        if len(nbrs) != 4 or len({cls[n] for n in nbrs}) != 4:
            continue
        heavy = [n for n in nbrs if mass.get(atype[n], 12.0) >= HYDROGEN_MASS]
        chosen = sorted(heavy if len(heavy) == 3 else nbrs, key=lambda n: cls[n])[:3]
        a, b, c = chosen
        (d,) = set(nbrs) - set(chosen)
        origin = pos[center]
        phis = [
            _dihedral(
                origin + image(pos[a] - origin),
                origin + image(pos[b] - origin),
                origin + image(pos[x] - origin),
                origin,
            )
            for x in (c, d)
        ]
        if min(abs(phi) for phi in phis) < 1e-3:
            continue
        for x, phi in zip((c, d), phis, strict=True):
            out.append((a, b, x, center, float(np.sign(phi)) * IDEAL_IMPROPER_DEG))
    return out


def write_stereo_files(directory: Path, system: System, k_kj: float) -> int:
    """Write ``stereo.lmp`` (fix lines) and ``stereo_release.lmp`` (unfix lines).

    Returns the number of restrained centers (two restraints each).
    """
    centers = find_stereocenters(system)
    k = units.energy(k_kj)
    lines = ["# stereocenter restraints (generated by backmap-prep, simulation.stereo_restraint)"]
    release = ["# release of the stereocenter restraints"]
    for n, (a, b, c, center, phi0) in enumerate(centers, 1):
        lines.append(
            f"fix stereo{n} all restrain dihedral {a} {b} {c} {center} {k:.6g} {k:.6g} "
            f"{restrain_target(phi0):.4f}"
        )
        release.append(f"unfix stereo{n}")
    (directory / "stereo.lmp").write_text("\n".join(lines) + "\n")
    (directory / "stereo_release.lmp").write_text("\n".join(release) + "\n")
    return len({center for _, _, _, center, _ in centers})
