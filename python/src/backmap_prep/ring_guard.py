"""Rings of the atomistic part and the backmapping-time guard that keeps atoms out of them.

At the start of the lambda ramp the bonded terms act but the atomistic nonbonded repulsion is
still scaled down, so a substituent can fold into its own ring and a bond ends up passing
through it (ring spearing). Once the full potential is on, the defect cannot undo itself. In the
epoxy example the CH2-NH2 arm of isophorone diamine does this
(``research/experiments/20261027_rim135-small-r-peaks.md``). ``simulation.ring_guard`` puts a
repulsive wall of radius rc around the centroid of every ring (``fix backmap/ringguard``) from the
first minimization until the hybrid structure is written.
"""

from __future__ import annotations

from collections import deque
from typing import TYPE_CHECKING

from . import units

if TYPE_CHECKING:
    from pathlib import Path

    from .builder import System

# Ring sizes guarded. Smaller rings (epoxide, cyclopropane) have their centroid closer than rc
# to their own substituents.
MIN_RING_SIZE = 5
MAX_RING_SIZE = 7


def _shortest_cycle_through(
    adjacency: dict[int, list[int]], u: int, v: int, max_size: int
) -> tuple[int, ...] | None:
    """Shortest path from u to v that does not use the bond u-v, as a ring of at most max_size."""
    parent = {u: u}
    queue: deque[tuple[int, int]] = deque([(u, 1)])
    while queue:
        node, size = queue.popleft()
        if size >= max_size:
            continue
        for nxt in adjacency[node]:
            if node == u and nxt == v:
                continue
            if nxt in parent:
                continue
            parent[nxt] = node
            if nxt == v:
                path = [v]
                while path[-1] != u:
                    path.append(parent[path[-1]])
                return tuple(reversed(path))
            queue.append((nxt, size + 1))
    return None


def _adjacency(system: System) -> dict[int, list[int]]:
    """Bonded neighbours of every atomistic atom (CG beads and bonds to them left out)."""
    at_ids = {a.atom_id for a in system.atoms if not a.is_cg}
    adjacency: dict[int, list[int]] = {i: [] for i in at_ids}
    for bond in system.bonds:
        if bond.i in at_ids and bond.j in at_ids:
            adjacency[bond.i].append(bond.j)
            adjacency[bond.j].append(bond.i)
    for nbrs in adjacency.values():
        nbrs.sort()
    return adjacency


def find_rings(system: System) -> list[tuple[int, ...]]:
    """Rings of 5 to 7 atomistic atoms: the shortest cycle through each bond, deduplicated.

    Atoms are listed in ring order starting from the lowest ID. Fused rings give one ring each
    (naphthalene: two six-membered rings); the envelope cycle of a fused system is not the
    shortest cycle through any bond and is not listed.
    """
    adjacency = _adjacency(system)
    seen: set[frozenset[int]] = set()
    rings: list[tuple[int, ...]] = []
    for u in sorted(adjacency):
        for v in adjacency[u]:
            if v < u:
                continue
            cycle = _shortest_cycle_through(adjacency, u, v, MAX_RING_SIZE)
            if cycle is None or len(cycle) < MIN_RING_SIZE:
                continue
            key = frozenset(cycle)
            if key in seen:
                continue
            seen.add(key)
            start = cycle.index(min(cycle))
            rings.append(cycle[start:] + cycle[:start])
    return sorted(rings)


def write_ring_guard_files(
    directory: Path, system: System, radius_nm: float, k_kj: float, bonded_nm: float = 0.12
) -> int:
    """Write ``rings.txt`` and the fix that applies the guard; return the number of rings.

    ``ring_guard.lmp`` defines ``fix backmap/ringguard`` on ``at_atoms`` with
    E = K (rc - d)^2 for an atom at distance d < rc from a ring centroid, K converted from
    kJ/mol/nm^2 and rc from nm. ``ring_guard_release.lmp`` removes it. Nothing is written when
    there are no rings.

    Each line of ``rings.txt`` is a ring, then ``|`` and the atoms bonded to it (its hydrogens
    and the first atom of each substituent). They are repelled only inside ``bonded_nm``: an
    axial hydrogen of cyclohexane sits about 0.18 nm from the centroid, while a substituent
    folded into its own ring (seen in the epoxy example) comes within 0.07-0.11 nm.
    """
    rings = find_rings(system)
    if not rings:
        return 0
    adjacency = _adjacency(system)
    k = units.energy(k_kj) / units.distance(1.0) ** 2
    rc = units.distance(radius_nm)
    rcb = units.distance(bonded_nm)
    table = [
        "# atom IDs of one ring per line, then | and the atoms bonded to it, which are repelled"
        " only inside the smaller radius (backmap-prep, simulation.ring_guard)"
    ]
    for ring in rings:
        bonded = sorted({b for a in ring for b in adjacency[a]} - set(ring))
        table.append(" ".join(map(str, ring)) + " | " + " ".join(map(str, bonded)))
    (directory / "rings.txt").write_text("\n".join(table) + "\n")
    (directory / "ring_guard.lmp").write_text(
        "# ring guard: repulsive wall around ring centroids (backmap-prep, simulation.ring_guard)\n"
        f"fix ring_guard at_atoms backmap/ringguard {k:.6g} {rc:.6g} {rcb:.6g} rings.txt\n"
        "fix_modify ring_guard energy yes\n"
        "# energy and number of atoms inside rc, every 1000 steps\n"
        'fix ring_guard_log all print 1000 "$(step) $(f_ring_guard[1]:%.6g) $(f_ring_guard[2]:%.0f)"'
        ' file ring_guard.log screen no title "# step energy contacts"\n'
    )
    (directory / "ring_guard_release.lmp").write_text(
        "# release of the ring guard\nunfix ring_guard_log\nunfix ring_guard\n"
    )
    return len(rings)
