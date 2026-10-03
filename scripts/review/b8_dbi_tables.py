"""Review item B8 (R1-3, co-author Z1): CG pair tables from the atomistic system by direct Boltzmann inversion.

A CG force field derived from the atomistic reference in one step, with no iterations:
V(r) = -kT ln g(r), where g(r) is the bead-centre radial distribution of the independent atomistic
melt, mapped to the CG beads (two united-atom sites per bead, mass-weighted centre).

Pairs the CG force field excludes (same chain and at most three bonds apart, as in the generated
``special_bonds lj 0 0 0``) are left out of g(r); g is normalised to the number of non-excluded pairs.
Where g(r) is small the potential is continued linearly into the core. Beyond the cut-off the
potential is shifted to zero. The tables replace ``table_A_A``, ``table_A_B``, ``table_B_B`` of an
example copy in the 7-column GROMACS layout of the IBI files (potential in column 6, force in column 7).
Bead types: A = the two end beads of a chain, B = the inner beads.

    uv run --with numpy --with MDAnalysis b8_dbi_tables.py --dir <example copy> \\
        --data dodecane_at_ref_final.data --traj traj_reference.dcd [--temperature 298] [--rc 1.4]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

KB_KJ = 0.0083144626  # kJ/mol/K
EXCLUDED_BONDS = 3
BIN_NM = 0.005
SITES_PER_BEAD = 2
G_CORE = 0.1  # g below this: the wall is continued linearly
PAIRS = {"A_A": ("A", "A"), "A_B": ("A", "B"), "B_B": ("B", "B")}


def bead_layout(universe) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per bead: first atom index, chain index, position in the chain, type (0 = A, 1 = B)."""
    atoms = universe.atoms
    mol = np.asarray(atoms.resids)
    n_atoms = len(atoms)
    order = np.argsort(atoms.ids, kind="stable")
    if not np.array_equal(order, np.arange(n_atoms)):
        raise ValueError("atoms are not in id order; sort the topology first")
    chain_ids, start = np.unique(mol, return_index=True)
    length = np.diff(np.append(start, n_atoms))
    if len(set(length)) != 1 or length[0] % SITES_PER_BEAD:
        raise ValueError(f"chains of unequal or odd length: {sorted(set(length))}")
    beads_per_chain = int(length[0] // SITES_PER_BEAD)
    n_chains = len(chain_ids)
    first = (start[:, None] + SITES_PER_BEAD * np.arange(beads_per_chain)[None, :]).ravel()
    chain = np.repeat(np.arange(n_chains), beads_per_chain)
    pos = np.tile(np.arange(beads_per_chain), n_chains)
    kind = np.where((pos == 0) | (pos == beads_per_chain - 1), 0, 1)
    return first, chain, pos, kind


def n_pairs(kind: np.ndarray, beads_per_chain: int, a: int, b: int) -> int:
    """Number of non-excluded unordered pairs of types a, b."""
    n_chains = len(kind) // beads_per_chain
    n_a, n_b = int((kind == a).sum()), int((kind == b).sum())
    total = n_a * (n_a - 1) // 2 if a == b else n_a * n_b
    pos = np.arange(beads_per_chain)
    k = np.where((pos == 0) | (pos == beads_per_chain - 1), 0, 1)
    excluded = 0
    for i in range(beads_per_chain):
        for j in range(i + 1, beads_per_chain):
            if j - i <= EXCLUDED_BONDS and sorted((k[i], k[j])) == sorted((a, b)):
                excluded += 1
    return total - excluded * n_chains


def bead_rdf(universe, rc_nm: float) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    from MDAnalysis.lib.distances import capped_distance

    first, chain, pos, kind = bead_layout(universe)
    beads_per_chain = int(pos.max()) + 1
    masses = np.asarray(universe.atoms.masses)
    weights = np.stack([masses[first + s] for s in range(SITES_PER_BEAD)], axis=1)
    edges = np.arange(0.0, rc_nm + BIN_NM, BIN_NM)
    counts = {name: np.zeros(len(edges) - 1) for name in PAIRS}
    pair_n = {
        name: n_pairs(kind, beads_per_chain, "AB".index(a), "AB".index(b))
        for name, (a, b) in PAIRS.items()
    }
    volume_sum, n_frames = 0.0, 0
    for ts in universe.trajectory:
        box = ts.dimensions
        x = np.asarray(ts.positions)
        sites = np.stack([x[first + s] for s in range(SITES_PER_BEAD)], axis=1)
        com = (sites * weights[..., None]).sum(axis=1) / weights.sum(axis=1)[:, None]
        pairs, dist = capped_distance(
            com, com, max_cutoff=rc_nm * 10.0, box=box, return_distances=True
        )
        keep = pairs[:, 0] < pairs[:, 1]
        pairs, dist = pairs[keep], dist[keep] / 10.0
        same = chain[pairs[:, 0]] == chain[pairs[:, 1]]
        near = np.abs(pos[pairs[:, 0]] - pos[pairs[:, 1]]) <= EXCLUDED_BONDS
        ok = ~(same & near)
        pairs, dist = pairs[ok], dist[ok]
        ka, kb = kind[pairs[:, 0]], kind[pairs[:, 1]]
        for name, (a, b) in PAIRS.items():
            ia, ib = "AB".index(a), "AB".index(b)
            sel = ((ka == ia) & (kb == ib)) | ((ka == ib) & (kb == ia))
            counts[name] += np.histogram(dist[sel], bins=edges)[0]
        volume_sum += box[0] * box[1] * box[2] / 1000.0  # nm^3
        n_frames += 1
    volume = volume_sum / n_frames
    shell = 4.0 * np.pi / 3.0 * (edges[1:] ** 3 - edges[:-1] ** 3)
    g = {name: counts[name] / n_frames / (pair_n[name] * shell / volume) for name in PAIRS}
    return 0.5 * (edges[1:] + edges[:-1]), g


def potential(r_bins: np.ndarray, g: np.ndarray, r_grid: np.ndarray, kt: float, rc: float):
    """V(r) = -kT ln g on the table grid; linear core where g < G_CORE; shifted to 0 at rc."""
    good = g >= G_CORE
    first = int(np.argmax(good))
    v_bins = np.full_like(g, np.nan)
    v_bins[good] = -kt * np.log(g[good])
    r0, r1 = r_bins[first], r_bins[first + 3]
    slope = (v_bins[first + 3] - v_bins[first]) / (r1 - r0)
    if slope >= 0:
        raise ValueError("potential does not rise towards small r at the core edge")
    v = np.interp(r_grid, r_bins[good], v_bins[good], left=np.nan, right=0.0)
    core = r_grid < r_bins[first]
    v[core] = v_bins[first] + slope * (r_grid[core] - r_bins[first])
    inside = r_grid <= rc
    v[inside] -= v[np.searchsorted(r_grid, rc)]
    v[~inside] = 0.0
    k = 5
    force = -np.gradient(v, r_grid)
    force = np.convolve(np.pad(force, k // 2, mode="edge"), np.ones(k) / k, mode="valid")
    return v, force, float(r_bins[first]), float(slope)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--data", type=Path, required=True, help="LAMMPS data file of the AT reference")
    ap.add_argument("--traj", type=Path, required=True, help="DCD trajectory of the AT reference")
    ap.add_argument("--temperature", type=float, default=298.0)
    ap.add_argument("--rc", type=float, default=1.4, help="CG cut-off, nm")
    args = ap.parse_args()

    from MDAnalysis import Universe

    universe = Universe(str(args.data), str(args.traj), format="DCD", topology_format="DATA")
    r_bins, g = bead_rdf(universe, args.rc)
    kt = KB_KJ * args.temperature
    summary = {}
    for name in PAIRS:
        path = args.dir / f"table_{name}.xvg"
        r_grid = np.loadtxt(path, comments=("#", "@"))[:, 0]
        v, f, r_core, slope = potential(r_bins, g[name], r_grid, kt, args.rc)
        out = np.zeros((len(r_grid), 7))
        out[:, 0], out[:, 5], out[:, 6] = r_grid, v, f
        np.savetxt(
            path,
            out,
            fmt="%.10e",
            header=f"B8 direct Boltzmann inversion of the bead RDF of the atomistic reference, "
            f"T {args.temperature} K, shifted at {args.rc} nm, linear core below {r_core:.3f} nm",
            comments="# ",
        )
        i_min = int(np.argmin(v))
        summary[name] = {
            "first_peak_nm": float(r_bins[int(np.argmax(g[name]))]),
            "g_max": float(g[name].max()),
            "v_min_kJmol": float(v[i_min]),
            "r_min_nm": float(r_grid[i_min]),
            "core_edge_nm": r_core,
            "core_slope_kJmol_nm": slope,
        }
        print(
            f"{name}: g_max {summary[name]['g_max']:.3f} at {summary[name]['first_peak_nm']:.3f} nm, "
            f"V_min {v[i_min]:.3f} kJ/mol at {r_grid[i_min]:.3f} nm, core edge {r_core:.3f} nm"
        )
    (args.dir / "b8_tables.json").write_text(
        json.dumps(
            {"variant": "dbi", "T": args.temperature, "rc_nm": args.rc, "pairs": summary}, indent=2
        )
        + "\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
