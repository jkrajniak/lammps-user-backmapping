"""Time-averaged RDF computation for the regenerated (pet_v3, noshake)
PET/Dacron production run. Adapted from pet/pet_rdf_multiframe.py for the
new atom-type numbering (this system has 18 total types: 6 CG + 12 AT,
water O=17/H=18) and the new dump file, which spans the whole run
(Phase 0-3, every 1000 steps) rather than a dedicated post-production
sampling continuation -- so this script filters to production-only frames
(step >= PROD_START) before accumulating histograms.

See research/experiments/20260726_pet-rdf-tier-c-validation.md for full
method notes (reference selection semantics, exclusion/ring-COM logic,
peak-finding convention) which this replicates unchanged.
"""
from __future__ import annotations
import sys
from collections import defaultdict, deque
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

DATA_FILE = Path("/home/sc/sc/pet_v3/pet_pipeline_check.data")
DUMP_FILE = Path("/home/sc/sc/pet_v3/dump.pet_bakery_faithful.custom")
OUT_DIR = Path("/home/sc/sc/pet_v3/rdf_out_bakery")

# Type numbering for this pipeline-regenerated system (verified against the
# Masses section: element identified by OPLS comment / mass, not assumed).
CARBON = {2, 9, 12}
OXYGEN = {4, 8, 10, 17}
HYDROGEN = {3, 5, 13, 15, 18}
RING_C_TYPE = 12  # opls_145, aromatic ring carbon

# Phase 3 (NVT production at lambda=1) runs from step 66000 to 76000 in this
# run (Phase 0 minimize ~50000 steps + Phase 1 6000 + Phase 2 ramp 10000).
# Only these frames reflect the final, fully-backmapped AT structure.
PROD_START = 12000

R_BIN = 0.05
R_MAX_SHORT = 10.0
R_MAX_LONG = 20.0


def parse_topology(path: Path):
    lines = path.read_text().splitlines()
    natoms = nbonds = 0
    i = 0
    while i < len(lines):
        s = lines[i].strip()
        if s.endswith("atoms"):
            natoms = int(s.split()[0])
        elif s.endswith("bonds"):
            nbonds = int(s.split()[0])
        elif s == "Atoms" or s.startswith("Atoms "):
            i += 2
            break
        i += 1
    ids = np.empty(natoms, dtype=np.int64)
    types = np.empty(natoms, dtype=np.int32)
    n = 0
    while n < natoms:
        parts = lines[i].split()
        i += 1
        if not parts:
            continue
        aid = int(parts[0])
        typ = int(parts[2])
        ids[aid - 1] = aid
        types[aid - 1] = typ
        n += 1
    bonds = np.empty((nbonds, 2), dtype=np.int64)
    while i < len(lines):
        s = lines[i].strip()
        if s == "Bonds" or s.startswith("Bonds "):
            i += 2
            break
        i += 1
    n = 0
    while n < nbonds:
        parts = lines[i].split()
        i += 1
        if not parts:
            continue
        bonds[n] = (int(parts[2]), int(parts[3]))
        n += 1
    return ids, types, bonds


def iter_frames(path: Path, step_min: int):
    lines = path.read_text().splitlines()
    i = 0
    n = len(lines)
    while i < n:
        if lines[i].strip() != "ITEM: TIMESTEP":
            i += 1
            continue
        i += 1
        step = int(lines[i].strip())
        i += 1
        i += 1  # ITEM: NUMBER OF ATOMS
        natoms = int(lines[i].strip())
        i += 1
        i += 1  # ITEM: BOX BOUNDS
        lo = np.empty(3)
        hi = np.empty(3)
        for d in range(3):
            a, b = lines[i].split()[:2]
            lo[d], hi[d] = float(a), float(b)
            i += 1
        header = lines[i].split()  # e.g. ITEM: ATOMS id mol type x y z f_bm
        cols = header[2:]
        ix, iy, iz = cols.index("x"), cols.index("y"), cols.index("z")
        i += 1
        if step < step_min:
            i += natoms
            continue
        ids = np.empty(natoms, dtype=np.int64)
        pos = np.empty((natoms, 3))
        for k in range(natoms):
            parts = lines[i].split()
            i += 1
            aid = int(parts[0])
            ids[k] = aid
            pos[k] = (float(parts[ix]), float(parts[iy]), float(parts[iz]))
        yield step, lo, hi, ids, pos


def build_exclusions(ids: np.ndarray, bonds: np.ndarray):
    adj: dict[int, set[int]] = defaultdict(set)
    for a1, a2 in bonds:
        adj[a1].add(a2)
        adj[a2].add(a1)
    excl: dict[int, set[int]] = {}
    for aid in ids:
        aid = int(aid)
        visited = {aid: 0}
        q = deque([aid])
        while q:
            cur = q.popleft()
            d = visited[cur]
            if d >= 3:
                continue
            for nb in adj.get(cur, ()):
                if nb not in visited:
                    visited[nb] = d + 1
                    q.append(nb)
        visited.pop(aid)
        excl[aid] = set(visited.keys())
    return excl


def find_rings(ids, types, bonds):
    ring_ids = set(int(x) for x in ids[types == RING_C_TYPE])
    adj = defaultdict(set)
    for a1, a2 in bonds:
        a1, a2 = int(a1), int(a2)
        if a1 in ring_ids and a2 in ring_ids:
            adj[a1].add(a2)
            adj[a2].add(a1)
    seen = set()
    rings = []
    for aid in ring_ids:
        if aid in seen:
            continue
        comp = set()
        q = deque([aid])
        seen.add(aid)
        while q:
            cur = q.popleft()
            comp.add(cur)
            for nb in adj[cur]:
                if nb not in seen:
                    seen.add(nb)
                    q.append(nb)
        if len(comp) != 6:
            print(f"WARNING: ring size {len(comp)} for {aid}", file=sys.stderr)
        rings.append(sorted(comp))
    return rings


def hist_cross(posA, posB, L, r_max, same_group, excl_pairs=None, idsA=None, idsB=None):
    treeA = cKDTree(posA, boxsize=L)
    if same_group:
        pairs = treeA.query_pairs(r_max, output_type="ndarray")
        d = posA[pairs[:, 0]] - posA[pairs[:, 1]]
        d -= L * np.round(d / L)
        dist = np.linalg.norm(d, axis=1)
    else:
        treeB = cKDTree(posB, boxsize=L)
        coo = treeA.sparse_distance_matrix(treeB, r_max, output_type="coo_matrix")
        rows, cols, dist = coo.row, coo.col, coo.data
        if excl_pairs is not None:
            keep = np.ones(len(dist), dtype=bool)
            for k, (r, c) in enumerate(zip(rows, cols)):
                a, b = int(idsA[r]), int(idsB[c])
                if b in excl_pairs.get(a, ()):
                    keep[k] = False
            dist = dist[keep]
    nbins = int(r_max / R_BIN)
    hist, edges = np.histogram(dist, bins=nbins, range=(0.0, r_max))
    return hist, edges


def first_peak(r, g, r_min=1.5):
    mask = r > r_min
    idx = np.argmax(g[mask])
    return float(r[mask][idx]), float(g[mask][idx])


def main() -> None:
    print(f"Parsing topology from {DATA_FILE} ...")
    ids, types, bonds = parse_topology(DATA_FILE)
    print(f"  {len(ids)} atoms, {len(bonds)} bonds")
    excl = build_exclusions(ids, bonds)
    rings = find_rings(ids, types, bonds)
    print(f"  {len(rings)} rings found")
    idx_by_id = {int(a): i for i, a in enumerate(ids)}
    idx_by_type = {t: np.where(types == t)[0] for t in set(types.tolist())}

    def group_idx(type_set):
        return np.concatenate([idx_by_type[t] for t in type_set if t in idx_by_type])

    c_idx = group_idx(CARBON)
    o_idx = group_idx(OXYGEN)
    h_idx = group_idx(HYDROGEN)
    n_c, n_o, n_h = len(c_idx), len(o_idx), len(h_idx)
    print(f"  C={n_c} O={n_o} H={n_h} (sum={n_c + n_o + n_h}, AT total={int((types >= 2).sum())})")
    n_ring = len(rings)
    ring_atom_idx = [[idx_by_id[a] for a in ring] for ring in rings]

    accum = {
        "C_O_excl": np.zeros(int(R_MAX_LONG / R_BIN)),
        "C_O": np.zeros(int(R_MAX_SHORT / R_BIN)),
        "C_H": np.zeros(int(R_MAX_SHORT / R_BIN)),
        "O_H": np.zeros(int(R_MAX_SHORT / R_BIN)),
        "ring_ring": np.zeros(int(R_MAX_LONG / R_BIN)),
    }
    n_frames = 0
    V_sum = 0.0
    steps_used = []
    for step, lo, hi, fids, fpos in iter_frames(DUMP_FILE, PROD_START):
        L = hi - lo
        pos_w = np.mod(fpos - lo, L)
        assert np.array_equal(fids, ids), "dump atom order mismatch"
        h, _ = hist_cross(pos_w[c_idx], pos_w[o_idx], L, R_MAX_LONG, False,
                           excl_pairs=excl, idsA=ids[c_idx], idsB=ids[o_idx])
        accum["C_O_excl"] += h
        h, _ = hist_cross(pos_w[c_idx], pos_w[o_idx], L, R_MAX_SHORT, False)
        accum["C_O"] += h
        h, _ = hist_cross(pos_w[c_idx], pos_w[h_idx], L, R_MAX_SHORT, False)
        accum["C_H"] += h
        h, _ = hist_cross(pos_w[o_idx], pos_w[h_idx], L, R_MAX_SHORT, False)
        accum["O_H"] += h
        ring_com = []
        for atoms in ring_atom_idx:
            ref = pos_w[atoms[0]]
            rel = pos_w[atoms] - ref
            rel -= L * np.round(rel / L)
            com = np.mod(ref + rel.mean(axis=0), L)
            ring_com.append(com)
        ring_com = np.array(ring_com)
        h, _ = hist_cross(ring_com, ring_com, L, R_MAX_LONG, True)
        accum["ring_ring"] += h
        n_frames += 1
        steps_used.append(step)
        V_sum += float(np.prod(L))
    print(f"  averaged over {n_frames} frames, steps={steps_used}")
    if n_frames == 0:
        print("ERROR: no production frames found", file=sys.stderr)
        sys.exit(1)
    V_avg = V_sum / n_frames
    results = {}
    for name, hist in accum.items():
        r_max = R_MAX_LONG if name in ("C_O_excl", "ring_ring") else R_MAX_SHORT
        nbins = len(hist)
        edges = np.linspace(0, r_max, nbins + 1)
        r = 0.5 * (edges[:-1] + edges[1:])
        shell_vol = 4.0 * np.pi * r**2 * R_BIN
        if name == "ring_ring":
            ideal = n_frames * n_ring * (n_ring - 1) / (2.0 * V_avg)
        else:
            nA, nB = {"C_O_excl": (n_c, n_o), "C_O": (n_c, n_o),
                      "C_H": (n_c, n_h), "O_H": (n_o, n_h)}[name]
            ideal = n_frames * nA * nB / V_avg
        with np.errstate(divide="ignore", invalid="ignore"):
            g = hist / (shell_vol * ideal)
        results[name] = (r, g)
        pr, pg = first_peak(r, g)
        print(f"{name:12s} first peak: r={pr:.3f} g={pg:.3f}")
    OUT_DIR.mkdir(exist_ok=True)
    for name, (r, g) in results.items():
        np.savetxt(OUT_DIR / f"rdf_{name}_backmap_avg.txt", np.column_stack([r, g]),
                   header=f"r(A) g(r)  [{n_frames}-frame average, steps={steps_used}]")
    print(f"Wrote {len(results)} averaged RDF files to {OUT_DIR}")


if __name__ == "__main__":
    main()
