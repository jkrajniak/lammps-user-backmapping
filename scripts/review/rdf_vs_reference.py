"""Tier C: element-group RDFs of a LAMMPS trajectory against published gmx rdf data.

Mirrors how the JCC 2018 reference RDFs were made (``compute_gmx_rdf.sh`` in
the data bundle): ``gmx rdf -ref 'name "C*"' -sel 'name "N*"'`` on the NVT
trajectory, all atom pairs (default) or without the topology exclusions
(``-excl``, files ``*_excl.xvg``), bin width 0.002 nm, normalized by the
selection's number density. Elements come from the masses in the LAMMPS data
file; GROMACS selected by the first letter of the atom name, which for these
force fields is the element.

    uv run --with numpy --with scipy --with matplotlib rdf_vs_reference.py \\
        --data rim135_at_final.data --dump dump.at_prod \\
        --ref-dir .../paper/rim135_small/aa/rdf --ref-glob 'rdf_s*_*_{pair}.xvg' \\
        --pair C:O --pair C:N --pair C:N:excl --out tierc_rim135

Writes ``<out>_<pair>.dat`` (r in nm, g), ``<out>_metrics.txt`` and, with
``--plot``, an overlay of the computed RDF and the reference mean and spread.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ELEMENT_MASSES = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "S": 32.06}


@dataclass
class Frame:
    step: int
    box: np.ndarray  # (3,) lengths, A
    ids: np.ndarray
    types: np.ndarray
    xyz: np.ndarray  # A


def element_of_mass(mass: float) -> str:
    element, ref = min(ELEMENT_MASSES.items(), key=lambda kv: abs(kv[1] - mass))
    if abs(ref - mass) > 0.1:
        raise ValueError(f"mass {mass} matches no element")
    return element


def read_data(path: Path) -> tuple[dict[int, str], list[tuple[int, int]]]:
    """Type -> element (from Masses) and the bond list (atom ids)."""
    lines = path.read_text().splitlines()
    elements: dict[int, str] = {}
    bonds: list[tuple[int, int]] = []
    section = None
    for line in lines:
        text = line.split("#")[0].strip()
        if not text:
            continue
        head = text.split()[0]
        if head in {
            "Masses",
            "Atoms",
            "Bonds",
            "Angles",
            "Dihedrals",
            "Impropers",
            "Velocities",
        } or (text.endswith("Coeffs")):
            section = head
            continue
        if section == "Masses":
            t, m = text.split()[:2]
            elements[int(t)] = element_of_mass(float(m))
        elif section == "Bonds":
            _, _, a, b = text.split()[:4]
            bonds.append((int(a), int(b)))
    return elements, bonds


def read_dump(path: Path, first_step: int, last_step: int | None) -> list[Frame]:
    frames: list[Frame] = []
    with path.open() as fh:
        while True:
            line = fh.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            step = int(fh.readline())
            fh.readline()
            n = int(fh.readline())
            fh.readline()
            box = np.array(
                [np.diff([float(v) for v in fh.readline().split()[:2]])[0] for _ in range(3)]
            )
            cols = fh.readline().split()[2:]
            rows = [fh.readline() for _ in range(n)]
            if step < first_step or (last_step is not None and step > last_step):
                continue
            data = np.loadtxt(rows, ndmin=2)
            c = {name: i for i, name in enumerate(cols)}
            order = np.argsort(data[:, c["id"]])
            data = data[order]
            frames.append(
                Frame(
                    step=step,
                    box=box,
                    ids=data[:, c["id"]].astype(int),
                    types=data[:, c["type"]].astype(int),
                    xyz=data[:, [c["x"], c["y"], c["z"]]],
                )
            )
    return frames


def excluded_pairs(bonds: list[tuple[int, int]], nrexcl: int = 3) -> set[tuple[int, int]]:
    """Atom-id pairs within nrexcl bonds (GROMACS topology exclusions)."""
    graph: dict[int, set[int]] = {}
    for a, b in bonds:
        graph.setdefault(a, set()).add(b)
        graph.setdefault(b, set()).add(a)
    pairs: set[tuple[int, int]] = set()
    for start in graph:
        seen = {start: 0}
        frontier = [start]
        for depth in range(1, nrexcl + 1):
            nxt = []
            for u in frontier:
                for v in graph[u]:
                    if v not in seen:
                        seen[v] = depth
                        nxt.append(v)
            frontier = nxt
        pairs |= {(min(start, v), max(start, v)) for v in seen if v != start}
    return pairs


def rdf(
    frames: list[Frame],
    elements: dict[int, str],
    ref_el: str,
    sel_el: str,
    rmax_nm: float,
    dr_nm: float,
    excl: set[tuple[int, int]] | None,
) -> tuple[np.ndarray, np.ndarray]:
    # gmx rdf convention: bin k is centred on k * dr
    centers = np.arange(0.0, rmax_nm + dr_nm / 2, dr_nm)
    edges = np.concatenate([[0.0], centers[:-1] + dr_nm / 2, [centers[-1] + dr_nm / 2]])
    counts = np.zeros(len(edges) - 1)
    norm = 0.0
    for fr in frames:
        el = np.array([elements[t] for t in fr.types])
        ref_idx = np.where(el == ref_el)[0]
        sel_idx = np.where(el == sel_el)[0]
        box_nm = fr.box / 10.0
        xyz = np.mod(fr.xyz / 10.0, box_nm)
        tree_sel = cKDTree(xyz[sel_idx], boxsize=box_nm)
        tree_ref = cKDTree(xyz[ref_idx], boxsize=box_nm)
        pairs = tree_ref.query_ball_tree(tree_sel, rmax_nm)
        ii = np.repeat(np.arange(len(ref_idx)), [len(p) for p in pairs])
        jj = (
            np.concatenate([np.asarray(p, dtype=int) for p in pairs])
            if len(ii)
            else np.array([], dtype=int)
        )
        a, b = ref_idx[ii], sel_idx[jj]
        keep = a != b
        if excl is not None:
            ia, ib = fr.ids[a], fr.ids[b]
            key = np.minimum(ia, ib) * (int(fr.ids.max()) + 1) + np.maximum(ia, ib)
            excl_keys = np.fromiter(
                (p * (int(fr.ids.max()) + 1) + q for p, q in excl), dtype=np.int64, count=len(excl)
            )
            keep &= ~np.isin(key, excl_keys)
        d = xyz[a[keep]] - xyz[b[keep]]
        d -= box_nm * np.round(d / box_nm)
        r = np.linalg.norm(d, axis=1)
        counts += np.histogram(r, bins=edges)[0]
        rho_sel = len(sel_idx) / np.prod(box_nm)
        norm += len(ref_idx) * rho_sel
    shell = 4.0 / 3.0 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    return centers, counts / (norm * shell)


def read_xvg(path: Path) -> tuple[np.ndarray, np.ndarray]:
    rows = [ln.split() for ln in path.read_text().splitlines() if ln and ln[0] not in "#@"]
    arr = np.array(rows, dtype=float)
    return arr[:, 0], arr[:, 1]


def first_peak(r: np.ndarray, g: np.ndarray, r_min: float) -> tuple[float, float]:
    mask = (r >= r_min) & np.isfinite(g)
    k = int(np.argmax(g[mask]))
    return float(r[mask][k]), float(g[mask][k])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--dump", type=Path, required=True)
    parser.add_argument("--first-step", type=int, default=0)
    parser.add_argument("--last-step", type=int, default=None)
    parser.add_argument("--ref-dir", type=Path, required=True)
    parser.add_argument("--ref-glob", default="rdf_s*_*_{pair}.xvg")
    parser.add_argument(
        "--pair", action="append", required=True, help="REF:SEL[:excl], e.g. C:N:excl"
    )
    parser.add_argument("--rmax", type=float, default=1.0, help="nm")
    parser.add_argument("--dr", type=float, default=0.002, help="nm (gmx rdf default)")
    parser.add_argument("--peak-rmin", type=float, default=0.2, help="nm, skip bonded peaks")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()

    elements, bonds = read_data(args.data)
    frames = read_dump(args.dump, args.first_step, args.last_step)
    if not frames:
        print("no frames in the window", file=sys.stderr)
        return 1
    excl = None
    report = [
        f"data {args.data}  dump {args.dump}  frames {len(frames)} "
        f"(steps {frames[0].step}-{frames[-1].step})",
        f"{'pair':10s} {'n_ref':>5s} {'peak r ours':>11s} {'ref mean':>9s} {'dr nm':>7s} "
        f"{'peak g ours':>11s} {'ref mean':>9s} {'ref min-max':>15s} {'L2(g)':>7s}",
    ]
    curves = {}
    for spec in args.pair:
        parts = spec.split(":")
        ref_el, sel_el, use_excl = parts[0], parts[1], len(parts) > 2 and parts[2] == "excl"
        if use_excl and excl is None:
            excl = excluded_pairs(bonds)
        r, g = rdf(frames, elements, ref_el, sel_el, args.rmax, args.dr, excl if use_excl else None)
        name = f"{ref_el}_{sel_el}" + ("_excl" if use_excl else "")
        np.savetxt(f"{args.out}_{name}.dat", np.column_stack([r, g]), header="r_nm g")
        refs = sorted(args.ref_dir.glob(args.ref_glob.format(pair=name)))
        if not refs:
            report.append(f"{name:10s} no reference files for {args.ref_glob.format(pair=name)}")
            continue
        ref_g = []
        for p in refs:
            rr, gg = read_xvg(p)
            ref_g.append(np.interp(r, rr, gg, left=np.nan, right=np.nan))
        ref_g = np.array(ref_g)
        mean = np.nanmean(ref_g, axis=0)
        pr, pg = first_peak(r, g, args.peak_rmin)
        rr_mean, rg_mean = first_peak(r, mean, args.peak_rmin)
        ref_peaks_g = [first_peak(r, x, args.peak_rmin)[1] for x in ref_g]
        valid = ~np.isnan(mean)
        l2 = float(
            np.sqrt(
                np.trapezoid((g[valid] - mean[valid]) ** 2, r[valid])
                / np.trapezoid(mean[valid] ** 2, r[valid])
            )
        )
        report.append(
            f"{name:10s} {len(refs):5d} {pr:11.3f} {rr_mean:9.3f} {pr - rr_mean:7.3f} "
            f"{pg:11.3f} {rg_mean:9.3f} {min(ref_peaks_g):7.3f}-{max(ref_peaks_g):7.3f} {l2:7.3f}"
        )
        curves[name] = (r, g, mean, np.nanmin(ref_g, axis=0), np.nanmax(ref_g, axis=0))
    text = "\n".join(report)
    Path(f"{args.out}_metrics.txt").write_text(text + "\n")
    print(text)
    if args.plot and curves:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, len(curves), figsize=(3.2 * len(curves), 2.8), squeeze=False)
        for ax, (name, (r, g, mean, lo, hi)) in zip(axes[0], curves.items(), strict=True):
            ax.fill_between(r, lo, hi, color="C1", alpha=0.25, lw=0, label="reference range")
            ax.plot(r, mean, "--", color="C1", lw=1.5, label="reference mean")
            ax.plot(r, g, color="C0", lw=1.5, label="backmapped")
            ax.set_title(name.replace("_", "\N{EN DASH}", 1).replace("_excl", " (excl.)"))
            ax.set_xlabel("r / nm")
        axes[0][0].set_ylabel("g(r)")
        axes[0][0].legend(fontsize=7, frameon=False)
        fig.tight_layout()
        fig.savefig(f"{args.out}_rdf.pdf")
    return 0


if __name__ == "__main__":
    sys.exit(main())
