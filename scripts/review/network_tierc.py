"""Tier C for the networks: backmapped RDFs against the published GROMACS reference RDFs.

Mirrors how the JCC 2018 reference RDFs were made (``compute_gmx_rdf.sh`` of the data bundle, ``gmx rdf``
with name-glob selections, all pairs or without the topology exclusions, residue centres of mass for the
rings) on the frames of the RDF window of the atomistic continuation, with the metrics and tolerances of
the melts: first-peak position (+-0.03 nm), first-peak height (+-20 %) and the RMS difference of g(r)
(0.2) against the mean of the reference runs, quoted with the standard error over blocks of frames.
The bonded range r < ``--rmin`` (0.2 nm) is left out of the peak search and of the RMS difference.

Inputs: the atomistic data file (masses, bonds), the merged dump of the window (``dump.at_prod``) and the
AT topology (``at_hyb_topol.top``: atom names and residues; its atom ids are those of the data file).

    uv run --with numpy --with scipy --with matplotlib network_tierc.py --system pet \\
        --data pet_at_final.data --dump dump.at_prod --top at_hyb_topol.top \\
        --ref-root <bundle>/paper-reverse-mapping-polymer-networks/paper --out tierc_pet --plot
"""

from __future__ import annotations

import argparse
import fnmatch
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rdf_vs_reference as rv

PEAK_POS_TOL, PEAK_HT_TOL, RMS_TOL = 0.03, 0.20, 0.20  # nm, relative, absolute


@dataclass(frozen=True)
class Pair:
    """One RDF: reference file name stem, groups, and whether topology exclusions apply."""

    name: str  # e.g. C_O, C_O_excl, ring_ring
    ref: str  # element symbol, or a COM group key
    sel: str
    excl: bool = False


# reference RDFs of the JCC 2018 bundle (paper/<dir>/rdf/rdf_s*_*_<name>.xvg) and the COM groups
# (residue-name glob : atom-name glob) of its compute_gmx_rdf.sh
SYSTEMS = {
    "pet": {
        "ref_dir": "dacron/rdf",
        "pairs": [
            Pair("C_O", "C", "O"),
            Pair("C_O_excl", "C", "O", True),
            Pair("C_H", "C", "H"),
            Pair("O_H", "O", "H"),
            Pair("ring_ring", "ring", "ring"),
        ],
        "com": {"ring": ("TER", "C[2-7]")},
    },
    "rim135": {
        "ref_dir": "rim135_small/aa/rdf",
        "pairs": [
            Pair("C_O", "C", "O"),
            Pair("C_N", "C", "N"),
            Pair("C_N_excl", "C", "N", True),
        ],
        "com": {},
    },
    "melamine_network": {
        "ref_dir": "mf3/rdf",
        "pairs": [
            Pair("C_O", "C", "O"),
            Pair("C_N", "C", "N"),
            Pair("C_H", "C", "H"),
            Pair("O_H", "O", "H"),
            Pair("N_H", "N", "H"),
            Pair("O_N", "O", "N"),
            Pair("ring_ring", "ring", "ring"),
        ],
        "com": {"ring": ("*", "[CN]1[1-3]")},
    },
}


@dataclass
class Topology:
    ids: np.ndarray
    names: list[str]
    resnames: list[str]
    resnrs: np.ndarray
    masses: np.ndarray


def read_topology(path: Path) -> Topology:
    """[ atoms ] of a GROMACS topology: id, atom name, residue name and number, mass."""
    rows = []
    in_atoms = False
    for line in path.read_text().splitlines():
        text = line.split(";")[0].strip()
        if not text:
            continue
        if text.startswith("["):
            in_atoms = text.strip("[] ").lower() == "atoms"
            continue
        parts = text.split()
        if in_atoms and len(parts) >= 8:
            rows.append((int(parts[0]), parts[4], parts[3], int(parts[2]), float(parts[7])))
    if not rows:
        raise ValueError(f"no [ atoms ] in {path}")
    ids, names, resnames, resnrs, masses = zip(*rows, strict=True)
    return Topology(np.array(ids), list(names), list(resnames), np.array(resnrs), np.array(masses))


def check_topology(top: Topology, elements: dict[int, str], frame: rv.Frame) -> None:
    """The topology's atom ids must be the data file's: same count, same element per atom."""
    if len(top.ids) != len(frame.ids) or not np.array_equal(top.ids, frame.ids):
        raise ValueError("topology and dump do not have the same atom ids")
    for t_mass, t_type in zip(top.masses, frame.types, strict=True):
        if rv.element_of_mass(float(t_mass)) != elements[int(t_type)]:
            raise ValueError("topology masses do not match the data file's elements")


def com_group(top: Topology, res_glob: str, name_glob: str) -> rv.Group:
    """Mass-weighted centre of the selected atoms of each residue (``gmx rdf -selrpos res_com``)."""
    keep = np.array(
        [
            fnmatch.fnmatchcase(n, name_glob) and fnmatch.fnmatchcase(r, res_glob)
            for n, r in zip(top.names, top.resnames, strict=True)
        ]
    )
    if not keep.any():
        raise ValueError(f"no atoms match residue {res_glob!r} and name {name_glob!r}")
    atom_ids = top.ids[keep]
    res = top.resnrs[keep]
    order = np.argsort(res, kind="stable")
    atom_ids, res, mass = atom_ids[order], res[order], top.masses[keep][order]
    starts = np.flatnonzero(np.r_[True, res[1:] != res[:-1]])
    sizes = np.diff(np.r_[starts, len(res)])
    residues = res[starts]

    def group(fr: rv.Frame) -> tuple[np.ndarray, np.ndarray]:
        box = fr.box / 10.0
        x = fr.xyz[np.searchsorted(fr.ids, atom_ids)] / 10.0
        d = x - np.repeat(x[starts], sizes, axis=0)
        d -= box * np.round(d / box)
        com = (
            x[starts]
            + np.add.reduceat(mass[:, None] * d, starts, axis=0)
            / np.add.reduceat(mass, starts)[:, None]
        )
        return residues, np.mod(com, box)

    return group


def first_peak(r: np.ndarray, g: np.ndarray, rmin: float) -> tuple[float, float]:
    mask = (r >= rmin) & np.isfinite(g)
    k = int(np.argmax(g[mask]))
    return float(r[mask][k]), float(g[mask][k])


def rms_difference(r: np.ndarray, g: np.ndarray, ref: np.ndarray, rmin: float) -> float:
    mask = (r >= rmin) & np.isfinite(ref)
    return float(np.sqrt(np.mean((g[mask] - ref[mask]) ** 2)))


def sem(values: list[float]) -> float:
    a = np.asarray(values)
    return float(a.std(ddof=1) / np.sqrt(len(a))) if len(a) > 1 else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--system", choices=sorted(SYSTEMS), required=True)
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--dump", type=Path, required=True)
    ap.add_argument("--top", type=Path, required=True)
    ap.add_argument("--ref-root", type=Path, required=True, help="the bundle's paper/ directory")
    ap.add_argument("--first-step", type=int, default=0, help="first step of the window")
    ap.add_argument("--blocks", type=int, default=5)
    ap.add_argument("--rmax", type=float, default=1.0, help="nm")
    ap.add_argument("--dr", type=float, default=0.002, help="nm (gmx rdf default)")
    ap.add_argument("--rmin", type=float, default=0.2, help="nm; bonded range left out")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()

    preset = SYSTEMS[args.system]
    elements, bonds = rv.read_data(args.data)
    frames = rv.read_dump(args.dump, args.first_step, None)
    if len(frames) < args.blocks:
        print(
            f"{len(frames)} frames in the window, fewer than {args.blocks} blocks", file=sys.stderr
        )
        return 1
    top = read_topology(args.top)
    check_topology(top, elements, frames[0])
    excl = rv.excluded_pairs(bonds)
    groups = {"com:" + k: com_group(top, *v) for k, v in preset["com"].items()}

    def group_of(key: str) -> rv.Group:
        return groups["com:" + key] if "com:" + key in groups else rv.element_group(elements, key)

    ref_dir = args.ref_root / preset["ref_dir"]
    lines = [
        f"system {args.system}  data {args.data.name}  frames {len(frames)} "
        f"(steps {frames[0].step}-{frames[-1].step}), {args.blocks} blocks, r >= {args.rmin} nm",
        f"{'pair':10s} {'runs':>4s} {'peak r nm':>17s} {'ref':>7s} {'peak g':>15s} {'ref':>7s} "
        f"{'ref min-max':>13s} {'rel err':>8s} {'RMS':>7s}  verdict",
    ]
    curves = {}
    n_pass = n_fail = 0
    for pair in preset["pairs"]:
        centers, counts, norm, shell = rv.rdf_per_frame(
            frames,
            group_of(pair.ref),
            group_of(pair.sel),
            args.rmax,
            args.dr,
            excl if pair.excl else None,
        )
        g = counts.sum(axis=0) / (norm.sum() * shell)
        np.savetxt(f"{args.out}_{pair.name}.dat", np.column_stack([centers, g]), header="r_nm g")
        refs = sorted(ref_dir.glob(f"rdf_s*_*_{pair.name}.xvg"))
        if not refs:
            lines.append(f"{pair.name:10s} no reference files in {ref_dir}")
            continue
        ref_g = np.array(
            [np.interp(centers, *rv.read_xvg(p), left=np.nan, right=np.nan) for p in refs]
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)  # bins outside every reference file
            mean = np.nanmean(ref_g, axis=0)
            ref_lo, ref_hi = np.nanmin(ref_g, axis=0), np.nanmax(ref_g, axis=0)
        # block statistics of the backmapped curve
        block_idx = np.array_split(np.arange(len(frames)), args.blocks)
        peaks, heights, rmss = [], [], []
        for idx in block_idx:
            gb = counts[idx].sum(axis=0) / (norm[idx].sum() * shell)
            pr, ph = first_peak(centers, gb, args.rmin)
            peaks.append(pr)
            heights.append(ph)
            rmss.append(rms_difference(centers, gb, mean, args.rmin))
        pos, ht = first_peak(centers, g, args.rmin)
        rpos, rht = first_peak(centers, mean, args.rmin)
        rmin_max = [first_peak(centers, x, args.rmin)[1] for x in ref_g]
        rel = abs(ht - rht) / rht
        rms = rms_difference(centers, g, mean, args.rmin)
        ok_pos, ok_ht, ok_rms = abs(pos - rpos) <= PEAK_POS_TOL, rel <= PEAK_HT_TOL, rms <= RMS_TOL
        verdict = (
            "PASS"
            if (ok_pos and ok_ht and ok_rms)
            else "FAIL "
            + "".join(
                f"{tag}"
                for tag, ok in (("pos ", ok_pos), ("height ", ok_ht), ("rms", ok_rms))
                if not ok
            )
        )
        n_pass += verdict == "PASS"
        n_fail += verdict != "PASS"
        lines.append(
            f"{pair.name:10s} {len(refs):4d} {pos:8.3f}+-{sem(peaks):.3f} {rpos:7.3f} "
            f"{ht:8.3f}+-{sem(heights):.3f} {rht:7.3f} {min(rmin_max):6.2f}-{max(rmin_max):<6.2f} "
            f"{rel:8.1%} {rms:7.3f}+-{sem(rmss):.3f}  {verdict}"
        )
        # block standard error of the curve
        gblocks = np.array([counts[i].sum(axis=0) / (norm[i].sum() * shell) for i in block_idx])
        curves[pair.name] = (
            centers,
            g,
            gblocks.std(axis=0, ddof=1) / np.sqrt(args.blocks),
            mean,
            ref_lo,
            ref_hi,
        )
    lines.append(f"{n_pass} pair(s) within all tolerances, {n_fail} outside")
    text = "\n".join(lines)
    Path(f"{args.out}_metrics.txt").write_text(text + "\n")
    print(text)
    if args.plot and curves:
        make_figure(curves, args.out, args.rmax)
    return 0


def make_figure(curves: dict, out: Path, rmax: float) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(curves)
    cols = min(3, n)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(6.8, 2.1 * rows), squeeze=False)
    blue, verm = "#0072B2", "#D55E00"
    for k, (name, (r, g, err, mean, lo, hi)) in enumerate(curves.items()):
        ax = axes[k // cols][k % cols]
        ax.fill_between(r, lo, hi, color=verm, alpha=0.2, lw=0)
        ax.plot(r, mean, color=verm, lw=1.2, ls="--", label="published reference")
        ax.fill_between(r, g - err, g + err, color=blue, alpha=0.3, lw=0)
        ax.plot(r, g, color=blue, lw=1.2, label="backmapped")
        pair = name.replace("_excl", "").replace("_", "\N{EN DASH}")
        ax.set_title(
            f"({'abcdefghi'[k]}) {pair}" + (" (excl.)" if name.endswith("_excl") else ""),
            fontsize=8,
            loc="left",
        )
        ax.set_xlim(0.15, rmax - 0.01)  # the last bin is half a bin wide
        ax.tick_params(labelsize=7)
        ax.set_xlabel("$r$ (nm)", fontsize=8)
        if k % cols == 0:
            ax.set_ylabel("$g(r)$", fontsize=8)
    for k in range(n, rows * cols):
        axes[k // cols][k % cols].set_axis_off()
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, fontsize=7, frameon=False)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(f"{out}_rdf.pdf")


if __name__ == "__main__":
    sys.exit(main())
