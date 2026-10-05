"""Figure: melamine network RDFs of the complete and of the reference topology against the published runs.

Reads the ``<prefix>_<pair>.dat`` curves written by network_tierc.py for two runs (the complete topology and the
topology without the angles and dihedrals of the cross-link bonds) and the published GROMACS reference RDFs, and
draws O-H, O-N, C-N and the ring centers of mass.

    uv run --with numpy --with matplotlib melamine_topology_figure.py \\
        --complete melnet/out --reftopo reftopo/out --ref-root <bundle>/paper --out fig6_melamine_network_rdf.pdf
"""

from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np

PAIRS = [
    ("O_H", "O$\\mathrm{-}$H"),
    ("O_N", "O$\\mathrm{-}$N"),
    ("C_N", "C$\\mathrm{-}$N"),
    ("ring_ring", "ring COM"),
]
RMIN = {"O_H": 0.2, "O_N": 0.2, "C_N": 0.2, "ring_ring": 0.35}
RMAX = {"O_H": 0.6, "O_N": 0.6, "C_N": 0.6, "ring_ring": 1.0}


def load(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    rows = [
        [float(x) for x in ln.split()[:2]]
        for ln in Path(path).read_text().splitlines()
        if ln.strip() and ln[0] not in "#@"
    ]
    a = np.array(rows)
    return a[:, 0], a[:, 1]


def reference(root: Path, pair: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    files = sorted(glob.glob(str(root / "mf3" / "rdf" / f"rdf_s?_*_{pair}.xvg")))
    curves = [load(f) for f in files]
    r = curves[0][0]
    g = np.array([c[1][: len(r)] for c in curves])
    return r, g.mean(axis=0), g.min(axis=0), g.max(axis=0)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument(
        "--complete", required=True, help="prefix of the network_tierc.py output, complete topology"
    )
    ap.add_argument(
        "--reftopo", required=True, help="prefix of the network_tierc.py output, reference topology"
    )
    ap.add_argument("--ref-root", type=Path, required=True, help="the bundle's paper/ directory")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    blue, verm, green = "#0072B2", "#D55E00", "#009E73"
    fig, axes = plt.subplots(2, 2, figsize=(6.8, 4.6))
    for k, (pair, title) in enumerate(PAIRS):
        ax = axes[k // 2][k % 2]
        rr, mean, lo, hi = reference(args.ref_root, pair)
        ax.fill_between(rr, lo, hi, color=verm, alpha=0.2, lw=0)
        ax.plot(rr, mean, color=verm, lw=1.2, ls="--", label="published reference")
        r1, g1 = load(f"{args.complete}_{pair}.dat")
        r2, g2 = load(f"{args.reftopo}_{pair}.dat")
        ax.plot(r1, g1, color=blue, lw=1.2, label="complete topology")
        ax.plot(r2, g2, color=green, lw=1.2, label="topology of the reference")
        lo_r, hi_r = RMIN[pair], RMAX[pair]
        ax.set_xlim(lo_r, hi_r)
        win = (r1 >= lo_r) & (r1 <= hi_r)
        ax.set_ylim(
            0, 1.12 * max(np.nanmax(g1[win]), np.nanmax(g2[: len(r2)][(r2 >= lo_r) & (r2 <= hi_r)]))
        )
        ax.set_title(f"({'abcd'[k]}) {title}", fontsize=8, loc="left")
        ax.tick_params(labelsize=7)
        ax.set_xlabel("$r$ (nm)", fontsize=8)
        if k % 2 == 0:
            ax.set_ylabel("$g(r)$", fontsize=8)
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=7, frameon=False)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(args.out)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
