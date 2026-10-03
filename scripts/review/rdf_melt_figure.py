"""Paper figure: backmapped vs independent all-atom reference RDFs of a melt, r in nm.

Each ``--row`` is one system: three panels (CH3-CH3, CH3-CH2, CH2-CH2 carbon pairs) with the
block-averaged g(r) of the backmapped and the reference run and +-1 standard-error bands over the
averaging blocks (the first block is dropped as further equilibration). The RDF files are the
``fix ave/time ... mode vector`` outputs of the Tier C runs (r in Angstrom, one block per 200 ps).

    uv run --with numpy --with matplotlib rdf_melt_figure.py \\
        --row "all-atom" pe_aa_backmap/rdf_backmap.dat pe_aa_reference/rdf_reference.dat \\
        --row "united-atom" pe4_backmap/rdf_backmap.dat pe4_reference/rdf_reference.dat \\
        --out fig3_pe_rdf.pdf [--skip-first 1] [--xmax 1.4]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

PAIRS = ("CH$_3$\u2013CH$_3$", "CH$_3$\u2013CH$_2$", "CH$_2$\u2013CH$_2$")
BLUE, VERMILLION = "#0072B2", "#D55E00"
ANGSTROM_TO_NM = 0.1


def parse_blocks(path: Path) -> tuple[np.ndarray, list[np.ndarray]]:
    """Bin centres (Angstrom) and, per pair, an (n_blocks, n_bins) array of g(r)."""
    blocks: list[np.ndarray] = []
    rows: list[list[float]] = []
    n_expected = 0
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) == 2 and len(rows) == n_expected:
            if rows:
                blocks.append(np.array(rows))
            rows, n_expected = [], int(parts[1])
            continue
        rows.append([float(x) for x in parts])
    if rows and len(rows) == n_expected:
        blocks.append(np.array(rows))
    if not blocks:
        raise ValueError(f"no RDF blocks in {path}")
    n_pairs = (blocks[0].shape[1] - 2) // 2
    gr = [np.stack([b[:, 2 + 2 * p] for b in blocks]) for p in range(n_pairs)]
    return blocks[0][:, 1], gr


def mean_sem(blocks: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = blocks.shape[0]
    sem = blocks.std(axis=0, ddof=1) / np.sqrt(n) if n > 1 else np.zeros(blocks.shape[1])
    return blocks.mean(axis=0), sem


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument(
        "--row", nargs=3, action="append", metavar=("LABEL", "BACKMAP", "REFERENCE"), required=True
    )
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--skip-first", type=int, default=1)
    ap.add_argument("--xmax", type=float, default=1.4, help="nm")
    ap.add_argument("--width", type=float, default=6.8, help="inches")
    args = ap.parse_args()

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_rows = len(args.row)
    fig, axes = plt.subplots(
        n_rows, 3, figsize=(args.width, 2.1 * n_rows), squeeze=False, sharex=True
    )
    letters = "abcdefghi"
    for i, (label, bm_path, ref_path) in enumerate(args.row):
        r_bm, g_bm = parse_blocks(Path(bm_path))
        r_ref, g_ref = parse_blocks(Path(ref_path))
        for j in range(3):
            ax = axes[i, j]
            bm_mean, bm_sem = mean_sem(g_bm[j][args.skip_first :])
            ref_mean, ref_sem = mean_sem(g_ref[j][args.skip_first :])
            x_bm, x_ref = r_bm * ANGSTROM_TO_NM, r_ref * ANGSTROM_TO_NM
            ax.fill_between(x_bm, bm_mean - bm_sem, bm_mean + bm_sem, color=BLUE, alpha=0.25, lw=0)
            ax.fill_between(
                x_ref, ref_mean - ref_sem, ref_mean + ref_sem, color=VERMILLION, alpha=0.25, lw=0
            )
            ax.plot(x_bm, bm_mean, color=BLUE, lw=1.4, label="backmapped")
            ax.plot(
                x_ref, ref_mean, color=VERMILLION, lw=1.4, ls="--", label="independent reference"
            )
            ax.set_xlim(0, args.xmax)
            ax.set_title(f"({letters[3 * i + j]}) {PAIRS[j]}, {label}", fontsize=8, loc="left")
            ax.tick_params(labelsize=7)
            if j == 0:
                ax.set_ylabel("$g(r)$", fontsize=8)
            if i == n_rows - 1:
                ax.set_xlabel("$r$ (nm)", fontsize=8)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, fontsize=7, frameon=False)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(args.out)
    print(f"wrote {args.out} ({n_rows} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
