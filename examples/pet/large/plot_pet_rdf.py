# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "matplotlib>=3.8",
#     "numpy>=1.26",
# ]
# ///

"""Plot the PET/Dacron Tier C RDF comparison (manuscript Figure 5).

Reads the backmapped curves produced by `pet_rdf_bakery.py` (archived under
`rdf_out_bakery/`) and the published Dacron reference `.xvg` curves from
paper-reverse-mapping-polymer-networks, and writes the five-panel overlay.

This script did not previously exist: the committed `fig5_pet_rdf.png` was
produced ad hoc and never archived, so the figure could not be regenerated
(vector output, restyling) without redoing the whole analysis. It is written
here to close that gap.

Reference curves are averaged over every seed and rate available for a pair
(4 seeds x 3 rates), matching `report_dacron_ref_peaks.py`.

Line widths, font sizes and vector output follow the conventions applied to
the dodecane (PR #7), melamine_network (PR #11) and rim135 (PR #18) figures.

Usage:
    uv run examples/pet/large/plot_pet_rdf.py \\
      --rdf-dir examples/pet/large/rdf_out_bakery \\
      --ref-dir ../paper-reverse-mapping-polymer-networks/paper/dacron/rdf \\
      --out examples/pet/large/rdf_comparison_bakery.pdf
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

DASH = "\u2013"  # en dash, matching the panel titles of the original figure

# (pair key, panel title, x limit in nm) -- the panel set and x ranges match
# the originally promoted figure.
PANELS: list[tuple[str, str, float]] = [
    ("C_O", f"C{DASH}O", 1.0),
    ("C_O_excl", f"C{DASH}O (excl. bonded)", 2.0),
    ("C_H", f"C{DASH}H", 1.0),
    ("O_H", f"O{DASH}H", 1.0),
    ("ring_ring", f"ring{DASH}ring (COM)", 2.0),
]


def parse_xvg(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read a GROMACS .xvg; returns r in nm and g(r)."""
    r, g = [], []
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith(("#", "@")):
            continue
        parts = s.split()
        if len(parts) >= 2:
            try:
                r.append(float(parts[0]))
                g.append(float(parts[1]))
            except ValueError:
                continue
    return np.array(r), np.array(g)


def load_reference(ref_dir: Path, pair: str) -> tuple[np.ndarray, np.ndarray] | None:
    """Mean reference curve over all seeds/rates for `pair`, r in nm."""
    curves = sorted(ref_dir.glob(f"rdf_s*_*_{pair}.xvg"))
    if not curves:
        return None
    r_ref: np.ndarray | None = None
    gs: list[np.ndarray] = []
    for c in curves:
        r, g = parse_xvg(c)
        if r_ref is None:
            r_ref = r
        if len(g) == len(r_ref):
            gs.append(g)
    if r_ref is None or not gs:
        return None
    return r_ref, np.mean(gs, axis=0)


def load_backmap(rdf_dir: Path, pair: str) -> tuple[np.ndarray, np.ndarray] | None:
    """Backmapped curve written by pet_rdf_bakery.py; r stored in A, returned in nm."""
    path = rdf_dir / f"rdf_{pair}_backmap_avg.txt"
    if not path.exists():
        return None
    data = np.loadtxt(path)
    return data[:, 0] / 10.0, data[:, 1]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rdf-dir", type=Path, required=True)
    parser.add_argument("--ref-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    for path in (args.rdf_dir, args.ref_dir):
        if not path.exists():
            print(f"error: missing input {path}", file=sys.stderr)
            return 2

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(18, 9))
    flat = axes.flatten()

    missing: list[str] = []
    for ax, (pair, title, x_max) in zip(flat, PANELS, strict=False):
        bm = load_backmap(args.rdf_dir, pair)
        ref = load_reference(args.ref_dir, pair)
        if bm is None or ref is None:
            missing.append(pair)
            ax.set_visible(False)
            continue
        ax.plot(bm[0], bm[1], color="C0", lw=2.4, label="backmapped")
        ax.plot(ref[0], ref[1], "--", color="C1", lw=2.4, label="reference")
        ax.set_xlim(0, x_max)
        ax.set_xlabel("r (nm)", fontsize=14)
        ax.set_ylabel("g(r)", fontsize=14)
        ax.set_title(title, fontsize=15)
        ax.tick_params(axis="both", labelsize=12)
        ax.legend(fontsize=12)

    # The panel set is odd-numbered; blank the unused sixth slot.
    for ax in flat[len(PANELS) :]:
        ax.set_visible(False)

    if missing:
        print(f"error: no data for {', '.join(missing)}", file=sys.stderr)
        return 1

    fig.tight_layout()
    # A .pdf path yields vector output; dpi only affects raster formats.
    fig.savefig(args.out)
    plt.close(fig)
    print(f"Wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
