"""Review item B6 (R1-5): bond-angle and dihedral distributions, backmapped vs reference.

Reads the ``fix ave/histo ... mode vector kind local`` files of the Tier C
productions (one block per Nfreq, e.g. 200 ps), normalizes each block to a
probability distribution and reports, per distribution:

- Jensen-Shannon divergence (base 2, 0 = identical, 1 = disjoint) of the
  block-averaged distributions, and its spread from pairing blocks;
- the same for two halves of the reference run (the sampling-noise floor);
- for dihedrals, the trans fraction (|phi| > 120 deg) with block SEM;
- mean value with block SEM.

    uv run --with numpy --with matplotlib b6_bonded_distributions.py \\
        --backmap angle_hist_backmap.dat --reference angle_hist_reference.dat \\
        --kind angle --out b6_pe_angle [--skip-first 1] [--plot]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np


def parse_blocks(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Bin centres and counts per block (n_blocks x n_bins)."""
    lines = [ln for ln in path.read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    centres: np.ndarray | None = None
    blocks: list[np.ndarray] = []
    i = 0
    while i < len(lines):
        head = lines[i].split()
        nbins = int(head[1])
        rows = np.array([ln.split() for ln in lines[i + 1 : i + 1 + nbins]], dtype=float)
        if centres is None:
            centres = rows[:, 1]
        blocks.append(rows[:, 2])
        i += 1 + nbins
    if centres is None:
        raise ValueError(f"no blocks in {path}")
    return centres, np.array(blocks)


def normalize(counts: np.ndarray) -> np.ndarray:
    total = counts.sum(axis=-1, keepdims=True)
    return counts / np.where(total > 0, total, 1.0)


def js_divergence(p: np.ndarray, q: np.ndarray) -> float:
    m = 0.5 * (p + q)

    def kl(a: np.ndarray, b: np.ndarray) -> float:
        mask = a > 0
        return float(np.sum(a[mask] * np.log2(a[mask] / b[mask])))

    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


def mean_sem(values: np.ndarray) -> tuple[float, float]:
    return float(values.mean()), float(values.std(ddof=1) / np.sqrt(len(values))) if len(
        values
    ) > 1 else 0.0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--backmap", type=Path, required=True)
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--kind", choices=["angle", "dihedral"], required=True)
    ap.add_argument("--skip-first", type=int, default=1, help="blocks dropped as equilibration")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()

    x_b, cb = parse_blocks(args.backmap)
    x_r, cr = parse_blocks(args.reference)
    if not np.allclose(x_b, x_r):
        print("bin centres differ between the two files", file=sys.stderr)
        return 1
    cb, cr = cb[args.skip_first :], cr[args.skip_first :]
    pb, pr = normalize(cb), normalize(cr)
    p_mean, r_mean = normalize(cb.sum(axis=0)), normalize(cr.sum(axis=0))

    js = js_divergence(p_mean, r_mean)
    n = min(len(pb), len(pr))
    js_blocks = np.array([js_divergence(pb[k], pr[k]) for k in range(n)])
    half = len(pr) // 2
    js_ref_halves = (
        js_divergence(normalize(cr[:half].sum(axis=0)), normalize(cr[half:].sum(axis=0)))
        if half
        else float("nan")
    )

    report = [
        f"{args.kind}: backmap {args.backmap.name} ({len(pb)} blocks), reference {args.reference.name} "
        f"({len(pr)} blocks), first {args.skip_first} dropped",
        f"JS(backmap, reference) block-averaged      {js:.3e}",
        f"JS per matched block pair, mean +- SEM     {js_blocks.mean():.3e} +- "
        f"{js_blocks.std(ddof=1) / np.sqrt(n) if n > 1 else 0.0:.1e}",
        f"JS(reference 1st half, 2nd half) (noise)   {js_ref_halves:.3e}",
    ]
    for name, p in (("backmap", pb), ("reference", pr)):
        means = (p * x_b).sum(axis=1)
        m, s = mean_sem(means)
        line = f"{name:10s} mean {m:8.3f} +- {s:.3f} deg"
        if args.kind == "dihedral":
            trans = p[:, np.abs(x_b) > 120.0].sum(axis=1)
            tm, ts = mean_sem(trans)
            line += f"   trans fraction {tm:.4f} +- {ts:.4f}"
        report.append(line)
    text = "\n".join(report)
    Path(f"{args.out}.txt").write_text(text + "\n")
    np.savetxt(
        f"{args.out}.dat",
        np.column_stack([x_b, p_mean, r_mean]),
        header="x_deg p_backmap p_reference",
    )
    print(text)

    if args.plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        width = x_b[1] - x_b[0]
        fig, ax = plt.subplots(figsize=(3.4, 2.6))
        ax.plot(x_b, r_mean / width, "--", color="C1", lw=1.5, label="reference")
        ax.plot(x_b, p_mean / width, color="C0", lw=1.2, label="backmapped")
        ax.set_xlabel(f"{'bond angle' if args.kind == 'angle' else 'dihedral'} / deg")
        ax.set_ylabel("probability density / deg$^{-1}$")
        ax.legend(frameon=False, fontsize=7)
        fig.tight_layout()
        fig.savefig(f"{args.out}.pdf")
    return 0


if __name__ == "__main__":
    sys.exit(main())
