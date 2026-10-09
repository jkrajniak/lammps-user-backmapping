"""Figure: strong-scaling speedup of the lambda-ramp block from the b2_scaling.sh logs.

uv run --with numpy --with matplotlib b2_figure.py --system "Dodecane" <b2-dir> --system ... --out fig_mpi_speedup.pdf
"""

from __future__ import annotations

import argparse
import statistics
from pathlib import Path

import b2_summary as b2


def speedups(d: Path) -> tuple[list[int], list[float], list[float], list[float]]:
    """Ranks, median speedup, and the speedup range from the slowest and the fastest repetition."""
    by_np: dict[int, list[float]] = {}
    for log in sorted(d.glob("log.b2_np*_r*")):
        parsed = b2.parse(log)
        if parsed is not None:
            by_np.setdefault(parsed[1], []).append(parsed[0])
    t1 = statistics.median(by_np[1])
    ranks = sorted(by_np)
    med = [t1 / statistics.median(by_np[n]) for n in ranks]
    lo = [t1 / max(by_np[n]) for n in ranks]
    hi = [t1 / min(by_np[n]) for n in ranks]
    return ranks, med, lo, hi


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--system", nargs=2, action="append", metavar=("LABEL", "DIR"), required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colours = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00"]
    markers = ["o", "s", "^", "D", "v"]
    fig, ax = plt.subplots(figsize=(3.6, 3.2))
    for k, (label, d) in enumerate(args.system):
        ranks, med, lo, hi = speedups(Path(d))
        ax.errorbar(
            ranks,
            med,
            yerr=[
                [m - low for m, low in zip(med, lo, strict=False)],
                [h - m for m, h in zip(med, hi, strict=False)],
            ],
            color=colours[k % 5],
            marker=markers[k % 5],
            ms=4,
            lw=1.0,
            capsize=2,
            label=label,
        )
    ax.plot([1, 8], [1, 8], color="#888888", lw=0.8, ls="--", label="ideal")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log", base=2)
    ax.set_xticks([1, 2, 4, 8], ["1", "2", "4", "8"])
    ax.set_yticks([1, 2, 4, 8], ["1", "2", "4", "8"])
    ax.set_xlabel("MPI ranks", fontsize=8)
    ax.set_ylabel("speedup $t_1/t_N$", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=6.5, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(args.out)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
