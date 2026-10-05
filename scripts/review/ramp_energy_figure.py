"""Figure: temperature of the atomistic sub-system and energy terms through the staged hybrid run of two systems.

Reads the CSV written by b3_ramp_energy.py for each system:

    uv run --with numpy --with matplotlib ramp_energy_figure.py \\
        --system "Dodecane melt" b3_dodecane.csv 6000 --system "Epoxy network" b3_rim135.csv 12953 --out fig_ramp_energy.pdf

The energies are given per atomistic site (the kinetic and potential energy of the hybrid system in the LAMMPS thermo output).
The shaded interval is the lambda ramp (0 < lambda < 1); the horizontal axis ends after the first 8 ps.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

TERMS = [
    ("E_bond", "bond", "#0072B2"),
    ("E_angle", "angle", "#E69F00"),
    ("E_dihed", "dihedral", "#009E73"),
    ("E_vdwl", "LJ", "#CC79A7"),
    ("PotEng", "total", "#000000"),
]


def read(path: Path) -> dict[str, np.ndarray]:
    with path.open() as fh:
        rows = [r for r in csv.DictReader(fh) if r["stage"] != "minimize"]
    cols = rows[0].keys()
    return {c: np.array([float(r[c]) if c != "stage" else 0.0 for r in rows]) for c in cols}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument(
        "--system", nargs=3, action="append", metavar=("TITLE", "CSV", "NATOMS"), required=True
    )
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--tmax", type=float, default=8.0, help="end of the time axis in ps")
    args = ap.parse_args()

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(args.system)
    fig, axes = plt.subplots(2, n, figsize=(3.4 * n, 4.6), sharex=True, squeeze=False)
    for k, (title, csv_path, natoms) in enumerate(args.system):
        d = read(Path(csv_path))
        t = d["time_fs"] / 1000.0
        lam = d["lambda"]
        ramp = (lam > 0) & (lam < 1)
        t0, t1 = t[ramp].min(), t[ramp].max()
        ax_t, ax_e = axes[0][k], axes[1][k]
        for ax in (ax_t, ax_e):
            ax.axvspan(t0, t1, color="#999999", alpha=0.18, lw=0)
            ax.set_xlim(0, args.tmax)
            ax.tick_params(labelsize=7)
        ax_t.plot(t, d["Temp"], color="#D55E00", lw=1.0)
        ax_t.set_yscale("log")
        ax_t.set_title(f"({'ab'[k]}) {title}", fontsize=8, loc="left")
        for key, label, colour in TERMS:
            ax_e.plot(t, d[key] / float(natoms), color=colour, lw=1.0, label=label)
        ax_e.set_xlabel("time (ps)", fontsize=8)
        if k == 0:
            ax_t.set_ylabel("$T$ of the atomistic sites (K)", fontsize=8)
            ax_e.set_ylabel("energy per atom (kcal mol$^{-1}$)", fontsize=8)
            ax_e.legend(fontsize=6.5, frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(args.out)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
