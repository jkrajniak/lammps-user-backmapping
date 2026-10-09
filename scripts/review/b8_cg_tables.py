"""Review item B8 (R1-3): degraded CG nonbonded tables for the dodecane melt.

Replaces the IBI pair tables (table_A_A/A_B/B_B.xvg) of an example copy by one
of two cruder models, keeping the r grid and the 7-column GROMACS layout
(r, f, -f', g, -g', h, -h'; the potential in h, as in the IBI files):

- ``lj``: 12-6 LJ fitted the way a user would without IBI: sigma = the zero
  crossing of the IBI potential, epsilon = the depth of its first minimum;
  potential-shifted to zero at the CG cutoff.
- ``wca``: the repulsive part of that LJ only (cut and shifted at 2^(1/6) sigma).

Below 0.8 sigma the potential is continued linearly with the slope at 0.8 sigma,
so the table stays finite at r = 0 like the IBI tables. The fitted parameters
are printed and written to ``b8_tables.json``.

    uv run --with numpy b8_cg_tables.py --dir <example copy> --variant lj|wca [--rc 1.4]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

PAIRS = ("A_A", "A_B", "B_B")


def ibi_fit(r: np.ndarray, u: np.ndarray) -> tuple[float, float]:
    """(sigma, epsilon) from the zero crossing and first minimum of an IBI potential."""
    inner = r > 0.25
    i_min = int(np.argmin(np.where(inner, u, np.inf)))
    cross = np.where((u[:-1] > 0) & (u[1:] <= 0) & (r[:-1] < r[i_min]))[0]
    if len(cross) == 0:
        raise ValueError("no zero crossing before the first minimum")
    k = int(cross[-1])
    sigma = r[k] + (r[k + 1] - r[k]) * u[k] / (u[k] - u[k + 1])
    return float(sigma), float(-u[i_min])


def model(r: np.ndarray, sigma: float, eps: float, variant: str, rc: float) -> np.ndarray:
    def lj(x: np.ndarray) -> np.ndarray:
        s6 = (sigma / x) ** 6
        return 4.0 * eps * (s6 * s6 - s6)

    def dlj(x: float) -> float:
        s6 = (sigma / x) ** 6
        return -24.0 * eps * (2.0 * s6 * s6 - s6) / x

    cut = 2.0 ** (1.0 / 6.0) * sigma if variant == "wca" else rc
    shift = float(lj(np.array([cut]))[0])
    r_core = 0.8 * sigma
    u = np.zeros_like(r)
    body = (r >= r_core) & (r < cut)
    u[body] = lj(r[body]) - shift
    core = r < r_core
    u_core = float(lj(np.array([r_core]))[0]) - shift
    u[core] = u_core + dlj(r_core) * (r[core] - r_core)
    return u


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--variant", choices=("lj", "wca"), required=True)
    ap.add_argument("--rc", type=float, default=1.4, help="CG cutoff, nm")
    args = ap.parse_args()

    params = {}
    for pair in PAIRS:
        path = args.dir / f"table_{pair}.xvg"
        data = np.loadtxt(path, comments=("#", "@"))
        r = data[:, 0]
        sigma, eps = ibi_fit(r, data[:, 5])
        u = model(r, sigma, eps, args.variant, args.rc)
        f = -np.gradient(u, r)
        out = np.zeros((len(r), 7))
        out[:, 0] = r
        out[:, 5] = u
        out[:, 6] = f
        header = (
            f"B8 {args.variant} replacement of the IBI table {path.name}: "
            f"sigma {sigma:.5f} nm, epsilon {eps:.5f} kJ/mol, "
            + (
                "cut and shifted at 2^(1/6) sigma"
                if args.variant == "wca"
                else f"shifted at {args.rc} nm"
            )
        )
        np.savetxt(path, out, fmt="%.10e", header=header, comments="# ")
        params[pair] = {"sigma_nm": sigma, "epsilon_kJmol": eps}
        print(f"{pair}: sigma {sigma:.4f} nm, epsilon {eps:.4f} kJ/mol ({args.variant})")
    (args.dir / "b8_tables.json").write_text(
        json.dumps({"variant": args.variant, "rc_nm": args.rc, "pairs": params}, indent=2) + "\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
