"""Fit the NVE total-energy drift of the B4 runs (log.b4_<case>).

For each case: linear fit of etotal (kcal/mol) against time over the NVE run
(the last thermo block), reported per AT atom per ns, and the RMS deviation
of etotal from the fit per AT atom. In the hybrid runs the AT atom count is
the size of group at_atoms (CG beads are slaved to their fragment COM).
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np


def last_thermo_block(text: str) -> tuple[list[str], np.ndarray]:
    blocks = re.findall(r"^\s+(Step .*?)\n(.*?)^Loop time", text, re.MULTILINE | re.DOTALL)
    if not blocks:
        raise ValueError("no thermo block")
    header, body = blocks[-1]
    cols = header.split()
    rows = [ln.split() for ln in body.splitlines() if ln.strip() and ln.split()[0].isdigit()]
    return cols, np.array(rows, dtype=float)


def n_at_atoms(text: str) -> int:
    match = re.search(r"^\s*(\d+) atoms in group at_atoms", text, re.MULTILINE)
    if match:
        return int(match.group(1))
    return int(re.search(r"^\s*(\d+) atoms$", text, re.MULTILINE).group(1))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("workdir", type=Path)
    args = parser.parse_args()
    print(
        f"{'case':22s} {'N_AT':>6s} {'ps':>6s} {'drift kcal/mol/atom/ns':>23s} {'rms kcal/mol/atom':>18s} {'<T> K':>7s}"
    )
    for log in sorted(args.workdir.glob("log.b4_*")):
        text = log.read_text()
        cols, data = last_thermo_block(text)
        c = {name: i for i, name in enumerate(cols)}
        t_fs = data[:, c["Time"]] - data[0, c["Time"]]
        e = data[:, c["TotEng"]]
        n = n_at_atoms(text)
        slope, icpt = np.polyfit(t_fs, e, 1)
        rms = float(np.sqrt(np.mean((e - (slope * t_fs + icpt)) ** 2)))
        print(
            f"{log.name.removeprefix('log.b4_'):22s} {n:6d} {t_fs[-1] / 1000:6.1f} "
            f"{slope * 1e6 / n:23.3e} {rms / n:18.3e} {data[:, c['c_tgrp']].mean():7.1f}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
