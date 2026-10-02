"""Review item B3 (R1-1): production-window energy terms, backmapped vs independent reference.

Reads the thermo output of two LAMMPS logs, takes the last ``run`` of each (the production
window, thermo every 5000 steps) and prints the per-atom energy terms, temperature, pressure and
density with 5-block standard errors, then the differences in units of the combined error.
The same script serves the Tier C logs (each system at its own density) and the matched-density
control written by ``b3_matched_density.sh`` (both at one density).

    python b3_energy_terms.py <log_backmap> <log_reference> <n_atoms>
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np

TERMS = ("Temp", "PotEng", "E_bond", "E_angle", "E_dihed", "E_vdwl", "E_coul", "Press", "Density")
N_BLOCKS = 5


def last_run(log: Path) -> tuple[list[str], np.ndarray]:
    """Thermo header and rows of the last ``run`` in a LAMMPS log."""
    header: list[str] = []
    runs: list[list[list[float]]] = []
    current: list[list[float]] | None = None
    for line in log.read_text().splitlines():
        words = line.split()
        if words and words[0] == "Step":
            header = words
            current = []
            runs.append(current)
        elif current is not None and len(words) == len(header) and re.fullmatch(r"-?\d+", words[0]):
            try:
                current.append([float(w) for w in words])
            except ValueError:
                current = None
        elif words and words[0] in ("Loop", "WARNING"):
            current = None
    if not runs or not runs[-1]:
        raise ValueError(f"no thermo rows in {log}")
    return header, np.array(runs[-1])


def block_mean_sem(x: np.ndarray) -> tuple[float, float]:
    blocks = np.array([b.mean() for b in np.array_split(x, N_BLOCKS)])
    return float(x.mean()), float(blocks.std(ddof=1) / np.sqrt(N_BLOCKS))


def summarize(log: Path, n_atoms: float) -> dict[str, tuple[float, float]]:
    header, rows = last_run(log)
    out = {}
    for term in TERMS:
        if term not in header:
            continue
        mean, sem = block_mean_sem(rows[:, header.index(term)])
        scale = 1.0 / n_atoms if term.startswith("E_") or term == "PotEng" else 1.0
        out[term] = (mean * scale, sem * scale)
    print(f"{log}: {len(rows)} rows, steps {int(rows[0, 0])}..{int(rows[-1, 0])}")
    return out


def main() -> int:
    if len(sys.argv) != 4:
        print(__doc__)
        return 1
    n_atoms = float(sys.argv[3])
    backmap = summarize(Path(sys.argv[1]), n_atoms)
    reference = summarize(Path(sys.argv[2]), n_atoms)
    print(f"\n{'term':9s} {'backmap':>20s} {'reference':>20s} {'difference':>22s}")
    for term, (a, ea) in backmap.items():
        b, eb = reference[term]
        err = float(np.hypot(ea, eb))
        tail = f"({(a - b) / err:+.1f} sigma)" if err > 1e-9 else "(fixed)"
        print(f"{term:9s} {a:12.5f} +-{ea:.5f} {b:12.5f} +-{eb:.5f} {a - b:+11.5f} {tail}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
