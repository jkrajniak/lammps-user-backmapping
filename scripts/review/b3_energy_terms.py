"""Review item B3 (R1-1): production-window energy terms, backmapped vs independent reference.

Reads the thermo output of two LAMMPS logs, takes the last ``run`` of each (the production
window, thermo every 5000 steps) and prints the per-atom energy terms, temperature, pressure and
density with 5-block standard errors, then the differences in units of the combined error.
The same script serves the Tier C logs (each system at its own density) and the matched-density
control written by ``b3_matched_density.sh`` (both at one density).

    python b3_energy_terms.py <log_backmap> <log_reference> <n_atoms> [--from-step N]

A log argument may be a glob pattern or a comma-separated list: the logs of the segments of a resumed
run (log.<name>.s1.lammps, ...). With --from-step N all rows with step > N from all runs of all
segments are used (one row per step); without it the last run of a single log.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np

TERMS = ("Temp", "PotEng", "E_bond", "E_angle", "E_dihed", "E_vdwl", "E_coul", "Press", "Density")
N_BLOCKS = 5


def all_runs(log: Path) -> tuple[list[str], list[np.ndarray]]:
    """Thermo header and the rows of every ``run`` in a LAMMPS log."""
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
    return header, [np.array(r) for r in runs if r]


def last_run(log: Path) -> tuple[list[str], np.ndarray]:
    """Thermo header and rows of the last ``run`` in a LAMMPS log."""
    header, runs = all_runs(log)
    if not runs:
        raise ValueError(f"no thermo rows in {log}")
    return header, runs[-1]


def rows_after(logs: list[Path], from_step: int) -> tuple[list[str], np.ndarray]:
    """Rows with step > from_step from all runs of all logs (the segments of a resumed run), one per
    step: a step that occurs in several segments is taken from the latest log."""
    merged: dict[int, np.ndarray] = {}
    header: list[str] = []
    for log in logs:
        head, runs = all_runs(log)
        header = header or head
        for run in runs:
            for row in run:
                if row[0] > from_step:
                    merged[int(row[0])] = row
    if not merged:
        raise ValueError(f"no thermo rows after step {from_step} in {[str(p) for p in logs]}")
    return header, np.array([merged[s] for s in sorted(merged)])


def block_mean_sem(x: np.ndarray) -> tuple[float, float]:
    blocks = np.array([b.mean() for b in np.array_split(x, N_BLOCKS)])
    return float(x.mean()), float(blocks.std(ddof=1) / np.sqrt(N_BLOCKS))


def expand(spec: str) -> list[Path]:
    """A log argument: one file, a comma-separated list, or a glob pattern; segments in numeric order."""
    paths: list[Path] = []
    for part in spec.split(","):
        matches = sorted(
            Path().glob(part), key=lambda p: [int(x) for x in re.findall(r"\d+", p.name)]
        )
        paths += matches or [Path(part)]
    return paths


def summarize(spec: str, n_atoms: float, from_step: int | None) -> dict[str, tuple[float, float]]:
    logs = expand(spec)
    if from_step is None:
        if len(logs) != 1:
            raise SystemExit("several logs (segments of a resumed run) need --from-step")
        header, rows = last_run(logs[0])
    else:
        header, rows = rows_after(logs, from_step)
    out = {}
    for term in TERMS:
        if term not in header:
            continue
        mean, sem = block_mean_sem(rows[:, header.index(term)])
        scale = 1.0 / n_atoms if term.startswith("E_") or term == "PotEng" else 1.0
        out[term] = (mean * scale, sem * scale)
    print(f"{spec}: {len(rows)} rows, steps {int(rows[0, 0])}..{int(rows[-1, 0])}")
    return out


def main() -> int:
    args = sys.argv[1:]
    from_step = None
    if "--from-step" in args:
        k = args.index("--from-step")
        from_step = int(args[k + 1])
        del args[k : k + 2]
    if len(args) != 3:
        print(__doc__)
        return 1
    n_atoms = float(args[2])
    backmap = summarize(args[0], n_atoms, from_step)
    reference = summarize(args[1], n_atoms, from_step)
    print(f"\n{'term':9s} {'backmap':>20s} {'reference':>20s} {'difference':>22s}")
    for term, (a, ea) in backmap.items():
        b, eb = reference[term]
        err = float(np.hypot(ea, eb))
        tail = f"({(a - b) / err:+.1f} sigma)" if err > 1e-9 else "(fixed)"
        print(f"{term:9s} {a:12.5f} +-{ea:.5f} {b:12.5f} +-{eb:.5f} {a - b:+11.5f} {tail}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
