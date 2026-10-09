"""Collect the strong-scaling logs of b2_scaling.sh: loop time of the measured block, speedup, ghost atoms.

    uv run b2_summary.py <outdir> [<outdir> ...]

Prints, for every directory (one system), the median loop time over the repetitions for each rank count, the speedup
t_1 / t_N, the parallel efficiency, the share of the loop time spent in "Comm", and the average number of ghost
atoms per rank of the measured run. The speedup uses the median over repetitions; the spread (min and max) is shown.
"""

from __future__ import annotations

import re
import statistics
import sys
from pathlib import Path

LOOP = re.compile(r"Loop time of ([0-9.eE+-]+) on (\d+) procs for (\d+) steps with (\d+) atoms")
COMM = re.compile(
    r"^Comm\s+\|\s+[0-9.eE+-]+\s+\|\s+[0-9.eE+-]+\s+\|\s+[0-9.eE+-]+\s+\|\s+[0-9.eE+-]+\s+\|\s+([0-9.]+)"
)
GHOST = re.compile(r"^Nghost:\s+([0-9.eE+-]+) ave")


def parse(log: Path) -> tuple[float, int, float | None, float | None] | None:
    """Loop time of the last run, number of ranks, Comm share (percent) and average ghosts of that run."""
    text = log.read_text().splitlines()
    idx = [i for i, ln in enumerate(text) if LOOP.search(ln)]
    if not idx:
        return None
    start = idx[-1]
    m = LOOP.search(text[start])
    assert m is not None
    comm = ghost = None
    for ln in text[start : start + 60]:
        if (c := COMM.match(ln)) is not None:
            comm = float(c.group(1))
        if (g := GHOST.match(ln)) is not None:
            ghost = float(g.group(1))
    return float(m.group(1)), int(m.group(2)), comm, ghost


def main() -> int:
    for d in sys.argv[1:]:
        by_np: dict[int, list[tuple[float, float | None, float | None]]] = {}
        for log in sorted(Path(d).glob("log.b2_np*_r*")):
            parsed = parse(log)
            if parsed is None:
                continue
            t, np_, comm, ghost = parsed
            by_np.setdefault(np_, []).append((t, comm, ghost))
        if not by_np:
            print(f"{d}: no finished logs")
            continue
        cutoff = Path(d, "b2_comm_cutoff.txt")
        print(
            f"== {d}  ({cutoff.read_text().strip() if cutoff.exists() else 'comm cutoff not recorded'})"
        )
        t1 = statistics.median(x[0] for x in by_np[min(by_np)]) if min(by_np) == 1 else None
        print(
            "ranks  reps  loop time (s) median [min..max]   speedup  efficiency  Comm %   ghosts/rank"
        )
        for np_ in sorted(by_np):
            ts = [x[0] for x in by_np[np_]]
            tm = statistics.median(ts)
            comms = [x[1] for x in by_np[np_] if x[1] is not None]
            ghosts = [x[2] for x in by_np[np_] if x[2] is not None]
            sp = f"{t1 / tm:6.2f}" if t1 else "   n/a"
            eff = f"{100 * t1 / tm / np_:6.0f} %" if t1 else "   n/a"
            print(
                f"{np_:5d}  {len(ts):4d}  {tm:10.2f} [{min(ts):.2f}..{max(ts):.2f}]   {sp}   {eff}   "
                f"{statistics.median(comms) if comms else float('nan'):6.1f}   "
                f"{statistics.median(ghosts) if ghosts else float('nan'):10.0f}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
