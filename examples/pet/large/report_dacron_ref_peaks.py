# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy"]
# ///
"""Report first-peak position/height of the published Dacron RDF references.

Reads the 60 RDF files in paper/dacron/rdf/ (4 seeds x 3 rates x 5 pairs),
averages g(r) across the 12 curves per pair, and reports the first peak
above r_min. These are the Tier C comparison targets for the backmapped
PET/Dacron system.

Pairs (GROMACS gmx rdf selections):
  ring_ring : resname TER && name C[2-7], residue-COM vs residue-COM
  C_O       : TER/DIO carbon - oxygen
  C_O_excl  : C-O excluding bonded (1-2) pairs
  C_H       : carbon - hydrogen
  O_H       : oxygen - hydrogen
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

PAIRS = ["ring_ring", "C_O", "C_O_excl", "C_H", "O_H"]
R_MIN = 1.5  # Å, skip the exclusion zone


def parse_xvg(path: Path) -> tuple[np.ndarray, np.ndarray]:
    r, g = [], []
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith(("#", "@")):
            continue
        parts = s.split()
        if len(parts) >= 2:
            try:
                r.append(float(parts[0]) * 10.0)  # nm -> Å
                g.append(float(parts[1]))
            except ValueError:
                continue
    return np.array(r), np.array(g)


def first_peak(r: np.ndarray, g: np.ndarray, r_min: float = R_MIN) -> tuple[float, float]:
    mask = r >= r_min
    if not mask.any():
        return 0.0, 0.0
    rc, gc = r[mask], g[mask]
    idx = int(np.argmax(gc))
    return float(rc[idx]), float(gc[idx])


def main(rdf_dir: Path) -> int:
    print("Published Dacron RDF reference peaks (avg of 4 seeds x 3 rates)")
    print(f"  source: {rdf_dir}")
    print(f"  {'pair':<12} {'peak r (Å)':<12} {'peak g(r)':<10} {'n curves'}")
    for pair in PAIRS:
        curves = sorted(rdf_dir.glob(f"rdf_s*_*_{pair}.xvg"))
        if not curves:
            print(f"  {pair:<12} (no files)")
            continue
        gs = []
        r_ref = None
        for c in curves:
            r, g = parse_xvg(c)
            if r_ref is None:
                r_ref = r
            if len(g) == len(r_ref):
                gs.append(g)
        g_mean = np.mean(gs, axis=0) if gs else r_ref
        pr, pg = first_peak(r_ref, g_mean)
        print(f"  {pair:<12} {pr:<12.3f} {pg:<10.3f} {len(gs)}")
    return 0


if __name__ == "__main__":
    rdf_dir = (
        Path(sys.argv[1])
        if len(sys.argv) > 1
        else Path(__file__).resolve().parents[4]
        / "paper-reverse-mapping-polymer-networks/paper/dacron/rdf"
    )
    sys.exit(main(rdf_dir))
