"""Review item B8: write an equilibrated CG LAMMPS data frame back into cg_conf.gro.

Keeps the .gro file's atom order, names and residue fields and replaces the
coordinates (wrapped into the box, nm) and the box by those of the data file.
The k-th atom of the data file in ID order is the k-th atom of the .gro file
(a CG frame cut out of a hybrid data file keeps the molecule and bead order
but not contiguous IDs). As a check, every .gro atom name must map to a single
data-file atom type. Velocities are written as zero; the backmapping protocol
assigns its own.

    python b8_data_to_gro.py cg_equil.data cg_conf.gro
"""

from __future__ import annotations

import sys
from pathlib import Path

ANGSTROM_TO_NM = 0.1


def read_data(
    path: Path,
) -> tuple[list[float], dict[int, tuple[float, float, float]], dict[int, int]]:
    lines = path.read_text().splitlines()
    lo_hi = {}
    for ln in lines:
        parts = ln.split()
        if len(parts) >= 4 and parts[2] in ("xlo", "ylo", "zlo"):
            lo_hi[parts[2][0]] = (float(parts[0]), float(parts[1]))
    start = next(i for i, ln in enumerate(lines) if ln.startswith("Atoms"))
    coords = {}
    types = {}
    for ln in lines[start + 2 :]:
        if not ln.strip():
            break
        p = ln.split()
        coords[int(p[0])] = (float(p[4]), float(p[5]), float(p[6]))
        types[int(p[0])] = int(p[2])
    lows = [lo_hi[a][0] for a in "xyz"]
    lengths = [lo_hi[a][1] - lo_hi[a][0] for a in "xyz"]
    wrapped = {
        i: tuple((c - lo) % length for c, lo, length in zip(x, lows, lengths, strict=True))
        for i, x in coords.items()
    }
    return lengths, wrapped, types


def main() -> int:
    data, gro = Path(sys.argv[1]), Path(sys.argv[2])
    lengths, coords, types = read_data(data)
    lines = gro.read_text().splitlines()
    n = int(lines[1])
    if n != len(coords):
        print(f"atom count differs: gro {n}, data {len(coords)}", file=sys.stderr)
        return 1
    ids = sorted(coords)
    name_type: dict[str, set[int]] = {}
    for atom_id, ln in zip(ids, lines[2 : 2 + n], strict=True):
        name_type.setdefault(ln[10:15].strip(), set()).add(types[atom_id])
    mixed = {name: t for name, t in name_type.items() if len(t) > 1}
    if mixed:
        print(f"order mismatch: gro names with several data types {mixed}", file=sys.stderr)
        return 1
    out = lines[:2]
    for atom_id, ln in zip(ids, lines[2 : 2 + n], strict=True):
        x, y, z = (c * ANGSTROM_TO_NM for c in coords[atom_id])
        out.append(f"{ln[:20]}{x:8.3f}{y:8.3f}{z:8.3f}{0.0:8.4f}{0.0:8.4f}{0.0:8.4f}")
    out.append("".join(f"{length * ANGSTROM_TO_NM:10.5f}" for length in lengths))
    gro.write_text("\n".join(out) + "\n")
    print(f"{gro}: {n} atoms from {data}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
