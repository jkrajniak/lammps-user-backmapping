"""Backmapped LAMMPS frame -> GROMACS .gro in the order of gromacs/topol_at.top.

The hybrid lists each residue as bead, fragment atoms, bead, ... (POPC atoms
grouped by bead; one W residue = one bead + 4 TIP3P waters named OWk/HWk1/HWk2).
The AT continuation (Slipids, PME) runs in GROMACS with the settings of the
independent AT reference, so the frame is reordered to POPC.itp order and the
waters are split into TIP3P molecules. Velocities are kept (A/fs -> nm/ps).

    uv run python to_gromacs.py popc_hybrid.data at_backmapped.gro
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _itp_order(path: Path) -> list[str]:
    names, sec = [], None
    for line in path.read_text().splitlines():
        m = re.match(r"\s*\[\s*(\w+)", line)
        if m:
            sec = m.group(1)
            continue
        t = line.split(";")[0].split()
        if sec == "atoms" and len(t) >= 5 and t[0].isdigit():
            names.append(t[4])
    return names


def _read_data(path: Path) -> tuple[list[float], dict[int, list[float]], dict[int, list[float]]]:
    box, xyz, vel, sec = [], {}, {}, None
    for line in path.read_text().splitlines():
        body = line.split("#")[0].strip()
        if not body:
            continue
        if body in ("Atoms", "Velocities") or line.startswith(("Atoms", "Velocities")):
            sec = body.split()[0]
            continue
        if re.match(r"^[A-Z][A-Za-z ]+$", body):
            sec = None
            continue
        t = body.split()
        if len(t) >= 4 and t[2] in ("xlo", "ylo", "zlo"):
            box.append(float(t[1]) - float(t[0]))
        elif sec == "Atoms":
            xyz[int(t[0])] = [float(v) for v in t[4:7]]
        elif sec == "Velocities":
            vel[int(t[0])] = [float(v) for v in t[1:4]]
    return box, xyz, vel


def main(frame: Path, out: Path) -> None:
    popc_order = _itp_order(HERE / "POPC.itp")
    gro = (HERE / "hyb_conf.gro").read_text().splitlines()
    n = int(gro[1])
    residues: dict[int, tuple[str, dict[str, int]]] = {}
    for atom_id, line in enumerate(gro[2 : 2 + n], start=1):
        resid, resname, name = int(line[0:5]), line[5:10].strip(), line[10:15].strip()
        residues.setdefault(resid, (resname, {}))[1][name] = atom_id
    box, xyz, vel = _read_data(frame)
    rows: list[tuple[str, str, int]] = []  # (resname, atom name, hybrid id)
    waters: list[tuple[str, str, int]] = []
    for resid in sorted(residues):
        resname, by_name = residues[resid]
        if resname == "POPC":
            rows += [("POPC", a, by_name[a]) for a in popc_order]
        elif resname == "W":
            for k in range(1, 5):
                waters += [
                    ("SOL", "OW", by_name[f"OW{k}"]),
                    ("SOL", "HW1", by_name[f"HW{k}1"]),
                    ("SOL", "HW2", by_name[f"HW{k}2"]),
                ]
        else:
            raise SystemExit(f"unexpected residue {resname}")
    rows += waters
    lines = [f"Backmapped MARTINI POPC from {frame.name}", f"{len(rows)}"]
    resnr = count = 0
    per_res = {"POPC": len(popc_order), "SOL": 3}
    for k, (resname, name, aid) in enumerate(rows):
        if k == 0 or rows[k - 1][0] != resname or count == per_res[resname]:
            resnr, count = resnr + 1, 0
        count += 1
        x = [c / 10.0 for c in xyz[aid]]
        v = [c * 100.0 for c in vel.get(aid, [0.0, 0.0, 0.0])]
        lines.append(
            f"{resnr % 100000:5d}{resname:<5s}{name:>5s}{(k + 1) % 100000:5d}"
            f"{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}{v[0]:8.4f}{v[1]:8.4f}{v[2]:8.4f}"
        )
    lines.append(" ".join(f"{b / 10.0:10.5f}" for b in box))
    out.write_text("\n".join(lines) + "\n")
    print(
        f"{len(rows)} atoms ({len(rows) - len(waters)} lipid, {len(waters) // 3} waters) -> {out}"
    )


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
