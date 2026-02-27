"""GROMACS .gro file parser."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class GroAtom:
    resid: int
    resname: str
    name: str
    index: int
    x: float  # nm
    y: float  # nm
    z: float  # nm
    vx: float = 0.0
    vy: float = 0.0
    vz: float = 0.0


@dataclass
class GroFile:
    title: str
    atoms: list[GroAtom]
    box: tuple[float, float, float]  # nm


def parse_gro(path: Path) -> GroFile:
    """Parse a GROMACS .gro coordinate file."""
    lines = Path(path).read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"Invalid .gro file: {path}")

    title = lines[0].strip()
    natoms = int(lines[1].strip())

    atoms: list[GroAtom] = []
    for i in range(2, 2 + natoms):
        line = lines[i]
        resid = int(line[0:5])
        resname = line[5:10].strip()
        name = line[10:15].strip()
        index = int(line[15:20])
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        vx = vy = vz = 0.0
        if len(line) >= 68:
            vx = float(line[44:52])
            vy = float(line[52:60])
            vz = float(line[60:68])
        atoms.append(GroAtom(resid, resname, name, index, x, y, z, vx, vy, vz))

    box_line = lines[2 + natoms].split()
    box = (float(box_line[0]), float(box_line[1]), float(box_line[2]))

    return GroFile(title=title, atoms=atoms, box=box)
