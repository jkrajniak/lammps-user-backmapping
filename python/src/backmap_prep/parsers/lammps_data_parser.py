"""LAMMPS data file parser for reading equilibrated CG configurations."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path


@dataclass
class LammpsDataAtom:
    atom_id: int
    mol_id: int
    type_id: int
    charge: float
    x: float
    y: float
    z: float
    ix: int = 0
    iy: int = 0
    iz: int = 0


@dataclass
class LammpsDataBond:
    bond_id: int
    type_id: int
    i: int
    j: int


@dataclass
class LammpsDataFile:
    atoms: list[LammpsDataAtom] = field(default_factory=list)
    bonds: list[LammpsDataBond] = field(default_factory=list)
    box: tuple[float, float, float] = (0.0, 0.0, 0.0)
    n_atom_types: int = 0
    n_bond_types: int = 0


def parse_lammps_data(path: Path) -> LammpsDataFile:
    """Parse a LAMMPS data file (atom_style full).

    Reads atom positions (with image flags if present), bonds, and box
    dimensions. Used to read equilibrated CG frames for the rebuild workflow.
    """
    text = path.read_text()
    lines = text.splitlines()

    result = LammpsDataFile()
    section: str | None = None
    box = [0.0, 0.0, 0.0]

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        if stripped.startswith("Atoms"):
            section = "atoms"
            continue
        if stripped.startswith("Bonds"):
            section = "bonds"
            continue
        if stripped.startswith("Velocities"):
            section = "velocities"
            continue
        if stripped.startswith("Masses"):
            section = "masses"
            continue
        if stripped.startswith("Angles"):
            section = "angles"
            continue

        if "xlo" in stripped:
            parts = stripped.split()
            box[0] = float(parts[1]) - float(parts[0])
        elif "ylo" in stripped:
            parts = stripped.split()
            box[1] = float(parts[1]) - float(parts[0])
        elif "zlo" in stripped:
            parts = stripped.split()
            box[2] = float(parts[1]) - float(parts[0])
        elif "atom types" in stripped:
            result.n_atom_types = int(stripped.split()[0])
        elif "bond types" in stripped:
            result.n_bond_types = int(stripped.split()[0])
        elif section == "atoms":
            parts = stripped.split()
            if len(parts) < 7:
                continue
            result.atoms.append(
                LammpsDataAtom(
                    atom_id=int(parts[0]),
                    mol_id=int(parts[1]),
                    type_id=int(parts[2]),
                    charge=float(parts[3]),
                    x=float(parts[4]),
                    y=float(parts[5]),
                    z=float(parts[6]),
                    ix=int(parts[7]) if len(parts) > 7 else 0,
                    iy=int(parts[8]) if len(parts) > 8 else 0,
                    iz=int(parts[9]) if len(parts) > 9 else 0,
                )
            )
        elif section == "bonds":
            parts = stripped.split()
            if len(parts) < 4:
                continue
            result.bonds.append(
                LammpsDataBond(
                    bond_id=int(parts[0]),
                    type_id=int(parts[1]),
                    i=int(parts[2]),
                    j=int(parts[3]),
                )
            )
        elif section in ("velocities", "masses", "angles"):
            pass

    result.box = (box[0], box[1], box[2])
    return result
