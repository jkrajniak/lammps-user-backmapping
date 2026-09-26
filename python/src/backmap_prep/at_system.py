"""Atomistic-only systems and force field derived from a hybrid system.

Tier C runs the backmapped melt as a plain AT system, next to an independent
AT reference built from the same molecule. Both must carry exactly the AT
force field of the hybrid, so this module derives it from the hybrid
``System`` (never restated by hand) with one fixed type numbering:

- ``<prefix>.at.ff.lmp``: plain LAMMPS styles (the backmap styles at full
  weight), 1-4 interactions through ``special_bonds`` (fudgeLJ, fudgeQQ).
- AT-only data files in that numbering, either extracted from a hybrid frame
  (CG beads removed) or built as an independent reference: copies of the AT
  template molecule, randomly rotated on a lattice in an expanded box.
"""

from __future__ import annotations

import math
import random
import re
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

from . import units
from .parsers import parse_gro

if TYPE_CHECKING:
    from collections.abc import Sequence
    from pathlib import Path

    from .builder import System
    from .schema import Settings

__all__ = [
    "AtTypeMaps",
    "at_type_maps",
    "write_at_forcefield",
    "write_at_from_hybrid_frame",
    "write_at_reference",
]

_BOND_ANGLE_STYLE = {"harmonic": "harmonic", "backmap/harmonic": "harmonic"}
_ANGLE_STYLE = {**_BOND_ANGLE_STYLE, "backmap/charmm": "charmm"}
_IMPROPER_STYLE = {"backmap/harmonic": "harmonic"}
_DIHEDRAL_STYLE = {
    "ryckaert": "ryckaert",
    "backmap/ryckaert": "ryckaert",
    "harmonic": "harmonic",
    "backmap/harmonic": "harmonic",
    "charmm": "charmm",
    "backmap/charmm": "charmm",
    "backmap/fourier": "fourier",
}
_DATA_SECTIONS = ("Masses", "Atoms", "Velocities", "Bonds", "Angles", "Dihedrals", "Impropers")


@dataclass
class _AtType:
    style: str
    params: tuple[float, ...]


@dataclass
class AtTypeMaps:
    """Hybrid type ID -> AT-only type ID, and the AT-only type definitions."""

    atom: dict[int, int] = field(default_factory=dict)
    bond: dict[int, int] = field(default_factory=dict)
    angle: dict[int, int] = field(default_factory=dict)
    dihedral: dict[int, int] = field(default_factory=dict)
    improper: dict[int, int] = field(default_factory=dict)
    bond_types: list[_AtType] = field(default_factory=list)
    angle_types: list[_AtType] = field(default_factory=list)
    dihedral_types: list[_AtType] = field(default_factory=list)
    improper_types: list[_AtType] = field(default_factory=list)


def _map_types(
    hybrid_types: Sequence[Any], styles: dict[str, str], kind: str
) -> tuple[dict[int, int], list[_AtType]]:
    mapping: dict[int, int] = {}
    unique: dict[tuple[str, tuple[float, ...]], int] = {}
    out: list[_AtType] = []
    for t in hybrid_types:
        if t.keyword == "cg":
            continue
        if t.style not in styles:
            raise ValueError(f"no AT-only equivalent for {kind} style {t.style!r}")
        key = (styles[t.style], tuple(round(p, 12) for p in t.params))
        if key not in unique:
            out.append(_AtType(style=key[0], params=tuple(t.params)))
            unique[key] = len(out)
        mapping[t.type_id] = unique[key]
    return mapping, out


def at_type_maps(system: System) -> AtTypeMaps:
    """Fixed AT-only numbering: AT atom types in hybrid order, bonded types by content."""
    maps = AtTypeMaps()
    for at in system.atom_types:
        if not at.is_cg:
            maps.atom[at.type_id] = len(maps.atom) + 1
    maps.bond, maps.bond_types = _map_types(system.bond_types, _BOND_ANGLE_STYLE, "bond")
    maps.angle, maps.angle_types = _map_types(system.angle_types, _ANGLE_STYLE, "angle")
    maps.dihedral, maps.dihedral_types = _map_types(
        system.dihedral_types, _DIHEDRAL_STYLE, "dihedral"
    )
    maps.improper, maps.improper_types = _map_types(
        system.improper_types, _IMPROPER_STYLE, "improper"
    )
    return maps


def _style_block(kind: str, types: list[_AtType]) -> list[str]:
    if not types:
        return []
    styles = sorted({t.style for t in types})
    hybrid = len(styles) > 1
    lines = [f"{kind}_style {'hybrid ' if hybrid else ''}{' '.join(styles)}"]
    for n, t in enumerate(types, start=1):
        p = list(t.params)
        if kind == "angle" and t.style == "charmm":
            values = " ".join(f"{v:.10g}" for v in p[:4])
        elif kind == "angle":
            values = f"{p[0]:.10g} {p[1]:.10g}"
        elif t.style == "fourier":
            m = int(p[0])
            values = f"{m} " + " ".join(
                f"{p[1 + 3 * i]:.10g} {int(p[2 + 3 * i])} {p[3 + 3 * i]:.10g}" for i in range(m)
            )
        elif t.style == "harmonic" and kind == "dihedral":
            values = f"{p[0]:.10g} {int(p[1])} {int(p[2])}"
        elif t.style == "charmm":
            # weight 0: 1-4 pairs come from special_bonds, not from the dihedral
            values = f"{p[0]:.10g} {int(p[1])} {round(p[2])} 0.0"
        else:
            values = " ".join(f"{v:.10g}" for v in p)
        lines.append(f"{kind}_coeff {n} {t.style + ' ' if hybrid else ''}{values}")
    return [*lines, ""]


def write_at_forcefield(
    system: System, settings: Settings, maps: AtTypeMaps, path: Path, data_name: str
) -> None:
    """Write the AT-only force field (``<prefix>.at.ff.lmp``)."""
    sim = settings.simulation
    lj_cut = units.distance(sim.lj_cutoff)
    coul_cut = units.distance(sim.coulomb_cutoff)
    lines = [
        f"# AT-only force field for {data_name} -- generated by backmap-prep, do not edit.",
        "# Same force field as the hybrid at lambda = 1; 1-4 terms via special_bonds.",
        "",
        f"pair_style lj/cut/coul/cut {lj_cut:.10g} {coul_cut:.10g}",
    ]
    for pt in system.pair_types:
        if pt.kind == "atomistic":
            i, j = sorted((maps.atom[pt.itype], maps.atom[pt.jtype]))
            lines.append(f"pair_coeff {i} {j} {pt.epsilon:.10g} {pt.sigma:.10g}")
    lines.append("")
    lines += _style_block("bond", maps.bond_types)
    lines += _style_block("angle", maps.angle_types)
    lines += _style_block("dihedral", maps.dihedral_types)
    lines += _style_block("improper", maps.improper_types)
    if system.has_cross_pairs:
        lj14, qq14 = system.fudge_lj, system.fudge_qq
    else:
        lj14 = qq14 = 0.0
    lines += [
        f"special_bonds lj 0.0 0.0 {lj14:.10g} coul 0.0 0.0 {qq14:.10g}",
        "",
        "neigh_modify delay 0 every 1 check yes",
        "",
    ]
    path.write_text("\n".join(lines))


def _read_sections(path: Path) -> tuple[list[tuple[float, float]], dict[str, list[list[str]]]]:
    bounds: list[tuple[float, float]] = []
    sections: dict[str, list[list[str]]] = {name: [] for name in _DATA_SECTIONS}
    current: str | None = None
    for line in path.read_text().splitlines():
        body = line.split("#")[0].strip()
        if not body:
            continue
        if body in _DATA_SECTIONS or re.match(r"^[A-Z][A-Za-z ]+$", body):
            current = body if body in _DATA_SECTIONS else None
            continue
        parts = body.split()
        if len(parts) >= 4 and parts[2] in ("xlo", "ylo", "zlo"):
            bounds.append((float(parts[0]), float(parts[1])))
        elif current:
            sections[current].append(parts)
    return bounds, sections


def _write_data(
    path: Path,
    title: str,
    bounds: list[tuple[float, float]],
    system: System,
    maps: AtTypeMaps,
    atoms: list[str],
    velocities: list[str],
    terms: dict[str, list[str]],
) -> None:
    masses = {maps.atom[at.type_id]: (at.mass, at.name) for at in system.atom_types if not at.is_cg}
    counts = [
        f"{len(atoms)} atoms",
        f"{len(terms['Bonds'])} bonds",
        f"{len(terms['Angles'])} angles",
        f"{len(terms['Dihedrals'])} dihedrals",
        f"{len(terms.get('Impropers', []))} impropers",
        "",
        f"{len(masses)} atom types",
        f"{len(maps.bond_types)} bond types",
        f"{len(maps.angle_types)} angle types",
        f"{len(maps.dihedral_types)} dihedral types",
        f"{len(maps.improper_types)} improper types",
        "",
    ]
    out = [title, "", *counts]
    for (lo, hi), axis in zip(bounds, "xyz", strict=True):
        out.append(f"{lo:.10g} {hi:.10g} {axis}lo {axis}hi")
    out += ["", "Masses", ""]
    out += [f"{t} {m:.10g} # {name}" for t, (m, name) in sorted(masses.items())]
    out += ["", "Atoms # full", "", *atoms]
    if velocities:
        out += ["", "Velocities", "", *velocities]
    for name in ("Bonds", "Angles", "Dihedrals", "Impropers"):
        if terms.get(name):
            out += ["", name, "", *terms[name]]
    path.write_text("\n".join(out) + "\n")


def write_at_from_hybrid_frame(system: System, maps: AtTypeMaps, frame: Path, out: Path) -> int:
    """Write the AT atoms of a hybrid frame (e.g. after the ramp) as an AT-only system.

    The frame must come from the hybrid built with the same settings (same
    type numbering). Returns the number of AT atoms.
    """
    bounds, sec = _read_sections(frame)
    new_id: dict[int, int] = {}
    atoms: list[str] = []
    for parts in sorted(sec["Atoms"], key=lambda p: int(p[0])):
        old_type = int(parts[2])
        if old_type not in maps.atom:
            continue
        new_id[int(parts[0])] = len(new_id) + 1
        flags = " ".join(parts[7:10]) if len(parts) >= 10 else "0 0 0"
        atoms.append(
            f"{new_id[int(parts[0])]} {parts[1]} {maps.atom[old_type]} {parts[3]} "
            f"{parts[4]} {parts[5]} {parts[6]} {flags}"
        )
    velocities = [
        f"{new_id[int(p[0])]} {p[1]} {p[2]} {p[3]}"
        for p in sec["Velocities"]
        if int(p[0]) in new_id
    ]
    type_map = {
        "Bonds": maps.bond,
        "Angles": maps.angle,
        "Dihedrals": maps.dihedral,
        "Impropers": maps.improper,
    }
    terms: dict[str, list[str]] = {}
    for name, tmap in type_map.items():
        rows: list[str] = []
        for p in sec[name]:
            ids = [int(v) for v in p[2:]]
            if all(i in new_id for i in ids) and int(p[1]) in tmap:
                rows.append(
                    f"{len(rows) + 1} {tmap[int(p[1])]} " + " ".join(str(new_id[i]) for i in ids)
                )
        terms[name] = rows
    _write_data(
        out, f"AT-only frame from {frame.name}", bounds, system, maps, atoms, velocities, terms
    )
    return len(atoms)


def _rotation(rng: random.Random) -> list[list[float]]:
    q = [rng.gauss(0.0, 1.0) for _ in range(4)]
    norm = math.sqrt(sum(v * v for v in q))
    w, x, y, z = (v / norm for v in q)
    return [
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)],
    ]


def write_at_reference(
    system: System,
    maps: AtTypeMaps,
    hybrid_gro: Path,
    template_gro: Path,
    n_mol: int,
    seed: int,
    out: Path,
) -> float:
    """Independent AT reference: ``n_mol`` template molecules on a lattice.

    Atom types, charges and bonded terms come from the first molecule of the
    hybrid (matched to the template by atom name); coordinates from the AT
    template. The box is expanded (dilute); compress it in the protocol.
    Returns the box length (A).
    """
    names = [(int(ln[0:5]), ln[10:15].strip()) for ln in hybrid_gro.read_text().splitlines()[2:-1]]
    first_res = names[0][0]
    by_id = {a.atom_id: a for a in system.atoms}
    mol_ids = [i + 1 for i, (res, _) in enumerate(names) if res == first_res]
    at_ids = [i for i in mol_ids if not by_id[i].is_cg]
    template = parse_gro(template_gro)
    tpl_index = {a.name: k for k, a in enumerate(template.atoms)}
    order = sorted(at_ids, key=lambda i: tpl_index[names[i - 1][1]])
    local = {aid: k + 1 for k, aid in enumerate(order)}
    coords = [
        [units.distance(c) for c in (a.x, a.y, a.z)]
        for a in sorted(template.atoms, key=lambda a: tpl_index[a.name])
    ]
    centre = [sum(c[k] for c in coords) / len(coords) for k in range(3)]
    coords = [[c[k] - centre[k] for k in range(3)] for c in coords]
    extent = 2.0 * max(math.sqrt(sum(v * v for v in c)) for c in coords)
    spacing = extent * 1.25
    n_side = math.ceil(n_mol ** (1.0 / 3.0))
    box = n_side * spacing
    rng = random.Random(seed)

    atoms: list[str] = []
    per_mol = len(order)
    for m in range(n_mol):
        ix, iy, iz = m % n_side, (m // n_side) % n_side, m // (n_side * n_side)
        origin = [(ix + 0.5) * spacing, (iy + 0.5) * spacing, (iz + 0.5) * spacing]
        rot = _rotation(rng)
        for k, aid in enumerate(order):
            c = coords[k]
            x = [origin[r] + sum(rot[r][s] * c[s] for s in range(3)) for r in range(3)]
            a = by_id[aid]
            atoms.append(
                f"{m * per_mol + k + 1} {m + 1} {maps.atom[a.type_id]} {a.charge:.10g} "
                f"{x[0]:.6f} {x[1]:.6f} {x[2]:.6f} 0 0 0"
            )
    template_terms = {
        "Bonds": [(maps.bond.get(b.type_id), [b.i, b.j]) for b in system.bonds],
        "Angles": [(maps.angle.get(t.type_id), [t.i, t.j, t.k]) for t in system.angles],
        "Dihedrals": [
            (maps.dihedral.get(d.type_id), [d.i, d.j, d.k, d.l]) for d in system.dihedrals
        ],
        "Impropers": [
            (maps.improper.get(d.type_id), [d.i, d.j, d.k, d.l]) for d in system.impropers
        ],
    }
    terms: dict[str, list[str]] = {}
    for name, rows in template_terms.items():
        mol_rows = [
            (t, [local[i] for i in ids]) for t, ids in rows if t and all(i in local for i in ids)
        ]
        out_rows: list[str] = []
        for m in range(n_mol):
            for t, ids in mol_rows:
                out_rows.append(
                    f"{len(out_rows) + 1} {t} " + " ".join(str(m * per_mol + i) for i in ids)
                )
        terms[name] = out_rows
    _write_data(
        out,
        f"Independent AT reference: {n_mol} molecules, seed {seed}",
        [(0.0, box)] * 3,
        system,
        maps,
        atoms,
        [],
        terms,
    )
    return float(box)
