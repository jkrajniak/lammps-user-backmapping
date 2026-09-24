"""Convert LAMMPS-format sources to GROMACS files for the hybrid builder.

The hybrid builder (vendored bakery) reads GROMACS ``.gro``/``.top``. A LAMMPS
CG data file or a LAMMPS AT fragment (data file plus bounded input script) is
written out as GROMACS files in ``<work_dir>/lammps_sources/`` and the settings
are pointed at them, so every input format goes through the same builder.

Force-field parameters are converted exactly (LAMMPS ``real`` units and
conventions to GROMACS): bond and angle K -> k = 2K, RB coefficients from the
IUPAC phi convention to GROMACS psi = phi - 180 ((-1)^n), kcal -> kJ, A -> nm.
Coordinates go through ``.gro`` at 0.001 nm, the precision bakery reads; they
only set the starting geometry.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from backmap_prep import units
from backmap_prep.parsers.lammps_data_parser import (
    parse_at_fragment,
    parse_cg_system,
    parse_lammps_data,
)

if TYPE_CHECKING:
    from pathlib import Path

    from backmap_prep.parsers.gro_parser import GroFile
    from backmap_prep.parsers.top_parser import Topology
    from backmap_prep.schema import CrossAngle, CrossBond, CrossDihedral, Settings

SOURCE_DIR = "lammps_sources"
KCAL_TO_KJ = 1.0 / units.KJ_TO_KCAL
ANG_TO_NM = 1.0 / units.NM_TO_ANGSTROM


def _num(value: float) -> str:
    return f"{value:.12g}"


def _write_gro(path: Path, gro: GroFile, resname: str, names: list[str]) -> None:
    lines = [gro.title or "converted from LAMMPS", f"{len(gro.atoms):5d}"]
    for atom, name in zip(gro.atoms, names, strict=True):
        lines.append(
            f"{atom.resid % 100000:5d}{resname[:5]:<5s}{name[:5]:>5s}{atom.index % 100000:5d}"
            f"{atom.x * ANG_TO_NM:8.3f}{atom.y * ANG_TO_NM:8.3f}{atom.z * ANG_TO_NM:8.3f}"
        )
    bx, by, bz = gro.box
    lines.append(f"{bx * ANG_TO_NM:10.5f}{by * ANG_TO_NM:10.5f}{bz * ANG_TO_NM:10.5f}")
    path.write_text("\n".join(lines) + "\n")


def _at_top_text(top: Topology, resname: str) -> str:
    """GROMACS topology for a LAMMPS AT fragment, parameters converted exactly."""
    (mol,) = top.molecule_types.values()
    out = ["[ defaults ]", "; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ", "1 2 no 1.0 1.0", ""]
    out += ["[ atomtypes ]", "; name mass charge ptype sigma(nm) epsilon(kJ/mol)"]
    for name, at in top.atom_types.items():
        out.append(
            f"{name} {_num(at.mass)} 0.0 A {_num(at.sigma * ANG_TO_NM)} "
            f"{_num(at.epsilon * KCAL_TO_KJ)}"
        )
    out += ["", "[ moleculetype ]", f"{resname} 3", "", "[ atoms ]"]
    for a in mol.atoms:
        out.append(
            f"{a.index} {a.type} 1 {resname} {a.name} {a.index} {_num(a.charge)} {_num(a.mass)}"
        )
    out += ["", "[ bonds ]"]
    for b in mol.bonds:
        r0, k = b.params
        out.append(
            f"{b.i} {b.j} 1 {_num(r0 * ANG_TO_NM)} {_num(2.0 * k * KCAL_TO_KJ / ANG_TO_NM**2)}"
        )
    out += ["", "[ angles ]"]
    for ang in mol.angles:
        theta0, k = ang.params
        out.append(f"{ang.i} {ang.j} {ang.k} 1 {_num(theta0)} {_num(2.0 * k * KCAL_TO_KJ)}")
    out += ["", "[ dihedrals ]"]
    for d in mol.dihedrals:
        rb = [((-1.0) ** n) * c * KCAL_TO_KJ for n, c in enumerate(d.params[:6])]
        out.append(f"{d.i} {d.j} {d.k} {d.atom_l} 3 " + " ".join(_num(c) for c in rb))
    out += ["", "[ system ]", resname, "", "[ molecules ]", f"{resname} 1", ""]
    return "\n".join(out)


def _bead_index(ref: str, bead_names: list[str]) -> int:
    """1-based index of ``MOL:BEAD`` (or ``BEAD``) in the template bead order."""
    return bead_names.index(ref.split(":")[-1]) + 1


def _cg_bonded_lines(
    entries: list[CrossBond] | list[CrossAngle] | list[CrossDihedral],
    members: str,
    bead_names: list[str],
) -> list[str]:
    lines: list[str] = []
    for entry in entries:
        if not entry.cg_bonded:
            continue
        if not entry.params:
            raise ValueError(f"cg_bonded cross interaction {entry} has no params")
        for group in getattr(entry, members):
            idx = [_bead_index(ref, bead_names) for ref in group]
            lines.append(" ".join(str(i) for i in idx) + f" {entry.params}")
    return lines


def _cg_top_text(
    top: Topology,
    settings: Settings,
    resname: str,
    bead_names: list[str],
    data_bonds: list[tuple[int, int]],
) -> str:
    (mol,) = top.molecule_types.values()
    out = ["[ defaults ]", "1 1 no 1.0 1.0", "", "[ atomtypes ]"]
    for name, at in top.atom_types.items():
        out.append(f"{name} {_num(at.mass)} 0.0 V 1.0 1.0")
    out += ["", "[ moleculetype ]", f"{resname} 3", "", "[ atoms ]"]
    for a, bead in zip(mol.atoms, bead_names, strict=True):
        out.append(f"{a.index} {a.type} 1 {resname} {bead} {a.index} 0.0 {_num(a.mass)}")
    ci = settings.cross_interactions
    bonds = _cg_bonded_lines(ci.bonds, "pairs", bead_names)
    known = {tuple(sorted(int(t) for t in line.split()[:2])) for line in bonds}
    for i, j in data_bonds:
        if tuple(sorted((i, j))) not in known:
            raise ValueError(
                f"CG bond {i}-{j} from the LAMMPS data file has no cg_bonded entry in "
                "cross_interactions.bonds (its table and parameters are needed)"
            )
    out += ["", "[ bonds ]", *bonds]
    out += ["", "[ angles ]", *_cg_bonded_lines(ci.angles, "triples", bead_names)]
    out += ["", "[ dihedrals ]", *_cg_bonded_lines(ci.dihedrals, "quadruples", bead_names)]
    n_mol = top.molecules[0][1]
    out += ["", "[ system ]", resname, "", "[ molecules ]", f"{resname} {n_mol}", ""]
    return "\n".join(out)


def materialize_lammps_sources(settings: Settings, work_dir: Path) -> Settings:
    """Return settings whose LAMMPS-format sources point at converted GROMACS files."""
    uses_lammps = settings.cg_system is not None and settings.cg_system.format == "lammps"
    uses_lammps |= any(mol.source.format == "lammps" for mol in settings.molecules)
    if not uses_lammps:
        return settings
    if len(settings.molecules) != 1:
        raise ValueError("LAMMPS-format sources are supported for one molecule type")
    mol = settings.molecules[0]
    out_dir = work_dir / SOURCE_DIR
    out_dir.mkdir(exist_ok=True)
    new = settings.model_copy(deep=True)
    new_mol = new.molecules[0]

    if mol.source.format == "lammps":
        assert mol.source.data is not None and mol.source.input_script is not None
        gro, top = parse_at_fragment(work_dir / mol.source.data, work_dir / mol.source.input_script)
        names = [a.name for a in gro.atoms]
        _write_gro(out_dir / "at_fragment.gro", gro, mol.name, names)
        (out_dir / "at_fragment.top").write_text(_at_top_text(top, mol.name))
        new_mol.source.format = "gromacs"
        new_mol.source.coordinates = f"{SOURCE_DIR}/at_fragment.gro"
        new_mol.source.topology = f"{SOURCE_DIR}/at_fragment.top"
        new_mol.source.data = None
        new_mol.source.input_script = None

    if settings.cg_system is not None and settings.cg_system.format == "lammps":
        assert settings.cg_system.data is not None and new.cg_system is not None
        data_path = work_dir / settings.cg_system.data
        gro, top = parse_cg_system(data_path)
        bead_names = [bead.name for bead in mol.beads]
        (template,) = top.molecule_types.values()
        if len(template.atoms) != len(bead_names):
            raise ValueError(
                f"{settings.cg_system.data}: {len(template.atoms)} CG atoms per molecule, "
                f"but molecules[0].beads defines {len(bead_names)}"
            )
        n_per = len(bead_names)
        names = [bead_names[i % n_per] for i in range(len(gro.atoms))]
        raw = parse_lammps_data(data_path)
        first_ids = {a.atom_id for a in sorted(raw.atoms, key=lambda a: a.atom_id)[:n_per]}
        offset = min(first_ids) - 1
        data_bonds = [
            (b.i - offset, b.j - offset) for b in raw.bonds if b.i in first_ids and b.j in first_ids
        ]
        _write_gro(out_dir / "cg_system.gro", gro, mol.name, names)
        (out_dir / "cg_system.top").write_text(
            _cg_top_text(top, settings, mol.name, bead_names, data_bonds)
        )
        new.cg_system.format = "gromacs"
        new.cg_system.coordinates = f"{SOURCE_DIR}/cg_system.gro"
        new.cg_system.topology = f"{SOURCE_DIR}/cg_system.top"
        new.cg_system.data = None
    return new
