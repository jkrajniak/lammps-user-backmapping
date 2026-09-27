"""Every backmap style's forces are minus the gradient of its own energy.

Central finite differences of the LAMMPS potential energy against the forces
LAMMPS reports, for each atom and coordinate of a small, distorted molecule.
The earlier checks compared energies (with GROMACS) but not forces, which let
two defects through: backmap/table angles and dihedrals used the table force
per degree where the force expression needs it per radian (angles 57x too
weak; dihedrals in addition missed the 1/sin(phi) of the dE/d(cos phi)
formulation, about 1e-4 of the correct force).

CG styles act on beads (cg_type) and are checked through fix backmap's
per-bead CG force (peratom full, columns 5-7, before redistribution); AT
styles act on AT atoms that belong to no bead, at lambda = 1 (weight 1).

Skipped unless ``BACKMAP_LMP`` points to a LAMMPS binary with the backmap package.
"""

from __future__ import annotations

import math
import os
import re
import subprocess
from pathlib import Path

import numpy as np
import pytest

LMP_ENV = "BACKMAP_LMP"
H = 1.0e-4  # A
# Pair cases: atoms at non-bonded distances (3.4-4.6 A), where the LJ tables are
# smooth enough for a finite-difference check of a linearly interpolated table.
XYZ_PAIR = [
    (10.0, 10.0, 10.0),
    (13.4, 11.2, 10.6),
    (11.2, 13.8, 11.1),
    (14.1, 14.6, 12.9),
    (10.9, 11.6, 13.8),
]
XYZ = [
    (10.0, 10.0, 10.0),
    (11.3, 10.4, 10.2),
    (12.0, 11.6, 10.9),
    (13.4, 11.9, 10.3),
    (12.2, 12.9, 11.8),
]


def _write_tables(d: Path) -> None:
    """Analytic tables in LAMMPS format (force = -dE/dx, angles per degree).

    Fine enough that linear interpolation of E and F separately stays within the
    tolerance of the finite-difference check.
    """
    with (d / "bond.table").open("w") as fh:
        n = 20001
        fh.write(f"B\nN {n}\n\n")
        for i in range(n):
            r = 0.5 + 3.0 * i / (n - 1)
            e = 50.0 * (r - 1.5) ** 2 + 3.0 * (r - 1.5) ** 3
            f = -(100.0 * (r - 1.5) + 9.0 * (r - 1.5) ** 2)
            fh.write(f"{i + 1} {r:.10f} {e:.12g} {f:.12g}\n")
    with (d / "angle.table").open("w") as fh:
        n = 1801
        fh.write(f"A\nN {n}\n\n")
        for i in range(n):
            t = 180.0 * i / (n - 1)
            e = 0.002 * (t - 110.0) ** 2 + 0.3 * math.cos(math.radians(2 * t))
            f = -(0.004 * (t - 110.0) - 0.6 * math.sin(math.radians(2 * t)) * math.pi / 180.0)
            fh.write(f"{i + 1} {t:.8f} {e:.12g} {f:.12g}\n")
    with (d / "dihedral.table").open("w") as fh:
        n = 3601
        fh.write(f"D\nN {n}\n\n")
        for i in range(n):
            p = -180.0 + 360.0 * i / (n - 1)
            pr = math.radians(p)
            e = 1.2 * (1 + math.cos(pr)) + 0.4 * math.sin(2 * pr) + 0.3 * math.cos(3 * pr - 0.4)
            dedp_rad = -1.2 * math.sin(pr) + 0.8 * math.cos(2 * pr) - 0.9 * math.sin(3 * pr - 0.4)
            f = -dedp_rad * math.pi / 180.0
            fh.write(f"{i + 1} {p:.8f} {e:.12g} {f:.12g}\n")
    with (d / "pair.table").open("w") as fh:
        n = 20000
        fh.write(f"P\nN {n} R 1.0 12.0\n\n")
        for i in range(n):
            r = 1.0 + 11.0 * i / (n - 1)
            e = 4 * 0.5 * ((3.2 / r) ** 12 - (3.2 / r) ** 6)
            f = 4 * 0.5 * (12 * 3.2**12 / r**13 - 6 * 3.2**6 / r**7)
            fh.write(f"{i + 1} {r:.10f} {e:.12g} {f:.12g}\n")


def _data(n_atoms: int, kinds: dict[str, list[tuple[int, ...]]], charges: bool, pair: bool) -> str:
    lines = ["fd test", "", f"{n_atoms} atoms", "2 atom types"]
    for name, terms in kinds.items():
        lines += [f"{len(terms)} {name}", f"1 {name[:-1]} types"]
    lines += [
        "",
        "0.0 40.0 xlo xhi",
        "0.0 40.0 ylo yhi",
        "0.0 40.0 zlo zhi",
        "",
        "Masses",
        "",
        "1 12.0",
        "2 12.0",
        "",
        "Atoms # full",
        "",
    ]
    xyz = XYZ_PAIR if pair else XYZ
    for k in range(n_atoms):
        q = (0.3 if k % 2 else -0.3) if charges else 0.0
        x, y, z = xyz[k]
        lines.append(f"{k + 1} {k + 1} 1 {q} {x} {y} {z}")
    for name, terms in kinds.items():
        lines += ["", name.capitalize(), ""]
        lines += [f"{i + 1} 1 " + " ".join(map(str, t)) for i, t in enumerate(terms)]
    return "\n".join(lines) + "\n"


CASES = {
    # name: (n_atoms, topology, style lines, cg?)
    "bond_table_cg": (
        2,
        {"bonds": [(1, 2)]},
        "bond_style backmap/table linear 20001\nbond_coeff 1 cg bond.table B",
        True,
    ),
    "angle_table_cg": (
        3,
        {"angles": [(1, 2, 3)]},
        "angle_style backmap/table linear 1801\nangle_coeff 1 cg angle.table A",
        True,
    ),
    "dihedral_table_cg": (
        4,
        {"dihedrals": [(1, 2, 3, 4)]},
        "dihedral_style backmap/table linear 3601\ndihedral_coeff 1 cg dihedral.table D",
        True,
    ),
    "pair_table_cg": (2, {}, "PAIR_CG", True),
    "bond_harmonic_at": (
        2,
        {"bonds": [(1, 2)]},
        "bond_style backmap/harmonic\nbond_coeff 1 at 300.0 1.2",
        False,
    ),
    "angle_harmonic_at": (
        3,
        {"angles": [(1, 2, 3)]},
        "angle_style backmap/harmonic\nangle_coeff 1 at 60.0 100.0",
        False,
    ),
    "angle_charmm_at": (
        3,
        {"angles": [(1, 2, 3)]},
        "angle_style backmap/charmm\nangle_coeff 1 at 60.0 100.0 20.0 2.1",
        False,
    ),
    "dihedral_ryckaert_at": (
        4,
        {"dihedrals": [(1, 2, 3, 4)]},
        "dihedral_style backmap/ryckaert\ndihedral_coeff 1 at 1.0 -0.8 0.5 0.7 -0.2 0.1",
        False,
    ),
    "dihedral_harmonic_at": (
        4,
        {"dihedrals": [(1, 2, 3, 4)]},
        "dihedral_style backmap/harmonic\ndihedral_coeff 1 at 1.5 -1 2",
        False,
    ),
    "dihedral_fourier_at": (
        4,
        {"dihedrals": [(1, 2, 3, 4)]},
        "dihedral_style backmap/fourier\ndihedral_coeff 1 at 2 1.2 1 30.0 0.5 3 0.0",
        False,
    ),
    "improper_harmonic_at": (
        4,
        {"impropers": [(1, 2, 3, 4)]},
        "improper_style backmap/harmonic\nimproper_coeff 1 at 40.0 5.0",
        False,
    ),
    "pair_ljcoul_at": (5, {}, "PAIR_AT", False),
}


def _input(case: str, displaced: tuple[int, int, float] | None) -> str:
    _, _, styles, cg = CASES[case]
    if styles == "PAIR_CG":
        styles = (
            "pair_style backmap 12.0 zero 12.0 table linear 20000\n"
            "pair_coeff 1 1 cg pair.table P\npair_coeff 1 2 none\npair_coeff 2 2 none"
        )
    elif styles == "PAIR_AT":
        styles = (
            "pair_style backmap 12.0 lj/cut/coul/cut 12.0 12.0 zero 12.0\n"
            "pair_coeff 1 1 atomistic 0.3 3.0\npair_coeff 1 2 none\npair_coeff 2 2 none"
        )
    else:
        styles = "pair_style zero 12.0\npair_coeff * *\n" + styles
    # CG cases: the atoms are beads (type 1), lambda 0 (CG weight 1).
    # AT cases: the atoms belong to no bead (the CG type 2 has no atoms),
    # lambda 1 (AT weight 1).
    lam, cg_type = ("0.0", "1") if cg else ("1.0", "2")
    move = ""
    if displaced is not None:
        atom, dim, h = displaced
        vec = ["0", "0", "0"]
        vec[dim] = repr(h)
        move = f"group d id {atom}\ndisplace_atoms d move {' '.join(vec)} units box\n"
    bm = f"fix bm all backmap cg_type {cg_type} alpha 0.0001 lambda0 {lam} peratom full"
    return f"""units real
atom_style full
boundary p p p
atom_modify map array
read_data d.data
{move}{styles}
special_bonds lj 0.0 0.0 0.0 coul 0.0 0.0 0.0
{bm}
dump d all custom 1 f.dump id fx fy fz{" f_bm[5] f_bm[6] f_bm[7]" if cg else ""}
dump_modify d sort id format float %.15g
thermo_style custom step pe
run 0
print "RESULT $(pe:%.15g)"
"""


def _run(lmp: str, d: Path, text: str) -> tuple[float, np.ndarray]:
    (d / "in.t").write_text(text)
    proc = subprocess.run(
        [lmp, "-in", "in.t", "-log", "log.t", "-screen", "none"],
        cwd=d,
        capture_output=True,
        text=True,
        check=False,
    )
    log = (d / "log.t").read_text() if (d / "log.t").exists() else ""
    assert proc.returncode == 0, (log + proc.stderr)[-2000:]
    energy = float(re.search(r"^RESULT (\S+)", log, re.MULTILINE).group(1))
    lines = (d / "f.dump").read_text().splitlines()
    head = next(i for i, ln in enumerate(lines) if ln.startswith("ITEM: ATOMS"))
    return energy, np.loadtxt(lines[head + 1 :], ndmin=2)


@pytest.mark.integration
@pytest.mark.parametrize("case", sorted(CASES))
def test_force_is_minus_energy_gradient(tmp_path: Path, case: str) -> None:
    lmp = os.environ.get(LMP_ENV)
    if not lmp or not Path(lmp).is_file():
        pytest.skip(f"set {LMP_ENV} to a LAMMPS binary built with the backmap package")
    n_atoms, topo, _, cg = CASES[case]
    _write_tables(tmp_path)
    (tmp_path / "d.data").write_text(
        _data(n_atoms, topo, charges=case == "pair_ljcoul_at", pair=case.startswith("pair"))
    )
    _, rows = _run(lmp, tmp_path, _input(case, None))
    forces = rows[:, 4:7] if cg else rows[:, 1:4]
    scale = max(np.abs(forces).max(), 1e-3)
    for atom in range(1, n_atoms + 1):
        for dim in range(3):
            e_plus, _ = _run(lmp, tmp_path, _input(case, (atom, dim, H)))
            e_minus, _ = _run(lmp, tmp_path, _input(case, (atom, dim, -H)))
            fd = -(e_plus - e_minus) / (2 * H)
            got = forces[atom - 1, dim]
            assert got == pytest.approx(fd, abs=2e-3 * scale), (case, atom, "xyz"[dim], got, fd)
