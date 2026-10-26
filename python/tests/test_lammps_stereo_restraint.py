"""LAMMPS regression test for the stereocenter restraint (``simulation.stereo_restraint``).

The restraints written by ``backmap_prep.stereo`` (``fix backmap/stereo``) must have zero energy
at the geometry they keep and a higher energy at the mirror image and when the hydrogen crosses
to the side of the heavy neighbors. The forces must not depend on how the atoms are distributed
over MPI ranks: ``fix restrain dihedral`` loses the forces it adds to ghost atoms.

Skipped unless ``BACKMAP_LMP`` points to a LAMMPS binary with the backmap package, and for the
parallel run ``mpirun``.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from backmap_prep.builder import AtomTypeInfo, LammpsAtom, LammpsBond, System
from backmap_prep.stereo import IDEAL_IMPROPER_DEG, write_stereo_files

LMP_ENV = "BACKMAP_LMP"
K_KJ = 83.68
K_KCAL = K_KJ / 4.184
TET = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float) / np.sqrt(3)
BOND = 1.09


def _lmp() -> Path:
    path = os.environ.get(LMP_ENV)
    if not path or not Path(path).is_file():
        pytest.skip(f"set {LMP_ENV} to a LAMMPS binary")
    return Path(path)


def _center(mirror: bool, crossed: bool = False) -> System:
    """A carbon with hydrogen, F, Cl and Br on the corners of an ideal tetrahedron.

    ``crossed`` moves the hydrogen to the side of the three heavy neighbors.
    """
    directions = TET * np.array([-1.0, 1.0, 1.0]) if mirror else TET.copy()
    if crossed:
        directions[0] *= -0.6
    system = System(box=(30.0, 30.0, 30.0))
    system.atom_types = [
        AtomTypeInfo(1, "C", 12.0, False),
        AtomTypeInfo(2, "H", 1.008, False),
        AtomTypeInfo(3, "F", 19.0, False),
        AtomTypeInfo(4, "Cl", 35.5, False),
        AtomTypeInfo(5, "Br", 79.9, False),
    ]
    system.atoms.append(LammpsAtom(1, 1, 1, 0.0, 15.0, 15.0, 15.0))
    for n, d in enumerate(directions, 2):
        system.atoms.append(LammpsAtom(n, 1, n, 0.0, *(15.0 + BOND * d)))
        system.bonds.append(LammpsBond(n - 1, 1, 1, n))
    return system


def _data(system: System) -> str:
    lines = ["stereo restraint test", "", f"{len(system.atoms)} atoms", "5 atom types", ""]
    lines += ["0.0 30.0 xlo xhi", "0.0 30.0 ylo yhi", "0.0 30.0 zlo zhi", "", "Masses", ""]
    lines += [f"{t.type_id} {t.mass}" for t in system.atom_types]
    lines += ["", "Atoms # atomic", ""]
    lines += [f"{a.atom_id} {a.type_id} {a.x:.10f} {a.y:.10f} {a.z:.10f}" for a in system.atoms]
    return "\n".join(lines) + "\n"


def _run(d: Path, restrained: System, evaluated: System, ranks: int = 1) -> str:
    d.mkdir()
    assert write_stereo_files(d, restrained, K_KJ) == 1
    (d / "test.data").write_text(_data(evaluated))
    (d / "in.test").write_text(
        f"units real\natom_style atomic\natom_modify map array\nboundary p p p\n"
        f"processors {ranks} 1 1\nread_data test.data\n"
        "pair_style zero 5.0\npair_coeff * *\ninclude stereo.lmp\n"
        "dump f all custom 1 forces.txt id fx fy fz\ndump_modify f sort id format float %.12e\n"
        "run 0\n"
        'print "RESULT $(f_stereo[1]:%.12e) $(f_stereo[2]:%.12e) $(pe:%.12e)"\n'
    )
    cmd = [str(_lmp()), "-in", "in.test", "-log", "log.test", "-screen", "none"]
    if ranks > 1:
        mpirun = shutil.which("mpirun")
        if not mpirun:
            pytest.skip("mpirun not found")
        cmd = [mpirun, "-np", str(ranks), *cmd]
    proc = subprocess.run(cmd, cwd=d, capture_output=True, text=True, check=False)
    log = (d / "log.test").read_text() if (d / "log.test").exists() else ""
    assert proc.returncode == 0, (log + proc.stdout + proc.stderr)[-3000:]
    return log


def _energies(d: Path, restrained: System, evaluated: System) -> tuple[float, float]:
    log = _run(d, restrained, evaluated)
    match = re.search(r"^RESULT (\S+) (\S+) (\S+)$", log, re.MULTILINE)
    assert match, log[-2000:]
    # fix_modify energy yes: the restraint energy is part of the potential energy
    assert float(match.group(3)) == pytest.approx(float(match.group(1)) + float(match.group(2)))
    return float(match.group(1)), float(match.group(2))


@pytest.mark.integration
@pytest.mark.parametrize("mirror", [False, True])
def test_restraint_is_zero_at_the_kept_geometry_and_raises_the_mirror_image(
    tmp_path: Path, mirror: bool
) -> None:
    kept, inverted = _center(mirror), _center(not mirror)
    assert _energies(tmp_path / "kept", kept, kept) == pytest.approx((0.0, 0.0), abs=1e-8)
    # the mirror image is at -phi: E = K[1 - cos(2 * 35.26 deg)] = 2K/3 for each
    expected = K_KCAL * (1.0 - np.cos(np.radians(2.0 * IDEAL_IMPROPER_DEG)))
    assert _energies(tmp_path / "inverted", kept, inverted) == pytest.approx(
        (expected, expected), rel=1e-4
    )


@pytest.mark.integration
@pytest.mark.parametrize("mirror", [False, True])
def test_hydrogen_crossing_to_the_heavy_neighbors_raises_the_second_restraint(
    tmp_path: Path, mirror: bool
) -> None:
    heavy, second = _energies(tmp_path / "crossed", _center(mirror), _center(mirror, crossed=True))
    assert heavy == pytest.approx(0.0, abs=1e-8)
    assert second > K_KCAL


@pytest.mark.integration
def test_forces_do_not_depend_on_the_mpi_decomposition(tmp_path: Path) -> None:
    """The center sits on the boundary of a two-rank split, so its restraints span both ranks."""
    kept, crossed = _center(False), _center(False, crossed=True)
    _run(tmp_path / "serial", kept, crossed)
    _run(tmp_path / "np2", kept, crossed, ranks=2)
    serial = np.loadtxt(tmp_path / "serial" / "forces.txt", skiprows=9)
    parallel = np.loadtxt(tmp_path / "np2" / "forces.txt", skiprows=9)
    assert np.abs(serial[:, 1:]).max() > 1.0  # the crossed hydrogen is pushed back
    np.testing.assert_allclose(parallel, serial, rtol=0, atol=1e-9)
