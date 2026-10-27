"""LAMMPS regression test for the ring guard (``simulation.ring_guard``, ``fix backmap/ringguard``).

An atom inside the wall radius rc of a ring centroid gets E = K (rc - d)^2 and a force along the
centroid-atom axis, the ring atoms share the reaction, ring atoms never repel their own ring, atoms
outside the fix group are ignored, and the forces do not depend on the MPI decomposition.

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
from backmap_prep.ring_guard import write_ring_guard_files

LMP_ENV = "BACKMAP_LMP"
RC_NM = 0.2
K_KJ = 41840.0
RC = 2.0  # A
K = 100.0  # kcal/mol/A^2
CENTER = np.array([15.0, 15.0, 15.0])
RING_RADIUS = 1.4


def _lmp() -> Path:
    path = os.environ.get(LMP_ENV)
    if not path or not Path(path).is_file():
        pytest.skip(f"set {LMP_ENV} to a LAMMPS binary")
    return Path(path)


def _ring_and_atom(offset: np.ndarray, guest_type: int = 1) -> System:
    """A planar six-membered ring (atoms 1-6) around CENTER and atom 7 at CENTER + offset.

    Type 1 atoms are in ``at_atoms``, type 2 is not.
    """
    system = System(box=(30.0, 30.0, 30.0))
    system.atom_types = [AtomTypeInfo(1, "C", 12.011, False), AtomTypeInfo(2, "X", 12.011, False)]
    for k in range(6):
        angle = np.pi / 3.0 * k
        pos = CENTER + RING_RADIUS * np.array([np.cos(angle), np.sin(angle), 0.0])
        system.atoms.append(LammpsAtom(k + 1, 1, 1, 0.0, *pos))
        system.bonds.append(LammpsBond(k + 1, 1, k + 1, (k + 1) % 6 + 1))
    system.atoms.append(LammpsAtom(7, 2, guest_type, 0.0, *(CENTER + offset)))
    return system


def _data(system: System) -> str:
    lines = ["ring guard test", "", f"{len(system.atoms)} atoms", "2 atom types", ""]
    lines += ["0.0 30.0 xlo xhi", "0.0 30.0 ylo yhi", "0.0 30.0 zlo zhi", "", "Masses", ""]
    lines += [f"{t.type_id} {t.mass}" for t in system.atom_types]
    lines += ["", "Atoms # atomic", ""]
    lines += [f"{a.atom_id} {a.type_id} {a.x:.10f} {a.y:.10f} {a.z:.10f}" for a in system.atoms]
    return "\n".join(lines) + "\n"


def _run(d: Path, system: System, ranks: int = 1) -> tuple[float, float, float, np.ndarray]:
    """Energy, contact count, pe and forces (id, fx, fy, fz) of one evaluation."""
    d.mkdir()
    assert write_ring_guard_files(d, system, RC_NM, K_KJ) == 1
    (d / "test.data").write_text(_data(system))
    (d / "in.test").write_text(
        f"units real\natom_style atomic\natom_modify map array\nboundary p p p\n"
        f"processors {ranks} 1 1\nread_data test.data\ngroup at_atoms type 1\n"
        "pair_style zero 5.0\npair_coeff * *\ninclude ring_guard.lmp\n"
        "dump f all custom 1 forces.txt id fx fy fz\ndump_modify f sort id format float %.12e\n"
        "run 0\n"
        'print "RESULT $(f_ring_guard[1]:%.12e) $(f_ring_guard[2]:%.1f) $(pe:%.12e)"\n'
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
    match = re.search(r"^RESULT (\S+) (\S+) (\S+)$", log, re.MULTILINE)
    assert match, log[-2000:]
    forces = np.loadtxt(d / "forces.txt", skiprows=9)
    return float(match.group(1)), float(match.group(2)), float(match.group(3)), forces


@pytest.mark.integration
def test_atom_outside_the_wall_and_the_ring_atoms_feel_nothing(tmp_path: Path) -> None:
    energy, contacts, pe, forces = _run(tmp_path / "out", _ring_and_atom(np.array([0, 0, 2.5])))
    # the ring atoms are 1.4 A from their own centroid, inside rc, and are not repelled
    assert (energy, contacts, pe) == (0.0, 0.0, 0.0)
    assert np.abs(forces[:, 1:]).max() == 0.0


@pytest.mark.integration
@pytest.mark.parametrize("height", [0.3, 1.0, 1.8])
def test_atom_inside_the_wall_is_pushed_out_along_the_axis(tmp_path: Path, height: float) -> None:
    energy, contacts, pe, forces = _run(
        tmp_path / "in", _ring_and_atom(np.array([0.0, 0.0, height]))
    )
    assert contacts == 1.0
    assert energy == pytest.approx(K * (RC - height) ** 2, rel=1e-10)
    # fix_modify energy yes: the guard energy is part of the potential energy
    assert pe == pytest.approx(energy, rel=1e-10)
    push = 2.0 * K * (RC - height)
    np.testing.assert_allclose(forces[6, 1:], [0.0, 0.0, push], atol=1e-9)
    np.testing.assert_allclose(forces[:6, 3], -push / 6.0, rtol=1e-10)
    np.testing.assert_allclose(forces[:, 1:].sum(axis=0), 0.0, atol=1e-9)


@pytest.mark.integration
def test_force_is_minus_the_energy_gradient(tmp_path: Path) -> None:
    offset = np.array([0.4, -0.3, 0.9])
    _, _, _, forces = _run(tmp_path / "f", _ring_and_atom(offset))
    h = 1e-5
    grad = []
    for k in range(3):
        step = np.zeros(3)
        step[k] = h
        e_plus = _run(tmp_path / f"p{k}", _ring_and_atom(offset + step))[0]
        e_minus = _run(tmp_path / f"m{k}", _ring_and_atom(offset - step))[0]
        grad.append((e_plus - e_minus) / (2 * h))
    np.testing.assert_allclose(forces[6, 1:], -np.array(grad), rtol=1e-5)


@pytest.mark.integration
def test_atoms_outside_the_group_are_ignored(tmp_path: Path) -> None:
    energy, contacts, _, forces = _run(
        tmp_path / "x", _ring_and_atom(np.array([0.0, 0.0, 0.5]), guest_type=2)
    )
    assert (energy, contacts) == (0.0, 0.0)
    assert np.abs(forces[:, 1:]).max() == 0.0


@pytest.mark.integration
def test_forces_do_not_depend_on_the_mpi_decomposition(tmp_path: Path) -> None:
    """The ring sits on the boundary of a two-rank split (x = 15), the guest on the other side
    of its anchor atom 1 (x = 16.4)."""
    system = _ring_and_atom(np.array([-0.6, 0.2, 0.7]))
    serial = _run(tmp_path / "serial", system)
    parallel = _run(tmp_path / "np2", system, ranks=2)
    assert serial[1] == parallel[1] == 1.0
    assert parallel[0] == pytest.approx(serial[0], rel=1e-12)
    assert np.abs(serial[3][:, 1:]).max() > 1.0
    np.testing.assert_allclose(parallel[3], serial[3], rtol=0, atol=1e-9)
