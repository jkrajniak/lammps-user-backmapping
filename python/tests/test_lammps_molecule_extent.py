"""fix backmap stops with a clear error when a molecule ID extends beyond the ghost region.

The bead map is built from the atoms one rank sees (local + ghost). With one molecule ID per
chain, a chain longer than the ghost region gives ranks different partial views: the run stopped
with a misleading divisibility error, or a partial view whose counts happened to divide gave a
wrong map, and the check could hang the ranks that passed it (code review F12). Molecules must
be complete in local + ghost on every rank that owns one of their atoms.

A straight chain of 12 beads (2 AT atoms each) along x on 4 ranks: one molecule ID per bead runs,
one ID for the whole chain stops with the extent error; serially both run.

Skipped unless ``BACKMAP_LMP`` points to a LAMMPS binary with the backmap package, and for the
parallel runs ``mpirun``.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

LMP_ENV = "BACKMAP_LMP"
N_BEADS = 12
SPACING = 3.2  # A between beads along x


def _data(per_bead: bool) -> str:
    atoms, bonds = [], []
    for b in range(N_BEADS):
        mol = b + 1 if per_bead else 1
        x = 2.0 + SPACING * b
        base = 3 * b
        atoms.append(f"{base + 1} {mol} 1 0.0 {x:.3f} 10.0 10.0")
        atoms.append(f"{base + 2} {mol} 2 0.0 {x - 0.75:.3f} 10.4 10.0")
        atoms.append(f"{base + 3} {mol} 2 0.0 {x + 0.75:.3f} 9.6 10.0")
        bonds.append((base + 2, base + 3))
        if b > 0:
            bonds.append((base - 0, base + 2))
    lines = ["molecule extent test", "", f"{3 * N_BEADS} atoms", "2 atom types"]
    lines += [f"{len(bonds)} bonds", "1 bond types", ""]
    lines += [f"0.0 {SPACING * N_BEADS + 2.0:.3f} xlo xhi", "0.0 20.0 ylo yhi", "0.0 20.0 zlo zhi"]
    lines += ["", "Masses", "", "1 28.0", "2 14.0", "", "Atoms # full", "", *atoms, "", "Bonds", ""]
    lines += [f"{k} 1 {i} {j}" for k, (i, j) in enumerate(bonds, 1)]
    return "\n".join(lines) + "\n"


def _run(tmp: Path, nprocs: int, per_bead: bool) -> tuple[int, str]:
    lmp = os.environ.get(LMP_ENV)
    if not lmp or not Path(lmp).is_file():
        pytest.skip(f"set {LMP_ENV} to a LAMMPS binary built with the backmap package")
    d = tmp / f"np{nprocs}_{'bead' if per_bead else 'chain'}"
    d.mkdir()
    (d / "d.data").write_text(_data(per_bead))
    (d / "in.t").write_text(
        f"""units real
atom_style full
boundary p p p
processors {nprocs} 1 1
read_data d.data
pair_style zero 4.0
pair_coeff * *
bond_style backmap/harmonic
bond_coeff 1 at 300.0 1.53
fix bm all backmap cg_type 1 alpha 0.0001 lambda0 0.5
run 0
"""
    )
    cmd = [lmp, "-in", "in.t", "-log", "log.t", "-screen", "none"]
    if nprocs > 1:
        mpirun = shutil.which("mpirun")
        if mpirun is None:
            pytest.skip("mpirun not available")
        cmd = [mpirun, "-np", str(nprocs), *cmd]
    proc = subprocess.run(cmd, cwd=d, capture_output=True, text=True, check=False, timeout=120)
    log = (d / "log.t").read_text() if (d / "log.t").exists() else ""
    return proc.returncode, log + proc.stdout + proc.stderr


@pytest.mark.integration
@pytest.mark.parametrize(("nprocs", "per_bead"), [(1, True), (1, False), (4, True)])
def test_complete_molecules_run(tmp_path: Path, nprocs: int, per_bead: bool) -> None:
    code, out = _run(tmp_path, nprocs, per_bead)
    assert code == 0, out[-2000:]


@pytest.mark.integration
def test_chain_molecule_beyond_the_ghost_region_is_an_error(tmp_path: Path) -> None:
    code, out = _run(tmp_path, 4, per_bead=False)
    assert code != 0
    assert "extends beyond the ghost region" in out, out[-2000:]
