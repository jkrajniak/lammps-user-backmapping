"""``dihedral_style backmap/table`` reads the table at the IUPAC dihedral, as ``dihedral_style table``.

GROMACS tabulated dihedrals and the native LAMMPS table style use the IUPAC sign; backmap/table
took the sign of dihedral_style opls (-phi), so an asymmetric table was mirrored between the
CG-only input (native table) and the hybrid stage (code review F10).

With an asymmetric table V(phi) = 1 - cos(phi - 60 deg) the energy and the forces of backmap/table
(atoms in no bead, lambda 1, weight 1) must equal those of dihedral_style table, for a geometry
and its mirror image, and the energy at the test points of the review (IUPAC phi = +90 deg) must be
V(+90).

Skipped unless ``BACKMAP_LMP`` points to a LAMMPS binary with the backmap and MOLECULE packages.
"""

from __future__ import annotations

import math
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

LMP_ENV = "BACKMAP_LMP"
N = 3600  # -180 to 179.9 deg: dihedral_style table needs a range below 360 deg
GEOMETRIES = {
    "chain": [(10.0, 10.0, 10.0), (11.3, 10.4, 10.2), (12.0, 11.6, 10.9), (13.4, 11.9, 10.3)],
    "mirror": [(10.0, 10.0, 10.0), (11.3, 10.4, 9.8), (12.0, 11.6, 9.1), (13.4, 11.9, 9.7)],
    "iupac_plus_90": [
        (11.0, 10.0, 10.0),
        (10.0, 10.0, 10.0),
        (10.0, 10.0, 11.0),
        (10.0, 11.0, 11.0),
    ],
}


def _v(phi_deg: float) -> tuple[float, float]:
    """V and -dV/dphi (per degree) of the asymmetric test table."""
    p = math.radians(phi_deg - 60.0)
    return 1.0 - math.cos(p), -math.sin(p) * math.pi / 180.0


def _table(d: Path) -> None:
    with (d / "dihedral.table").open("w") as fh:
        fh.write(f"D\nN {N}\n\n")
        for i in range(N):
            phi = -180.0 + 360.0 * i / N
            e, f = _v(phi)
            fh.write(f"{i + 1} {phi:.8f} {e:.15g} {f:.15g}\n")


def _run(tmp: Path, geometry: str, style: str) -> tuple[float, np.ndarray]:
    lmp = os.environ.get(LMP_ENV)
    if not lmp or not Path(lmp).is_file():
        pytest.skip(f"set {LMP_ENV} to a LAMMPS binary built with the backmap package")
    d = tmp / f"{geometry}_{style.split()[0].replace('/', '_')}"
    d.mkdir()
    _table(d)
    atoms = "\n".join(
        f"{k} {k} 1 0.0 {x} {y} {z}" for k, (x, y, z) in enumerate(GEOMETRIES[geometry], 1)
    )
    (d / "d.data").write_text(
        "dihedral table sign\n\n4 atoms\n2 atom types\n1 dihedrals\n1 dihedral types\n\n"
        "0.0 30.0 xlo xhi\n0.0 30.0 ylo yhi\n0.0 30.0 zlo zhi\n\nMasses\n\n1 12.0\n2 12.0\n\n"
        f"Atoms # full\n\n{atoms}\n\nDihedrals\n\n1 1 1 2 3 4\n"
    )
    coeff = "1 at dihedral.table D" if style.startswith("backmap") else "1 dihedral.table D"
    (d / "in.t").write_text(
        f"""units real
atom_style full
boundary p p p
atom_modify map array
read_data d.data
pair_style zero 8.0
pair_coeff * *
dihedral_style {style}
dihedral_coeff {coeff}
fix bm all backmap cg_type 2 alpha 0.0001 lambda0 1.0
dump d all custom 1 f.dump id fx fy fz
dump_modify d sort id format float %.15g
thermo_style custom step edihed
run 0
print "RESULT $(edihed:%.15g)"
"""
    )
    proc = subprocess.run(
        [lmp, "-in", "in.t", "-log", "log.t", "-screen", "none"],
        cwd=d,
        capture_output=True,
        text=True,
        check=False,
    )
    log = (d / "log.t").read_text() if (d / "log.t").exists() else ""
    assert proc.returncode == 0, (log + proc.stderr)[-2000:]
    energy = float(next(ln.split()[1] for ln in log.splitlines() if ln.startswith("RESULT ")))
    lines = (d / "f.dump").read_text().splitlines()
    head = next(i for i, ln in enumerate(lines) if ln.startswith("ITEM: ATOMS"))
    return energy, np.loadtxt(lines[head + 1 :], ndmin=2)[:, 1:4]


@pytest.mark.integration
@pytest.mark.parametrize("geometry", sorted(GEOMETRIES))
def test_backmap_table_equals_native_table(tmp_path: Path, geometry: str) -> None:
    e_native, f_native = _run(tmp_path, geometry, f"table linear {N}")
    e_backmap, f_backmap = _run(tmp_path, geometry, f"backmap/table linear {N}")
    assert e_backmap == pytest.approx(e_native, rel=1e-6, abs=1e-9)
    scale = max(np.abs(f_native).max(), 1e-3)
    np.testing.assert_allclose(f_backmap, f_native, rtol=0, atol=1e-4 * scale)


@pytest.mark.integration
def test_energy_is_the_table_at_the_iupac_angle(tmp_path: Path) -> None:
    energy, _ = _run(tmp_path, "iupac_plus_90", f"backmap/table linear {N}")
    assert energy == pytest.approx(_v(90.0)[0], abs=1e-6)
    assert energy != pytest.approx(_v(-90.0)[0], abs=1e-2)
