"""Check rdf_vs_reference.py against the real ``gmx rdf`` on the same frames.

The published reference RDFs of the networks were made with ``gmx rdf -ref 'name "C*"' -sel 'name "N*"'``
(all pairs, bin 0.002 nm). This converts a LAMMPS data file and dump into a full-precision GROMACS
trajectory (TRR, atoms named after their element), runs ``gmx rdf`` for each requested pair, computes the
same RDF with rdf_vs_reference.rdf, and prints the largest difference.

    uv run --with numpy --with scipy --with MDAnalysis validate_rdf_vs_gmx.py \\
        --data rim135_at_final.data --dump dump.at_prod --gmx /path/to/gmx --pair C:N --pair C:O --pair O:H
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rdf_vs_reference as ref


def write_gromacs_files(data: Path, dump: Path, workdir: Path) -> tuple[Path, Path, int]:
    """Structure (gro) and trajectory (trr) of the dump frames, atoms named C1, N2, ... by element."""
    import MDAnalysis as mda  # noqa: N813

    elements, _ = ref.read_data(data)
    frames = ref.read_dump(dump, 0, None)
    n = len(frames[0].ids)
    names = [f"{elements[t]}{i + 1}" for i, t in enumerate(frames[0].types)]
    u = mda.Universe.empty(
        n, n_residues=1, atom_resindex=np.zeros(n, dtype=int), trajectory=True, n_frames=len(frames)
    )
    u.add_TopologyAttr("name", names)
    u.add_TopologyAttr("resname", ["MOL"])
    u.add_TopologyAttr("resid", [1])
    gro, trr = workdir / "ref.gro", workdir / "traj.trr"
    with mda.Writer(str(trr), n) as w:
        for fr in frames:
            u.dimensions = [*fr.box, 90.0, 90.0, 90.0]
            u.atoms.positions = np.mod(fr.xyz, fr.box)
            w.write(u.atoms)
    u.dimensions = [*frames[0].box, 90.0, 90.0, 90.0]
    u.atoms.positions = np.mod(frames[0].xyz, frames[0].box)
    u.atoms.write(str(gro))
    return gro, trr, len(frames)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--dump", type=Path, required=True)
    ap.add_argument("--gmx", required=True)
    ap.add_argument("--pair", action="append", required=True, help="REF:SEL, e.g. C:N")
    ap.add_argument("--rmax", type=float, default=1.0)
    ap.add_argument("--dr", type=float, default=0.002)
    args = ap.parse_args()

    elements, _ = ref.read_data(args.data)
    frames = ref.read_dump(args.dump, 0, None)
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        gro, trr, n_frames = write_gromacs_files(args.data, args.dump, work)
        print(f"{n_frames} frames, {len(frames[0].ids)} atoms")
        import MDAnalysis as mda  # noqa: N813

        # the exact single-precision coordinates GROMACS reads: no unit conversion (nm), scaled to
        # Angstrom in double precision
        u = mda.Universe(str(gro), str(trr), convert_units=False)
        trr_frames = [
            ref.Frame(
                fr.step,
                np.array(ts.dimensions[:3], dtype=float) * 10.0,
                fr.ids,
                fr.types,
                ts.positions.astype(float) * 10.0,
            )
            for fr, ts in zip(frames, u.trajectory, strict=True)
        ]
        worst = 0.0
        for spec in args.pair:
            a, b = spec.split(":")[:2]
            out = work / f"rdf_{a}_{b}.xvg"
            subprocess.run(
                [
                    args.gmx,
                    "rdf",
                    "-f",
                    str(trr),
                    "-s",
                    str(gro),
                    "-ref",
                    f'name "{a}*"',
                    "-sel",
                    f'name "{b}*"',
                    "-bin",
                    str(args.dr),
                    "-rmax",
                    str(args.rmax),
                    "-o",
                    str(out),
                ],
                check=True,
                capture_output=True,
                text=True,
                cwd=work,
            )
            r_gmx, g_gmx = ref.read_xvg(out)
            r, g = ref.rdf(frames, elements, a, b, args.rmax, args.dr, None)
            # the same frames as GROMACS read them (single-precision TRR coordinates)
            _, g2 = ref.rdf(trr_frames, elements, a, b, args.rmax, args.dr, None)
            n = min(len(r), len(r_gmx))
            dr_axis = float(np.abs(r[:n] - r_gmx[:n]).max())
            diff = np.abs(g[:n] - g_gmx[:n])
            diff_trr = np.abs(g2[:n] - g_gmx[:n])
            worst = max(worst, float(diff_trr.max()))
            peak_gmx = g_gmx[:n].max()
            print(
                f"{a}-{b}: bins {len(r)} vs {len(r_gmx)}, max |r - r_gmx| {dr_axis:.1e} nm; "
                f"peak gmx {peak_gmx:.4f}; max |g - g_gmx|: dump coordinates {diff.max():.2e}, "
                f"TRR coordinates {diff_trr.max():.2e} (gmx writes 3 decimals)"
            )
    print(f"largest difference over all pairs {worst:.2e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
