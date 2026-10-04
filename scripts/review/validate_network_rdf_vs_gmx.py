"""Check the network RDF tool (residue centres of mass, topology exclusions) against ``gmx rdf``.

One configuration of a network (the atomistic data file) is written as a GROMACS structure with the
atom, residue names and numbers of the AT topology, a run input file is made with ``grompp`` (the
topology's force-field include is pointed at ``--ff-include``), and the RDFs the reference used are
computed by ``gmx rdf`` and by network_tierc / rdf_vs_reference on the same coordinates:

- an element pair with exclusions:   ``-excl`` (``C:N:excl``)
- an element pair without exclusions: ``--plain`` (the reference's O-H, O-N, ...)
- a residue centre-of-mass pair:     ``-selrpos res_com -seltype res_com`` (``ring``)

    uv run --with numpy --with scipy --with MDAnalysis validate_network_rdf_vs_gmx.py \\
        --data melamine_network_at.data --top at_hyb_topol.top --gmx /path/to/gmx \\
        --ff-include /path/to/oplsaa.ff/forcefield.itp --excl C:N --com '*:[CN]1[1-3]'
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import network_tierc as nt
import rdf_vs_reference as rv


def frame_from_data(data: Path) -> rv.Frame:
    """The configuration of an LAMMPS data file as a frame (wrapped coordinates)."""
    lines = data.read_text().splitlines()
    box = []
    for ln in lines:
        w = ln.split()
        if len(w) >= 4 and w[2] in ("xlo", "ylo", "zlo"):
            box.append((float(w[0]), float(w[1])))
    lo = np.array([b[0] for b in box])
    length = np.array([b[1] - b[0] for b in box])
    start = next(i for i, ln in enumerate(lines) if ln.startswith("Atoms")) + 2
    rows = []
    for ln in lines[start:]:
        if not ln.strip():
            break
        rows.append(ln.split("#")[0].split())
    arr = np.array([[float(x) for x in r[:7]] for r in rows])
    arr = arr[np.argsort(arr[:, 0])]
    xyz = np.mod(arr[:, 4:7] - lo, length)
    return rv.Frame(0, length, arr[:, 0].astype(int), arr[:, 2].astype(int), xyz)


def write_gro(top: nt.Topology, frame: rv.Frame, path: Path) -> None:
    """GROMACS structure with 6 decimals (nm); the reader takes the precision from the field width."""
    n = len(top.ids)
    out = ["network configuration", f"{n}"]
    for i in range(n):
        x, y, z = frame.xyz[i] / 10.0
        out.append(
            f"{top.resnrs[i] % 100000:5d}{top.resnames[i]:<5s}{top.names[i]:>5s}"
            f"{(i + 1) % 100000:5d}{x:12.6f}{y:12.6f}{z:12.6f}"
        )
    out.append(" ".join(f"{b / 10.0:.6f}" for b in frame.box))
    path.write_text("\n".join(out) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--top", type=Path, required=True)
    ap.add_argument("--gmx", required=True)
    ap.add_argument("--ff-include", required=True, help="local oplsaa.ff/forcefield.itp")
    ap.add_argument("--excl", action="append", default=[], help="REF:SEL element pair, e.g. C:N")
    ap.add_argument(
        "--plain", action="append", default=[], help="REF:SEL element pair, no exclusions"
    )
    ap.add_argument("--com", action="append", default=[], help="RESGLOB:NAMEGLOB, e.g. TER:C[2-7]")
    ap.add_argument("--rmax", type=float, default=1.0)
    ap.add_argument("--dr", type=float, default=0.002)
    args = ap.parse_args()

    elements, bonds = rv.read_data(args.data)
    frame = frame_from_data(args.data)
    # both tools get the coordinates as written to the structure file (1e-6 nm)
    frame.xyz = np.round(frame.xyz / 10.0, 6) * 10.0
    frame.box = np.round(frame.box / 10.0, 6) * 10.0
    top = nt.read_topology(args.top)
    nt.check_topology(top, elements, frame)
    excl = rv.excluded_pairs(bonds)
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        gro = work / "conf.gro"
        write_gro(top, frame, gro)
        topfile = work / "topol.top"
        text = args.top.read_text()
        text = re.sub(
            r"^#include\s+\S+", f'#include "{args.ff_include}"', text, count=1, flags=re.M
        )
        topfile.write_text(text)
        (work / "x.mdp").write_text("integrator = md\nnsteps = 0\ncutoff-scheme = Verlet\n")
        subprocess.run(
            [
                args.gmx,
                "grompp",
                "-f",
                "x.mdp",
                "-c",
                "conf.gro",
                "-p",
                "topol.top",
                "-o",
                "t.tpr",
                "-maxwarn",
                "50",
            ],
            check=True,
            capture_output=True,
            text=True,
            cwd=work,
        )

        def gmx_rdf(ref_sel: str, sel_sel: str, *extra: str) -> tuple[np.ndarray, np.ndarray]:
            subprocess.run(
                [
                    args.gmx,
                    "rdf",
                    "-f",
                    "conf.gro",
                    "-s",
                    "t.tpr",
                    "-ref",
                    ref_sel,
                    "-sel",
                    sel_sel,
                    "-bin",
                    str(args.dr),
                    "-rmax",
                    str(args.rmax),
                    "-o",
                    "o.xvg",
                    *extra,
                ],
                check=True,
                capture_output=True,
                text=True,
                cwd=work,
            )
            return rv.read_xvg(work / "o.xvg")

        def report(label: str, r_gmx, g_gmx, r, g) -> None:
            n = min(len(r), len(r_gmx))
            diff = np.abs(g[:n] - g_gmx[:n])
            mask = r[:n] >= 0.2
            print(
                f"{label}: peak gmx {g_gmx[:n][mask].max():.4f} ours {g[:n][mask].max():.4f}; "
                f"max |g - g_gmx| {diff.max():.2e}; for r >= 0.2 nm {diff[mask].max():.2e}; "
                f"bins {len(r)} vs {len(r_gmx)}"
            )

        for spec in args.excl:
            a, b = spec.split(":")
            r_gmx, g_gmx = gmx_rdf(f'name "{a}*"', f'name "{b}*"', "-excl")
            r, g = rv.rdf([frame], elements, a, b, args.rmax, args.dr, excl)
            report(f"{a}-{b} excl", r_gmx, g_gmx, r, g)
            _, g_all = rv.rdf([frame], elements, a, b, args.rmax, args.dr, None)
            print(
                f"   (without exclusions the bonded peak is {g_all.max():.2f}; with {g.max():.2f})"
            )
        for spec in args.plain:
            a, b = spec.split(":")
            r_gmx, g_gmx = gmx_rdf(f'name "{a}*"', f'name "{b}*"')
            r, g = rv.rdf([frame], elements, a, b, args.rmax, args.dr, None)
            report(f"{a}-{b} plain", r_gmx, g_gmx, r, g)
        for spec in args.com:
            res_glob, name_glob = spec.split(":")
            sel = f'name "{name_glob}"' + ("" if res_glob == "*" else f' && resname "{res_glob}"')
            r_gmx, g_gmx = gmx_rdf(sel, sel, "-selrpos", "res_com", "-seltype", "res_com")
            group = nt.com_group(top, res_glob, name_glob)
            _, counts, norm, shell = rv.rdf_per_frame(
                [frame], group, group, args.rmax, args.dr, None
            )
            report(
                f"COM {spec}",
                r_gmx,
                g_gmx,
                np.arange(counts.shape[1]) * args.dr,
                counts.sum(axis=0) / (norm.sum() * shell),
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
