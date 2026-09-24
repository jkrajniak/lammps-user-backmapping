#!/usr/bin/env python3
"""Check a generated hybrid system against GROMACS on its AT source topology.

At lambda = 1 the atomistic part of a hybrid system must carry exactly the force
field of the AT source topology. This script evaluates both on the same
coordinates and compares per-term energies:

* LAMMPS: the generated input's force-field block, ``fix backmap lambda0 1.0``
  with the ramp inactive, ``run 0``.
* GROMACS (double precision): the AT source ``.top`` (one molecule type,
  repeated), the AT atoms of the same data file in template order, ``nsteps 0``.

Needs the compute VM (LAMMPS with the backmap package, ``gmx_d``). Standard
library only.

Usage:
  gmx_energy_check.py --dir EXAMPLE_DIR --data pe.data --input in.pe \\
      --gro hyb_conf.gro --top topol_aa.top --lmp LMP --gmx GMX_D [--rtol 1e-5]
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

KCAL_TO_KJ = 4.184
TERMS = {
    # LAMMPS thermo keyword -> GROMACS energy term(s), summed
    "ebond": ["Bond"],
    "eangle": ["Angle"],
    "edihed": ["Ryckaert-Bell.", "Proper Dih.", "Improper Dih."],
    "evdwl": ["LJ (SR)"],
}


def read_data(
    path: Path,
) -> tuple[list[float], dict[int, tuple[int, int, list[float]]], dict[int, float]]:
    box: list[float] = []
    atoms: dict[int, tuple[int, int, list[float]]] = {}
    masses: dict[int, float] = {}
    section = None
    for line in path.read_text().splitlines():
        s = line.split("#")[0].strip()
        if not s:
            continue
        if re.match(r"^[A-Z][A-Za-z ]+$", s):
            section = s
            continue
        p = s.split()
        if len(p) >= 4 and p[2] in ("xlo", "ylo", "zlo"):
            box.append(float(p[1]) - float(p[0]))
        elif section == "Masses":
            masses[int(p[0])] = float(p[1])
        elif section == "Atoms":
            x = [float(v) for v in p[4:7]]
            img = [int(v) for v in p[7:10]] if len(p) >= 10 else [0, 0, 0]
            unwrapped = [x[k] + img[k] * box[k] for k in range(3)]
            atoms[int(p[0])] = (int(p[1]), int(p[2]), unwrapped)
    return box, atoms, masses


def read_gro_names(path: Path) -> list[tuple[int, str]]:
    """(residue number, atom name) per atom; the residue groups one AT molecule."""
    lines = path.read_text().splitlines()
    n = int(lines[1])
    return [(int(lines[2 + i][0:5]), lines[2 + i][10:15].strip()) for i in range(n)]


def template_order(top: Path) -> tuple[str, list[str]]:
    """Molecule type name and atom names, in order, of the first moleculetype."""
    name = ""
    names: list[str] = []
    section = None
    for line in top.read_text().splitlines():
        s = line.split(";")[0].strip()
        if not s:
            continue
        m = re.match(r"^\[\s*(\w+)\s*\]$", s)
        if m:
            section = m.group(1)
            if section == "moleculetype" and names:
                break
            continue
        if section == "moleculetype" and not name:
            name = s.split()[0]
        elif section == "atoms":
            names.append(s.split()[4])
    return name, names


def at_coordinates(
    data: Path, gro: Path, top: Path, at_types: set[int]
) -> tuple[list[float], list[list[float]], int]:
    box, atoms, _ = read_data(data)
    gro_names = read_gro_names(gro)
    if len(gro_names) != len(atoms):
        raise SystemExit(f"{gro}: {len(gro_names)} atoms, {data}: {len(atoms)}")
    _, tpl = template_order(top)
    index = {n: i for i, n in enumerate(tpl)}
    by_mol: dict[int, list[tuple[int, list[float]]]] = {}
    # Group by the .gro residue, not the data-file molecule ID: the hybrid
    # gives each bead and its fragment its own molecule ID for fix backmap.
    for aid in sorted(atoms):
        _, typ, x = atoms[aid]
        if typ not in at_types:
            continue
        resid, name = gro_names[aid - 1]
        by_mol.setdefault(resid, []).append((index[name], x))
    coords: list[list[float]] = []
    for mol in sorted(by_mol):
        entries = sorted(by_mol[mol])
        if [e[0] for e in entries] != list(range(len(tpl))):
            raise SystemExit(f"molecule {mol}: AT atoms do not cover the template once")
        coords.extend(e[1] for e in entries)
    return box, coords, len(by_mol)


def write_g96(
    path: Path, box: list[float], coords: list[list[float]], resname: str, names: list[str]
) -> None:
    n_tpl = len(names)
    out = ["TITLE", "gmx_energy_check", "END", "POSITION"]
    for i, x in enumerate(coords):
        res = i // n_tpl + 1
        out.append(
            f"{res:5d} {resname[:5]:5s} {names[i % n_tpl][:5]:5s}{i + 1:7d}"
            f"{x[0] / 10:15.9f}{x[1] / 10:15.9f}{x[2] / 10:15.9f}"
        )
    out += ["END", "BOX", f"{box[0] / 10:15.9f}{box[1] / 10:15.9f}{box[2] / 10:15.9f}", "END"]
    path.write_text("\n".join(out) + "\n")


def gromacs_energies(
    gmx: str, top: Path, g96: Path, n_mol: int, cutoff_nm: float, work: Path
) -> dict[str, float]:
    text = top.read_text()
    text = re.sub(
        r"(\[\s*molecules\s*\][^\[]*?\n\s*)(\S+)(\s+)\d+",
        lambda m: f"{m.group(1)}{m.group(2)}{m.group(3)}{n_mol}",
        text,
        count=1,
    )
    (work / "topol.top").write_text(text)
    (work / "run.mdp").write_text(
        "\n".join(
            [
                "integrator = md",
                "nsteps = 0",
                "cutoff-scheme = Verlet",
                "verlet-buffer-tolerance = -1",
                f"rlist = {cutoff_nm}",
                "vdwtype = Cut-off",
                "vdw-modifier = None",
                f"rvdw = {cutoff_nm}",
                "coulombtype = Cut-off",
                f"rcoulomb = {cutoff_nm}",
                "DispCorr = no",
                "pbc = xyz",
                "nstcalcenergy = 1",
                "nstenergy = 1",
                "constraints = none",
                "",
            ]
        )
    )

    def run(*args: str, input: str | None = None) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [gmx, *args], cwd=work, capture_output=True, text=True, check=True, input=input
        )

    run(
        "grompp",
        "-f",
        "run.mdp",
        "-c",
        g96.name,
        "-p",
        "topol.top",
        "-o",
        "run.tpr",
        "-maxwarn",
        "5",
    )
    run("mdrun", "-s", "run.tpr", "-deffnm", "run", "-nt", "1")
    wanted = sorted({t for terms in TERMS.values() for t in terms})
    available = run("dump", "-e", "run.edr").stdout
    present = [t for t in wanted if re.search(rf"(^|\s){re.escape(t)}\s", available, re.MULTILINE)]
    run(
        "energy",
        "-f",
        "run.edr",
        "-dp",
        "-o",
        "energy.xvg",
        input="\n".join(t.replace(" ", "-") for t in present) + "\n0\n",
    )
    xvg = (work / "energy.xvg").read_text().splitlines()
    legends = [m.group(1) for line in xvg if (m := re.match(r'@ s\d+ legend "(.*)"', line))]
    first = next(line for line in xvg if line and line[0] not in "#@")
    values = [float(v) for v in first.split()[1:]]
    energies: dict[str, float] = dict(zip(legends, values, strict=False))
    missing = [t for t in present if t not in energies]
    if missing:
        raise SystemExit(f"GROMACS terms not read back: {missing}")
    return energies


def lammps_energies(lmp: str, input_path: Path, data: Path, work: Path) -> dict[str, float]:
    lines = []
    for line in input_path.read_text().splitlines():
        if line.startswith("fix bm "):
            line = re.sub(r"lambda0\s+\S+", "lambda0 1.0", line)
            lines.append(line)
            break
        if line.split()[:1] in (["dump"], ["dump_modify"], ["run"], ["minimize"]):
            continue
        lines.append(line.replace(data.name, str(data.resolve())) if "read_data" in line else line)
    lines += [
        "thermo_style custom step ebond eangle edihed evdwl ecoul",
        "run 0",
        'print "RESULT ebond=$(ebond:%.12e) eangle=$(eangle:%.12e) edihed=$(edihed:%.12e) evdwl=$(evdwl:%.12e) ecoul=$(ecoul:%.12e)"',
    ]
    (input_path.parent / "in.check").write_text("\n".join(lines) + "\n")
    subprocess.run(
        [lmp, "-in", "in.check", "-log", "log.check", "-screen", "none"],
        cwd=input_path.parent,
        check=True,
    )
    log = (input_path.parent / "log.check").read_text()
    m = re.search(r"^RESULT (.*)$", log, re.MULTILINE)
    if not m:
        raise SystemExit(log[-2000:])
    return {k: float(v) for k, v in (kv.split("=") for kv in m.group(1).split())}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--data", required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--gro", required=True)
    ap.add_argument("--top", required=True)
    ap.add_argument("--lmp", required=True)
    ap.add_argument("--gmx", required=True)
    ap.add_argument("--cutoff", type=float, default=1.4, help="LJ cutoff, nm")
    ap.add_argument("--rtol", type=float, default=1e-5)
    args = ap.parse_args()

    d = args.dir.resolve()
    data, inp, gro, top = d / args.data, d / args.input, d / args.gro, d / args.top
    _, _, masses = read_data(data)
    at_types = {t for t, m in masses.items() if m < 20.0}
    box, coords, n_mol = at_coordinates(data, gro, top, at_types)
    resname, names = template_order(top)

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        write_g96(work / "conf.g96", box, coords, resname, names)
        gmx = gromacs_energies(args.gmx, top, work / "conf.g96", n_mol, args.cutoff, work)
        lmp = lammps_energies(args.lmp, inp, data, work)

    failed = False
    print(f"{'term':8s} {'LAMMPS (kJ/mol)':>18s} {'GROMACS (kJ/mol)':>18s} {'rel diff':>10s}")
    for key, gmx_terms in TERMS.items():
        ref = sum(gmx.get(t, 0.0) for t in gmx_terms)
        val = lmp[key] * KCAL_TO_KJ
        rel = abs(val - ref) / max(abs(ref), 1e-12)
        ok = rel <= args.rtol or abs(val - ref) < 1e-6
        failed |= not ok
        print(f"{key:8s} {val:18.6f} {ref:18.6f} {rel:10.2e} {'ok' if ok else 'FAIL'}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
