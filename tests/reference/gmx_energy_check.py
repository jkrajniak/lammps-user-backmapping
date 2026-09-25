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
# GROMACS terms with no comparable LAMMPS value (plain cut-off Coulomb is not
# available in GROMACS with the Verlet scheme); printed, not compared.
REPORTED = ["Coulomb (SR)"]
TERMS = {
    # LAMMPS thermo keyword -> GROMACS energy term(s), summed
    "ebond": ["Bond"],
    "eangle": ["Angle"],
    "edihed": ["Ryckaert-Bell.", "Proper Dih.", "Improper Dih."],
    "evdwl": ["LJ (SR)"],
    "lj14": ["LJ-14"],
    "coul14": ["Coulomb-14"],
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


def complete_pairs(text: str) -> str:
    """Rewrite each molecule type's [ pairs ] as the full topological 1-4 set.

    Independent of backmap-prep's implementation. A molecule type without a
    [ pairs ] section keeps none (e.g. united-atom alkanes). Otherwise every
    atom pair exactly three bonds apart (not closer) appears once: an explicit
    line with parameters is kept, a missing pair gets ``i j 1`` (gen-pairs),
    duplicates and non-1-4 lines are dropped.
    """
    blocks = re.split(r"(?=^\s*\[\s*moleculetype\s*\])", text, flags=re.MULTILINE)
    out = []
    for block in blocks:
        if not re.search(r"^\s*\[\s*pairs\s*\]", block, re.MULTILINE):
            out.append(block)
            continue
        bonds: list[tuple[int, int]] = []
        explicit: dict[tuple[int, int], str] = {}
        section = None
        kept: list[str] = []
        for line in block.splitlines():
            m = re.match(r"^\s*\[\s*(\w+)\s*\]", line)
            if m:
                section = m.group(1)
            body = line.split(";")[0].split()
            if section == "bonds" and len(body) >= 2 and body[0].isdigit():
                bonds.append((int(body[0]), int(body[1])))
            if section == "pairs" and not m:
                if len(body) >= 4 and body[0].isdigit():
                    key = (min(int(body[0]), int(body[1])), max(int(body[0]), int(body[1])))
                    explicit.setdefault(key, " ".join(body))
                continue
            kept.append(line)
        nb: dict[int, set[int]] = {}
        for i, j in bonds:
            nb.setdefault(i, set()).add(j)
            nb.setdefault(j, set()).add(i)
        pairs = set()
        for a in nb:
            first = nb[a]
            second = {c for b in first for c in nb[b]} - {a}
            third = {c for b in second for c in nb[b]}
            for c in third - second - first - {a}:
                pairs.add((min(a, c), max(a, c)))
        lines_out = []
        for line in kept:
            lines_out.append(line)
            if re.match(r"^\s*\[\s*pairs\s*\]", line):
                lines_out += [explicit.get(p, f"{p[0]} {p[1]} 1") for p in sorted(pairs)]
        out.append("\n".join(lines_out) + "\n")
    return "".join(out)


def hybrid_to_at_top(hyb_top: Path, out: Path) -> tuple[list[int], str]:
    """Standard GROMACS topology of the AT part of a bakery hybrid topology.

    Keeps atoms whose type is not a virtual (CG) type, renumbered in order;
    merges [ x ] and [ cross_x ] for bonds, angles, dihedrals and pairs, keeping
    terms whose atoms are all AT. The force-field #include is reduced to
    ``oplsaa.ff/forcefield.itp`` (resolved with --include-dir). Returns the
    hybrid indices of the AT atoms (1-based, in order) and the molecule name.
    """
    sections: dict[str, list[str]] = {}
    order: list[str] = []
    includes: list[str] = []
    current = None
    for line in hyb_top.read_text().splitlines():
        stripped = line.split(";")[0].strip()
        if stripped.startswith("#include"):
            name = stripped.split()[1].strip('"')
            includes.append(
                f'#include "oplsaa.ff/{Path(name).name}"' if "oplsaa.ff" in name else stripped
            )
            continue
        m = re.match(r"^\[\s*(\w+)\s*\]", stripped)
        if m:
            current = m.group(1)
            if current not in sections:
                sections[current] = []
                order.append(current)
            continue
        if current and stripped:
            sections[current].append(stripped)
    virtual = {
        ln.split()[0]
        for ln in sections.get("atomtypes", [])
        if len(ln.split()) >= 4 and ln.split()[3] == "V"
    }
    at_index: dict[int, int] = {}
    atoms_out: list[str] = []
    for ln in sections["atoms"]:
        parts = ln.split()
        if parts[1] in virtual:
            continue
        at_index[int(parts[0])] = len(at_index) + 1
        parts[0] = str(at_index[int(parts[0])])
        parts[5] = parts[0]
        atoms_out.append(" ".join(parts))
    n_atoms = {"bonds": 2, "angles": 3, "dihedrals": 4, "pairs": 2}
    body: dict[str, list[str]] = {}
    for name, n in n_atoms.items():
        rows = []
        for ln in sections.get(name, []) + sections.get(f"cross_{name}", []):
            parts = ln.split()
            ids = [int(v) for v in parts[:n]]
            if all(i in at_index for i in ids):
                rows.append(" ".join([*(str(at_index[i]) for i in ids), *parts[n:]]))
        body[name] = rows
    lines = [*includes, ""]
    type_width = {"bondtypes": 2, "angletypes": 3, "dihedraltypes": 4, "pairtypes": 2}
    for extra, width in type_width.items():
        rows = [
            ln for ln in sections.get(extra, []) if not virtual.intersection(ln.split()[:width])
        ]
        if rows:
            lines += [f"[ {extra} ]", *rows, ""]
    at_types = [ln for ln in sections.get("atomtypes", []) if ln.split()[0] not in virtual]
    if at_types:
        lines += ["[ atomtypes ]", *at_types, ""]
    lines += ["[ moleculetype ]", "ATSYS 3", "", "[ atoms ]", *atoms_out, ""]
    for name in n_atoms:
        lines += [f"[ {name} ]", *body[name], ""]
    lines += ["[ system ]", "AT part of hybrid", "", "[ molecules ]", "ATSYS 1", ""]
    out.write_text("\n".join(lines))
    return sorted(at_index, key=at_index.__getitem__), "ATSYS"


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


def _zero_top_charges(text: str) -> str:
    """Set the charge column of every [ atoms ] line (and atomtypes) to zero."""
    out = []
    section = None
    for line in text.splitlines():
        m = re.match(r"^\s*\[\s*(\w+)\s*\]", line)
        if m:
            section = m.group(1)
        elif section == "atoms" and line.split(";")[0].strip():
            parts = line.split(";")[0].split()
            if len(parts) >= 7:
                parts[6] = "0.0"
                line = " ".join(parts)
        out.append(line)
    return "\n".join(out) + "\n"


def gromacs_energies(
    gmx: str,
    top: Path,
    g96: Path,
    n_mol: int,
    cutoff_nm: float,
    work: Path,
    zero_charges: bool = False,
    include_dir: Path | None = None,
) -> dict[str, float]:
    text = top.read_text()
    text = re.sub(
        r"(\[\s*molecules\s*\][^\[]*?\n\s*)(\S+)(\s+)\d+",
        lambda m: f"{m.group(1)}{m.group(2)}{m.group(3)}{n_mol}",
        text,
        count=1,
    )
    if zero_charges:
        text = _zero_top_charges(text)
    text = complete_pairs(text)
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
                # Networks span the box (a molecule bonded to its own image);
                # without this GROMACS mishandles their exclusions.
                "periodic-molecules = yes",
                "nstcalcenergy = 1",
                "nstenergy = 1",
                "constraints = none",
                *([f"include = -I{include_dir}"] if include_dir else []),
                "",
            ]
        )
    )

    def run(*args: str, input: str | None = None) -> subprocess.CompletedProcess[str]:
        proc = subprocess.run(
            [gmx, *args], cwd=work, capture_output=True, text=True, check=False, input=input
        )
        if proc.returncode != 0:
            raise SystemExit(f"gmx {args[0]} failed:\n{proc.stdout[-2000:]}\n{proc.stderr[-4000:]}")
        return proc

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
    wanted = sorted({t for terms in TERMS.values() for t in terms} | set(REPORTED))
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


def _expanded_lines(path: Path) -> list[str]:
    """Input lines with ``include`` files expanded in place (relative to the input)."""
    out: list[str] = []
    for line in path.read_text().splitlines():
        parts = line.split()
        if parts[:1] == ["include"] and len(parts) >= 2:
            out.extend(_expanded_lines(path.parent / parts[1]))
        else:
            out.append(line)
    return out


def _ff_header(input_path: Path, data: Path) -> list[str]:
    """The generated input up to fix backmap, reading ``data``, at lambda = 1."""
    lines = []
    text = _expanded_lines(input_path)
    pairs_fix = [line for line in text if re.match(r"^fix\s+\S+\s+\S+\s+backmap/pairs\b", line)]
    for line in text:
        if line.startswith("fix bm "):
            line = re.sub(r"lambda0\s+\S+", "lambda0 1.0", line)
            lines.append(line)
            lines.extend(pairs_fix)
            break
        if line.split()[:1] in (["dump"], ["dump_modify"], ["run"], ["minimize"]):
            continue
        lines.append(f"read_data {data.resolve()}" if line.startswith("read_data") else line)
    if not any(line.startswith("fix bm ") for line in lines):
        raise SystemExit(f"{input_path}: no 'fix bm' found (after expanding includes)")
    return lines


def at_lj_cutoff_nm(input_path: Path) -> float:
    """LJ cutoff of the AT sub-style of ``pair_style backmap`` (nm).

    ``pair_style backmap <cut> <at-style> <lj-cut> [<coul-cut>] <cg-cut> <cg-style> ...``;
    the global and CG cutoffs can differ from the AT LJ cutoff.
    """
    for line in _expanded_lines(input_path):
        parts = line.split()
        if parts[:2] == ["pair_style", "backmap"] and len(parts) >= 5:
            return float(parts[4]) / 10.0
    raise SystemExit(f"{input_path}: no 'pair_style backmap' found")


def relax(lmp: str, input_path: Path, data: Path, steps: int) -> Path:
    """Minimize at lambda = 1 so the comparison frame has no overlapping atoms.

    The unrelaxed hybrid places fragments independently, so neighbouring
    fragments overlap; near-degenerate dihedrals and r -> 0 pairs are then
    evaluated differently by any two codes. Returns the relaxed data file.
    """
    out = input_path.parent / "relaxed_check.data"
    lines = [
        *_ff_header(input_path, data),
        "fix_modify bm active no",
        # Relax without charges: at lambda = 1 in an overlapping frame,
        # opposite charges (OPLS hydroxyl H has no LJ core) collapse onto
        # each other. The comparison itself uses the original charges.
        "fix qsave all store/state 0 q",
        "set group all charge 0.0",
        f"minimize 0.0 1.0e-6 {steps} {10 * steps}",
        "variable qrest atom f_qsave",
        "set group all charge v_qrest",
        "unfix qsave",
        f"write_data {out} nocoeff",
    ]
    (input_path.parent / "in.relax_check").write_text("\n".join(lines) + "\n")
    subprocess.run(
        [lmp, "-in", "in.relax_check", "-log", "log.relax_check", "-screen", "none"],
        cwd=input_path.parent,
        check=True,
    )
    return out


def lammps_energies(
    lmp: str, input_path: Path, data: Path, work: Path, zero_charges: bool = False
) -> dict[str, float]:
    lines = _ff_header(input_path, data)
    if zero_charges:
        lines.append("set group all charge 0.0")
    pairs_id = next(
        (
            line.split()[1]
            for line in lines
            if re.match(r"^fix\s+\S+\s+\S+\s+backmap/pairs\b", line)
        ),
        None,
    )
    lj14 = f"$(f_{pairs_id}[1]:%.12e)" if pairs_id else "0.0"
    coul14 = f"$(f_{pairs_id}[2]:%.12e)" if pairs_id else "0.0"
    lines += [
        "thermo_style custom step ebond eangle edihed evdwl ecoul",
        "run 0",
        'print "RESULT ebond=$(ebond:%.12e) eangle=$(eangle:%.12e) edihed=$(edihed:%.12e) '
        f'evdwl=$(evdwl:%.12e) ecoul=$(ecoul:%.12e) lj14={lj14} coul14={coul14}"',
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


def at_only_energies(
    lmp: str, at_data: Path, at_ff: Path, workdir: Path, zero_charges: bool
) -> dict[str, float]:
    """LAMMPS run 0 of an AT-only data file with its generated force field."""
    lines = [
        "units real",
        "atom_style full",
        "boundary p p p",
        f"read_data {at_data.resolve()}",
        f"include {at_ff.resolve()}",
        *(["set group all charge 0.0"] if zero_charges else []),
        "thermo_style custom step ebond eangle edihed evdwl ecoul",
        "run 0",
        'print "RESULT ebond=$(ebond:%.12e) eangle=$(eangle:%.12e) edihed=$(edihed:%.12e) '
        'evdwl=$(evdwl:%.12e) ecoul=$(ecoul:%.12e)"',
    ]
    (workdir / "in.at_check").write_text("\n".join(lines) + "\n")
    subprocess.run(
        [lmp, "-in", "in.at_check", "-log", "log.at_check", "-screen", "none"],
        cwd=workdir,
        check=True,
    )
    log = (workdir / "log.at_check").read_text()
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
    ap.add_argument(
        "--cutoff",
        type=float,
        default=None,
        help="LJ cutoff, nm (default: the AT LJ cutoff of pair_style backmap in the input)",
    )
    ap.add_argument("--rtol", type=float, default=1e-6)
    ap.add_argument(
        "--atol",
        type=float,
        default=1e-4,
        help="absolute tolerance, kJ/mol (data-file coordinates carry 1e-6 A)",
    )
    ap.add_argument(
        "--include-dir",
        type=Path,
        default=None,
        help="directory holding e.g. oplsaa.ff/ for #include in the AT topology",
    )
    ap.add_argument(
        "--zero-charges",
        action="store_true",
        help="compare with all charges zero on both sides (pair_style backmap reports "
        "LJ + cut-off Coulomb as evdwl; GROMACS has no plain cut-off Coulomb)",
    )
    ap.add_argument(
        "--relax",
        type=int,
        default=0,
        help="minimize this many steps at lambda = 1 first and compare on that frame",
    )
    ap.add_argument(
        "--hybrid-top",
        default=None,
        help="compare against the AT part of this bakery hybrid topology instead of "
        "the AT source topology (networks); --top/--gro are then ignored",
    )
    ap.add_argument(
        "--at-data",
        type=Path,
        default=None,
        help="AT-only data file with the same coordinates (backmap-prep at-system --from)",
    )
    ap.add_argument("--at-ff", type=Path, default=None, help="Its <prefix>.at.ff.lmp")
    args = ap.parse_args()

    d = args.dir.resolve()
    data, inp, gro, top = d / args.data, d / args.input, d / args.gro, d / args.top
    if args.relax:
        data = relax(args.lmp, inp, data, args.relax)
    if args.hybrid_top:
        hyb_top = d / args.hybrid_top
        top = d / "hybrid_at_check.top"
        at_ids, resname = hybrid_to_at_top(hyb_top, top)
        box, atoms_all, _ = read_data(data)
        coords = [atoms_all[i][2] for i in at_ids]
        n_mol = 1
        names = template_order(top)[1]
    else:
        _, _, masses = read_data(data)
        at_types = {t for t, m in masses.items() if m < 20.0}
        box, coords, n_mol = at_coordinates(data, gro, top, at_types)
        resname, names = template_order(top)

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        write_g96(work / "conf.g96", box, coords, resname, names)
        gmx = gromacs_energies(
            args.gmx,
            top,
            work / "conf.g96",
            n_mol,
            args.cutoff if args.cutoff is not None else at_lj_cutoff_nm(inp),
            work,
            args.zero_charges,
            args.include_dir.resolve() if args.include_dir else None,
        )
        lmp = lammps_energies(args.lmp, inp, data, work, args.zero_charges)

    failed = False
    print(f"{'term':8s} {'LAMMPS (kJ/mol)':>18s} {'GROMACS (kJ/mol)':>18s} {'rel diff':>10s}")
    charged = abs(gmx.get("Coulomb (SR)", 0.0)) > 0.0
    for key, gmx_terms in TERMS.items():
        if key == "evdwl" and charged:
            print(
                f"{key:8s} skipped: charges present (evdwl holds LJ + cut-off Coulomb); "
                "rerun with --zero-charges"
            )
            continue
        ref = sum(gmx.get(t, 0.0) for t in gmx_terms)
        val = lmp[key] * KCAL_TO_KJ
        rel = abs(val - ref) / max(abs(ref), 1e-12)
        ok = rel <= args.rtol or abs(val - ref) <= args.atol
        failed |= not ok
        print(f"{key:8s} {val:18.6f} {ref:18.6f} {rel:10.2e} {'ok' if ok else 'FAIL'}")
    if args.at_data is not None and args.at_ff is not None:
        at = at_only_energies(args.lmp, args.at_data, args.at_ff, d, args.zero_charges)
        # special_bonds puts the 1-4 terms into evdwl/ecoul of the AT-only run
        at_terms = {
            "ebond": ["Bond"],
            "eangle": ["Angle"],
            "edihed": ["Ryckaert-Bell.", "Proper Dih.", "Improper Dih."],
            "evdwl": ["LJ (SR)", "LJ-14"],
        }
        print("AT-only force field (at-system data + .at.ff.lmp):")
        for key, gmx_terms in at_terms.items():
            if key == "evdwl" and charged:
                print(f"  {key:8s} skipped: charges present; rerun with --zero-charges")
                continue
            ref = sum(gmx.get(t, 0.0) for t in gmx_terms)
            val = at[key] * KCAL_TO_KJ
            rel = abs(val - ref) / max(abs(ref), 1e-12)
            ok = rel <= args.rtol or abs(val - ref) <= args.atol
            failed |= not ok
            print(f"  {key:8s} {val:18.6f} {ref:18.6f} {rel:10.2e} {'ok' if ok else 'FAIL'}")
    for term in REPORTED:
        if term in gmx:
            print(f"{term:12s} GROMACS {gmx[term]:14.6f} kJ/mol (no LAMMPS thermo counterpart)")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
