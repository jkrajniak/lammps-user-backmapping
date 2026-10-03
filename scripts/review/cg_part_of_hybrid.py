"""Pure-CG force field from a generated hybrid force field (<prefix>.ff.lmp).

For timing the CG model on the hybrid's own CG system (review item B7): the
input reads the hybrid data file, deletes the AT atoms, and includes the file
written here. CG terms become the corresponding stock LAMMPS styles (tables,
harmonic, ...); AT types get the zero styles; CG exclusions follow cg_special
when the hybrid sets it.

    python cg_part_of_hybrid.py pe.ff.lmp > cg.ff.lmp

Stock dihedral_style table needs a range below 360 degrees; the -180..180
CG dihedral tables are rewritten next to the originals without the duplicate
180 degree point (<file>.stock).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

PLAIN = {
    "backmap/table": "table",
    "backmap/harmonic": "harmonic",
    "backmap/ryckaert": "ryckaert",
    "backmap/charmm": "charmm",
    "backmap/fourier": "fourier",
}


def _stock_dihedral_table(path: Path, keyword: str) -> str:
    """Copy of a dihedral table section without a final point at +180 deg."""
    lines = path.read_text().splitlines()
    k = next(i for i, ln in enumerate(lines) if ln.strip() == keyword)
    head = lines[k + 1].split()
    n = int(head[1])
    rows = [ln for ln in lines[k + 2 :] if ln.strip()][:n]
    if (
        abs(float(rows[-1].split()[1]) - 180.0) < 1e-9
        and abs(float(rows[0].split()[1]) + 180.0) < 1e-9
    ):
        rows = rows[:-1]
    out = path.with_name(path.name + ".stock")
    body = [f"{i + 1} " + " ".join(r.split()[1:]) for i, r in enumerate(rows)]
    extra = " ".join(head[2:])
    out.write_text(f"{keyword}\nN {len(body)} {extra}\n\n" + "\n".join(body) + "\n")
    return out.name


def convert(ff: str) -> str:
    out: list[str] = []
    pair = re.search(r"^pair_style backmap (\S+) .*?table linear (\d+)(.*)$", ff, re.MULTILINE)
    if pair is None:
        raise ValueError("no 'pair_style backmap ... table linear N' line")
    cut, npts, tail = pair.group(1), pair.group(2), pair.group(3)
    out.append(f"pair_style hybrid table linear {npts} zero {cut}")
    out.append("pair_coeff * * zero")
    for i, j, rest in re.findall(r"^pair_coeff (\d+) (\d+) cg (.*)$", ff, re.MULTILINE):
        out.append(f"pair_coeff {i} {j} table {rest}")
    for kind in ("bond", "angle", "dihedral", "improper"):
        style_line = re.search(rf"^{kind}_style (.*)$", ff, re.MULTILINE)
        if style_line is None:
            continue
        cg_lines = []
        styles: dict[str, str] = {}
        for tid, body in re.findall(rf"^{kind}_coeff (\d+) (.*)$", ff, re.MULTILINE):
            words = body.split()
            if words[0].startswith("backmap/"):  # hybrid form: style keyword args
                style, key, args = words[0], words[1], words[2:]
            else:  # single style: keyword args
                style, key, args = style_line.group(1).split()[0], words[0], words[1:]
            if key != "cg":
                continue
            plain = PLAIN.get(style, style)
            spec = re.search(rf"{re.escape(style)}( linear \d+)?", style_line.group(1))
            styles[plain] = plain + (spec.group(1) if spec and spec.group(1) else "")
            if kind == "dihedral" and plain == "table":
                args = [_stock_dihedral_table(Path(args[0]), args[1]), args[1]]
            cg_lines.append((tid, plain, " ".join(args)))
        if not cg_lines:
            out.append(f"{kind}_style zero")
            out.append(f"{kind}_coeff *")
            continue
        n_types = len(re.findall(rf"^{kind}_coeff \d+ ", ff, re.MULTILINE))
        use_zero = n_types > len(cg_lines)  # AT types of this kind exist
        if len(styles) == 1 and not use_zero:
            out.append(f"{kind}_style {next(iter(styles.values()))}")
            out += [f"{kind}_coeff {tid} {args}" for tid, _, args in cg_lines]
            continue
        out.append(f"{kind}_style hybrid {' '.join(styles.values())}{' zero' if use_zero else ''}")
        if use_zero:
            out.append(f"{kind}_coeff * zero")
        out += [f"{kind}_coeff {tid} {plain} {args}" for tid, plain, args in cg_lines]
    special = re.search(r"cg_special ([\d.]+) ([\d.]+) ([\d.]+)", tail)
    w = (
        special.groups()
        if special
        else re.search(r"^special_bonds lj (\S+) (\S+) (\S+)", ff, re.MULTILINE).groups()
    )
    out.append(f"special_bonds lj {' '.join(w)} coul {' '.join(w)}")
    comm = re.search(r"^comm_modify .*$", ff, re.MULTILINE)
    if comm:
        out.append(comm.group(0))
    out.append("neigh_modify delay 0 every 1 check yes")
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    sys.stdout.write(convert(Path(sys.argv[1]).read_text()))
