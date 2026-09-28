"""Renders of the POPC bilayer before, during and after backmapping (paper figure C9).

Reads the per-atom lambda dump of the hybrid run (``dump.backmap``: id mol type
x y z f_bm) and the hybrid data file (for the type names in Masses), and
renders side views of a slab of the bilayer with OVITO's software (Tachyon)
renderer, same camera for every panel:

- lambda = 0: CG beads only (MARTINI 3), coloured by chemical group;
- intermediate lambda: CG beads fading (transparency = lambda) over the
  atomistic fragments;
- lambda = 1: atoms only, element colours.

Water (CG W beads, TIP3P oxygens) is drawn small and translucent so the bilayer
stays visible; TIP3P hydrogens are left out. All panels share one camera.

    uv run --with ovito popc_render.py --dump dump.backmap --data popc.data \\
        --lambdas 0 0.5 1 --out popc_transition --slab 25
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

# Okabe-Ito colours (colourblind-safe); a bead and the atoms it becomes share
# a hue: choline / N blue, phosphate / P orange, glycerol-ester / O vermillion,
# tails / C grey.
BLUE, ORANGE, VERMILLION = (0.0, 0.447, 0.698), (0.902, 0.624, 0.0), (0.835, 0.369, 0.0)
GREEN, SKY = (0.0, 0.620, 0.451), (0.337, 0.706, 0.914)
CG_COLOURS = {
    "Q1": BLUE,  # choline
    "Q5": ORANGE,  # phosphate
    "SN4a": VERMILLION,  # glycerol / ester
    "N4a": VERMILLION,
    "C1": (0.70, 0.70, 0.70),  # saturated tail beads
    "C4h": GREEN,  # unsaturated tail bead
}
ELEMENT_COLOURS = {
    "H": (0.96, 0.96, 0.96),
    "C": (0.42, 0.42, 0.42),
    "N": BLUE,
    "O": VERMILLION,
    "P": ORANGE,
}
WATER_COLOUR = SKY
ELEMENT_RADII = {"H": 0.55, "C": 0.85, "N": 0.85, "O": 0.85, "P": 1.0}
CG_RADIUS = 2.35  # A, MARTINI 3 regular bead (sigma 0.47 nm / 2)


def type_table(data: Path) -> dict[int, tuple[str, bool, float]]:
    """type id -> (name, is_cg, mass) from the Masses comments of the data file."""
    table = {}
    section = False
    for line in data.read_text(encoding="utf-8").splitlines():
        if line.startswith("Masses"):
            section = True
            continue
        if section and line.strip() and not line[0].isdigit() and not line.startswith(" "):
            break
        m = re.match(r"\s*(\d+)\s+([\d.]+)\s*#\s*(\S+)(\s*\(CG\))?", line)
        if section and m:
            table[int(m.group(1))] = (m.group(3), bool(m.group(4)), float(m.group(2)))
    return table


def element_of(mass: float) -> str:
    for el, ref in (("H", 1.008), ("C", 12.011), ("N", 14.007), ("O", 15.999), ("P", 30.974)):
        if abs(mass - ref) < 0.1:
            return el
    return "C"


def frames_by_lambda(dump: Path) -> list[tuple[int, int, float]]:
    """(frame index, timestep, lambda) for every frame of the dump."""
    out = []
    lines = dump.read_text(encoding="utf-8").splitlines()
    i, k = 0, 0
    while i < len(lines):
        if lines[i].startswith("ITEM: TIMESTEP"):
            step = int(lines[i + 1])
            n = int(lines[i + 3])
            first = lines[i + 9].split()
            out.append((k, step, float(first[-1])))
            i += 9 + n
            k += 1
        else:
            i += 1
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dump", type=Path, required=True)
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--lambdas", type=float, nargs="+", default=[0.0, 0.5, 1.0])
    ap.add_argument("--slab", type=float, default=25.0, help="slab thickness along y, A")
    ap.add_argument("--size", type=int, nargs=2, default=[1600, 1600])
    ap.add_argument(
        "--gro",
        type=Path,
        default=None,
        help="optional last panel: a GROMACS frame of the continued AT system",
    )
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    from ovito.io import import_file
    from ovito.modifiers import DeleteSelectedModifier, ExpressionSelectionModifier
    from ovito.vis import TachyonRenderer, Viewport

    types = type_table(args.data)
    hw = {t for t, (name, _, _) in types.items() if name == "HW"}
    frames = frames_by_lambda(args.dump)
    print("frames:", len(frames), "lambda range", frames[0][2], frames[-1][2])

    pipeline = import_file(
        str(args.dump),
        columns=[
            "Particle Identifier",
            "Molecule Identifier",
            "Particle Type",
            "Position.X",
            "Position.Y",
            "Position.Z",
            "lambda",
        ],
    )

    def colour(frame_index, data):
        types_arr = np.asarray(data.particles.particle_types)
        lam = float(np.asarray(data.particles["lambda"])[0])
        n = len(types_arr)
        colours = np.zeros((n, 3))
        radii = np.zeros(n)
        transp = np.zeros(n)
        for t, (name, is_cg, mass) in types.items():
            sel = types_arr == t
            if name == "W":
                colours[sel] = WATER_COLOUR
                radii[sel] = CG_RADIUS
                transp[sel] = 0.9 + 0.1 * lam
                continue
            if name == "OW":
                colours[sel] = WATER_COLOUR
                radii[sel] = 1.2
                transp[sel] = 1.0 - 0.1 * lam
                continue
            if is_cg:
                colours[sel] = CG_COLOURS.get(name, (0.7, 0.7, 0.7))
                radii[sel] = CG_RADIUS * (0.8 if name in ("SN4a",) else 1.0)
                transp[sel] = lam  # beads fade out as lambda goes to 1
            else:
                el = element_of(mass)
                colours[sel] = ELEMENT_COLOURS[el]
                radii[sel] = ELEMENT_RADII[el]
                transp[sel] = 0.0 if lam > 0.0 else 1.0  # atoms appear with the ramp
        data.particles_.create_property("Color", data=colours)
        data.particles_.create_property("Radius", data=radii)
        data.particles_.create_property("Transparency", data=transp)

    cy = None
    wexpr = " || ".join(f"ParticleType == {t}" for t in sorted(hw))
    pipeline.modifiers.append(ExpressionSelectionModifier(expression=wexpr))
    pipeline.modifiers.append(DeleteSelectedModifier())
    pipeline.modifiers.append(colour)
    pipeline.add_to_scene()
    pipeline.source.data.cell.vis.enabled = False
    vp = None

    for lam in args.lambdas:
        k, step, got = min(frames, key=lambda f: abs(f[2] - lam))
        data = pipeline.compute(k)
        pos = np.asarray(data.particles.positions)
        if cy is None:
            cy = float(np.median(pos[:, 1]))
        # keep a slab along y (the viewing direction)
        slab = pipeline.modifiers
        sel = ExpressionSelectionModifier(expression=f"abs(Position.Y - {cy}) > {args.slab / 2}")
        slab.append(sel)
        slab.append(DeleteSelectedModifier())
        if vp is None:
            vp = Viewport(type=Viewport.Type.Front, fov=40.0)
            vp.zoom_all(size=tuple(args.size))
        out = f"{args.out}_lambda{got:.2f}.png"
        vp.render_image(
            filename=out,
            size=tuple(args.size),
            frame=k,
            background=(1.0, 1.0, 1.0),
            renderer=TachyonRenderer(ambient_occlusion=True),
        )
        slab.pop()
        slab.pop()
        print(f"lambda {got:.3f} (step {step}) -> {out}")
        p_type = next(t for t, (name, _, _) in types.items() if name == "PL")
        tt = np.asarray(data.particles.particle_types)
        mid_z = float(np.median(pos[tt == p_type, 2]))
        centre_xy = np.diag(np.asarray(data.cell[:, :3]))[:2] / 2.0 + np.asarray(data.cell[:2, 3])

    if args.gro is not None:
        from ovito.modifiers import AffineTransformationModifier, WrapPeriodicImagesModifier

        pipeline.remove_from_scene()
        gro = import_file(str(args.gro))
        g0 = gro.compute()
        gpos = np.asarray(g0.particles.positions)
        is_p = np.asarray(
            [
                g0.particles.particle_types.type_by_id(t).name == "P"
                for t in g0.particles.particle_types
            ]
        )
        g_centre = np.diag(np.asarray(g0.cell[:, :3]))[:2] / 2.0 + np.asarray(g0.cell[:2, 3])
        shift = [
            centre_xy[0] - g_centre[0],
            centre_xy[1] - g_centre[1],
            mid_z - float(np.median(gpos[is_p, 2])),
        ]
        gro.modifiers.append(
            AffineTransformationModifier(
                transformation=[[1, 0, 0, shift[0]], [0, 1, 0, shift[1]], [0, 0, 1, shift[2]]],
                operate_on={"particles", "cell"},
            )
        )
        gro.modifiers.append(WrapPeriodicImagesModifier())
        gro.modifiers.append(
            ExpressionSelectionModifier(expression='ResidueType == "SOL" && ParticleType == "H"')
        )
        gro.modifiers.append(DeleteSelectedModifier())

        def style_gro(frame_index, data):
            elements = [
                data.particles.particle_types.type_by_id(t).name
                for t in data.particles.particle_types
            ]
            water = np.asarray(data.particles["Residue Type"])
            n = len(elements)
            colours = np.array([ELEMENT_COLOURS.get(e, (0.5, 0.5, 0.5)) for e in elements])
            radii = np.array([ELEMENT_RADII.get(e, 0.85) for e in elements])
            transp = np.zeros(n)
            sol_id = next(t.id for t in data.particles["Residue Type"].types if t.name == "SOL")
            sol = water == sol_id
            colours[sol] = WATER_COLOUR
            radii[sol] = 1.2
            transp[sol] = 0.9
            data.particles_.create_property("Color", data=colours)
            data.particles_.create_property("Radius", data=radii)
            data.particles_.create_property("Transparency", data=transp)

        gro.modifiers.append(style_gro)
        gro.modifiers.append(
            ExpressionSelectionModifier(expression=f"abs(Position.Y - {cy}) > {args.slab / 2}")
        )
        gro.modifiers.append(DeleteSelectedModifier())
        gro.add_to_scene()
        gro.source.data.cell.vis.enabled = False
        out = f"{args.out}_gro.png"
        vp.render_image(
            filename=out,
            size=tuple(args.size),
            background=(1.0, 1.0, 1.0),
            renderer=TachyonRenderer(ambient_occlusion=True),
        )
        print(f"{args.gro.name} -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
