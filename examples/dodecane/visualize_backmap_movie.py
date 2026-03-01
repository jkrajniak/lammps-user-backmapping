# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "numpy",
#   "matplotlib",
# ]
# ///
"""Generate a movie of the 10-chain dodecane backmapping process.

Reads dump.movie.lammpstrj (custom dump with: id mol type x y z f_bm)
and produces an MP4 movie showing:
  - CG beads as large translucent spheres that fade as lambda -> 1
  - AT atoms as smaller spheres that appear as lambda -> 1
  - Bonds drawn within each molecule
  - Lambda value and timestep shown as text overlay

Usage:
    uv run visualize_backmap_movie.py [--input dump.movie.lammpstrj] [--output backmap_movie.mp4]
    uv run visualize_backmap_movie.py --png-sequence   # produce PNG frames instead of MP4
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

BEADS_PER_CHAIN = 6
ATOMS_PER_BEAD = 3  # 1 CG + 2 AT

TYPE_COLORS = {
    1: "#42A5F5",  # CG type A
    2: "#1E88E5",  # CG type B
    3: "#FF5722",  # AT CH3
    4: "#FF9800",  # AT CH2
}

CG_ALPHA = 0.35
CG_SIZE = {1: 350, 2: 300}
AT_SIZE = {3: 70, 4: 55}

CHAIN_PALETTE = [
    "#E53935",
    "#D81B60",
    "#8E24AA",
    "#5E35B1",
    "#3949AB",
    "#1E88E5",
    "#039BE5",
    "#00ACC1",
    "#00897B",
    "#43A047",
]


@dataclass
class Frame:
    timestep: int
    ids: np.ndarray
    mols: np.ndarray
    types: np.ndarray
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    lam: np.ndarray
    box: tuple[tuple[float, float], ...] = field(default_factory=tuple)


def parse_dump(path: Path) -> list[Frame]:
    """Parse LAMMPS custom dump: id mol type x y z f_bm."""
    frames: list[Frame] = []
    with open(path) as fh:
        while True:
            line = fh.readline()
            if not line:
                break
            if "ITEM: TIMESTEP" not in line:
                continue

            timestep = int(fh.readline().strip())
            fh.readline()  # ITEM: NUMBER OF ATOMS
            natoms = int(fh.readline().strip())
            fh.readline()  # ITEM: BOX BOUNDS

            box = []
            for _ in range(3):
                lo, hi = fh.readline().split()
                box.append((float(lo), float(hi)))

            fh.readline()  # ITEM: ATOMS ...

            ids = np.zeros(natoms, dtype=int)
            mols = np.zeros(natoms, dtype=int)
            types = np.zeros(natoms, dtype=int)
            x = np.zeros(natoms)
            y = np.zeros(natoms)
            z = np.zeros(natoms)
            lam = np.zeros(natoms)

            for i in range(natoms):
                tokens = fh.readline().split()
                ids[i] = int(tokens[0])
                mols[i] = int(tokens[1])
                types[i] = int(tokens[2])
                x[i] = float(tokens[3])
                y[i] = float(tokens[4])
                z[i] = float(tokens[5])
                lam[i] = float(tokens[6])

            order = np.argsort(ids)
            frames.append(
                Frame(
                    timestep=timestep,
                    ids=ids[order],
                    mols=mols[order],
                    types=types[order],
                    x=x[order],
                    y=y[order],
                    z=z[order],
                    lam=lam[order],
                    box=tuple(box),
                )
            )
    return frames


@dataclass
class ChainBonds:
    chain_id: int
    cg_bonds: list[tuple[int, int]]
    at_bonds: list[tuple[int, int]]


def build_chain_bonds(frame: Frame) -> list[ChainBonds]:
    """Reconstruct per-chain bonds from per-bead molecule IDs.

    Mol IDs are per-bead (1 CG + 2 AT each). Consecutive mol IDs within
    groups of BEADS_PER_CHAIN form a chain.
    """
    # Group atoms by mol ID
    mol_atoms: dict[int, list[tuple[int, int, int]]] = defaultdict(list)
    for idx in range(len(frame.ids)):
        mol_atoms[frame.mols[idx]].append((frame.ids[idx], frame.types[idx], idx))

    sorted_mols = sorted(mol_atoms.keys())
    nchains = len(sorted_mols) // BEADS_PER_CHAIN

    chains: list[ChainBonds] = []
    for c in range(nchains):
        bead_mols = sorted_mols[c * BEADS_PER_CHAIN : (c + 1) * BEADS_PER_CHAIN]

        cg_indices = []
        at_indices = []
        for mol_id in bead_mols:
            for aid, atype, idx in sorted(mol_atoms[mol_id]):
                if atype in (1, 2):
                    cg_indices.append(idx)
                else:
                    at_indices.append(idx)

        cg_bonds = [
            (cg_indices[i], cg_indices[i + 1]) for i in range(len(cg_indices) - 1)
        ]
        at_bonds = [
            (at_indices[i], at_indices[i + 1]) for i in range(len(at_indices) - 1)
        ]

        chains.append(ChainBonds(chain_id=c, cg_bonds=cg_bonds, at_bonds=at_bonds))

    return chains


def draw_frame(ax: Axes3D, frame: Frame, chains: list[ChainBonds]) -> None:
    ax.cla()

    mean_lam = frame.lam.mean()
    box = frame.box
    nchains = len(chains)

    # --- CG beads ---
    for atom_type in [1, 2]:
        mask = frame.types == atom_type
        if not mask.any():
            continue

        cg_fade = max(0.05, 1.0 - mean_lam)
        base_color = matplotlib.colors.to_rgba(TYPE_COLORS[atom_type])
        n = mask.sum()
        colors = np.array([(*base_color[:3], CG_ALPHA * cg_fade)] * n)

        ax.scatter(
            frame.x[mask],
            frame.y[mask],
            frame.z[mask],
            c=colors,
            s=CG_SIZE[atom_type] * cg_fade,
            edgecolors=(*base_color[:3], 0.5 * cg_fade),
            linewidths=0.6,
            depthshade=True,
            zorder=1,
        )

    # --- AT atoms (opacity grows with lambda) ---
    for atom_type in [3, 4]:
        mask = frame.types == atom_type
        if not mask.any():
            continue

        lam_vals = frame.lam[mask]
        alphas = np.clip(lam_vals, 0.05, 1.0)

        base_color = matplotlib.colors.to_rgba(TYPE_COLORS[atom_type])
        colors = np.array([(*base_color[:3], a) for a in alphas])
        sizes = AT_SIZE[atom_type] * np.clip(alphas, 0.3, 1.0)

        ax.scatter(
            frame.x[mask],
            frame.y[mask],
            frame.z[mask],
            c=colors,
            s=sizes,
            edgecolors="none",
            depthshade=True,
            zorder=2,
        )

    # --- Bonds per chain ---
    cg_fade = max(0.05, 1.0 - mean_lam)
    at_alpha_bond = max(0.1, mean_lam)

    for chain in chains:
        chain_rgb = matplotlib.colors.to_rgb(
            CHAIN_PALETTE[chain.chain_id % len(CHAIN_PALETTE)]
        )

        for i, j in chain.cg_bonds:
            ax.plot(
                [frame.x[i], frame.x[j]],
                [frame.y[i], frame.y[j]],
                [frame.z[i], frame.z[j]],
                color=(*matplotlib.colors.to_rgb("#1E88E5"), CG_ALPHA * cg_fade),
                linewidth=2.0 * cg_fade,
                zorder=0,
            )

        for i, j in chain.at_bonds:
            ax.plot(
                [frame.x[i], frame.x[j]],
                [frame.y[i], frame.y[j]],
                [frame.z[i], frame.z[j]],
                color=(*chain_rgb, at_alpha_bond),
                linewidth=1.2 * max(0.3, at_alpha_bond),
                zorder=1,
            )

    if box:
        ax.set_xlim(box[0][0], box[0][1])
        ax.set_ylim(box[1][0], box[1][1])
        ax.set_zlim(box[2][0], box[2][1])

    ax.set_xlabel("x (A)", fontsize=8)
    ax.set_ylabel("y (A)", fontsize=8)
    ax.set_zlabel("z (A)", fontsize=8)
    ax.tick_params(labelsize=6)

    phase = (
        "CG equilibration"
        if mean_lam < 0.01
        else ("Backmapping" if mean_lam < 0.99 else "AT production")
    )
    ax.set_title(
        f"Dodecane backmapping ({nchains} chains) — step {frame.timestep}\n"
        f"lambda = {mean_lam:.3f}  |  {phase}",
        fontsize=11,
        fontweight="bold",
    )

    ax.view_init(elev=20, azim=(30 + frame.timestep * 0.03) % 360)


def make_movie(
    frames: list[Frame],
    output: Path,
    fps: int = 15,
    dpi: int = 150,
) -> None:
    chains = build_chain_bonds(frames[0])

    fig = plt.figure(figsize=(10, 8), facecolor="white")
    ax = fig.add_subplot(111, projection="3d")

    draw_frame(ax, frames[0], chains)
    fig.tight_layout()

    def update(i: int):
        draw_frame(ax, frames[i], chains)

    anim = animation.FuncAnimation(
        fig, update, frames=len(frames), interval=1000 // fps
    )

    if output.suffix == ".gif":
        writer = animation.PillowWriter(fps=fps)
    else:
        writer = animation.FFMpegWriter(fps=fps, codec="libx264", bitrate=3000)

    print(f"Rendering {len(frames)} frames to {output} ...")
    anim.save(str(output), writer=writer, dpi=dpi)
    print(f"Done: {output}")
    plt.close(fig)


def save_png_sequence(
    frames: list[Frame],
    output_dir: Path,
    dpi: int = 150,
) -> None:
    chains = build_chain_bonds(frames[0])

    output_dir.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(10, 8), facecolor="white")
    ax = fig.add_subplot(111, projection="3d")

    for i, frame in enumerate(frames):
        draw_frame(ax, frame, chains)
        fig.tight_layout()
        out = output_dir / f"frame_{i:04d}.png"
        fig.savefig(out, dpi=dpi, bbox_inches="tight", facecolor="white")
        if (i + 1) % 20 == 0:
            print(f"  saved frame {i + 1}/{len(frames)}")

    plt.close(fig)
    print(f"PNG sequence saved to {output_dir}/")
    print(
        f"Combine with ffmpeg: ffmpeg -framerate 15 -i {output_dir}/frame_%04d.png "
        f"-c:v libx264 -pix_fmt yuv420p backmap_movie.mp4"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Visualize 10-chain dodecane backmapping trajectory as movie"
    )
    parser.add_argument(
        "--input", "-i", default="dump.movie.lammpstrj", help="LAMMPS custom dump file"
    )
    parser.add_argument(
        "--output",
        "-o",
        default="backmap_movie.mp4",
        help="Output movie file (.mp4 or .gif)",
    )
    parser.add_argument("--fps", type=int, default=15, help="Frames per second")
    parser.add_argument("--dpi", type=int, default=150, help="Resolution")
    parser.add_argument(
        "--png-sequence", action="store_true", help="Save PNG frames instead of movie"
    )
    args = parser.parse_args()

    dump_path = Path(args.input)
    if not dump_path.exists():
        print(
            f"ERROR: {dump_path} not found. Run LAMMPS first:\n"
            f"  lmp -in in.dodecane_movie"
        )
        raise SystemExit(1)

    frames = parse_dump(dump_path)
    print(f"Loaded {len(frames)} frames from {dump_path}")

    if not frames:
        print("No frames found.")
        raise SystemExit(1)

    nmols = len(set(frames[0].mols))
    lam_first = frames[0].lam.mean()
    lam_last = frames[-1].lam.mean()
    print(f"Molecules: {nmols}, lambda range: {lam_first:.4f} -> {lam_last:.4f}")

    if args.png_sequence:
        save_png_sequence(frames, Path("frames"), dpi=args.dpi)
    else:
        make_movie(frames, Path(args.output), fps=args.fps, dpi=args.dpi)


if __name__ == "__main__":
    main()
