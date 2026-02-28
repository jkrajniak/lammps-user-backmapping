# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "numpy",
#   "matplotlib",
# ]
# ///
"""Generate a movie of the backmapping process from LAMMPS custom dump.

Reads dump.movie.lammpstrj (custom dump with: id mol type x y z f_bm)
and produces an MP4 movie showing:
  - CG beads as large translucent spheres that fade as lambda → 1
  - AT atoms as smaller spheres that appear as lambda → 1
  - Bonds drawn between connected atoms
  - Lambda value and timestep shown as text overlay

Usage:
    uv run visualize_backmap_movie.py [--input dump.movie.lammpstrj] [--output backmap_movie.mp4]
    uv run visualize_backmap_movie.py --png-sequence   # produce PNG frames instead of MP4
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


@dataclass
class Frame:
    timestep: int
    ids: np.ndarray
    types: np.ndarray
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    lam: np.ndarray


# CG beads (types 1,2) mapped to AT atoms (types 3,4)
# Molecule 1: CG atoms 1-6, AT atoms 7-18
# Each CG bead maps to 2 AT atoms:
#   CG 1 (A1) → AT 7,8   (C1,C2)
#   CG 2 (B1) → AT 9,10  (C3,C4)
#   CG 3 (B2) → AT 11,12 (C5,C6)
#   CG 4 (B3) → AT 13,14 (C7,C8)
#   CG 5 (B4) → AT 15,16 (C9,C10)
#   CG 6 (A2) → AT 17,18 (C11,C12)

BONDS_AT = [
    (7, 8),
    (8, 9),
    (9, 10),
    (10, 11),
    (11, 12),
    (12, 13),
    (13, 14),
    (14, 15),
    (15, 16),
    (16, 17),
    (17, 18),
]

BONDS_CG = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)]

TYPE_COLORS = {
    1: "#42A5F5",  # CG type A — bright blue
    2: "#1E88E5",  # CG type B — blue
    3: "#FF5722",  # AT CH3 — orange-red
    4: "#FF9800",  # AT CH2 — amber
}

CG_ALPHA = 0.45
CG_SIZE = {1: 500, 2: 450}
AT_SIZE = {3: 100, 4: 80}

TYPE_LABELS = {1: "CG-A", 2: "CG-B", 3: "CH₃", 4: "CH₂"}


def parse_dump(path: Path) -> list[Frame]:
    """Parse LAMMPS custom dump: id mol type x y z f_bm."""
    frames: list[Frame] = []
    with open(path) as fh:
        while True:
            line = fh.readline()
            if not line:
                break
            if "ITEM: TIMESTEP" in line:
                timestep = int(fh.readline().strip())
                fh.readline()  # ITEM: NUMBER OF ATOMS
                natoms = int(fh.readline().strip())
                fh.readline()  # ITEM: BOX BOUNDS
                for _ in range(3):
                    fh.readline()
                fh.readline()  # ITEM: ATOMS ...

                ids = np.zeros(natoms, dtype=int)
                types = np.zeros(natoms, dtype=int)
                x = np.zeros(natoms)
                y = np.zeros(natoms)
                z = np.zeros(natoms)
                lam = np.zeros(natoms)

                for i in range(natoms):
                    tokens = fh.readline().split()
                    ids[i] = int(tokens[0])
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
                        types=types[order],
                        x=x[order],
                        y=y[order],
                        z=z[order],
                        lam=lam[order],
                    )
                )
    return frames


def compute_frame_limits(frame: Frame, pad: float = 5.0) -> tuple[tuple, tuple, tuple]:
    """Compute view limits centered on the molecule in this frame."""
    cx, cy, cz = frame.x.mean(), frame.y.mean(), frame.z.mean()
    span = max(np.ptp(frame.x), np.ptp(frame.y), np.ptp(frame.z)) / 2.0 + pad
    return (
        (cx - span, cx + span),
        (cy - span, cy + span),
        (cz - span, cz + span),
    )


def draw_frame(ax: Axes3D, frame: Frame) -> None:
    ax.cla()

    xlim, ylim, zlim = compute_frame_limits(frame)
    id_to_idx = {aid: i for i, aid in enumerate(frame.ids)}
    mean_lam = frame.lam.mean()

    # --- Draw CG beads first (always large, semi-transparent) ---
    for atom_type in [1, 2]:
        mask = frame.types == atom_type
        if not mask.any():
            continue

        base_color = matplotlib.colors.to_rgba(TYPE_COLORS[atom_type])
        n = mask.sum()
        colors = np.array([(*base_color[:3], CG_ALPHA)] * n)

        ax.scatter(
            frame.x[mask],
            frame.y[mask],
            frame.z[mask],
            c=colors,
            s=CG_SIZE[atom_type],
            edgecolors=(*base_color[:3], 0.7),
            linewidths=0.8,
            depthshade=True,
            zorder=1,
        )

    # --- Draw AT atoms (opacity grows with lambda) ---
    for atom_type in [3, 4]:
        mask = frame.types == atom_type
        if not mask.any():
            continue

        lam_vals = frame.lam[mask]
        alphas = np.clip(lam_vals, 0.5, 1.0)

        base_color = matplotlib.colors.to_rgba(TYPE_COLORS[atom_type])
        colors = np.array([(*base_color[:3], a) for a in alphas])
        sizes = AT_SIZE[atom_type] * np.clip(alphas, 0.5, 1.0)

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

    # --- CG bonds: always visible, semi-transparent ---
    for i_id, j_id in BONDS_CG:
        if i_id not in id_to_idx or j_id not in id_to_idx:
            continue
        i, j = id_to_idx[i_id], id_to_idx[j_id]
        ax.plot(
            [frame.x[i], frame.x[j]],
            [frame.y[i], frame.y[j]],
            [frame.z[i], frame.z[j]],
            color=(*matplotlib.colors.to_rgb("#1E88E5"), CG_ALPHA),
            linewidth=2.5,
            zorder=0,
        )

    # --- AT bonds: grow in with lambda ---
    for i_id, j_id in BONDS_AT:
        if i_id not in id_to_idx or j_id not in id_to_idx:
            continue
        i, j = id_to_idx[i_id], id_to_idx[j_id]
        alpha_bond = max(0.5, mean_lam)
        ax.plot(
            [frame.x[i], frame.x[j]],
            [frame.y[i], frame.y[j]],
            [frame.z[i], frame.z[j]],
            color=(*matplotlib.colors.to_rgb("#E65100"), alpha_bond),
            linewidth=1.8 * max(0.3, alpha_bond),
            zorder=1,
        )

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_zlim(*zlim)
    ax.set_xlabel("x (Å)", fontsize=9)
    ax.set_ylabel("y (Å)", fontsize=9)
    ax.set_zlabel("z (Å)", fontsize=9)
    ax.tick_params(labelsize=7)

    phase = (
        "CG equilibration"
        if mean_lam < 0.01
        else ("Backmapping" if mean_lam < 0.99 else "AT production")
    )
    ax.set_title(
        f"Dodecane backmapping — step {frame.timestep}\nλ = {mean_lam:.3f}  |  {phase}",
        fontsize=12,
        fontweight="bold",
    )

    ax.view_init(elev=25, azim=(30 + frame.timestep * 0.04) % 360)


def make_movie(
    frames: list[Frame],
    output: Path,
    fps: int = 15,
    dpi: int = 150,
) -> None:
    fig = plt.figure(figsize=(10, 8), facecolor="white")
    ax = fig.add_subplot(111, projection="3d")

    draw_frame(ax, frames[0])
    fig.tight_layout()

    def update(i: int):
        draw_frame(ax, frames[i])

    anim = animation.FuncAnimation(
        fig, update, frames=len(frames), interval=1000 // fps
    )

    if output.suffix == ".gif":
        writer = animation.PillowWriter(fps=fps)
    else:
        writer = animation.FFMpegWriter(fps=fps, codec="libx264", bitrate=2000)

    print(f"Rendering {len(frames)} frames to {output} ...")
    anim.save(str(output), writer=writer, dpi=dpi)
    print(f"Done: {output}")
    plt.close(fig)


def save_png_sequence(
    frames: list[Frame],
    output_dir: Path,
    dpi: int = 150,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(10, 8), facecolor="white")
    ax = fig.add_subplot(111, projection="3d")

    for i, frame in enumerate(frames):
        draw_frame(ax, frame)
        fig.tight_layout()
        out = output_dir / f"frame_{i:04d}.png"
        fig.savefig(out, dpi=dpi, bbox_inches="tight", facecolor="white")
        if (i + 1) % 20 == 0:
            print(f"  saved frame {i + 1}/{len(frames)}")

    plt.close(fig)
    print(f"PNG sequence saved to {output_dir}/")
    print(
        f"Combine with ffmpeg: ffmpeg -framerate 15 -i {output_dir}/frame_%04d.png -c:v libx264 -pix_fmt yuv420p backmap_movie.mp4"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Visualize backmapping trajectory as movie"
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
            f"ERROR: {dump_path} not found. Run LAMMPS first:\n  lmp -in in.dodecane_single_chain"
        )
        raise SystemExit(1)

    frames = parse_dump(dump_path)
    print(f"Loaded {len(frames)} frames from {dump_path}")

    if not frames:
        print("No frames found.")
        raise SystemExit(1)

    lam_first = frames[0].lam.mean()
    lam_last = frames[-1].lam.mean()
    print(f"Lambda range: {lam_first:.4f} → {lam_last:.4f}")

    if args.png_sequence:
        save_png_sequence(frames, Path("frames"), dpi=args.dpi)
    else:
        make_movie(frames, Path(args.output), fps=args.fps, dpi=args.dpi)


if __name__ == "__main__":
    main()
