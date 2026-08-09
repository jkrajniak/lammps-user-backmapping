# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "matplotlib"]
# ///
"""Render four static snapshots of the real single-chain dodecane backmapping
trajectory (CG-only, hybrid-build, mid-ramp, AT-only) for use as Fig. 1 insets
in the CPC manuscript, replacing the paper's earlier abstract dot icons with
an actual coarse-grained chain and its atomistic counterpart.

Reuses the parsing/coloring conventions of visualize_backmap_movie.py in the
same directory. Fixes the camera angle across all four snapshots (unlike the
movie, which slowly rotates) so the four panels read as one continuous story.
Frame indices were chosen from a scan of end-to-end/contour-length ratios so
each snapshot shows a visibly bent (not near-straight) real conformation --
see the module docstring notes below for the specific frames.

Usage:
  uv run render_fig1_snapshots.py [--output-dir DIR]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import (
    Axes3D,  # noqa: TC002 (needed at runtime to register the 3d projection)
)
from visualize_backmap_movie import (
    BONDS_AT,
    BONDS_CG,
    Frame,
    compute_frame_limits,
    parse_dump,
)

VIEW_ELEV = 18
VIEW_AZIM = 52

# Overrides visualize_backmap_movie.TYPE_COLORS (blue/orange, tuned for the
# standalone movie) with the paper's own CG/AT legend colors, so these
# snapshots match Fig. 1's legend swatches and its CG->AT lambda-ramp
# gradient bar instead of clashing with the LAMMPS/backmap-prep box-border
# colors (which are also blue/orange in generate_fig1_workflow.py).
PAPER_C_CG = "#8172B3"
PAPER_C_AT = "#55A868"
TYPE_COLORS = {
    1: PAPER_C_CG,  # CG type A
    2: "#6C5B9E",  # CG type B (darker shade for bead separation)
    3: PAPER_C_AT,  # AT CH3
    4: "#7FB96A",  # AT CH2 (lighter shade for atom-type separation)
}

# frame_idx chosen by scanning end-to-end/contour ratio across the trajectory
# for a visibly bent (not near-straight) real conformation; all are genuine
# simulation snapshots, not synthetic geometry.
STAGES = {
    "cg": {"frame_idx": 240, "show_cg": True, "show_at": False, "at_min_alpha": 0.0},
    "hybrid": {"frame_idx": 240, "show_cg": True, "show_at": True, "at_min_alpha": 0.35},
    "ramp": {"frame_idx": 450, "show_cg": True, "show_at": True, "at_min_alpha": 0.5},
    "at": {"frame_idx": 560, "show_cg": False, "show_at": True, "at_min_alpha": 1.0},
}


def draw_snapshot(
    ax: Axes3D,
    frame: Frame,
    show_cg: bool,
    show_at: bool,
    at_min_alpha: float,
) -> None:
    ax.cla()
    xlim, ylim, zlim = compute_frame_limits(frame, pad=2.0)
    id_to_idx = {aid: i for i, aid in enumerate(frame.ids)}
    mean_lam = frame.lam.mean()

    if show_cg:
        for atom_type in [1, 2]:
            mask = frame.types == atom_type
            if not mask.any():
                continue
            base_color = matplotlib.colors.to_rgba(TYPE_COLORS[atom_type])
            n = mask.sum()
            colors = np.array([(*base_color[:3], 0.55)] * n)
            ax.scatter(
                frame.x[mask],
                frame.y[mask],
                frame.z[mask],
                c=colors,
                s=650,
                edgecolors=(*base_color[:3], 0.85),
                linewidths=1.2,
                depthshade=True,
                zorder=1,
            )
        for i_id, j_id in BONDS_CG:
            if i_id not in id_to_idx or j_id not in id_to_idx:
                continue
            i, j = id_to_idx[i_id], id_to_idx[j_id]
            ax.plot(
                [frame.x[i], frame.x[j]],
                [frame.y[i], frame.y[j]],
                [frame.z[i], frame.z[j]],
                color=(*matplotlib.colors.to_rgb(PAPER_C_CG), 0.6),
                linewidth=4.0,
                zorder=0,
            )

    if show_at:
        for atom_type in [3, 4]:
            mask = frame.types == atom_type
            if not mask.any():
                continue
            base_color = matplotlib.colors.to_rgba(TYPE_COLORS[atom_type])
            n = mask.sum()
            alpha = max(at_min_alpha, min(1.0, mean_lam))
            colors = np.array([(*base_color[:3], alpha)] * n)
            ax.scatter(
                frame.x[mask],
                frame.y[mask],
                frame.z[mask],
                c=colors,
                s=140,
                edgecolors="none",
                depthshade=True,
                zorder=2,
            )
        for i_id, j_id in BONDS_AT:
            if i_id not in id_to_idx or j_id not in id_to_idx:
                continue
            i, j = id_to_idx[i_id], id_to_idx[j_id]
            alpha_bond = max(at_min_alpha, min(1.0, mean_lam))
            ax.plot(
                [frame.x[i], frame.x[j]],
                [frame.y[i], frame.y[j]],
                [frame.z[i], frame.z[j]],
                color=(*matplotlib.colors.to_rgb(PAPER_C_AT), alpha_bond),
                linewidth=2.6 * max(0.4, alpha_bond),
                zorder=1,
            )

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_zlim(*zlim)
    ax.set_axis_off()
    ax.view_init(elev=VIEW_ELEV, azim=VIEW_AZIM)
    ax.set_box_aspect((1, 1, 1))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory to write fig1_snapshot_{cg,hybrid,ramp,at}.png",
    )
    args = parser.parse_args()

    dump_path = Path(__file__).resolve().parent / "dump.movie.lammpstrj"
    frames = parse_dump(dump_path)
    print(f"Loaded {len(frames)} frames from {dump_path}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, cfg in STAGES.items():
        fig = plt.figure(figsize=(3, 3), facecolor="none")
        ax = fig.add_subplot(111, projection="3d")
        ax.patch.set_alpha(0.0)
        frame = frames[cfg["frame_idx"]]
        draw_snapshot(
            ax,
            frame,
            show_cg=cfg["show_cg"],
            show_at=cfg["show_at"],
            at_min_alpha=cfg["at_min_alpha"],
        )
        out = args.output_dir / f"fig1_snapshot_{name}.png"
        fig.savefig(out, dpi=300, bbox_inches="tight", pad_inches=0.02, transparent=True)
        plt.close(fig)
        print(
            f"Wrote {out} (frame {cfg['frame_idx']}, step {frame.timestep}, "
            f"mean_lam={frame.lam.mean():.3f})"
        )


if __name__ == "__main__":
    main()
