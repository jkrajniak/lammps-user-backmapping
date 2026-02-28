# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "ovito>=3.10",
# ]
# ///
"""Render backmapping movie using OVITO Python API.

OVITO handles lambda-based coloring natively via the custom property f_bm.
Produces high-quality ray-traced frames or an MP4 movie.

Usage:
    uv run visualize_backmap_ovito.py [--input dump.movie.lammpstrj] [--output backmap_movie.mp4]
    uv run visualize_backmap_ovito.py --frames-only   # PNG frames only
"""

from __future__ import annotations

import argparse
from pathlib import Path

from ovito.io import import_file
from ovito.modifiers import (
    ColorCodingModifier,
    ComputePropertyModifier,
)
from ovito.vis import ParticlesVis, Viewport


def setup_pipeline(dump_path: str):
    pipeline = import_file(dump_path, multiple_frames=True)

    # Color AT atoms by lambda: transparent (low lambda) → opaque (high lambda)
    # Use f_bm as the lambda property
    pipeline.modifiers.append(
        ColorCodingModifier(
            property="f_bm",
            start_value=0.0,
            end_value=1.0,
            gradient=ColorCodingModifier.Viridis(),
        )
    )

    # Set particle radii based on type
    pipeline.modifiers.append(
        ComputePropertyModifier(
            output_property="Radius",
            expressions=[
                "ParticleType==1 ? 1.8 : (ParticleType==2 ? 1.6 : (ParticleType==3 ? 0.9 : 0.7))"
            ],
        )
    )

    # Transparency: CG atoms fade out (alpha = 1 - lambda), AT atoms fade in (alpha = lambda)
    pipeline.modifiers.append(
        ComputePropertyModifier(
            output_property="Transparency",
            expressions=["ParticleType<=2 ? f_bm : (1.0 - f_bm)"],
        )
    )

    vis = pipeline.source.data.particles.vis
    if isinstance(vis, ParticlesVis):
        vis.radius = 0.5
        vis.shape = ParticlesVis.Shape.Sphere

    return pipeline


def render_movie(
    pipeline, output: Path, fps: int = 15, res: tuple[int, int] = (1280, 960)
):
    pipeline.add_to_scene()

    vp = Viewport(type=Viewport.Type.Perspective)
    vp.zoom_all()
    vp.camera_dir = (1, -0.5, -0.3)

    num_frames = pipeline.source.num_frames
    print(f"Rendering {num_frames} frames at {res[0]}x{res[1]}...")

    if output.suffix in (".mp4", ".avi"):
        vp.render_anim(
            filename=str(output),
            size=res,
            fps=fps,
            every_nth=1,
        )
    else:
        output.mkdir(parents=True, exist_ok=True)
        for frame_idx in range(num_frames):
            fname = output / f"frame_{frame_idx:04d}.png"
            vp.render_image(
                filename=str(fname),
                size=res,
                frame=frame_idx,
            )
            if (frame_idx + 1) % 20 == 0:
                print(f"  Frame {frame_idx + 1}/{num_frames}")

    pipeline.remove_from_scene()
    print(f"Done: {output}")


def main() -> None:
    parser = argparse.ArgumentParser(description="OVITO backmapping movie renderer")
    parser.add_argument("--input", "-i", default="dump.movie.lammpstrj")
    parser.add_argument("--output", "-o", default="backmap_movie.mp4")
    parser.add_argument("--fps", type=int, default=15)
    parser.add_argument("--width", type=int, default=1280)
    parser.add_argument("--height", type=int, default=960)
    parser.add_argument(
        "--frames-only", action="store_true", help="Save PNG frames to directory"
    )
    args = parser.parse_args()

    dump_path = Path(args.input)
    if not dump_path.exists():
        print(f"ERROR: {dump_path} not found. Run LAMMPS first.")
        raise SystemExit(1)

    pipeline = setup_pipeline(str(dump_path))

    output = Path("frames") if args.frames_only else Path(args.output)
    render_movie(pipeline, output, fps=args.fps, res=(args.width, args.height))


if __name__ == "__main__":
    main()
