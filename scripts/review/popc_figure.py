"""Assemble the POPC backmapping transition figure from the panel renders.

Panels (same camera, rendered by popc_render.py): CG (lambda = 0), mid-ramp,
end of the ramp (lambda = 1), and the atomistic bilayer after the continued
production. All panels are cropped to one common bounding box, so they keep a
common scale. Legend: MARTINI 3 bead groups and atom elements (a bead and the
atoms it becomes share a colour), water.

    uv run --with matplotlib --with numpy --with pillow popc_figure.py \\
        --panels popc_lambda0.00.png popc_lambda0.49.png popc_lambda1.00.png popc_gro.png \\
        --labels "CG, $\\lambda$ = 0" "$\\lambda$ = 0.49" "$\\lambda$ = 1" "AT, 10 ns later" \\
        --out popc_transition
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

BLUE, ORANGE, VERMILLION = "#0072B2", "#E69F00", "#D55E00"
GREEN, SKY = "#009E73", "#56B4E9"


def common_bbox(images: list[np.ndarray], pad: int = 12) -> tuple[int, int, int, int]:
    rows, cols = [], []
    for img in images:
        ink = np.any(img[..., :3] < 0.97, axis=2)
        r, c = np.where(ink)
        rows += [r.min(), r.max()]
        cols += [c.min(), c.max()]
    h, w = images[0].shape[:2]
    return (
        max(min(rows) - pad, 0),
        min(max(rows) + pad, h),
        max(min(cols) - pad, 0),
        min(max(cols) + pad, w),
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--panels", nargs="+", type=Path, required=True)
    ap.add_argument("--labels", nargs="+", required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--width", type=float, default=7.0, help="inches (double column)")
    args = ap.parse_args()
    if len(args.panels) != len(args.labels):
        print("one label per panel", file=sys.stderr)
        return 1

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from PIL import Image

    images = [np.asarray(Image.open(p).convert("RGB"), dtype=float) / 255.0 for p in args.panels]
    r0, r1, c0, c1 = common_bbox(images)
    crops = [img[r0:r1, c0:c1] for img in images]
    aspect = (r1 - r0) / (c1 - c0)
    n = len(crops)
    panel_w = args.width / n
    fig = plt.figure(figsize=(args.width, panel_w * aspect + 0.75))
    gs = fig.add_gridspec(2, n, height_ratios=[panel_w * aspect, 0.55], hspace=0.05, wspace=0.03)
    for k, (img, label) in enumerate(zip(crops, args.labels, strict=True)):
        ax = fig.add_subplot(gs[0, k])
        ax.imshow(img, interpolation="lanczos")
        ax.set_axis_off()
        ax.set_title(f"({'abcdefgh'[k]}) {label}", fontsize=8, loc="left", pad=3)
    beads = [
        Patch(color=BLUE, label="choline (Q1) / N"),
        Patch(color=ORANGE, label="phosphate (Q5) / P"),
        Patch(color=VERMILLION, label="glycerol, ester (N4a, SN4a) / O"),
        Patch(color="#B3B3B3", label="tail (C1) / C"),
        Patch(color=GREEN, label="unsaturated tail (C4h)"),
        Patch(color="#F5F5F5", ec="#999999", lw=0.5, label="H"),
        Patch(color=SKY, alpha=0.5, label="water (W / TIP3P O)"),
    ]
    lax = fig.add_subplot(gs[1, :])
    lax.set_axis_off()
    lax.legend(
        handles=beads,
        loc="center",
        ncol=4,
        fontsize=7,
        frameon=False,
        title="MARTINI 3 bead / atom",
        title_fontsize=7,
        handlelength=1.2,
        columnspacing=1.2,
    )
    for ext in ("pdf", "png"):
        fig.savefig(f"{args.out}.{ext}", dpi=400, bbox_inches="tight")
    meta = {
        "command": " ".join(sys.argv),
        "panels": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in args.panels},
        "generated_at": datetime.now(timezone.utc).isoformat(),
    }
    Path(f"{args.out}.meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"wrote {args.out}.pdf/.png ({n} panels, crop rows {r0}-{r1}, cols {c0}-{c1})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
