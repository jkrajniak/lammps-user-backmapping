# Single Dodecane Chain — Backmapping Verification & Movie

Visual verification of the CG → AT backmapping process on one dodecane molecule (C₁₂H₂₆).

## Structure

| CG bead | Type | AT atoms | AT types |
|---------|------|----------|----------|
| A1      | 1    | C1, C2   | CH₃, CH₂ |
| B1      | 2    | C3, C4   | CH₂, CH₂ |
| B2      | 2    | C5, C6   | CH₂, CH₂ |
| B3      | 2    | C7, C8   | CH₂, CH₂ |
| B4      | 2    | C9, C10  | CH₂, CH₂ |
| A2      | 1    | C11, C12 | CH₂, CH₃ |

## Quick Start

```bash
# 1. Run the single-chain LAMMPS simulation
lmp -in in.dodecane_single_chain

# 2a. Generate movie with matplotlib (no extra tools needed)
uv run visualize_backmap_movie.py

# 2b. Or with OVITO (higher quality, needs ovito package)
uv run visualize_backmap_ovito.py

# 2c. Or open interactively in VMD
vmd -e visualize_backmap.vmd
```

## Simulation Phases

| Phase | Steps | Lambda | Dump interval |
|-------|-------|--------|---------------|
| CG equilibration | 0–500 | 0 (frozen) | every 10 |
| Backmapping | 500–1500 | 0 → 1 (α=0.001) | every 10 |
| AT production | 1500–2000 | 1 (frozen) | every 10 |

Total: 200 frames covering the full backmapping transition.

## Visualization Scripts

### `visualize_backmap_movie.py` (matplotlib)
- Self-contained, only needs numpy + matplotlib
- 3D scatter plot with rotating camera
- CG beads fade out (opacity ∝ 1−λ), AT atoms fade in (opacity ∝ λ)
- Output: `backmap_movie.mp4` or PNG sequence

```bash
uv run visualize_backmap_movie.py --output backmap.gif     # animated GIF
uv run visualize_backmap_movie.py --png-sequence            # PNG frames
uv run visualize_backmap_movie.py --fps 30 --dpi 200        # high quality
```

### `visualize_backmap_ovito.py` (OVITO)
- Uses OVITO Python API for ray-traced rendering
- Color-codes atoms by lambda (Viridis colormap)
- CG atoms become transparent as λ → 1
- Output: MP4 or PNG frames

### `visualize_backmap.vmd` (VMD)
- Interactive: scrub through frames with the VMD slider
- CG beads shown as large blue spheres, AT atoms as smaller orange/yellow
- `render_movie backmap_frame` to render TGA frames

## What to Look For

1. **Phase 1** (CG equil.): Only CG beads visible, AT atoms at zero opacity
2. **Phase 2** (Backmapping):
   - CG beads gradually fade while AT atoms emerge
   - AT atoms initially overlap CG bead positions (COM constraint)
   - Bond lengths and angles relax to equilibrium values
   - No unphysical atom overlaps or explosions
3. **Phase 3** (AT production): Only AT atoms visible, correct dodecane chain topology
