# 10-Chain Dodecane — Backmapping with Movie

Backmapping of 10 dodecane (C₁₂H₂₆) molecules from CG to atomistic resolution,
with trajectory visualization.

## System

- **10 chains**, 180 atoms total (60 CG beads + 120 AT atoms)
- Box: 57.5 x 57.5 x 57.5 A (periodic)
- CG model: 6 beads per chain (A1, B1-B4, A2), tabulated IBI interactions
- AT model: GROMOS united-atom (CH₃ / CH₂)

Per-chain mapping (same as single_chain):

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
# 1. Run LAMMPS (movie input with frequent dumps)
lmp -in in.dodecane_movie

# 2. Generate movie
uv run visualize_backmap_movie.py
```

## Input Files

| File | Description |
|------|-------------|
| `in.dodecane` | Production input (dump every 1000 steps, 30k total) |
| `in.dodecane_movie` | Movie input (dump every 50 steps, 6k total, ~120 frames) |
| `dodecane.data` | LAMMPS data file (10 chains, 180 atoms) |
| `settings.yaml` | backmap-prep configuration |
| `table_b1.table` | Tabulated CG bond potential |
| `visualize_backmap_movie.py` | matplotlib movie generator |

## Simulation Phases (movie input)

| Phase | Steps | Lambda | Dump interval |
|-------|-------|--------|---------------|
| CG equilibration | 0-3000 | 0 (frozen) | every 50 |
| Backmapping | 3000-5000 | 0 → 1 (alpha=0.0005) | every 50 |
| AT production | 5000-6000 | 1 (frozen) | every 50 |

Total: ~120 frames covering the full backmapping transition.

## Visualization

```bash
# MP4 movie (default)
uv run visualize_backmap_movie.py

# Animated GIF
uv run visualize_backmap_movie.py --output backmap.gif

# PNG frame sequence
uv run visualize_backmap_movie.py --png-sequence

# High quality
uv run visualize_backmap_movie.py --fps 30 --dpi 200
```

Features:
- CG beads fade out as lambda → 1, AT atoms fade in
- AT bonds colored per molecule (10 distinct colors)
- Rotating camera view of the full simulation box
- Phase and lambda annotation overlay

## What to Look For

1. **Phase 1** (CG equil.): Large blue CG beads visible, AT atoms transparent
2. **Phase 2** (Backmapping):
   - CG beads shrink and fade while AT atoms emerge
   - AT atoms initially overlap CG bead COM positions
   - Bond lengths and angles relax to equilibrium
   - 10 separate chains should remain intact (no cross-linking)
3. **Phase 3** (AT production): Only AT atoms visible, 10 correct dodecane chains
