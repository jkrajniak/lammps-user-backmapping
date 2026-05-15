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

## RDF comparison (backmapped vs reference AT)

From this directory, after `backmap-prep` (if you regenerate inputs):

```bash
lmp -in in.dodecane
uv run extract_at_frame.py dodecane_hybrid.data dodecane_at.data --cg-types 1 2
lmp -in in.dodecane_at
lmp -in in.dodecane_at_ref    # independent reference, same box and 10 molecules
uv run compare_rdf.py --backmap rdf_backmap.dat --reference rdf_reference.dat --plot rdf_comparison.png
```

`in.dodecane_at` uses NVT equilibration at 298 K before the RDF window. The
reference run is longer; with only 10 molecules, `compare_rdf.py` tolerances may
still fail on RMSD or CH₃–CH₃ peak position — increase system size or sampling
for tighter agreement.

## Input Files

| File | Description |
|------|-------------|
| `in.dodecane` | Backmapping only: λ 0→1, `write_data dodecane_hybrid.data` |
| `in.dodecane_at` | Pure AT run after extraction; writes `rdf_backmap.dat` |
| `in.dodecane_at_ref` | Independent AT reference; writes `rdf_reference.dat` |
| `extract_at_frame.py` | Strip CG atoms from hybrid data → `dodecane_at.data` |
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

## Medium-scale (250 molecules)

`n250/` builds a **250-molecule** melt (subset of the bakery `large/cg_conf.gro`)
with `topol_aa_250.top` / `topol_cg_250.top`. See
[`n250/README.md`](n250/README.md) and `./n250/prepare_inputs.sh`.

## Large-scale variant

A larger system is available in `large/`. Topologies and tables in `large/` are the same as in the root (backmap-prep compatible); only the CG coordinates differ. To use bakery’s full-size coordinates:

1. Copy `cg_conf.gro` from `bakery/examples/dodecane/` into `examples/dodecane/large/` (or use the one already there if you copied from bakery).
2. From repo root: `uv run backmap-prep examples/dodecane/large/settings.yaml`
3. Run LAMMPS from `examples/dodecane/large/`: `lmp -in in.dodecane`

The large-scale input (`large/in.dodecane`) uses the robust multi-phase
protocol: energy minimisation with CG frozen, `nve/limit` relaxation at
dt = 0.01 fs, `nve/limit` lambda ramp at dt = 0.10 fs, then gradual NVT
equilibration (0.25 → 0.50 → 1.00 fs). This prevents "Bond atoms missing"
errors that occur with aggressive timesteps in large systems.

See `large/README.md` for details and [Large-scale examples](../../docs/large-scale-examples.md) in the docs.
