# Molecule-aware PBC unwrap for network LAMMPS export

## Why

Rim135 Tier B LAMMPS runs fail at ~1000 steps with `Bond atoms missing` and startup
warnings about **114 Å** bond communication estimates. Root cause: hybrid export
unwraps the full bond graph globally, then zeroes image flags while leaving **95 Å**
Euclidean bonded spans. Min-image bond length passes validation (< 20 Å) but LAMMPS
bonded ghost communication requires short **Cartesian** separations with `ix=iy=iz=0`.

## What changes

- Unwrap **per `mol_id`** (CG bead + mapped AT fragment), then translate whole
  molecules to shorten inter-mol (cross) bonds.
- Gate export on **max Euclidean bond length** ≤ comm cutoff budget, not min-image alone.
- Keep default export path: unwrapped Cartesian coords, `ix=iy=iz=0` (no image flags).

## Non-goals

- Image-flag-based export as default (fragile for network topologies).
- PR3 dynamics protocol changes (equilibration staging, CG ramp).
