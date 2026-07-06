# Tier B bakery protocol alignment for rim135 (PR4a)

## Why

PR3 soft-start (0.5 fs, α=0.00005) improved stability marginally but did not match the
bakery/E++ epoxy AA reference protocol. The legacy reference uses cap_force capping,
γ=15 thermostat, Maxwell velocity initialization, and 1 fs backmapping timestep with
α=0.0001 (10k ramp steps).

## What changes

- Wire `cap_force` from settings into `fix cap all backmap/capforce` (unit-converted).
- Emit `velocity all create` before integration fixes.
- Honor `thermostat_gamma: 15.0` for Langevin damping (≈66.7 fs).
- Revert rim135 `settings.v2.yaml` to bakery-aligned dt/alpha/temperature knobs.
- Integration tests assert PR4 protocol gates in generated `in.rim135`.

## Non-goals

- LAMMPS C++ `fix backmap/capforce` implementation (PR4b).
- CG-only thermostat mode or multi-stage alpha schedules.
- Tier C RDF validation.
