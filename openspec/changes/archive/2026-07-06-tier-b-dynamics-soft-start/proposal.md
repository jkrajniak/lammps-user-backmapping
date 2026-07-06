# Tier B dynamics soft-start for rim135

## Why

After PR1 (PBC export), rim135 LAMMPS still crashes at ~1000 λ-ramp steps with runaway
temperature and `Bond atoms missing` (3763–12623). Static topology and min-image bond
lengths are acceptable; the failure is **dynamics** during the λ ramp from a stiff hybrid
start state.

## What changes

- Enable **gentler λ ramp** via `alpha: 0.00005` and `timestep_backmapping: 0.0005` ps in
  `examples/epoxy/settings.v2.yaml`.
- Keep `equilibration_steps: 0` (CG pre-equilibrated; optional equil heated AT in smoke).
- Integration test asserts soft-start timestep and ramp length in generated `in.rim135`.

## Non-goals

- New writer phases or schema fields (existing `equilibration_steps` / `timestep_backmapping`).
- CG-only thermostat mode or multi-stage alpha schedules.
- Tier C RDF validation.
