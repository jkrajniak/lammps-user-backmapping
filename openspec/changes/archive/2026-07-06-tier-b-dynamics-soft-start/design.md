# Design: rim135 Tier B soft-start

## Protocol (settings.v2.yaml)

| Knob | Before | After | Rationale |
|------|--------|-------|-----------|
| `equilibration_steps` | 0 | 0 | CG melt pre-equilibrated; in-hybrid λ=0 equil **heated** AT to ~526 K in smoke |
| `timestep_backmapping` | 0.001 ps | 0.0005 ps | Halve dt during ramp |
| `alpha` | 0.0001 | 0.00005 | Doubles ramp length (20k steps) |
| `production_steps` | 10000 | 10000 | Unchanged |

## Writer behavior (no code change)

Existing `write_lammps_input` emits backmapping at `timestep_backmapping`, then optional production at `timestep`.

## Tier B pass criteria

1. ≥5000 λ-ramp steps without crash
2. PE finite in first 1000 ramp steps (no 1e4+ spike)
3. No `Bond atoms missing`
4. T ≈ 300 K (order-of-magnitude)

## Verification ladder

1. `pytest` — input structure gates
2. Local truncated smoke (≥6000 total dynamics steps)
3. VM full `in.rim135` when local passes
