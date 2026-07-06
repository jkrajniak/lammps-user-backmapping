# Design: rim135 Tier B bakery protocol (PR4a)

## Protocol (settings.v2.yaml)

| Knob | PR3 (soft-start) | PR4 (bakery) | Rationale |
|------|------------------|--------------|-----------|
| `alpha` | 0.00005 | 0.0001 | Bakery/E++ reference (10k ramp steps) |
| `timestep_backmapping` | 0.0005 ps | 0.001 ps | Revert to 1 fs during ramp |
| `timestep` | 0.001 ps | 0.001 ps | Unchanged |
| `temperature` | 300.0 K | 298.0 K | Bakery AA reference |
| `thermostat_gamma` | 0.5 (default) | 15.0 | Bakery Langevin γ (ps⁻¹) |
| `cap_force` | null | 50000.0 | Bakery force cap kJ/(mol·nm) |
| `equilibration_steps` | 0 | 0 | CG pre-equilibrated |
| `production_steps` | 10000 | 10000 | Unchanged |

## CapForce unit conversion

Bakery stores cap_force in GROMACS units: **kJ/(mol·nm)**.

LAMMPS `real` units expect **kcal/(mol·Å)**.

```
fmax = cap_force × units.FORCE
     = cap_force × (KJ_TO_KCAL / NM_TO_ANGSTROM)
     = 50000 × 0.239006 / 10
     = 1195.03 kcal/(mol·Å)
```

Implementation: `units.force(sim.cap_force)` in `_write_cap_force`.

## Writer ordering (_write_setup)

1. Groups + `neigh_modify`
2. **`velocity all create`** (new)
3. Integration fixes (NVE + Langevin)
4. `fix bm all backmap`
5. `fix pairs` (if cross pairs)
6. **`fix cap all backmap/capforce`** (new, when cap_force set)
7. `compute at_temp`, `timestep`, thermo, dump

## Langevin damping

`thermostat_gamma` is in ps⁻¹. Writer converts to LAMMPS damp parameter (fs):

```
damp_fs = units.time(1.0 / thermostat_gamma)
        = 1000 / 15 ≈ 66.7 fs
```

## Tier B pass criteria (PR4 protocol)

1. Generated input includes velocity create, capforce fix, langevin 298 K, damp ≈66.7 fs
2. Backmapping at `timestep 1.00`, `run 10000`
3. LAMMPS smoke (PR4b+) completes ramp without crash

## Verification ladder

1. `pytest` — input structure gates (this PR)
2. Local truncated smoke after PR4b C++ lands
3. VM full `in.rim135` when local passes
