## MODIFIED Requirements

### Requirement: Simulation parameters section

The `simulation` section SHALL configure backmapping parameters:

```yaml
simulation:
  alpha: 0.001
  initial_resolution: 0.0
  nonuniform_lambda: false
  timestep: 0.001
  timestep_backmapping: 0.001
  equilibration_steps: 10000
  production_steps: 10000
  temperature: 423.0
  thermostat: langevin
  thermostat_gamma: 0.5
  thermostat_target: atomistic
  lj_cutoff: 1.2
  cg_cutoff: 1.4
  coulomb_cutoff: 0.9
  table_groups: [WCG]
  exclusion_nrexcl: 3
  energy_interval: 1000
  trajectory_interval: 1000
  rng_seed: -1
  restart_interval: null      # NEW: steps between restart checkpoints (null = disabled)
```

When `restart_interval` is set to an integer N, the generator SHALL emit
LAMMPS `restart N file1 file2` commands in the generated input scripts and
generate per-phase scripts for restart support. When `null` (default), no
restart commands are generated and only the master input script is produced.

All parameters in GROMACS units. The generator converts to LAMMPS `real` units
in the output.

#### Scenario: Default parameters

- **WHEN** the `simulation` section omits optional fields
- **THEN** the generator SHALL use the documented default values, including
  `restart_interval: null`

#### Scenario: Parameter validation

- **WHEN** `alpha` is negative or `temperature` is zero
- **THEN** the generator SHALL abort with a validation error

#### Scenario: Restart interval set

- **WHEN** `restart_interval: 5000` is specified
- **THEN** the generated input scripts SHALL contain `restart 5000` commands
  and per-phase scripts SHALL be generated

#### Scenario: Restart interval validation

- **WHEN** `restart_interval` is set to a non-positive integer
- **THEN** the generator SHALL abort with a validation error
