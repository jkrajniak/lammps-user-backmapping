## MODIFIED Requirements

### Requirement: Simulation section — dynamics parameters

The generator SHALL honor `simulation.equilibration_steps`, `simulation.timestep_backmapping`,
`simulation.alpha`, `simulation.thermostat_gamma`, and `simulation.cap_force` when writing
LAMMPS input scripts for network hybrids.

For rim135 settings v2, the recommended Tier B bakery-aligned protocol SHALL use:

- `equilibration_steps: 0` (CG pre-equilibrated)
- `timestep_backmapping: 0.001` (ps, 1 fs)
- `alpha: 0.0001` (10k ramp steps)
- `temperature: 298.0` (K)
- `thermostat_gamma: 15.0` (ps⁻¹)
- `cap_force: 50000.0` (kJ/(mol·nm))

#### Scenario: Bakery-aligned rim135 input

- **WHEN** `settings.v2.yaml` is built with PR4 dynamics knobs
- **THEN** `in.rim135` SHALL use `timestep 1.00` and `run 10000` for the λ ramp

## ADDED Requirements

### Requirement: Initial velocity creation

The generator SHALL emit a Maxwell-Boltzmann velocity initialization before integration
fixes when writing the setup block.

#### Scenario: Velocity create at simulation temperature

- **WHEN** any backmapping input is generated
- **THEN** the setup block SHALL include `velocity all create <T> <seed> dist gaussian mom yes rot yes`
  after group/neigh_modify and before integration fixes

### Requirement: Cap force fix line

When `simulation.cap_force` is set to a positive value (kJ/(mol·nm)), the generator
SHALL emit a dedicated cap-force fix after `fix backmap` (and optional `fix pairs`).

#### Scenario: cap_force 50000 kJ/(mol·nm)

- **WHEN** `cap_force: 50000.0` is set
- **THEN** the input script SHALL include `fix cap all backmap/capforce 1195.0300`
  (converted via `units.force()` to kcal/(mol·Å))

#### Scenario: cap_force unset or zero

- **WHEN** `cap_force` is null or ≤ 0
- **THEN** the input script SHALL NOT include `fix cap all backmap/capforce`

### Requirement: Thermostat gamma for Langevin

The generator SHALL convert `simulation.thermostat_gamma` (ps⁻¹) to LAMMPS Langevin
damping (fs) as `units.time(1.0 / thermostat_gamma)`.

#### Scenario: gamma 15 ps⁻¹

- **WHEN** `thermostat_gamma: 15.0` and thermostat is langevin
- **THEN** the langevin fix SHALL use damp ≈ 66.7 fs
