## MODIFIED Requirements

### Requirement: Simulation section — dynamics parameters

The generator SHALL honor `simulation.equilibration_steps`, `simulation.timestep_backmapping`,
and `simulation.alpha` when writing LAMMPS input scripts for network hybrids.

For rim135 settings v2, the recommended Tier B soft-start protocol SHALL use:

- `equilibration_steps: 0` (CG pre-equilibrated)
- `timestep_backmapping: 0.0005` (ps)
- `alpha: 0.00005` (20k ramp steps)

#### Scenario: Soft-start rim135 input

- **WHEN** `settings.v2.yaml` is built with PR3 dynamics knobs
- **THEN** `in.rim135` SHALL use `timestep 0.50` and `run 20000` for the λ ramp
