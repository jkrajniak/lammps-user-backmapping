## ADDED Requirements

### Requirement: Rim135 Tier B soft-start input gates

The rim135 integration test SHALL verify that settings v2 generates a soft-start
LAMMPS input suitable for Tier B dynamics tuning.

#### Scenario: Gentler ramp parameters

- **WHEN** `test_rim135_build_v2_lammps_smoke` runs
- **THEN** `in.rim135` SHALL use `timestep 0.50` during backmapping and `run 20000` for the λ ramp

#### Scenario: No in-hybrid equilibration for pre-equilibrated CG

- **WHEN** `test_rim135_build_v2_lammps_smoke` runs
- **THEN** `in.rim135` SHALL NOT contain `fix_modify bm active no`
