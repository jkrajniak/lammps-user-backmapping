## ADDED Requirements

### Requirement: Rim135 Tier B bakery protocol input gates (PR4)

The rim135 integration test SHALL verify that settings v2 generates a bakery-aligned
LAMMPS input suitable for Tier B dynamics tuning.

#### Scenario: Velocity initialization

- **WHEN** `test_rim135_build_v2_lammps_smoke` runs
- **THEN** `in.rim135` SHALL contain `velocity all create`

#### Scenario: Cap force fix

- **WHEN** `test_rim135_build_v2_lammps_smoke` runs
- **THEN** `in.rim135` SHALL contain `fix cap all backmap/capforce`

#### Scenario: Langevin at 298 K with bakery gamma

- **WHEN** `test_rim135_build_v2_lammps_smoke` runs
- **THEN** `in.rim135` SHALL contain `langevin 298.0` with damp between 66.0 and 67.0 fs

#### Scenario: Bakery dt and ramp length

- **WHEN** `test_rim135_build_v2_lammps_smoke` runs
- **THEN** `in.rim135` SHALL use `timestep 1.00` during backmapping and `run 10000` for the λ ramp

#### Scenario: No in-hybrid equilibration for pre-equilibrated CG

- **WHEN** `test_rim135_build_v2_lammps_smoke` runs
- **THEN** `in.rim135` SHALL NOT contain `fix_modify bm active no`
