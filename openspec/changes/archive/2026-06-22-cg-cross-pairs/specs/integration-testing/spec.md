## ADDED Requirements

### Requirement: Rim135 cross-pair integration

#### Scenario: pairs.dat emitted

- **WHEN** `test_rim135_build_v2_lammps_smoke` runs
- **THEN** `pairs.dat` SHALL exist with > 0 pairs

#### Scenario: fix backmap/pairs in input

- **WHEN** `in.rim135` is generated
- **THEN** it SHALL contain `fix backmap/pairs`
