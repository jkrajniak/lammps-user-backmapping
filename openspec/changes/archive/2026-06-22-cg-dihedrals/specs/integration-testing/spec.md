## ADDED Requirements

### Requirement: Rim135 dihedral integration

#### Scenario: Non-zero dihedral count

- **WHEN** `test_rim135_build_v2_lammps_smoke` runs
- **THEN** generated data SHALL NOT contain `0 dihedrals`

#### Scenario: Hybrid dihedral styles in input

- **WHEN** `in.rim135` is generated
- **THEN** it SHALL contain `dihedral_style hybrid` with `ryckaert` and/or `backmap/ryckaert`
