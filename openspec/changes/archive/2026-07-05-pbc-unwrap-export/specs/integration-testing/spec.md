## ADDED Requirements

### Requirement: Rim135 min-image bond gate

#### Scenario: LAMMPS data bond spans after v2 build

- **WHEN** `test_rim135_build_v2_lammps_smoke` runs
- **THEN** max min-image bonded distance in the generated system SHALL be < 20 Å
