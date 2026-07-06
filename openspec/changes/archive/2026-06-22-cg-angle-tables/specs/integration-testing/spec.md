## ADDED Requirements

### Requirement: Rim135 CG angle table integration

The rim135 integration test SHALL verify that func-8 CG cross-angles are wired
to tabulated potentials.

#### Scenario: Generated angle table files

- **WHEN** `test_rim135_build_v2_lammps_smoke` runs successfully
- **THEN** `table_a1.table` and `table_a2.table` SHALL exist in the output directory

#### Scenario: No placeholder CG angle coefficients

- **WHEN** `in.rim135` is generated from settings v2
- **THEN** it SHALL contain `angle_coeff` lines with `backmap/table cg`
- **AND** it SHALL NOT contain `angle_coeff … cg 0.000000 0.0000`

#### Scenario: Distinct tablenr types

- **WHEN** rim135 hybrid TOP has func-8 angles with tablenr 1 and 2
- **THEN** the generated input SHALL reference both `table_a1.table` and
  `table_a2.table` in angle_coeff lines
