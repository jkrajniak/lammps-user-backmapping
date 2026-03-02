## MODIFIED Requirements

### Requirement: Integration test coverage for migrated examples

The integration test suite SHALL validate the migrated bakery examples in addition to the existing dodecane test. Each migrated example SHALL have a test that:

1. Runs `backmap-prep settings.yaml` to generate LAMMPS data and input files
2. Validates the generated data file structure (atom counts, bond counts, molecule IDs)
3. Optionally runs a short LAMMPS simulation to verify the input script is functional

#### Scenario: PE example generates valid output
- **WHEN** the CI pipeline runs `backmap-prep` on `examples/pe/settings.yaml`
- **THEN** the generated `pe.data` SHALL have the correct number of atoms, bonds, angles, and molecule IDs for the PE system

#### Scenario: PE4 example generates valid output
- **WHEN** the CI pipeline runs `backmap-prep` on `examples/pe4/settings.yaml`
- **THEN** the generated `pe4.data` SHALL have the correct number of atoms per molecule (25 CG + 100 AT = 125) and correct bond topology

#### Scenario: PE-10 example generates valid output
- **WHEN** the CI pipeline runs `backmap-prep` on `examples/pe_10/settings.yaml`
- **THEN** the generated `pe_10.data` SHALL have the correct atom counts reflecting the high AT-to-CG ratio

#### Scenario: PE-AA example generates valid output
- **WHEN** the CI pipeline runs `backmap-prep` on `examples/pe_aa/settings.yaml`
- **THEN** the generated `pe_aa.data` SHALL include hydrogen atoms with correct charges and OPLS/AA types

#### Scenario: Melamine example generates valid output
- **WHEN** the CI pipeline runs `backmap-prep` on `examples/melamine/settings.yaml`
- **THEN** the generated `melamine.data` SHALL have triangular CG bonding (3 CG bonds per molecule) and correct atom counts (30 atoms per molecule)

#### Scenario: All examples produce consistent output
- **WHEN** `backmap-prep` is run on any migrated example twice with the same settings
- **THEN** the output files SHALL be byte-identical (deterministic generation)
