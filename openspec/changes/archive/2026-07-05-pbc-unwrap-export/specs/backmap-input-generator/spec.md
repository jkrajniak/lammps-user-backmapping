## ADDED Requirements

### Requirement: Network hybrid molecule-aware unwrap

For network hybrid systems with cross bonded interactions, the generator SHALL
prepare coordinates by unwrapping each `mol_id` fragment independently, then
translating whole molecules to shorten inter-molecular bonds, before writing
LAMMPS data with `ix=iy=iz=0`.

#### Scenario: Cross-mol bond Cartesian span

- **WHEN** a cross bond connects atoms in different `mol_id` groups
- **THEN** export preparation SHALL translate entire molecule fragments and assign
  shared image flags per `mol_id` so min-image bonded distance stays below the
  communication cutoff budget

#### Scenario: Min-image bond validation

- **WHEN** `validate_bond_geometry` runs after coordinate preparation
- **THEN** it SHALL fail if any bonded pair exceeds `max(20 Å, 2 × lj/cg cutoff)`
  in minimum-image distance (using unwrapped coords and image flags when present)
