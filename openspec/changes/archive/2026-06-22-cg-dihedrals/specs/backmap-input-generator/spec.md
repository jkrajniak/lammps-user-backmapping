## ADDED Requirements

### Requirement: Hybrid TOP dihedral export

The network builder SHALL parse `[ dihedrals ]` and `[ cross_dihedrals ]`, resolve func-3
RB coefficients from dihedraltypes, and emit LAMMPS dihedral sections.

#### Scenario: Rim135 dihedral count

- **WHEN** `backmap-prep build examples/epoxy/settings.v2.yaml` runs
- **THEN** `rim135.data` SHALL contain approximately 33,421 dihedrals

#### Scenario: Missing dihedraltypes

- **WHEN** a func-3 dihedral has no inline params and no dihedraltypes match
- **THEN** the generator SHALL fail with a clear error naming the atom types
