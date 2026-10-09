## ADDED Requirements

### Requirement: OPLS dihedral backmap style

The package SHALL provide `dihedral_style backmap/opls` with coefficients
`at|cg K1 K2 K3 K4` and energy
`E = w [K1/2 (1 + cos phi) + K2/2 (1 - cos 2phi) + K3/2 (1 + cos 3phi) + K4/2 (1 - cos 4phi)]`,
where w is the `compute_weight3` weight of the dihedral.

#### Scenario: Equals the stock style at full weight
- **WHEN** an intra-bead `at` dihedral is evaluated at any lambda
- **THEN** its energy and forces SHALL equal those of LAMMPS `dihedral_style opls`
  with the same K1..K4

#### Scenario: Weighted across beads
- **WHEN** lambda_global = 0.5 and the dihedral is `at` across beads
- **THEN** energy and forces SHALL be half the stock values
