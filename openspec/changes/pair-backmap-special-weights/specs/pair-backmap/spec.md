## ADDED Requirements

### Requirement: Special-bond factors

`pair_style backmap` SHALL evaluate every listed pair with its `special_bonds`
LJ and Coulomb factors, treating factors below 1e-30 as an exclusion.

#### Scenario: Scaled 1-4 pair
- **WHEN** `special_bonds lj 0 0 0.5 coul 0 0 0.5` and a 1-4 AT pair across beads at lambda = 1
- **THEN** its energy SHALL be half the sub-style energy

### Requirement: CG-specific special weights

With `cg_special w12 w13 w14`, CG-CG pairs SHALL use these weights instead of
the `special_bonds` factors.

#### Scenario: MARTINI exclusions with an AT force field at nrexcl 3
- **WHEN** `special_bonds lj 0 1e-100 1e-100` and `cg_special 0 1 1`
- **THEN** CG 1-3 and 1-4 pairs SHALL be evaluated at full weight and AT 1-3 and 1-4 pairs SHALL be excluded
