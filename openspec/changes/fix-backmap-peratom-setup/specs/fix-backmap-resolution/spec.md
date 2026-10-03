## ADDED Requirements

### Requirement: CG force distribution at setup

`fix backmap` SHALL distribute the CG bead forces to the mapped AT atoms after
the force evaluation that opens every run, as it does after every step.

#### Scenario: run 0 at lambda = 0
- **WHEN** two CG beads interact through a CG pair potential and `run 0` is executed
- **THEN** each AT atom SHALL carry m_i / M_bead times its bead's force and the beads SHALL carry none

### Requirement: Per-atom force and COM output

With `peratom full`, `fix backmap` SHALL provide a 7-column per-atom array:
lambda, the atom's bead position, and the CG force share the atom received
(for a bead: its CG force before redistribution).
