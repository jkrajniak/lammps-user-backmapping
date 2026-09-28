## MODIFIED Requirements

### Requirement: Native LAMMPS AT fragment dihedral styles

The native AT-fragment script parser SHALL accept `dihedral_style ryckaert`
(coefficients used as given) and `dihedral_style opls` (K1..K4 converted to the
six Ryckaert coefficients C0 = K1/2 + K2 + K3/2, C1 = K1/2 - 3K3/2,
C2 = -K2 + 4K4, C3 = 2K3, C4 = -4K4, C5 = 0), and SHALL reject any other
dihedral style with a named error.

#### Scenario: OPLS fragment
- **WHEN** a fragment script declares `dihedral_style opls` and `dihedral_coeff 1 K1 K2 K3 K4`
- **THEN** the parsed dihedral coefficients SHALL be the converted Ryckaert six,
  and their energy SHALL equal the OPLS energy at every angle
