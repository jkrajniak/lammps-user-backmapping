## ADDED Requirements

### Requirement: Single hybrid engine

The generator SHALL build every hybrid system, linear or network, through one
engine. The settings key `prep.engine` SHALL NOT select a code path; if present
it SHALL be ignored with a deprecation warning.

#### Scenario: Linear melt and network use the same path
- **WHEN** `backmap-prep build` runs on a linear melt (e.g. PE) and on a network
  (e.g. rim135)
- **THEN** both SHALL be built by the same hybrid builder and written by the
  same writers

### Requirement: Fragment placement on the CG bead

Each AT fragment SHALL be translated so that its mass-weighted centre of mass
coincides with the position of its CG bead.

#### Scenario: Every bead sits on its fragment COM
- **WHEN** a hybrid data file is generated
- **THEN** for every CG bead, the distance between the bead and the
  mass-weighted COM of its AT atoms (bead map as used by `fix backmap`) SHALL
  be below 0.01 A

### Requirement: Complete CG bonded terms

The data file SHALL contain every CG bond, angle and dihedral defined by the CG
topology and by `cross_interactions` entries with `cg_bonded: true`.

#### Scenario: CG angles present for a linear melt
- **WHEN** the PE example is generated
- **THEN** the data file SHALL contain one CG angle per CG bead triple and one
  CG dihedral per CG bead quadruple of each chain

### Requirement: Force field matches the GROMACS source

At lambda = 1 the per-term energies (bond, angle, dihedral, LJ) of the
generated LAMMPS system SHALL match GROMACS energies of the AT source topology
for the same coordinates within 1e-5 relative.

#### Scenario: GROMACS rerun agrees with LAMMPS run 0
- **WHEN** the reference check runs on a regenerated example
- **THEN** every per-term energy SHALL agree within tolerance

## MODIFIED Requirements

### Requirement: Feature phasing

The Phase 1 and Phase 4 features (linear molecules, LAMMPS-native CG system and
AT fragment inputs) SHALL be provided by the single engine. LAMMPS-format inputs
SHALL be converted to GROMACS-equivalent inputs before the hybrid build, with an
exact unit round trip.
