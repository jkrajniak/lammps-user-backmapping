## ADDED Requirements

### Requirement: Crosslinked melamine example imports bakery's exact network

The `examples/melamine_network/` directory SHALL contain a melamine-formaldehyde
backmapping example using the same 500-molecule MF system as
`examples/melamine/`, but crosslinked to match bakery's reference network
exactly: 675 inter-molecular covalent bonds across 1500 CG beads (the same
pattern found in `network_backmapping/mf/backmapping/cg_topol.top`), imported
as a fixed, static topology -- not regenerated via live reaction dynamics.

#### Scenario: Crosslink count matches bakery's reference
- **WHEN** the generated hybrid topology for `examples/melamine_network/` is
  inspected
- **THEN** it SHALL contain exactly 675 inter-molecular AT bonds (plus their
  associated crosslink-site angles/dihedrals), matching bakery's `cg_topol.top`
  `; chem` bond count

#### Scenario: Unreacted arms match the uncrosslinked example
- **WHEN** a CG bead in `examples/melamine_network/` is not part of a crosslink
- **THEN** its AT fragment SHALL be identical (atom count, charges, bonded
  terms) to the corresponding bead in `examples/melamine/`

### Requirement: Crosslink sites have real force-field parameters

The generator SHALL NOT leave the C-O-C angle or C-C-O-C dihedral terms at
crosslink sites unparameterized. Real OPLS-AA parameters SHALL be sourced and
applied, resolving the gap present in bakery's own reference
(`missing_definitions.txt`, `; chem MISSING params type: A-A`).

#### Scenario: No missing definitions for generated crosslink terms
- **WHEN** `examples/melamine_network/`'s hybrid topology is generated
- **THEN** every angle/dihedral type actually present in the generated system
  SHALL have a defined force constant (no unresolved "MISSING params" entries
  for types used in this system)

### Requirement: Crosslink bonds are always-on, not lambda-weighted

Crosslink AT bonds/angles/dihedrals SHALL be written using plain always-on
bond/angle/dihedral styles (e.g. `harmonic`), not the lambda-weighted
`backmap/harmonic` inter-bead styles used for the CG/AT resolution ramp --
they represent real, permanent covalent chemistry, independent of the
current resolution state.

#### Scenario: Crosslink bond present at lambda=0
- **WHEN** the generated `in.*` LAMMPS input is inspected at any point in the
  eq/ramp/production sequence, including lambda=0
- **THEN** crosslink bond/angle/dihedral terms SHALL contribute their full,
  unweighted force -- independent of the current lambda value

### Requirement: Existing uncrosslinked melamine example is unaffected

Adding `examples/melamine_network/` SHALL NOT modify `examples/melamine/`,
its `settings.yaml`, its generated files, or
`openspec/specs/example-melamine/spec.md`.

#### Scenario: Uncrosslinked example still builds unchanged
- **WHEN** `backmap-prep build examples/melamine/large/settings.yaml` is run
  after this change is merged
- **THEN** it SHALL produce byte-identical output to before this change
