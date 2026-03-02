## MODIFIED Requirements

### Requirement: Feature phasing

The generator SHALL be implemented in phases. The YAML schema supports all features from Phase 1 through Phase 4, but deferred features SHALL produce a clear "not yet implemented" error if encountered in the settings file.

**Phase 1 (MVP)** — Linear molecules:
- YAML parsing with all five sections
- Single molecule type, simple `atoms` list per bead
- Cross bonds and angles from `cross_interactions`
- GROMACS source files (`.gro`/`.top`/`.itp`)
- Tabulated potential conversion (`.xvg` to `.table`)
- Unit conversion (GROMACS to LAMMPS real)
- Full LAMMPS `.data` and `.in` generation with `backmap/*` styles
- Validation: dodecane, PE (polyethylene), PE4, PE-10, PE-AA, melamine

**Phase 2** — Extended bonded interactions:
- Cross dihedrals and 1-4 pairs
- Tabulated bond/dihedral potentials
- Two-phase backmapping (`two_phase: true`)
- Energy minimization
- Force capping

**Phase 3** — Reactive networks:
- Degree-dependent bead definitions (`atoms_by_degree`)
- Active sites and bond formation logic
- Charge management (equilibration, transfers, charge_map, type_map)
- Atom removal on bond formation (`remove`)
- Predefined active sites file
- Restricted cross-bond patterns
- Multiple molecule types in one system
- Validation: epoxy networks (RIM135)

**Phase 4** — Format flexibility:
- Non-GROMACS source formats (PDB, LAMMPS data, XYZ)
- Non-GROMACS CG input formats
- Native LAMMPS `.table` format as input (skip conversion)

#### Scenario: Deferred feature used
- **WHEN** the YAML contains `atoms_by_degree` (Phase 3 feature) in Phase 1
- **THEN** the generator SHALL abort with: "Feature 'atoms_by_degree' is not yet implemented (planned for Phase 3)"

#### Scenario: Phase 2 dihedral feature used
- **WHEN** the YAML contains `cross_interactions.dihedrals` entries
- **THEN** the generator SHALL process them and produce dihedral sections in the output files

#### Scenario: Phase 3 feature used before implementation
- **WHEN** the YAML contains `atoms_by_degree` (Phase 3 feature) before Phase 3 is implemented
- **THEN** the generator SHALL abort with: "Feature 'atoms_by_degree' is not yet implemented (planned for Phase 3)"

#### Scenario: Phase 4 format before implementation
- **WHEN** `cg_system.format: pdb` is specified before Phase 4 is implemented
- **THEN** the generator SHALL abort with: "Format 'pdb' is not yet supported (planned for Phase 4)"

#### Scenario: PE validation in Phase 1
- **WHEN** `backmap-prep` processes the PE example `settings.yaml` with 50 beads/chain, 2 atoms/bead, tabulated CG bonds and angles
- **THEN** the output SHALL be a valid LAMMPS data file and input script that runs without errors

#### Scenario: Melamine validation in Phase 1
- **WHEN** `backmap-prep` processes the melamine example `settings.yaml` with 3 beads/molecule in a triangular topology
- **THEN** the output SHALL be a valid LAMMPS data file with correct triangular bonding and input script
