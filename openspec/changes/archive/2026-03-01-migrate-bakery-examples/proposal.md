## Why

The bakery project has a rich library of validated backmapping examples (PE variants, dodecane, melamine, epoxy networks) built for ESPResSo++. Dodecane has already been migrated to the new LAMMPS-based `backmap-prep` workflow. The remaining examples need migration to serve as validation targets, user tutorials, and regression tests for the LAMMPS backmapping package. Migrating them now exercises the tool across diverse chemistries and mapping resolutions, exposing gaps in the YAML schema and `backmap-prep` generator before users encounter them.

## What Changes

- **Add PE example** (`examples/pe/`): Polyethylene with 2 united-atoms per CG bead (75 chains, 50 beads/chain). Tabulated CG bonds, angles, and dihedrals. Exercises the standard linear polymer workflow with multiple CG bead types (A terminal, B backbone).
- **Add PE4 example** (`examples/pe4/`): Polyethylene with 4 united-atoms per CG bead (75 chains, 25 beads/chain). Coarser mapping variant validating that `backmap-prep` handles different mapping ratios.
- **Add PE-10 example** (`examples/pe_10/`): Polyethylene with ~30-36 all-atoms per CG bead (75 chains, 10 beads/chain). Tests very coarse mapping with OPLS/AA force field and many atoms per bead.
- **Add PE-AA example** (`examples/pe_aa/`): Polyethylene all-atom OPLS/AA with explicit hydrogens (6-7 atoms per bead, 50 beads/chain). Tests the full all-atom workflow including hydrogen handling.
- **Add melamine example** (`examples/melamine/`): Melamine-formaldehyde small molecule (3 CG beads per molecule, 500 molecules). Simple network-forming molecule — tests triangular CG topology.
- **Copy source GROMACS files** (GRO, TOP, XVG tables) from bakery into each example directory, adapting naming conventions.
- **Create `settings.yaml`** for each example defining molecule mapping, CG system, cross-interactions, and simulation parameters in the new YAML format.
- **Create LAMMPS input scripts** (`in.<system>`) as reference/template for each system.
- **Add per-example `README.md`** describing the system, source, and how to run it.
- **Add verification scripts** for validating output (molecule integrity, energy convergence).

## Capabilities

### New Capabilities
- `example-pe-systems`: Migration of polyethylene examples (pe, pe4, pe_10, pe_aa) from bakery XML workflow to YAML-driven backmap-prep. Covers linear polymer chains with varying CG-to-AT mapping ratios (2:1, 4:1, 10:1) and force field types (OPLS-UA, OPLS-AA).
- `example-melamine`: Migration of melamine-formaldehyde example from bakery. Covers small network-forming molecules with triangular CG topology and single bead type.

### Modified Capabilities
- `backmap-input-generator`: The settings schema may need minor extensions to handle dihedral interactions and multi-type CG systems encountered in PE and melamine examples. Specifically, the existing `cross_interactions.dihedrals` support and tabulated potential conversion must be validated against real data.
- `integration-testing`: New examples serve as additional integration test targets beyond dodecane.

## Impact

- **`examples/` directory**: Five new example subdirectories with settings, input files, and documentation.
- **`python/` (backmap-prep)**: May require bug fixes or minor schema extensions if edge cases are found during migration. No API changes expected.
- **Documentation**: New examples should be referenced in `docs/getting-started.md` and tutorial pages.
- **CI**: New examples can be added to integration test matrix.
- **Dependencies**: No new dependencies. All examples use existing GROMACS-format inputs and the `backmap-prep` tool.
