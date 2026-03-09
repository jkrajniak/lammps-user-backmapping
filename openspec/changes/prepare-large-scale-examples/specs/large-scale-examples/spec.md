## ADDED Requirements

### Requirement: Large-scale layout per example

Each existing example directory (dodecane, pe, pe4, pe_10, pe_aa, melamine) SHALL support a large-scale variant under a `large/` subdirectory. The `large/` subdirectory SHALL contain or reference inputs and generated outputs for production-scale system sizes (e.g. 75 chains for PE systems, 500 molecules for melamine), aligned with bakery’s reference sizes where applicable.

#### Scenario: Large-scale subdirectory exists for PE
- **WHEN** the user inspects `examples/pe/`
- **THEN** a `large/` subdirectory SHALL exist and SHALL contain (or document how to obtain) large-scale CG coordinates, topology, and LAMMPS input/data files

#### Scenario: Large-scale subdirectory exists for melamine
- **WHEN** the user inspects `examples/melamine/`
- **THEN** a `large/` subdirectory SHALL exist and SHALL contain or document inputs for ~500 molecules

### Requirement: Source of large-scale inputs from bakery

Large-scale inputs (GRO, TOP, XVG tables, single-molecule references) SHALL be sourced from the bakery project. Each example’s README or a central “Large-scale examples” doc SHALL state that large-scale inputs come from bakery and SHALL document the path or procedure to copy or refresh them.

#### Scenario: Documentation references bakery
- **WHEN** a user reads the large-scale example documentation
- **THEN** it SHALL indicate that inputs are obtained from bakery and SHALL describe how to get or update them

### Requirement: backmap-prep and LAMMPS workflow for large-scale

Large-scale variants SHALL use the same `backmap-prep` settings schema and LAMMPS workflow as the small-scale examples. Where paths differ (e.g. `cg_conf.gro` in `large/`), a `settings.yaml` in `large/` or equivalent configuration SHALL be provided so that `backmap-prep` generates LAMMPS input and data files for the large system.

#### Scenario: backmap-prep runs on large-scale settings
- **WHEN** the user runs `backmap-prep` with the large-scale settings for an example (e.g. `examples/pe/large/settings.yaml`)
- **THEN** it SHALL produce the corresponding `.data` and `in.<name>` files without errors, subject to valid bakery-sourced inputs being present

#### Scenario: LAMMPS runs with large-scale data
- **WHEN** the user runs LAMMPS with the generated large-scale input script and data file
- **THEN** the run SHALL complete without fatal errors for at least one documented large-scale example (e.g. dodecane or pe)

### Requirement: Documentation for large-scale examples

The repository SHALL document large-scale examples so users can run them. Documentation SHALL include: (1) where large-scale inputs come from (bakery), (2) layout of `large/` (or equivalent), (3) how to run `backmap-prep` and LAMMPS for the large variant. The main README and/or docs SHALL be updated to mention large-scale examples and link to the detailed instructions.

#### Scenario: README or docs mention large-scale
- **WHEN** a user reads the top-level README or the “Large-scale examples” doc
- **THEN** they SHALL find a reference to large-scale variants and how to obtain and run them

### Requirement: Optional CI or script for one large-scale run

An optional integration test, CI job, or script SHALL exist that runs `backmap-prep` and LAMMPS on at least one large-scale example (e.g. dodecane or pe). The test MAY be time-limited or run only on demand/schedule so that normal CI remains fast.

#### Scenario: Optional large-scale validation exists
- **WHEN** the optional large-scale validation is executed (manually or by CI)
- **THEN** it SHALL run `backmap-prep` and LAMMPS for at least one large-scale example and SHALL report success or failure
