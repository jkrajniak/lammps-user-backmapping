## ADDED Requirements

### Requirement: Documentation site builds and deploys via GitHub Pages

The project SHALL provide an MkDocs-based documentation site that builds from Markdown sources in the `docs/` directory and deploys to GitHub Pages automatically on push to `main`.

#### Scenario: Local documentation build

- **WHEN** a developer runs `make docs` from the repository root
- **THEN** MkDocs SHALL build the site locally and serve it for preview at `http://localhost:8000`

#### Scenario: CI deployment on push to main

- **WHEN** a commit is pushed to the `main` branch
- **THEN** the GitHub Actions workflow SHALL build the documentation and deploy it to the `gh-pages` branch

#### Scenario: Build failure on broken links

- **WHEN** the documentation contains broken internal links
- **THEN** the MkDocs build SHALL fail with an error identifying the broken link

### Requirement: Settings reference documents all YAML fields

The documentation site SHALL include a settings reference page (`settings-reference.md`) that documents every field in the YAML settings file with its type, default value, accepted values, and description.

#### Scenario: Top-level sections documented

- **WHEN** a user opens the settings reference page
- **THEN** they SHALL see documentation for all top-level sections: `molecules`, `cg_system`, `cross_interactions`, `simulation`, and `output`

#### Scenario: Each field has type and default

- **WHEN** a user looks up any settings field (e.g., `simulation.alpha`)
- **THEN** the documentation SHALL show the field's type (e.g., `float`), default value (e.g., `0.001`), and a description of its purpose

#### Scenario: Nested structures documented

- **WHEN** a user looks up a nested structure (e.g., `molecules[].beads[].atoms`)
- **THEN** the documentation SHALL show the full path, type, and description with an example

### Requirement: Tutorial for setting up a new system

The documentation site SHALL include a tutorial page (`tutorial-new-system.md`) that walks users through setting up backmapping for a new molecular system from scratch.

#### Scenario: Tutorial covers end-to-end workflow

- **WHEN** a user follows the tutorial
- **THEN** the tutorial SHALL cover: preparing CG and AT input files, writing the settings YAML, running `backmap-prep`, and executing the LAMMPS simulation

#### Scenario: Tutorial includes example snippets

- **WHEN** a user reads a tutorial section
- **THEN** each step SHALL include code/config snippets that the user can adapt for their own system

### Requirement: Theoretical background page

The documentation site SHALL include a theory page (`theory.md`) that explains the backmapping method.

#### Scenario: Lambda ramp explained

- **WHEN** a user reads the theory page
- **THEN** they SHALL find an explanation of the lambda parameter, how it ramps from 0 (CG) to 1 (AT), and the alpha parameter that controls the ramp rate

#### Scenario: Force weighting explained

- **WHEN** a user reads the theory page
- **THEN** they SHALL find an explanation of how bonded and non-bonded forces are weighted by lambda during the transition

### Requirement: Getting started quickstart guide

The documentation site SHALL include a getting-started page (`getting-started.md`) that provides minimal steps to install the package and run the example.

#### Scenario: Quickstart installation

- **WHEN** a new user reads the getting-started page
- **THEN** they SHALL find installation instructions for both the C++ LAMMPS package and the Python CLI tool

#### Scenario: Quickstart first run

- **WHEN** a new user follows the quickstart
- **THEN** they SHALL be able to run the dodecane example and see a successful backmapping simulation complete

### Requirement: LAMMPS component documentation

The documentation site SHALL include pages documenting each LAMMPS style provided by the package: `fix backmap`, `pair_style backmap`, `bond_style backmap/harmonic`, `bond_style backmap/table`, and `angle_style backmap/harmonic`.

#### Scenario: Fix backmap syntax documented

- **WHEN** a user opens the fix-backmap documentation page
- **THEN** they SHALL see the LAMMPS command syntax, all arguments with descriptions, and usage examples

#### Scenario: Pair and bonded style syntax documented

- **WHEN** a user opens a bonded or pair style page
- **THEN** they SHALL see the LAMMPS `pair_coeff` / `bond_coeff` / `angle_coeff` syntax, coefficient descriptions, and how lambda weighting applies

### Requirement: CLI reference documentation

The documentation site SHALL include a page documenting the `backmap-prep` CLI tool with all options and flags.

#### Scenario: CLI help documented

- **WHEN** a user opens the CLI reference page
- **THEN** they SHALL see the command syntax, all available flags, and example invocations
