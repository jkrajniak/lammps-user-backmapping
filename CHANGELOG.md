# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Conventional Commits](https://www.conventionalcommits.org/).

## [Unreleased]

### Added

- `apb` keyword in `fix backmap` for non-uniform atoms-per-bead mapping
  (`apb T1:N1 T2:N2 ...`). Required for systems where different CG bead
  types contain different numbers of AT atoms (e.g. all-atom PE with
  7-atom end beads and 6-atom interior beads).
- `atom->add_callback(Atom::GROW)` and `atom->add_callback(Atom::RESTART)` in
  `fix backmap` to dynamically resize per-atom arrays when LAMMPS reallocates
  atom storage, preventing heap-use-after-free segfaults.

### Fixed

- All PE example input scripts now use `cg_type 1 2` (both CG bead types)
  instead of `cg_type 1`, which caused incorrect bead-to-atom mapping.

- `fix backmap` now accepts multiple CG atom types via `cg_type T1 T2 ...`
  syntax, enabling correct bead-to-atom mapping in systems with more than one
  CG bead type.
- `fix backmap` COM tracking uses round-based PBC wrapping, correctly handling
  ghost atoms that may be more than one box length from the local CG bead.
- `fix backmap` bead map is rebuilt on every neighbor list rebuild (via
  `pre_force()` callback), preventing stale local indices after LAMMPS atom
  sorting.
- `pair_style backmap` defers AT pair interactions until both atoms reach
  &lambda; > 0.1, preventing LJ singularities from initial inter-molecular
  overlaps.
- Python writer emits `cg_type T1 T2 ...` with all CG type IDs and wraps
  atom coordinates into [0, L) in the LAMMPS data file.

### Changed

- Large-scale example input scripts use a robust multi-phase simulation
  protocol (minimise → `nve/limit` relaxation → `nve/limit` lambda ramp →
  gradual NVT equilibration) instead of aggressive single-phase NVE + Langevin,
  preventing "Bond atoms missing" errors in production-size systems.
- Documentation: README, `docs/index.md`, `docs/theory.md`, and `AGENTS.md` now
  tie the package motivation to **migrating from ESPResSo++** toward **LAMMPS**,
  and cite the **complex polymer / network** reverse-mapping paper with
  **December 2017** online publication date (DOI 10.1002/jcc.25129; print *J.
  Comput. Chem.* 2018).

- **backmap-prep** LAMMPS setup: emit `neigh_modify delay 0 every 1 check yes`
  after group definitions; use the **backmapping timestep** for the first
  dynamics segment (including λ-frozen equilibration when
  `equilibration_steps` > 0) instead of the larger production timestep, so
  hybrid relaxation does not start at an overly aggressive Δt.
- **examples/dodecane/n250** `settings.yaml`: `timestep_backmapping: 0.00025` ps
  and `equilibration_steps: 4000` so the default generated input completes the
  λ ramp without early missing-bond errors on the 250-molecule melt.
- **examples/dodecane/n250** RDF validation: `in.dodecane_n250_at`,
  `in.dodecane_n250_at_ref`, and `compare_rdf_n250.sh` (runs
  `../large/compare_rdf.py` strict + relaxed tolerances for small-N noise).

- **backmap-prep** default `simulation.production_steps` is now `0`: generated
  backmapping inputs end with `write_data …_hybrid.data` after the λ ramp only;
  atomistic production (e.g. RDF) is intended as a separate run on the extracted
  AT system. Set `production_steps` > 0 only to append an optional post-backmap
  segment in the same file.
- `scripts/run-backmap.sh`: supports single-segment restart layouts (only
  `*.phase1`, backmapping-only) in addition to two- and three-segment flows.
- Packaging and docs: Python project metadata (`pyproject.toml`, `uv.lock`) is now at repository root so `uv sync` from root installs `backmap-prep`; documentation and examples now use root-level `uv run backmap-prep ...` commands.
- README: repository layout section now lists all example directories and the `scripts/` validation script.

### Added

- **250-molecule dodecane** example layout: `examples/dodecane/n250/` with
  `build_cg_conf.py` (subset of `large/cg_conf.gro`), matching `topol_*_250.top`
  files, `settings.yaml`, and `prepare_inputs.sh`. Suited to laptops (hybrid
  ~4500 atoms; serial LAMMPS RSS typically low compared to 32 GB RAM).

- Small dodecane example: `examples/dodecane/in.dodecane_at` and
  `extract_at_frame.py` for the hybrid → AT → RDF workflow; README documents
  the `compare_rdf.py` sequence.

- Restart/checkpoint support for preemptible cloud instances: new
  `restart_interval` setting in `simulation` generates `write_restart`
  commands and per-phase input scripts. Includes `scripts/run-backmap.sh`
  entrypoint that detects restart files and resumes from the correct phase.
- Google Cloud Batch example (`examples/cloud-batch/`) with spot VM job
  template and setup instructions.
- "Running on Cloud / HPC" documentation page (`docs/cloud-hpc.md`) covering
  restart configuration, Cloud Batch, and Slurm/Apptainer workflows.
- Dockerfile for building LAMMPS with the backmapping package in a container
  (multi-stage build, configurable LAMMPS version via `LAMMPS_VERSION` build
  arg, includes `backmap-prep` CLI). Supports conversion to
  Singularity/Apptainer for HPC clusters.
- Docker documentation page (`docs/docker.md`) with build, run, MPI, and HPC
  conversion instructions.

- Large-scale example variants in `examples/<name>/large/` for dodecane, pe, pe4, pe_10, pe_aa, and melamine. Inputs are sourced from the [bakery](https://github.com/bakery-cg2at/bakery) project; each `large/` contains a README and `settings.yaml` so `backmap-prep` can generate LAMMPS data and input files. See [Large-scale examples](docs/large-scale-examples.md) and the main README.

- C++ LAMMPS styles: `fix backmap`, `pair_style backmap`,
  `bond_style backmap/harmonic`, `bond_style backmap/table`,
  `angle_style backmap/harmonic` for time-dependent CG-to-AT backmapping.
- CMake build integration (`src/CMakeLists.txt`) and legacy `Install.sh`.
- Python CLI tool `backmap-prep` that generates LAMMPS data files, input
  scripts, and interaction tables from a YAML settings file.
- Pydantic v2 settings schema with validation and deferred-feature guards.
- GROMACS topology and coordinate parsers.
- LAMMPS data and input file writers.
- XVG-to-LAMMPS table converter.
- Dodecane example (`examples/dodecane/`) demonstrating full backmapping
  workflow with tabulated CG interactions.
- Makefile with convenience targets for install, lint, format, typecheck,
  test, and pre-commit.
- Pre-commit configuration with ruff and mypy hooks.
- Project constitution and OpenSpec change tracking (`openspec/`).
- Documentation requirements: `CHANGELOG.md` and `README.md` must be kept
  up to date with every change.
