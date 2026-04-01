# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Conventional Commits](https://www.conventionalcommits.org/).

## [Unreleased]

### Changed

- Packaging and docs: Python project metadata (`pyproject.toml`, `uv.lock`) is now at repository root so `uv sync` from root installs `backmap-prep`; documentation and examples now use root-level `uv run backmap-prep ...` commands.
- README: repository layout section now lists all example directories and the `scripts/` validation script.

### Added

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
