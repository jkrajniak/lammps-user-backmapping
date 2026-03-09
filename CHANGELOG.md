# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Conventional Commits](https://www.conventionalcommits.org/).

## [Unreleased]

### Changed

- Documentation: getting-started now uses correct `backmap-prep` invocation from `python/`; index quick links and large-scale examples page mention the validation script; examples README documents large-scale variants and `scripts/validate-large-scale-prep.sh`.
- README: repository layout section now lists all example directories and the `scripts/` validation script.

### Added

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
