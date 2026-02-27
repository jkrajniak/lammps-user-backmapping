# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Conventional Commits](https://www.conventionalcommits.org/).

## [Unreleased]

### Added

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
