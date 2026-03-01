# LAMMPS Backmapping Package

Time-dependent backmapping (reverse mapping) from coarse-grained to atomistic
resolution, implemented as a LAMMPS user package.

The method ramps a per-atom resolution parameter **lambda** from 0 (pure CG) to
1 (pure AT) uniformly across the simulation box, gradually restoring atomistic
detail.

> Krajniak et al., "Generic Adaptive Resolution Method for Reverse Mapping of
> Polymers from Coarse-Grained to Atomistic Descriptions",
> *J. Chem. Theory Comput.* 2016, 12, 5549–5562.
> [DOI: 10.1021/acs.jctc.6b00595](https://doi.org/10.1021/acs.jctc.6b00595)

**[Full Documentation](https://jkrajniak.github.io/lammps-user-backmapping/)**
— settings reference, tutorials, theory, and LAMMPS component docs.

## Components

| Component | Description |
|-----------|-------------|
| `fix backmap` | Lambda ramp, CG-AT mapping, COM tracking, CG force distribution |
| `pair_style backmap` | Lambda-weighted non-bonded pair forces |
| `bond_style backmap/harmonic` | Lambda-weighted harmonic cross-CG bond forces |
| `bond_style backmap/table` | Lambda-weighted tabulated cross-CG bond forces |
| `angle_style backmap/harmonic` | Lambda-weighted harmonic cross-CG angle forces |

## Repository Layout

```
src/                        C++ LAMMPS styles (fix, pair, bond, angle)
python/
  src/backmap_prep/         Python package source (backmap-prep CLI)
  tests/                    pytest unit tests
  pyproject.toml            Python project metadata & tool config
examples/
  dodecane/                 Dodecane backmapping example
openspec/                   Specifications and change tracking
Makefile                    Top-level convenience targets
```

## Prerequisites

- **LAMMPS** source tree (C++17)
- **CMake** ≥ 3.16
- **Python** ≥ 3.10
- **uv** (Python dependency manager)

## Installation

### C++ Package (LAMMPS Styles)

Use the install script to copy files and register the package (auto-detects
CMake vs traditional make):

```bash
./install.sh /path/to/lammps           # install
./install.sh --uninstall /path/to/lammps   # remove
```

Then rebuild LAMMPS:

```bash
# CMake (recommended)
cd build
cmake -D PKG_BACKMAP=yes /path/to/lammps/cmake
cmake --build .

# Traditional make
cd /path/to/lammps/src
make yes-backmap
make mpi
```

### Python CLI (`backmap-prep`)

```bash
make install-dev     # install with dev dependencies
make install-hooks   # set up pre-commit hooks
```

Or manually:

```bash
cd python && uv sync --extra dev
```

This installs the `backmap-prep` command-line tool.

## Usage

`backmap-prep` reads a YAML settings file that describes the CG-to-AT mapping
and generates the LAMMPS data file, input script, and interaction tables needed
for a backmapping simulation.

```bash
backmap-prep settings.yaml
backmap-prep settings.yaml --output-prefix mysystem
```

See `examples/dodecane/` for a complete working example with a dodecane system
(6 CG beads mapped to 12 united-atom carbons).

### Settings File

The YAML settings file defines:

- **molecules** — CG bead definitions and their constituent AT atoms
- **cg_system** — paths to CG coordinate and topology files (GROMACS format)
- **cross_interactions** — cross-CG bonds, angles, and dihedrals with parameters
- **simulation** — backmapping parameters (alpha, timestep, temperature, cutoffs, …)
- **output** — output prefix and format

## Development

All Python tooling is managed with `uv`. Available Makefile targets:

```bash
make help           # show all targets
make lint           # run ruff linter
make format         # run ruff formatter
make typecheck      # run mypy (strict mode)
make test           # run pytest
make test-cov       # run pytest with coverage
make pre-commit     # run pre-commit hooks on staged files
make pre-commit-all # run pre-commit hooks on all files
make clean          # remove caches and build artifacts
```

## License

GPL-3.0-or-later

## Author

Jakub Krajniak (jkrajniak@gmail.com)
