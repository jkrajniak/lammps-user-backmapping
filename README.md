# LAMMPS Backmapping Package

Time-dependent backmapping (reverse mapping) from coarse-grained to atomistic
resolution, implemented as a LAMMPS user package.

The method ramps a per-atom resolution parameter **lambda** from 0 (pure CG) to
1 (pure AT) uniformly across the simulation box, gradually restoring atomistic
detail.

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
pyproject.toml              Python project metadata & tool config
examples/
  dodecane/, pe/, pe4/,      Backmapping examples (each has large/ for
  pe_10/, pe_aa/, melamine/ production-scale variants)
scripts/
  validate-large-scale-prep.sh   Optional: run backmap-prep on one large example
openspec/                   Specifications and change tracking
Makefile                    Top-level convenience targets
```

## Prerequisites

- **LAMMPS** ≥ `stable_22Jul2025` source tree (C++17) — this package uses the
  `Domain::minimum_image(FLERR, ...)` API introduced in that release
- **CMake** ≥ 3.16
- **Python** ≥ 3.10
- **uv** (Python dependency manager)

## Installation

This repository is now structured as a **Python-first distribution**: the main
installable artifact is `backmap-prep`. The LAMMPS C++ package is an optional
engine extension that you install into a LAMMPS source tree when you want to
run simulations.

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

### Docker (Recommended for HPC)

Build a container with LAMMPS + backmapping package — no manual compilation:

```bash
docker build -t lammps-backmap .
docker run --rm -v "$(pwd)":/work lammps-backmap lmp -in in.backmap
```

Override the LAMMPS version with `--build-arg LAMMPS_VERSION=<tag>`. Convert to
Singularity/Apptainer for HPC clusters with `apptainer build lammps-backmap.sif
docker-daemon://lammps-backmap:latest`. See the
[Docker docs](https://jkrajniak.github.io/lammps-user-backmapping/docker/) for
full details.

For cloud spot instances with automatic restart on preemption, add
`restart_interval: 5000` to your `settings.yaml` and use the included
`run-backmap.sh` entrypoint. See
[Running on Cloud / HPC](https://jkrajniak.github.io/lammps-user-backmapping/cloud-hpc/)
for Cloud Batch and Slurm examples.

### Python CLI (`backmap-prep`)

```bash
make install-dev     # install with dev dependencies
make install-hooks   # set up pre-commit hooks
```

Or manually:

```bash
uv sync --extra dev
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
(6 CG beads mapped to 12 united-atom carbons). Larger-scale variants (e.g. 75-chain
PE, 500-molecule melamine) are in each example’s `large/` subdirectory; see
[Large-scale examples](https://jkrajniak.github.io/lammps-user-backmapping/large-scale-examples/) in the docs.

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

## How to Cite

If you use this package in your research, please cite:

> Krajniak, Pandiyan, Nies, Samaey, "Generic Adaptive Resolution
> Method for Reverse Mapping of Polymers from Coarse-Grained to Atomistic
> Descriptions", *J. Chem. Theory Comput.* 2016.
> [DOI: 10.1021/acs.jctc.6b00595](https://doi.org/10.1021/acs.jctc.6b00595)

If you are reverse mapping complex polymer structures (e.g. polymer networks),
please also cite:

> Krajniak, Zhang, Pandiyan, Nies, Samaey, "Reverse Mapping Method for
> Complex Polymer Systems", *J. Comput. Chem.* 2018.
> [DOI: 10.1002/jcc.25129](https://doi.org/10.1002/jcc.25129)

```bibtex
@article{krajniak2016generic,
  title={Generic Adaptive Resolution Method for Reverse Mapping of Polymers from Coarse-Grained to Atomistic Descriptions},
  author={Krajniak, Jakub and Pandiyan, Sudharsan and Nies, Erik and Samaey, Giovanni},
  journal={Journal of Chemical Theory and Computation},
  year={2016},
  publisher={American Chemical Society}
}

@article{krajniak2018reverse,
  title={Reverse mapping method for complex polymer systems},
  author={Krajniak, Jakub and Zhang, Zidan and Pandiyan, Sudharsan and Nies, Eric and Samaey, Giovanni},
  journal={Journal of Computational Chemistry},
  volume={39},
  number={11},
  pages={648--664},
  year={2018}
}
```

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

GPL-3.0-or-later — see [LICENSE](LICENSE) for the full text.

## Authors

- Jakub Krajniak (jkrajniak@gmail.com)
- Zidan Zhang
