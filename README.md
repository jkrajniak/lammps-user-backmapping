# LAMMPS Backmapping Package

Adaptive-resolution backmapping (reverse mapping) from coarse-grained to
atomistic resolution, implemented as a LAMMPS user package.

The method ramps a per-atom resolution parameter **lambda** from 0 (pure CG) to
1 (pure AT) uniformly across the simulation box, gradually restoring atomistic
detail.

**Goal:** bring the methodology that was originally implemented with **custom
bonded extensions in ESPResSo++** into **standard LAMMPS** plus a small
**`backmap-prep`** toolchain, so the same physics can be run with familiar
inputs, restarts, and HPC/cloud workflows. The **linear-melt** examples in this
repository are the natural first validation target; **reaction-formed
topologies** (e.g. **epoxy and melamine networks**, hyperbranched polymers) were
demonstrated at the fine-grained level in the **complex-system reverse mapping**
publication cited below (first published online **December 2017**).

**[Full Documentation](https://jkrajniak.github.io/lammps-user-backmapping/)**
— settings reference, tutorials, theory, and LAMMPS component docs.

## Components


| Component                            | Description                                                     |
| ------------------------------------ | --------------------------------------------------------------- |
| `fix backmap`                        | Lambda ramp, CG-AT mapping, COM tracking, CG force distribution |
| `fix backmap/pairs`                  | Lambda-weighted explicit 1-4 pairs                              |
| `fix backmap/capforce`               | Per-atom force cap                                              |
| `fix backmap/stereo`                 | Stereocenter restraints, correct across MPI ranks               |
| `fix backmap/ringguard`              | Repulsive wall at ring centroids against ring spearing          |
| `pair_style backmap`                 | Lambda-weighted non-bonded pair forces (AT and CG sub-styles)   |
| `bond_style backmap/harmonic`        | Lambda-weighted harmonic bonds                                  |
| `bond_style backmap/gromos`          | Lambda-weighted GROMOS quartic bonds                            |
| `bond_style backmap/table`           | Lambda-weighted tabulated bonds                                 |
| `angle_style backmap/harmonic`       | Lambda-weighted harmonic angles                                 |
| `angle_style backmap/charmm`         | Lambda-weighted harmonic + Urey-Bradley angles                  |
| `angle_style backmap/cosine/squared` | Lambda-weighted GROMOS cosine-harmonic angles                   |
| `angle_style backmap/table`          | Lambda-weighted tabulated angles                                |
| `dihedral_style backmap/ryckaert`    | Lambda-weighted Ryckaert-Bellemans dihedrals                    |
| `dihedral_style backmap/harmonic`    | Lambda-weighted harmonic dihedrals                              |
| `dihedral_style backmap/fourier`     | Lambda-weighted multi-term periodic dihedrals                   |
| `dihedral_style backmap/opls`        | Lambda-weighted OPLS dihedrals                                  |
| `dihedral_style backmap/table`       | Lambda-weighted tabulated dihedrals                             |
| `improper_style backmap/harmonic`    | Lambda-weighted harmonic impropers                              |


## Repository Layout

```
src/                        C++ LAMMPS styles (fix, pair, bond, angle, dihedral, improper)
python/
  src/backmap_prep/         Python package source (backmap-prep CLI)
  tests/                    pytest unit tests
pyproject.toml              Python project metadata & tool config
examples/
  dodecane/, pe/, pe4/,      Backmapping examples (each has large/ for
  pe_10/, pe_aa/, melamine/ production-scale variants; dodecane also has
  n250/ for a 250-molecule laptop-friendly melt)
scripts/
  validate-large-scale-prep.sh   Optional: run backmap-prep on one large example
openspec/                   Specifications and change tracking
Makefile                    Top-level convenience targets
```

## Prerequisites

- **LAMMPS** ≥ `stable_22Jul2025` source tree (C++17), which introduced the
`Domain::minimum_image(FLERR, ...)` API this package uses. CI builds against
`stable_22Jul2025_update2` and `stable_30Sep2026_update1`.
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
Singularity/Apptainer for HPC clusters with `apptainer build lammps-backmap.sif docker-daemon://lammps-backmap:latest`. See the
[Docker docs](https://jkrajniak.github.io/lammps-user-backmapping/docker/) for
full details.

For cloud spot instances with automatic restart on preemption, add
`restart_interval: 5000` (with `protocol: standard`) to your `settings.yaml` and use the included
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

`build` also writes `<prefix>.ff.lmp` and `<prefix>.backmap.lmp` (include
them in your own protocols instead of restating coefficients) and
`<prefix>.at.ff.lmp` for AT-only runs; `backmap-prep at-system --from
<frame>` / `--reference N` writes the matching AT-only data files.

Linear melts and cross-linked networks go through the same hybrid builder:
each AT fragment is placed with its mass-weighted centre of mass on its CG
bead, and outputs are written next to the settings file.

See `examples/dodecane/` for a complete working example with a dodecane system
(6 CG beads mapped to 12 united-atom carbons). Larger-scale variants (e.g. 75-chain
PE, 500-molecule melamine) are in each example’s `large/` subdirectory; see
[Large-scale examples](https://jkrajniak.github.io/lammps-user-backmapping/large-scale-examples/) in the docs.

Large-scale systems use a robust multi-phase protocol (energy minimisation,
`nve/limit` relaxation, `nve/limit` lambda ramp, gradual NVT equilibration)
to avoid instabilities from initial AT overlaps. See
[Theory — Simulation Protocol](https://jkrajniak.github.io/lammps-user-backmapping/theory/#simulation-protocol)
for details.

### Settings File

The YAML settings file defines:

- **prep** — input directories (`data_dir`, `tables_dir`, `forcefield_dir`); only read
- **molecules** — CG bead definitions and their constituent AT atoms (`name` must equal the CG residue name); each molecule's AT reference is GROMACS `.gro`/`.top` (default) or a native LAMMPS `data` file + input script
- **cg_system** — paths to CG coordinate and topology files (GROMACS `.gro`/`.top`, or a native LAMMPS `data` file)
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

LAMMPS regression tests (`python/tests/test_lammps_*.py`) run only when
`BACKMAP_LMP` points to a LAMMPS binary built with the package:
`BACKMAP_LMP=/path/to/lmp uv run pytest python/tests/test_lammps_energy.py python/tests/test_lammps_write_data.py`.

## How to Cite

If you use this package in your research, please cite:

> Krajniak, Pandiyan, Nies, Samaey, "Generic Adaptive Resolution
> Method for Reverse Mapping of Polymers from Coarse-Grained to Atomistic
> Descriptions", *J. Chem. Theory Comput.* 2016.
> [DOI: 10.1021/acs.jctc.6b00595](https://doi.org/10.1021/acs.jctc.6b00595)

If you are reverse mapping **complex polymer structures** (connectivity created
or evolved at the CG level, e.g. **epoxy networks**, **trimethylol melamine
networks**, hyperbranched polymers, or **PET**-style step-growth products),
please also cite:

> Krajniak, Zhang, Pandiyan, Nies, Samaey, "Reverse Mapping Method for
> Complex Polymer Systems", *J. Comput. Chem.* **39**, 648--664 (**2018**);
> first published online **6 December 2017**.
> [DOI: 10.1002/jcc.25129](https://doi.org/10.1002/jcc.25129)

```bibtex
@article{krajniak2016generic,
  title={Generic Adaptive Resolution Method for Reverse Mapping of Polymers from Coarse-Grained to Atomistic Descriptions},
  author={Krajniak, Jakub and Pandiyan, Sudharsan and Nies, Erik and Samaey, Giovanni},
  journal={Journal of Chemical Theory and Computation},
  volume={12},
  number={11},
  pages={5549--5562},
  year={2016},
  doi={10.1021/acs.jctc.6b00595},
  publisher={American Chemical Society}
}

@article{krajniak2018reverse,
  title={Reverse mapping method for complex polymer systems},
  author={Krajniak, Jakub and Zhang, Zidan and Pandiyan, Sudharsan and Nies, Eric and Samaey, Giovanni},
  journal={Journal of Computational Chemistry},
  volume={39},
  number={11},
  pages={648--664},
  year={2018},
  doi={10.1002/jcc.25129},
  note={First published online 6 December 2017}
}
```

### Applications (optional)

The method has been applied in the studies below. If your system is close to
one of them, you may also cite the relevant paper.

If you are backmapping **charged polymers** (e.g. **polymeric ionic liquids**
of different chain lengths), you may also cite:

> Zhang, Wheatle, Krajniak, Keith, Ganesan, "Ion Mobilities, Transference
> Numbers, and Inverse Haven Ratios of Polymeric Ionic Liquids",
> *ACS Macro Lett.* **9**, 84--89 (**2020**).
> [DOI: 10.1021/acsmacrolett.9b00908](https://doi.org/10.1021/acsmacrolett.9b00908)

If you are backmapping **microphase-separated block copolymers** (morphology
equilibrated at the CG level, atomistic detail added by backmapping), you may
also cite:

> Zhang, Krajniak, Keith, Ganesan, "Mechanisms of Ion Transport in Block
> Copolymeric Polymerized Ionic Liquids", *ACS Macro Lett.* **8**, 1096--1101
> (**2019**).
> [DOI: 10.1021/acsmacrolett.9b00478](https://doi.org/10.1021/acsmacrolett.9b00478)

> Zhang, Krajniak, Ganesan, "A Multiscale Simulation Study of Influence of
> Morphology on Ion Transport in Block Copolymeric Ionic Liquids",
> *Macromolecules* **54**, 4997--5010 (**2021**).
> [DOI: 10.1021/acs.macromol.1c00025](https://doi.org/10.1021/acs.macromol.1c00025)

> Zhang, Sass, Krajniak, Ganesan, "Ion Correlations and Partial Ionicities in
> the Lamellar Phases of Block Copolymeric Ionic Liquids", *ACS Macro Lett.*
> **11**, 1265--1271 (**2022**).
> [DOI: 10.1021/acsmacrolett.2c00401](https://doi.org/10.1021/acsmacrolett.2c00401)

> Zhang, Krajniak, Sass, Sachar, Marioni, Duncan, Ganesan, "Conductivity and
> Transference Numbers in Lithium Salt-Doped Block Copolymeric Ionic Liquid
> Electrolytes", *Macromolecules* **56**, 9750--9765 (**2023**).
> [DOI: 10.1021/acs.macromol.3c01791](https://doi.org/10.1021/acs.macromol.3c01791)

If you are backmapping **topologies formed by polymerization at the CG level**
(e.g. **hyperbranched polymers**), you may also cite:

> Zhang, Krajniak, Samaey, Nies, "A Parallel Multiscale Simulation Framework
> for Complex Polymerization: AB2-Type Monomer Hyperbranched Polymerization as
> an Example", *Adv. Theory Simul.* **2**, 1800102 (**2019**).
> [DOI: 10.1002/adts.201800102](https://doi.org/10.1002/adts.201800102)

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

GPL-3.0-or-later — see [LICENSE](LICENSE) for the full text.

## Authors

- Jakub Krajniak ([jkrajniak@gmail.com](mailto:jkrajniak@gmail.com))
- Zidan Zhang ([zhangzidan@gmail.com](mailto:zhangzidan@gmail.com))
