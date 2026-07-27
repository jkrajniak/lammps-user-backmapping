## ADDED Requirements

### Requirement: Multi-stage Dockerfile builds LAMMPS with backmap package

The repository SHALL contain a `Dockerfile` at the project root that builds
LAMMPS from source with the backmap package enabled (`PKG_BACKMAP=yes`) using a
multi-stage build. The builder stage SHALL compile LAMMPS with MPI support. The
runtime stage SHALL contain only the compiled `lmp` binary, MPI runtime
libraries, and the `backmap-prep` Python CLI.

#### Scenario: Default build produces working LAMMPS binary

- **WHEN** a user runs `docker build -t lammps-backmap .` from the repository root
- **THEN** the resulting image contains a `lmp` binary at `/usr/local/bin/lmp`
  that lists `BACKMAP` in its package output (`lmp -h | grep BACKMAP`)

#### Scenario: Backmap styles are registered

- **WHEN** the built image runs `lmp -h`
- **THEN** the output includes `backmap/harmonic` (bond and angle) and
  `backmap/table` styles

### Requirement: LAMMPS version is configurable via build argument

The Dockerfile SHALL accept a `LAMMPS_VERSION` build argument that specifies
the LAMMPS Git tag to build against. The default value SHALL match the version
used in CI (`stable_22Jul2025_update2`).

#### Scenario: Override LAMMPS version

- **WHEN** a user runs `docker build --build-arg LAMMPS_VERSION=stable_29Aug2024_update1 -t lammps-backmap .`
- **THEN** the image is built against the specified LAMMPS version

### Requirement: Image includes backmap-prep Python CLI

The runtime image SHALL include Python 3 and the `backmap-prep` CLI installed
and available on the `PATH`.

#### Scenario: backmap-prep is available

- **WHEN** a user runs `docker run lammps-backmap backmap-prep --help`
- **THEN** the command succeeds and displays help output

### Requirement: Working directory convention

The image SHALL set `WORKDIR /work` so that users can bind-mount their
simulation directories to `/work` for running simulations.

#### Scenario: Run simulation with volume mount

- **WHEN** a user runs `docker run -v $(pwd):/work lammps-backmap lmp -in in.backmap`
- **THEN** LAMMPS reads input from the mounted directory and writes output there

### Requirement: Docker build context is filtered

The repository SHALL contain a `.dockerignore` file that excludes unnecessary
files (`.git`, `build/`, `lammps/`, `__pycache__/`, `.cursor/`, `openspec/`,
etc.) from the Docker build context.

#### Scenario: Build context excludes large directories

- **WHEN** the Docker build runs
- **THEN** the `.git` directory, any `build/` directories, and cached files are
  excluded from the context sent to the Docker daemon

### Requirement: Docker usage is documented

The documentation site (`docs/`) SHALL include a page describing how to build
and run the Docker image, including volume mounts, MPI usage, and
Singularity/Apptainer conversion. The top-level README SHALL mention Docker as
an installation option.

#### Scenario: Docs site includes Docker page

- **WHEN** a user visits the documentation site
- **THEN** there is a "Docker" page with build, run, and HPC instructions

#### Scenario: README mentions Docker

- **WHEN** a user reads the repository README
- **THEN** there is a section or link pointing to Docker-based usage
