## Why

There is no containerised way to build and run LAMMPS with the backmapping
package. Testing on HPC clusters or remote machines currently requires manually
cloning LAMMPS, running `install.sh`, installing dependencies, and building from
source. A Dockerfile based on an existing LAMMPS image would let users build a
ready-to-run container in one step and ship it to any machine with Docker or
Singularity/Apptainer support.

## What Changes

- Add a `Dockerfile` that starts from an official LAMMPS Docker image, copies
  the backmapping package sources, runs `install.sh`, and recompiles LAMMPS with
  `PKG_BACKMAP=yes`.
- Add a `.dockerignore` to keep the build context lean.
- Add a `docker-compose.yml` (optional convenience) for local build/run.
- Document the Docker workflow in `docs/` and update the README.

## Capabilities

### New Capabilities

- `docker-build`: Dockerfile and supporting files for building a LAMMPS+backmap
  container image, including multi-stage build, example mount points, and
  documentation.

### Modified Capabilities

_(none — this is purely additive infrastructure)_

## Impact

- **New files**: `Dockerfile`, `.dockerignore`, optionally `docker-compose.yml`.
- **Docs**: New "Docker" section in the documentation site and README.
- **CI**: No CI changes in this initial scope (can be added later).
- **Dependencies**: Relies on the official `lammps/lammps` Docker images on
  Docker Hub; no new Python or C++ dependencies.
