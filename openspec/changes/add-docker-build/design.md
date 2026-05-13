## Context

The backmapping package currently requires manual compilation against a LAMMPS
source tree. The CI workflow (`cpp-ci.yml`) already demonstrates the full build
sequence: clone LAMMPS at a pinned tag, run `install.sh`, cmake with
`PKG_BACKMAP=yes`, and build. A Dockerfile can encapsulate this same sequence
into a portable, reproducible image that works on HPC clusters via
Docker/Singularity/Apptainer.

The official `lammps/lammps` Docker Hub images haven't been updated since 2022
and ship outdated LAMMPS versions, so they cannot serve as a reliable base.

## Goals / Non-Goals

**Goals:**

- Single-command build of a container image with LAMMPS + backmap package.
- MPI-enabled build (OpenMPI) for production-scale simulations.
- Pin LAMMPS version to match CI (`stable_22Jul2025_update2`) by default, with
  easy override via build arg.
- Include the `backmap-prep` Python CLI in the same image.
- Minimal final image size via multi-stage build.
- Compatibility with Singularity/Apptainer (common on HPC).

**Non-Goals:**

- GPU/CUDA support (can be added later as a separate Dockerfile variant).
- CI integration (running tests inside Docker) — out of scope for this change.
- Publishing to Docker Hub or any registry.
- Supporting traditional-make builds inside the container (CMake only).

## Decisions

### 1. Build from `ubuntu:22.04` instead of official LAMMPS image

**Rationale**: The official `lammps/lammps` images are 4+ years stale and ship
LAMMPS versions that predate the APIs this package targets. Building from a
clean Ubuntu base gives full control over compiler, MPI, and LAMMPS version.

**Alternative considered**: Multi-FROM with `lammps/lammps` as base — rejected
because the outdated LAMMPS version would require a full rebuild anyway.

### 2. Multi-stage build

**Rationale**: The build stage needs cmake, g++, MPI dev headers, and the full
LAMMPS source tree (~2 GB). The runtime stage only needs the compiled `lmp`
binary, MPI runtime libraries, and Python. A two-stage build keeps the final
image under ~500 MB instead of ~3 GB.

- **Stage 1 (builder)**: Ubuntu 22.04 + build tools + LAMMPS source + compile.
- **Stage 2 (runtime)**: Ubuntu 22.04 + OpenMPI runtime + Python 3 + `lmp`
  binary + `backmap-prep` CLI.

### 3. LAMMPS version as a build argument

**Rationale**: Defaults to the same tag used in CI
(`stable_22Jul2025_update2`), but users can override with
`--build-arg LAMMPS_VERSION=<tag>` for forward compatibility.

### 4. Include `backmap-prep` Python CLI via `uv pip install`

**Rationale**: The Python tooling is the companion to the C++ package. Users
preparing input files on an HPC system benefit from having both tools in a
single container. Install from the local `python/` directory using pip (inside
the container, uv is not required for a simple install).

### 5. Working directory and volume mount convention

The image sets `WORKDIR /work` and documents that users should bind-mount their
simulation directory there:

```
docker run -v $(pwd):/work lammps-backmap mpirun -np 4 lmp -in in.backmap
```

## Risks / Trade-offs

- **[Large image build time]** → Building LAMMPS from source takes 5-15 min.
  Mitigated by multi-stage caching (the builder stage is cached across rebuilds
  if LAMMPS version and package sources haven't changed).
- **[MPI version mismatch on HPC]** → When converting to Singularity, the host
  MPI must be ABI-compatible with the container MPI. Mitigated by using
  OpenMPI, which is the most common HPC MPI, and documenting the bind approach.
- **[Stale LAMMPS version]** → The default pinned version will age. Mitigated
  by the `LAMMPS_VERSION` build arg and keeping the default in sync with CI.
