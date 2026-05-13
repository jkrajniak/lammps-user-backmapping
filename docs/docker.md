# Docker

Build and run LAMMPS with the backmapping package in a container — no manual
compilation required.

## Building the Image

From the repository root:

```bash
docker build -t lammps-backmap .
```

The default build uses LAMMPS `stable_22Jul2025_update2`. To build against a
different LAMMPS version:

```bash
docker build --build-arg LAMMPS_VERSION=stable_29Aug2024_update1 \
    -t lammps-backmap .
```

!!! note
    The first build compiles LAMMPS from source and takes 5–15 minutes
    depending on the machine. Subsequent rebuilds use Docker layer caching
    and are much faster unless the LAMMPS version or package sources change.

## Running Simulations

Mount your simulation directory to `/work` inside the container:

```bash
docker run --rm -v "$(pwd)":/work lammps-backmap lmp -in in.backmap
```

### MPI (multiple processes)

```bash
docker run --rm -v "$(pwd)":/work lammps-backmap \
    mpirun --allow-run-as-root -np 4 lmp -in in.backmap
```

### Using `backmap-prep`

The `backmap-prep` CLI is included in the image:

```bash
docker run --rm -v "$(pwd)":/work lammps-backmap \
    backmap-prep settings.yaml
```

## What's Included

| Component | Location |
|-----------|----------|
| `lmp` binary (MPI-enabled) | `/usr/local/bin/lmp` |
| `backmap-prep` CLI | on `PATH` (pip-installed) |
| OpenMPI runtime | system packages |
| Python 3 | system packages |

The image uses a multi-stage build: the final runtime image contains only the
compiled binary, MPI runtime, and Python — no compiler toolchain or LAMMPS
source tree.

## HPC: Singularity / Apptainer

Most HPC clusters don't allow Docker but support
[Apptainer](https://apptainer.org/) (formerly Singularity). Convert the Docker
image to a SIF file:

```bash
# Build the Docker image first, then convert
docker build -t lammps-backmap .
docker save lammps-backmap -o lammps-backmap.tar
apptainer build lammps-backmap.sif docker-archive://lammps-backmap.tar
```

Or build directly from a Docker daemon (if available on the build machine):

```bash
apptainer build lammps-backmap.sif docker-daemon://lammps-backmap:latest
```

Then on the cluster:

```bash
apptainer exec lammps-backmap.sif lmp -in in.backmap
```

!!! warning "MPI binding"
    For multi-node MPI runs, the host MPI must be ABI-compatible with the
    container's OpenMPI. Use the `--bind` flag to mount the host MPI into the
    container, or rebuild the image with the same MPI version as your cluster.

    ```bash
    mpirun -np 16 apptainer exec --bind /opt/openmpi:/opt/openmpi \
        lammps-backmap.sif lmp -in in.backmap
    ```

## Verifying the Build

Check that the backmap package is loaded:

```bash
docker run --rm lammps-backmap lmp -h | grep BACKMAP
```

List registered backmap styles:

```bash
docker run --rm lammps-backmap lmp -h | grep backmap
```

Expected output includes `backmap/harmonic` (bond and angle) and
`backmap/table`.
