## ADDED Requirements

### Requirement: Per-phase input scripts with shared setup

When `restart_interval` is set, `backmap-prep` SHALL generate per-phase input
scripts in addition to the master `in.<prefix>`:

- `in.<prefix>.setup` — shared pair/bond/angle/fix/group definitions.
- `in.<prefix>.phase1` — CG equilibration (includes setup, runs phase 1).
- `in.<prefix>.phase2` — backmapping ramp (reads restart, includes setup, runs
  phase 2).
- `in.<prefix>.phase3` — AT production (reads restart, includes setup, runs
  phase 3).

#### Scenario: Per-phase scripts generated when restart enabled

- **WHEN** `restart_interval: 5000` is set in the YAML settings
- **THEN** `backmap-prep` generates `in.backmap.setup`, `in.backmap.phase1`,
  `in.backmap.phase2`, `in.backmap.phase3` alongside the master `in.backmap`

#### Scenario: Per-phase scripts not generated when restart disabled

- **WHEN** `restart_interval` is not set (default `None`)
- **THEN** `backmap-prep` generates only the master `in.backmap` (backward
  compatible)

### Requirement: Periodic write_restart commands in generated scripts

When `restart_interval` is set, the generated input scripts SHALL include
LAMMPS `restart N file1 file2` commands that write alternating restart files
during each phase, plus a `write_restart` at the end of each phase.

#### Scenario: Restart commands in phase scripts

- **WHEN** `restart_interval: 5000` is set
- **THEN** each phase script contains `restart 5000 restart.backmap restart.backmap2`
  and a `write_restart restart.backmap` after the `run` command

#### Scenario: Phase sentinel files

- **WHEN** a phase completes
- **THEN** the script writes a sentinel file (e.g., `shell echo done > phase_1.done`)
  to indicate that the phase finished successfully

### Requirement: Restart-aware entrypoint script

The repository SHALL contain a `scripts/run-backmap.sh` shell script that:

1. Checks for existing restart files (`restart.backmap`, `restart.backmap2`).
2. If no restart file exists, runs `lmp -in in.backmap` (fresh start via master
   script).
3. If a restart file exists, checks sentinel files (`phase_1.done`,
   `phase_2.done`) to determine the current phase, then runs the appropriate
   `in.<prefix>.phaseN` script.
4. Passes through MPI arguments (number of processes).

#### Scenario: Fresh start (no restart file)

- **WHEN** `run-backmap.sh` is invoked and no `restart.backmap*` files exist
- **THEN** it runs `lmp -in in.backmap` (or with `mpirun` if `-np` is provided)

#### Scenario: Resume after preemption in phase 2

- **WHEN** `run-backmap.sh` is invoked, `restart.backmap` exists, and
  `phase_1.done` exists but `phase_2.done` does not
- **THEN** it runs `lmp -in in.backmap.phase2` which reads the restart file
  and continues the backmapping ramp

#### Scenario: Resume after preemption in phase 3

- **WHEN** `run-backmap.sh` is invoked, `restart.backmap` exists, and both
  `phase_1.done` and `phase_2.done` exist but `phase_3.done` does not
- **THEN** it runs `lmp -in in.backmap.phase3`

#### Scenario: All phases complete

- **WHEN** `run-backmap.sh` is invoked and `phase_3.done` exists
- **THEN** it prints "All phases complete" and exits 0

### Requirement: Docker image includes entrypoint script

The `Dockerfile` SHALL copy `scripts/run-backmap.sh` into the runtime image at
`/usr/local/bin/run-backmap.sh` so it is available on `PATH`.

#### Scenario: Entrypoint available in container

- **WHEN** a user runs `docker run lammps-backmap run-backmap.sh --help`
- **THEN** the script displays usage information

### Requirement: Cloud Batch job example

The repository SHALL contain an `examples/cloud-batch/` directory with:

- `job.json` — a Google Cloud Batch job template configured for spot VMs with
  `maxRetryCount` for automatic restart on preemption, a GCS volume mount for
  input/output, and the `run-backmap.sh` entrypoint.
- `README.md` — instructions for pushing the Docker image to Artifact Registry,
  configuring the GCS bucket, and submitting the job.

#### Scenario: Cloud Batch example exists

- **WHEN** a user navigates to `examples/cloud-batch/`
- **THEN** they find `job.json` and `README.md` with complete setup instructions

### Requirement: Cloud/HPC documentation page

The documentation site SHALL include a "Running on Cloud / HPC" page
(`docs/cloud-hpc.md`) covering restart configuration, the entrypoint script,
Cloud Batch setup, and Singularity/Apptainer on traditional HPC clusters.

#### Scenario: Docs site includes cloud/HPC page

- **WHEN** a user visits the documentation site
- **THEN** there is a "Running on Cloud / HPC" page with restart, Cloud Batch,
  and HPC instructions
