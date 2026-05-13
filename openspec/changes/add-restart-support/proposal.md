## Why

Running backmapping simulations on cloud spot instances (e.g., Google Cloud
Batch) requires the ability to checkpoint and resume — spot VMs can be preempted
at any time. The C++ styles already implement full LAMMPS restart support
(`write_restart`/`read_restart` for fix, pair, bond, and angle), but the
generated input scripts contain no `write_restart` commands, and there is no
entrypoint logic to detect a restart file and resume from it.

## What Changes

- Add `restart_interval` setting to the simulation parameters schema so
  `backmap-prep` can emit periodic `write_restart` commands in the generated
  input script.
- Modify the LAMMPS input writer to emit `write_restart` at the end of each
  phase and optionally every N steps within a phase.
- Add a restart-aware entrypoint shell script (`scripts/run-backmap.sh`) that
  detects an existing restart file and either starts fresh or resumes from the
  correct phase.
- Add a Cloud Batch job configuration example with documentation.

## Capabilities

### New Capabilities

- `restart-checkpointing`: Restart-aware input script generation, entrypoint
  script for preemptible environments, and Cloud Batch job example.

### Modified Capabilities

- `backmap-input-generator`: Adding `restart_interval` to simulation parameters
  and `write_restart`/`read_restart` commands to the generated input script.

## Impact

- **Python code**: `schema.py` (new field), `writers.py` (restart commands in
  generated script).
- **New files**: `scripts/run-backmap.sh`, `examples/cloud-batch/` job config.
- **Docker**: Dockerfile updated to include the entrypoint script.
- **Docs**: New "Running on Cloud / HPC" documentation page.
