## 1. Schema & Input Writer

- [x] 1.1 Add `restart_interval: int | None = None` field to `SimulationParams` in `schema.py` with a validator rejecting non-positive values
- [x] 1.2 Extract shared setup block from `write_lammps_input` into a helper `_write_setup(f, system, settings)` in `writers.py`
- [x] 1.3 Add `_write_phase_script(path, phase, system, settings, data_filename)` that generates a per-phase input script using `read_restart` + `include in.<prefix>.setup`
- [x] 1.4 Modify `write_lammps_input` to emit `restart N restart.backmap restart.backmap2` and `write_restart` + sentinel commands when `restart_interval` is set
- [x] 1.5 Generate per-phase scripts (`in.<prefix>.setup`, `.phase1`, `.phase2`, `.phase3`) from `write_lammps_input` when `restart_interval` is set
- [x] 1.6 Add unit tests for restart-related generation (restart commands present/absent, per-phase file generation, sentinel commands)

## 2. Entrypoint Script

- [x] 2.1 Create `scripts/run-backmap.sh` — detects restart files, checks phase sentinels, runs correct LAMMPS input script with optional MPI
- [x] 2.2 Update `Dockerfile` to copy `scripts/run-backmap.sh` into `/usr/local/bin/` in the runtime stage

## 3. Cloud Batch Example

- [x] 3.1 Create `examples/cloud-batch/job.json` — Cloud Batch job template with spot VM, GCS mount, `maxRetryCount`, and `run-backmap.sh` entrypoint
- [x] 3.2 Create `examples/cloud-batch/README.md` — setup instructions (Artifact Registry push, GCS bucket, job submission)

## 4. Documentation

- [x] 4.1 Add `docs/cloud-hpc.md` page covering restart configuration, entrypoint script, Cloud Batch setup, and traditional HPC (Singularity/Apptainer)
- [x] 4.2 Add "Running on Cloud / HPC" entry to `mkdocs.yml` navigation
- [x] 4.3 Update `README.md` to mention restart support and cloud deployment
- [x] 4.4 Add entries to `CHANGELOG.md` under `[Unreleased]` → `Added`
