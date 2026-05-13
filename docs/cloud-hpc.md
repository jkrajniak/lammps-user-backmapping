# Running on Cloud / HPC

Run backmapping simulations on cloud spot instances or HPC clusters with
checkpoint/restart support for fault tolerance.

## Enabling Restart Checkpoints

Add `restart_interval` to your `settings.yaml`:

```yaml
simulation:
  restart_interval: 5000   # write checkpoint every 5000 steps
  alpha: 0.001
  # ... other settings
```

When set, `backmap-prep` generates:

| File | Purpose |
|------|---------|
| `in.backmap` | Master script (includes restart commands) |
| `in.backmap.setup` | Shared style/coefficient definitions |
| `in.backmap.phase1` | First segment: CG equilibration if `equilibration_steps > 0`, else backmapping λ 0→1 |
| `in.backmap.phase2` | Only if `equilibration_steps > 0`: backmapping ramp; **or** if `equilibration_steps == 0` and `production_steps > 0`: optional AT production |
| `in.backmap.phase3` | Only if `equilibration_steps > 0` and `production_steps > 0`: optional AT production |

With `equilibration_steps: 0` and `production_steps: 0` (default), only `phase1`
is generated for per-phase restarts (backmapping only, writes `*_hybrid.data`).

Without `restart_interval`, only the master `in.backmap` is generated (backward
compatible).

### How checkpointing works

1. LAMMPS writes alternating restart files (`restart.backmap`,
   `restart.backmap2`) every N steps — if one is corrupted mid-write, the other
   is still valid.
2. At each phase boundary, a clean `write_restart` and a sentinel file
   (`phase_1.done`, `phase_2.done`, `phase_3.done`) are written.
3. On restart, the entrypoint script checks sentinels to determine which phase
   to resume.

## Entrypoint Script

The `run-backmap.sh` script handles restart logic automatically:

```bash
# Fresh start
run-backmap.sh -np 4 -in in.backmap

# After preemption — detects restart files and resumes
run-backmap.sh -np 4 -in in.backmap
```

The script:

1. Checks for `restart.backmap` / `restart.backmap2`.
2. If no restart file: runs `in.backmap` (fresh start).
3. If restart exists: reads sentinel files to determine the phase, runs the
   appropriate `in.backmap.phaseN`.

## Google Cloud Batch (Spot VMs)

Cloud Batch runs containerised jobs on spot (preemptible) VMs with automatic
retry on preemption.

### Quick start

```bash
# 1. Build and push image
docker build -t lammps-backmap .
docker tag lammps-backmap REGION-docker.pkg.dev/PROJECT/REPO/lammps-backmap:latest
docker push REGION-docker.pkg.dev/PROJECT/REPO/lammps-backmap:latest

# 2. Upload simulation data to GCS
gsutil -m cp *.data in.* *.table gs://BUCKET/sim-data/

# 3. Submit job
gcloud batch jobs submit my-backmap-job \
    --location=us-central1 \
    --config=examples/cloud-batch/job.json
```

See [`examples/cloud-batch/`](https://github.com/jkrajniak/lammps-user-backmapping/tree/main/examples/cloud-batch)
for the full job template and setup instructions.

### Key settings in `job.json`

| Setting | Value | Purpose |
|---------|-------|---------|
| `provisioningModel` | `SPOT` | Use preemptible VMs (up to 90% cheaper) |
| `maxRetryCount` | `3` | Retry on preemption |
| `maxRunDuration` | `14400s` | 4-hour timeout per attempt |
| GCS volume mount | `/work` | Persistent storage survives preemption |

## Traditional HPC (Singularity / Apptainer)

For clusters without Docker support, convert to a SIF image:

```bash
docker save lammps-backmap -o lammps-backmap.tar
apptainer build lammps-backmap.sif docker-archive://lammps-backmap.tar
```

### Running with restart on HPC

Copy `scripts/run-backmap.sh` to your working directory along with the
generated input files, then use it in your job script:

```bash
#!/bin/bash
#SBATCH --job-name=backmap
#SBATCH --ntasks=64
#SBATCH --time=04:00:00
#SBATCH --requeue

srun apptainer exec lammps-backmap.sif \
    run-backmap.sh -np 64 -in in.backmap
```

With `--requeue`, Slurm resubmits the job if it's preempted. On restart,
`run-backmap.sh` detects the checkpoint and resumes from the correct phase.

!!! warning "MPI compatibility"
    For multi-node runs, the host MPI must be ABI-compatible with the
    container's OpenMPI. See the [Docker page](docker.md) for binding
    instructions.
