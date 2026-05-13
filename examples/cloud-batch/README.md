# Running on Google Cloud Batch (Spot Instances)

This example demonstrates running a backmapping simulation on Cloud Batch with
spot (preemptible) VMs that automatically restart on preemption.

## Prerequisites

- Google Cloud project with Batch API enabled
- `gcloud` CLI installed and configured
- Docker image pushed to Artifact Registry

## Setup

### 1. Build and push the Docker image

```bash
# From the repository root
docker build -t lammps-backmap .

# Tag and push to Artifact Registry
export PROJECT_ID=your-project-id
export REGION=us-central1
export REPO=your-repo

docker tag lammps-backmap ${REGION}-docker.pkg.dev/${PROJECT_ID}/${REPO}/lammps-backmap:latest
docker push ${REGION}-docker.pkg.dev/${PROJECT_ID}/${REPO}/lammps-backmap:latest
```

### 2. Prepare simulation data in GCS

```bash
export BUCKET=your-bucket-name

# Generate input files locally
cd examples/dodecane
backmap-prep settings.yaml

# Upload to GCS (include restart_interval in settings.yaml!)
gsutil -m cp *.data in.* *.table gs://${BUCKET}/simulation-data/
```

Make sure your `settings.yaml` includes `restart_interval` to enable
checkpointing:

```yaml
simulation:
  restart_interval: 5000
  # ... other settings
```

### 3. Submit the job

Edit `job.json` — replace `PROJECT_ID`, `REGION`, `REPO`, and `BUCKET_NAME`
with your values.

```bash
gcloud batch jobs submit backmap-run-1 \
    --location=us-central1 \
    --config=job.json
```

### 4. Monitor

```bash
gcloud batch jobs describe backmap-run-1 --location=us-central1
gcloud batch jobs list --location=us-central1
```

## How restart works

1. `run-backmap.sh` is the container entrypoint.
2. LAMMPS writes `restart.backmap` / `restart.backmap2` every N steps and at
   each phase boundary.
3. Sentinel files (`phase_1.done`, `phase_2.done`, `phase_3.done`) mark
   completed phases.
4. If the VM is preempted, Cloud Batch retries the task (up to
   `maxRetryCount` times).
5. On retry, `run-backmap.sh` finds the restart file and sentinels on the GCS
   mount and resumes from the correct phase.

## Customizing

| Field in `job.json` | Description |
|---------------------|-------------|
| `machineType` | VM size (e.g., `c2-standard-8` for 4 physical cores) |
| `provisioningModel` | `SPOT` for preemptible, `STANDARD` for on-demand |
| `maxRetryCount` | Number of retries on preemption (3 is reasonable) |
| `maxRunDuration` | Wall-clock timeout per attempt |
| `-np 4` | MPI processes — match to available vCPUs |
