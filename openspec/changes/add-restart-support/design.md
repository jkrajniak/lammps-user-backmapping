## Context

The backmapping simulation has a 3-phase sequential structure:

1. **CG equilibration** — lambda fixed at 0, CG forces only
2. **Backmapping** — lambda ramps 0→1 over `1/alpha` steps
3. **AT production** — lambda = 1, full atomistic dynamics

All C++ styles (`fix_backmap`, `pair_backmap`, bond and angle styles) already
implement LAMMPS restart support — per-atom lambda is packed/unpacked, pair and
bonded coefficients are serialized. But the generated input scripts contain no
`write_restart` commands and no logic to resume from a restart file.

For preemptible cloud instances (Google Cloud Batch spot VMs, AWS Spot, etc.),
the simulation must be able to checkpoint at regular intervals and resume from
the most recent checkpoint after preemption.

## Goals / Non-Goals

**Goals:**

- Generated input scripts include periodic `write_restart` commands.
- A shell entrypoint script detects existing restart files and resumes from the
  correct simulation phase.
- A Cloud Batch example demonstrates the full workflow (GCS data, spot VMs).
- Existing workflows are unaffected — restart is opt-in via a new setting.

**Non-Goals:**

- Modifying the C++ restart implementation (already complete).
- Auto-splitting into 3 separate jobs (too complex for this change).
- Supporting LAMMPS `read_restart` with different MPI decompositions
  (out of scope — same `-np` assumed).

## Decisions

### 1. Use LAMMPS variable-based phase tracking

LAMMPS has no built-in "which phase am I in" state that survives a restart.
To track the current phase, the input script uses a LAMMPS `variable` that is
saved in the restart file via `fix store/state` or, simpler, a **file-based
sentinel** approach:

- At the start of each phase, the script writes a sentinel file
  (`phase_1.done`, `phase_2.done`).
- The entrypoint script checks which sentinels exist to determine the resume
  point.
- On resume, the entrypoint calls LAMMPS with `read_restart` and a
  phase-specific continuation script.

**Rationale**: File sentinels are simpler and more debuggable than trying to
encode phase state inside LAMMPS variables. They also work with any job
scheduler.

**Alternative considered**: Single monolithic script with LAMMPS `if` variables
— rejected because LAMMPS scripting is limited and the resulting input file
would be fragile and hard to read.

### 2. Generate per-phase input scripts

Instead of one monolithic `in.backmap`, `backmap-prep` generates:

- `in.backmap` — master script (runs all 3 phases sequentially, default).
- `in.backmap.phase1` — CG equilibration only.
- `in.backmap.phase2` — backmapping only (can `read_restart`).
- `in.backmap.phase3` — AT production only (can `read_restart`).
- `in.backmap.restart` — restart dispatcher (determines which phase to resume).

The master `in.backmap` remains for non-restart use (backward compatible). The
per-phase scripts share the same pair/bond/angle/fix setup via an
`in.backmap.setup` include file.

**Rationale**: LAMMPS `include` and `read_restart` work well with this pattern.
Each phase script is self-contained and testable. The restart dispatcher is a
thin LAMMPS script that checks for the restart file and includes the right
phase.

### 3. `restart_interval` setting with smart defaults

Add `restart_interval: int | None` to `SimulationParams`:

- `None` (default): no `write_restart` commands — backward compatible.
- An integer N: emit `restart N restart.backmap restart.backmap2` in each phase
  (LAMMPS alternating restart files for crash safety).

When `restart_interval` is set, also emit `write_restart` at the end of each
phase to guarantee a clean checkpoint between phases.

### 4. Entrypoint script: `scripts/run-backmap.sh`

A bash script that:

1. Checks for `restart.backmap` or `restart.backmap2` files.
2. If no restart: runs `lmp -in in.backmap` (fresh start).
3. If restart exists: determines the phase from sentinel files, runs the
   appropriate `in.backmap.phaseN` with `read_restart`.

This script is copied into the Docker image and documented for HPC use.

### 5. Cloud Batch example

A minimal `examples/cloud-batch/` directory with:

- `job.json` — Cloud Batch job template (spot VM, GCS mount, restart on
  preemption via `maxRetryCount`).
- `README.md` — Setup instructions (push image to Artifact Registry, configure
  GCS bucket, submit job).

## Risks / Trade-offs

- **[Restart file size]** → LAMMPS restart files scale with atom count. For
  large systems (millions of atoms), files can be hundreds of MB. Mitigated by
  using LAMMPS's alternating restart mechanism (two files, always one valid).
- **[Phase boundary edge case]** → If preempted exactly at a phase boundary
  before the sentinel is written, the restart may re-run a few steps of the
  previous phase. Acceptable since LAMMPS phases are idempotent at the boundary.
- **[GCS mount latency]** → Writing restart files to GCS-FUSE can be slow for
  large files. Document that users should write to local SSD and sync to GCS
  periodically if this is an issue.
