# Resumable runs on spot or preemptible machines

The long atomistic stages of the large examples (days on a workstation) can be interrupted and
continued without losing more than one checkpoint interval. This page describes the pieces and
how to use them.

## What is resumable

| Run | Input | Checkpoint |
|-----|-------|------------|
| Tier C continuation of a network or of the melamine melt | `examples/common/in.at_reference_protocol` | every `ckpt` steps (default 50 000 = 50 ps) |
| Tier C continuation of the dodecane, pe4 and pe_aa melts | `examples/<example>/large/in.<example>_at` | every `ckpt` steps (default 200 000, one RDF block) |
| Matched-density energy control | `scripts/review/b3_matched_density.sh` | every `CKPT` steps (default 50 000) |
| Dodecane degraded-CG experiment (B8) | `scripts/review/b8_degraded_cg.sh` | stage markers, atomistic stage as above |

The hybrid (backmapping) stage, which takes minutes to an hour, is not restartable: after an
interruption it runs again from its start. Reference runs (independent atomistic melts) run to
completion or run again.

## How an input resumes

Every stage runs to an **absolute end step** (`run N upto`) and a restart file is written every
`ckpt` steps (`rst.s<seg>.a`, `rst.s<seg>.b`, alternating). To continue, an input is started with

```
-var resume 1 -var resume_file <restart file> -var seg <segment number>
```

It reads the restart file, defines the style and force-field settings again, and jumps to the
stage that contains the current step. Output files of a segment carry the segment number
(`dump.at_prod.s2`, `rdf_backmap.s2.dat`, ...) so that a continuation never overwrites earlier data.

## The wrapper

`scripts/review/resumable_lmp.sh` does all of this:

```bash
LMP=/path/to/lmp NP=8 scripts/review/resumable_lmp.sh at -in in.protocol -log log.run.lammps -var prefix pet
```

Run the same command again after an interruption. It picks the newest valid restart file (every
candidate is read by LAMMPS to get its step, so a half-written file is ignored; the two newest are
kept), numbers the segment, writes one log per segment (`log.run.s2.lammps`), and merges the
per-segment dump, `fix ave/time` and `fix ave/histo` files at the end
(`scripts/review/merge_segments.py`: sorted by timestep, a step that occurs twice is taken from the
latest segment). After two failures from the same checkpoint it gives up (`.at.fail`), so that a
deterministic crash does not loop. The step of a restart is independent of the number of MPI ranks:
a continuation may use a different machine and rank count.

The example scripts `run_tier_bc.sh` use it when `RESUMABLE=<path to resumable_lmp.sh>` is set, and
write a `.done.<step>` marker after each finished step.

To analyse energies from the logs of a resumed run, give `scripts/review/b3_energy_terms.py` the
glob of the segment logs and `--from-step` (the end of the equilibration).

## A queue for a spot machine

```bash
scripts/review/queue_runner.sh <queue-dir>
```

`<queue-dir>/queue.txt` has one job per line, `name|source-dir|command`; `queue.env` holds `LMP`,
`NP`, `PREP`, `AT_IN`, ... Starting the runner again (after a reboot, on a new machine, or by hand)
is always safe: only one runner per queue directory, finished jobs are skipped, a job that failed is
given up (`.job.failed`), and an unfinished job continues from its newest checkpoint. `SLOTS` jobs
run at once.

For a machine that can disappear, copy the queue directory to persistent storage while it runs:

```bash
SYNC_DEST=user@host:/path SYNC_EVERY=600 scripts/review/queue_runner.sh <queue-dir>
```

(`scripts/review/ckpt_sync.sh`: restart files, logs, markers and results; trajectories are left
out; nothing is ever deleted at the destination). On the replacement machine:

```bash
rsync -a user@host:/path/ <queue-dir>/ && scripts/review/queue_runner.sh <queue-dir>
```

A systemd unit that starts the runner at boot is in `scripts/review/lub-queue.service`.

## What a continuation changes

A restart is exact for the state LAMMPS keeps in the restart file (positions, velocities, box,
Nose-Hoover thermostat and barostat variables). Two small differences remain:

- the PPPM grid is chosen again from the current box, so energies right after a restart differ in
  the fifth digit at the 1e-5 accuracy setting;
- the annealing ramp of the RIM135 protocol runs in chunks of `ckpt` steps with a fresh thermostat per
  chunk, also without any interruption.

In tests (kills with SIGKILL during the stages) the thermodynamics of a resumed run agree with the
uninterrupted run (temperature 1e-6, potential energy 4e-5, volume 1e-6 relative; positions diverge
as in any two chaotic runs), and for a small melt without PPPM the RDF and histogram blocks were
identical.
