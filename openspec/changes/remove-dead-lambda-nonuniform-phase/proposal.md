## Why

`fix backmap` evolved from an earlier, broader AdResS-style implementation
(spatially-varying per-particle resolution, ported from ESPResSo++'s
`VerletListAdressInteractionTemplate`) into a narrower-scoped backmapping
tool. A prior change (`2026-07-27-fix-lambda-weighting-model`, archived)
already replaced per-particle `λ_i × λ_j` weighting with a single global
`lambda_global` scalar + per-atom CG-bead membership (`atom2cg`). But three
pieces of the old AdResS-era design were left behind and are now purely
inert:

1. A per-atom `lambda[]` array and `nonuniform yes/no` fix keyword
   (staggered per-atom ramp onset) — zero external callers of the per-atom
   accessor anywhere in `src/`; every interaction style reads only
   `lambda_global` + `atom2cg`. `nonuniform` is `no` in every one of the
   ~45 checked-in example scripts, so the per-atom array has always been
   numerically identical to `lambda_global` for every atom, in every
   validated run.
2. A `phase` fix-keyword reserved in the argument parser but with **no
   handling branch at all** — a user typing `fix bm ... phase 1` today
   gets a hard LAMMPS error. `compute_weight3()`'s `phase` parameter was
   never passed explicitly by any of the 14 production call sites, only
   by 2 unit tests.
3. Stale documentation and specs (`pair-backmap`, `bonded-backmap-styles`)
   that still describe the pre-2026-07-27 `λ_i × λ_j` model, and a
   `fix-backmap-resolution` "Simulation phases" requirement describing a
   two-phase `fix_modify bm phase 1/2` state machine that was never
   implemented.

This is confusing to readers of both the code and the companion CPC paper,
which currently has to explain "atoms carry a per-atom resolution value,
but weights use a single global λ... not a per-particle product" —
confusing precisely because the per-atom value is real state but
functionally inert. Repo owner requested removal, and explicitly
confirmed folding the `phase` cleanup into the same change.

This is a pure refactor: no interaction style's computed weight changes
for any existing validated example (proven by the `nonuniform`-always-`no`
argument above), and `phase` was already unreachable from valid input.

## What Changes

- **BREAKING (LAMMPS input-script syntax)**: `fix backmap` no longer
  accepts the `nonuniform yes/no` keyword. Every checked-in example script
  that specified it (~45 files, always `nonuniform no`) has the token
  removed — a pure syntax fix, zero physics change. The already-unreachable
  `phase` keyword is also removed (a user could never legally invoke it
  before; nothing changes for any working script).
- Remove the per-atom `lambda[]` array's independent state: rename to
  `lambda_display`, a pure mirror of `lambda_global` refreshed
  unconditionally every `end_of_step()`, existing solely to back `f_bm`
  dump/thermo output (`compute_scalar()`, `vector_atom`). No interaction
  style ever reads it.
- Drop `pack_exchange`/`unpack_exchange` and all restart machinery
  (`pack_restart`/`unpack_restart`/`size_restart`/`maxsize_restart`,
  `restart_peratom`) — nothing atom-identity-specific remains to migrate
  or restart; the mirror self-heals from `lambda_global` every step.
  `lambda_global` itself was never restart-persisted before this change
  either — that pre-existing (and previously undocumented) gap is now
  documented instead of silently masked by unrelated per-atom restart code.
- Shrink `pack_forward_comm`/`unpack_forward_comm`'s ghost buffer from 5 to
  4 doubles/atom (drops the lambda slot — nothing ever read a ghost atom's
  lambda).
- Remove `compute_weight3()`'s `phase` parameter (always called at its
  default in production) and simplify `compute_scalar()` to
  `return lambda_global;` (provably equal to the old per-atom mean for
  every real run).
- Remove `nonuniform_lambda` from `backmap-prep`'s `SimulationParams`
  schema and stop emitting the `nonuniform` token from the generated
  `fix backmap` line.
- Rewrite `pair-backmap` and `bonded-backmap-styles` specs' weighting
  requirements to describe the current three-case `compute_weight3()`
  model instead of the superseded `λ_i × λ_j` product. Remove
  `fix-backmap-resolution`'s "Simulation phases" requirement (never
  implemented). Correct `backmap-input-generator`'s two-phase-backmapping
  scenarios to reflect that `simulation.two_phase` is rejected at
  validation time (pre-existing, unrelated to this change) rather than
  targeting the now-removed `phase` keyword.
- Correct `docs/components/fix-backmap.md`'s Restart section, which
  claimed "seamless continuation" — already false before this change
  (`lambda_global` was never restart-persisted) and now cleanly false
  once the per-atom restart code is gone.

## Capabilities

### Modified Capabilities

- `fix-backmap-resolution`: `nonuniform` and `phase` keywords removed from
  command syntax; "Simulation phases" requirement removed; lambda-ramp and
  lambda-output requirements rewritten around a single `lambda_global`
  with no per-atom state.
- `pair-backmap`: weighting requirements rewritten from the superseded
  `λ_i × λ_j` product model to the current three-case `compute_weight3()`
  model (already implemented since 2026-07-27; the spec had not been
  updated to match).
- `bonded-backmap-styles`: same weighting-model correction as
  `pair-backmap`, across bond/angle/dihedral `backmap/*` styles and
  `fix backmap/pairs`; `phase`-parameter requirements removed.
- `backmap-input-generator`: `nonuniform_lambda` removed from the
  `simulation` schema and fix-line syntax template; two-phase-backmapping
  scenarios corrected to describe actual (rejected) behavior instead of
  targeting the removed `phase` keyword.

## Impact

- **C++ source** (`lammps-user-backmapping/src/`): `fix_backmap.h`/`.cpp`,
  `backmap_lambda.h`, `backmap_lambda_weights.h`.
- **C++ unit tests**: `tests/unit/test_backmap_lambda.cpp` (drop the
  phase-1-override test).
- **Python** (`python/src/backmap_prep/`): `schema.py`, `writers.py`.
- **Examples**: ~45 `.in.*` files (token removed from the `fix bm ...`
  line), 14 `settings.yaml` files (`nonuniform_lambda` line removed).
- **Docs**: `docs/components/fix-backmap.md`, `docs/settings-reference.md`,
  `docs/tutorial-new-system.md`.
- **CHANGELOG.md**: new entry under `## [Unreleased]`.
- No breaking change to force/energy output for any existing validated
  example — every real run always had `nonuniform no`, and `phase` was
  already unreachable.
- **Paper follow-up** (separate, out of scope here): the companion CPC
  paper's §3.1 can drop its "per-atom resolution value... not a
  per-particle product" explanatory aside once this ships.
