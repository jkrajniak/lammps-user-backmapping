## 1. Isolate the work

- [x] 1.1 Create a fresh git worktree `.worktrees/remove-dead-lambda-nonuniform-phase`
  on branch `refactor/remove-dead-lambda-nonuniform-phase` off
  `secure-upstream/main`, separate from the dirty, unrelated
  `fix/pair-table-core-extension-knots` worktree.

## 2. C++ core (`src/`)

- [x] 2.1 `fix_backmap.h`: remove `int nonuniform;`; rename
  `double *lambda;` → `double *lambda_display;`; remove
  `pack_exchange`/`unpack_exchange`/`pack_restart`/`unpack_restart`/
  `size_restart`/`maxsize_restart` declarations.
- [x] 2.2 `fix_backmap.cpp`: drop `"nonuniform"` and `"phase"` from
  `KNOWN_KEYWORDS`; remove the `nonuniform` keyword-parsing branch and the
  `RanMars`-based staggering block; `comm_forward = 4` (was 5); drop
  `restart_peratom`/`Atom::RESTART` callback registration and its
  destructor teardown.
- [x] 2.3 `end_of_step()`: decouple the mirror refresh from `ramp_active`
  (unconditional every call); gate only the `lambda_global` increment.
- [x] 2.4 `compute_scalar()`: simplify to `return lambda_global;`.
- [x] 2.5 `extract()`: drop the `"lambda"` case.
- [x] 2.6 `grow_arrays()`/`copy_arrays()`: rename `lambda`→`lambda_display`
  throughout.
- [x] 2.7 Delete `pack_exchange`/`unpack_exchange`/`pack_restart`/
  `unpack_restart`/`size_restart`/`maxsize_restart` method bodies.
- [x] 2.8 `pack_forward_comm`/`unpack_forward_comm`: drop the lambda slot,
  shrink to 4 doubles/ghost atom.
- [x] 2.9 Fix the stale `pre_force()`/`post_force()` comments describing
  the forward-comm call as being about lambda.
- [x] 2.10 `backmap_lambda.h`: delete `extract_lambda()`.
- [x] 2.11 `backmap_lambda_weights.h`: remove `compute_weight3()`'s `phase`
  parameter; fix `clamp_lambda()`'s comment.
- [x] 2.12 `tests/unit/test_backmap_lambda.cpp`: delete the phase-1-override
  test.

## 3. Python (`python/src/backmap_prep/`)

- [x] 3.1 `schema.py`: delete `nonuniform_lambda` from `SimulationParams`.
- [x] 3.2 `writers.py`: delete the `nonuniform` local and its fragment from
  the generated `fix_line`.

## 4. Local verification

- [x] 4.1 `uv run pytest` (python/): 276 passed, 7 skipped (pre-existing
  skips), 0 failures.
- [x] 4.2 `cmake --build tests/unit/build && ctest --test-dir
  tests/unit/build`: 20/20 passed (was 22 before dropping the phase test).

## 5. Examples (~45 files)

- [x] 5.1 Spot-checked regenerating `examples/dodecane/in.dodecane` via
  `backmap-prep build`: produced a large unrelated diff (halved bond/angle
  K constants, added `comm_modify`/`neigh_modify`/velocity lines, changed
  `thermo_style`) — the checked-in examples are not kept in sync with the
  tool and reflect other already-landed fixes never reapplied here.
  Reverted; regeneration is unsafe without a full revalidation pass,
  out of scope for this change.
- [x] 5.2 Mechanically stripped the ` nonuniform (yes|no)` token from all
  45 affected `.in.*`/`in.*` files instead — confirmed exactly one
  occurrence per file, one line changed per file, zero other content
  change.
- [x] 5.3 Removed the stale `nonuniform_lambda:` line from all 14 affected
  `settings.yaml` files (safe: `SimulationParams` has no `model_config`
  override, pydantic's default `extra='ignore'` applies).

## 6. Docs

- [x] 6.1 `docs/components/fix-backmap.md`: dropped `[nonuniform yes/no]`
  from syntax/keyword docs; deleted the `### nonuniform (optional)`
  section and the `extract("lambda", dim)` bullet; corrected the Restart
  section's "seamless continuation" claim.
- [x] 6.2 `docs/settings-reference.md`: deleted
  `#### simulation.nonuniform_lambda`.
- [x] 6.3 `docs/tutorial-new-system.md`: deleted the `nonuniform_lambda:
  false` example line.

## 7. openspec specs

- [x] 7.1 `fix-backmap-resolution/spec.md`: rewrote the lambda-ramp and
  lambda-output requirements around `lambda_global`; removed the
  "Simulation phases" requirement (never implemented); corrected the fix
  command syntax requirement.
- [x] 7.2 `pair-backmap/spec.md`: rewrote the weighting requirements from
  the superseded `λ_i × λ_j` product to the current three-case
  `compute_weight3()` model.
- [x] 7.3 `bonded-backmap-styles/spec.md`: same weighting-model correction
  across bond/angle/dihedral `backmap/*` styles and `fix backmap/pairs`;
  removed `phase`-parameter requirements; restored three "Restart file
  support" scenarios (bond/dihedral coefficient restart — a separate,
  unrelated, still-accurate LAMMPS mechanism) that were dropped in an
  earlier draft of this rewrite and needed putting back.
- [x] 7.4 `backmap-input-generator/spec.md`: removed `nonuniform_lambda`
  from the example YAML and fix-line syntax template; corrected the
  two-phase-backmapping scenarios to describe the actual (validation-time
  rejection) behavior instead of targeting the removed `phase` keyword,
  without touching the separately-scoped `two_phase` schema guard itself.
- [x] 7.5 `pe-examples-findings.md`: updated the one example fix line.
- [x] 7.6 Confirmed via `openspec validate` that this change's delta specs
  introduce no new validation errors beyond a pre-existing, repo-wide
  "missing ## Purpose section" issue that predates this change and affects
  most specs in `openspec/specs/` (including many never touched here) —
  out of scope to fix as part of this change.

## 8. CHANGELOG

- [x] 8.1 Added an entry under `## [Unreleased]` documenting the removal,
  cross-referencing the earlier "Lambda-weighting formula" fix this
  follows up.

## 9. VM verification (pending, requires `sc@<vm_host>` — never build/run
   LAMMPS locally per this workspace's standing rule)

- [ ] 9.1 Rebuild the LAMMPS package with the modified C++.
- [ ] 9.2 Rerun `examples/dodecane/` and confirm numerically identical
  output to a pre-change baseline.
- [ ] 9.3 Rerun `examples/pe/large/in.pe_robust` and confirm the `f_bm`
  thermo column is identical before/after.
- [ ] 9.4 Rerun the MPI serial-vs-4-rank parity test if available (the
  `comm_forward` shrink touches wire-level MPI behavior).
