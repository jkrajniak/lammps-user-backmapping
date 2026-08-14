## Context

See `proposal.md` for the full rationale. This section records the
verification trail and the specific design decisions made while scoping
the removal, since several of them are not obvious from the diff alone.

## Verification that this is a pure refactor

- `grep`'d every `.cpp`/`.h` in `src/` for callers of
  `BackmapLambda::extract_lambda()` (the per-atom accessor): zero, outside
  its own definition. All 11 interaction style source files
  (`pair_backmap.cpp`, the six `bond/angle/dihedral_backmap_*.cpp` files,
  `fix_backmap_pairs.cpp`, `fix_backmap_capforce.cpp`) call only
  `extract_lambda_global()` + `extract_atom2cg()`.
- `nonuniform` is `no` in every one of the ~45 checked-in example scripts
  (confirmed by grep across `examples/`). Since the constructor
  initializes `lambda[i] = lambda0` uniformly when `nonuniform` is unset,
  and `end_of_step()` increments `lambda[i]` and `lambda_global` by the
  same `alpha` in lockstep, `lambda[i] == lambda_global` for every atom,
  on every step, in every real run — provably, not just empirically. This
  is why `compute_scalar()`'s old per-atom mean and the new
  `return lambda_global;` are exactly equal for every existing case.
- `phase` has no handling branch in the constructor's keyword parser
  (confirmed: only one occurrence of the string `"phase"` in
  `fix_backmap.cpp`, in `KNOWN_KEYWORDS`, with no matching
  `else if (strcmp(arg[iarg], "phase") == 0)`). A user specifying it today
  falls through to `error->all(FLERR, "Illegal fix backmap argument: {}",
  arg[iarg])`. `compute_weight3()`'s `phase` parameter defaults to `2` and
  is never passed a different value by any of the 14 production call
  sites — only 2 unit tests exercised `phase=1`.
- No C++ unit test or Python test references `nonuniform`, `vector_atom`,
  `compute_scalar`, or the per-atom lambda array before this change.

## Design decision: keep a display mirror, don't remove `f_bm` output

`examples/pe/large/in.pe_robust` uses `thermo_style custom ... f_bm` +
`thermo_modify colname f_bm lambda` (the fix's scalar output), and
`examples/pe/in.pe` / `examples/pe4/in.pe4` use
`dump traj all custom 1000 dump.backmap id mol type x y z f_bm` (the fix's
per-atom vector output). Removing `vector_atom`/`compute_scalar` entirely
would break these real, checked-in scripts' column layout. Instead:

- `compute_scalar()` → `return lambda_global;` directly (see the
  equivalence proof above — behavior-preserving for every real case, and
  now cheaper: no more MPI-reduced mean).
- The per-atom array (renamed `lambda_display`) stays, but purely as a
  broadcast mirror: `end_of_step()` fills every local slot with the
  current `lambda_global`, **unconditionally**, not gated on
  `ramp_active`. This is what makes dropping `pack_exchange` safe — the
  mirror self-heals after any atom migration (reneighboring, `fix
  balance`, etc.) instead of relying on exchange to carry per-atom
  identity across ranks. Gating the mirror refresh on `ramp_active` (as a
  naive port of the old `if (ramp_active) { lambda[i] += alpha; ... }`
  loop would do) would leave stale/uninitialized values in newly-migrated
  slots during a frozen ramp (e.g. the documented equilibration phase,
  `fix_modify bm active no`) until the ramp reactivates.
- Ghost-atom forwarding of the mirror is dropped entirely (shrinking
  `pack_forward_comm`'s buffer from 5 to 4 doubles/atom) since dump/thermo
  only ever read local per-atom vectors, never ghost ones — nothing on the
  receiving end of that slot was ever reading it.

## Design decision: drop restart machinery entirely, don't try to preserve it

`lambda_global` — the only value any interaction style weights by — was
**not** restart-persisted before this change either (`restart_global = 0`
throughout; only the per-atom `lambda[]` had `pack_restart`/
`unpack_restart`). So restart continuation never actually resumed the ramp
physics; it only carried over a per-atom value nothing read. Removing the
restart callbacks doesn't regress anything real — it removes dead
machinery that was masking an existing gap. `docs/components/fix-backmap.md`
previously claimed "seamless continuation," which was already inaccurate;
this change corrects that claim rather than just removing the
now-doubly-false version of it.

## Design decision: `phase` cleanup scope (LAMMPS side only)

The Python `backmap-prep` schema has a separate, legitimately-deferred
`simulation.two_phase` field, explicitly guarded at validation time
("Feature 'two_phase' backmapping is not yet implemented (planned for
Phase 2)") — this is intentional forward-looking schema design, not dead
code, and is **out of scope** for this change. What *is* in scope: the
`backmap-input-generator` spec's scenarios describing what the generator
would emit if `two_phase` were implemented (`fix backmap ... phase 1`,
`fix_modify bm phase 2`) targeted the now-removed C++ `phase` keyword —
those scenarios are corrected to describe the actual (rejected) behavior,
without touching the `two_phase` schema field, its guard, or its test.

## Scenario, not requirement: dead code found but left alone

`compute_weight3()`'s CG-case branch keyed on `phase` was the only
`phase`-dependent code path; removing it is in scope. The rest of the file
(`same_bead()` overloads, `clamp_lambda()`, `is_almost_zero()`) is
untouched — only `clamp_lambda()`'s comment is corrected (it attributed
negative values to "nonuniform init," which is no longer a real source;
negative `lambda_global` can still occur from a negative `lambda0`).
