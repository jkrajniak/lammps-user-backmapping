## Context

`fix backmap` drives a time-dependent CG→AT resolution change via a lambda
ramp. Every interaction style that needs lambda-dependent weighting
(`pair_backmap`, the six `bond/angle/dihedral_backmap_*` styles, and
`fix_backmap_pairs`) currently calls one shared helper:

```cpp
// backmap_lambda.h (current)
inline double compute_weight(double lambda_i, double lambda_j, bool is_cg) {
  double w = lambda_i * lambda_j;
  return is_cg ? (1.0 - w) : w;
}
```

`lambda_i`/`lambda_j` come from `fix_backmap`'s **per-atom** `lambda[]` array
(uniform across atoms today since `nonuniform no` is the default and no
in-flight examples use `nonuniform yes` for these interactions). The
intended physics (confirmed with the user) is:

1. **CG-CG** (intramolecular CG network): `w = 1 − λ` (global λ, not
   per-atom product).
2. **AT-AT intra-bead** (both atoms in the same CG bead): `w = 1` once
   `λ > 0` — full strength immediately, never scaled down. Today this case
   does not exist as a distinct code path: `pair_backmap.cpp` treats *all*
   AT-AT pairs (same-bead or not) as one category weighted by `λ_i × λ_j`.
3. **AT-AT inter-bead** (atoms in different CG beads): `w = λ` (global,
   linear).

`fix_backmap.h` already builds `atom2cg[]` (per-atom bead index) for its own
internal COM/force-distribution bookkeeping, but does not expose it to
styles. Bonded styles (`bond/angle/dihedral_backmap_*`) largely sidestep the
intra/inter-bead question already: intra-bead bonded terms are emitted by
`backmap-prep` as plain, unweighted LAMMPS styles (e.g. `bond_coeff 1
harmonic ...`) inside a `hybrid` style, while only genuinely cross-bead
("cross-CG") bonded terms get the `backmap/*` weighted style with an `at`/`cg`
tag. `pair_style backmap`, by contrast, classifies purely by **atom type**
(AT vs CG), so it has no way today to tell an intra-bead AT-AT pair from an
inter-bead one — this is the one gap with real physical consequences,
motivating the `atom2cg` same-bead check at the style level. We add the same
check to the six bonded `backmap/*` styles and `fix_backmap/pairs` too, both
for consistency (the plan the user reviewed and approved) and because
individual cross-bead angles/dihedrals can span atoms where not all pairs
within the same interaction are guaranteed to be inter-bead in every possible
future topology.

This repo has no C++ unit tests. `cpp-ci.yml` only builds LAMMPS+the package
and checks `clang-format --dry-run` against an *undeclared* fallback style
(no `.clang-format` file exists). Full-scale Tier B validation runs
(Dodecane, PE, rim135, melamine, PET) take minutes to hours each; a change to
core weighting math must be checked by a fast, deterministic test before any
of those runs.

## Goals / Non-Goals

**Goals:**

- Replace the per-particle `λ_i × λ_j` weighting formula with a 3-way
  formula driven only by a single global λ scalar and same-bead membership.
- Expose `lambda_global` and `atom2cg` from `fix_backmap` via `extract()`.
- Update all 8 call sites consistently; remove the old formula entirely (no
  dual code paths).
- Add fast, standalone (no LAMMPS build) unit tests for the new formula and
  the same-bead helper, runnable in CI and locally in seconds.
- Wire the new tests into `cpp-ci.yml` and `make test-cpp`; fix the
  currently-implicit `.clang-format` fallback.
- Regress-test the four already-validated systems (Dodecane, PE, rim135,
  melamine) and re-check the PET CG-angle-table anomaly under the corrected
  model.

**Non-Goals:**

- Two-phase mode (`fix_modify bm phase 1|2`) semantics beyond a straight
  formula substitution. The existing `bonded-backmap-styles` spec documents
  a Phase-1-forces-CG-to-full-strength override; this change preserves that
  override's *behavior* (same phase semantics) but re-expresses it in terms
  of `(same_bead, is_cg, lambda_global, phase)` instead of
  `(lambda_i, lambda_j, is_cg, phase)`. No example in this repo currently
  exercises two-phase mode with non-trivial per-atom lambda variance, so
  this is a mechanical carry-over, not a re-design.
- `nonuniform yes` (staggered per-molecule λ0) is unaffected in how it's
  *set*; it only stops being read by the weighting formula (which now reads
  `lambda_global` instead of per-atom `lambda[]`). `nonuniform` still varies
  the per-atom `lambda[]` array used for thermo/dump/restart output; whether
  that combination is physically meaningful is out of scope here — no
  validated example uses `nonuniform yes` today.
- A LAMMPS-native, full-build `unittest/force-styles`-style YAML integration
  test that drives real per-atom forces through `fix backmap` + a style.
  Flagged as a deferred follow-up (see Risks) — valuable, but it needs the
  full LAMMPS build either way, so it does not serve the "fast pre-simulation
  gate" goal and would slow this change down.
- Re-architecting the CG-angle-table anomaly (`table_a0.table`, angle type
  36) found during PET debugging. That is data/table content, not a
  weighting-formula bug; this change only re-verifies whether it persists
  once the formula fix lands.

## Decisions

### D1: Drop the per-particle product entirely, keep a single new helper

**Decision:** Remove `compute_weight(lambda_i, lambda_j, is_cg)` from
`backmap_lambda.h` outright and replace every call site with
`compute_weight3(same_bead, is_cg, lambda_global)`:

```cpp
inline double compute_weight3(bool same_bead, bool is_cg, double lambda_global) {
  if (is_cg) return 1.0 - lambda_global;
  if (same_bead) return (lambda_global > 0.0) ? 1.0 : 0.0;
  return lambda_global;
}
```

**Why over alternatives:**
- *Keep both functions, add `compute_weight3` alongside* — rejected per
  explicit user direction ("you can also drop the formula... just use
  global lambda value"). Keeping a dead per-particle formula around also
  invites accidental reuse.
- *Overload `compute_weight()` with a different signature* — rejected;
  same name with different semantics (global vs per-particle) is a latent
  bug magnet. A new name (`compute_weight3`) makes the 3-way behavior and
  the signature change explicit at every call site.

### D2: Same-bead determined from `atom2cg`, not from a topology tag

**Decision:** Add a small header-only helper,
`same_bead(const int *atom2cg, int i1, int i2, ..., iN)`, that returns true
iff all given local atom indices map to the same bead index in `atom2cg`.
Pair/bond styles pass 2 indices, angle styles 3, dihedral styles 4.

**Why over alternatives:**
- *Encode intra/inter-bead as a third `bond_coeff`/`pair_coeff` tag
  (`intra`/`inter`/`cg`) set at generation time* — rejected for pairs:
  `pair_style backmap` operates per atom-type-pair, not per specific atom
  pair, so it cannot statically know at `pair_coeff` time whether a given
  runtime i-j pair is same-bead or not (that depends on which molecule/bead
  the atoms belong to, only known at `compute()` time from `atom2cg`).
  Using the same runtime mechanism for the six bonded styles keeps one
  code path instead of two.
- *Compute same-bead from `atom->molecule` equality* — insufficient: a
  single molecule can contain many CG beads (see `fix-backmap-resolution`'s
  existing per-bead mapping), so two AT atoms in the same molecule can still
  be in different beads.

### D3: Expose `lambda_global` and `atom2cg` via `extract()`, keep per-atom `lambda[]` for other uses

**Decision:** Add a `double lambda_global` member to `FixBackmap`, updated
identically to the per-atom ramp (`lambda_global = min(1.0, lambda_global +
alpha)`), and extend `extract()`:

```cpp
void *FixBackmap::extract(const char *str, int &dim) {
  if (strcmp(str, "lambda") == 0) { dim = 1; return static_cast<void*>(lambda); }
  if (strcmp(str, "atom2cg") == 0) { dim = 1; return static_cast<void*>(atom2cg); }
  if (strcmp(str, "lambda_global") == 0) { dim = 0; return static_cast<void*>(&lambda_global); }
  return nullptr;
}
```

The per-atom `lambda[]` array is **not removed** — it remains the backing
store for thermo (`f_bm`), dump output, and restart, and continues to support
`nonuniform yes` staggering for those purposes. It is simply no longer read
by the weighting formula in styles.

**Why over alternatives:**
- *Derive `lambda_global` on the fly from `lambda[0]` in each style* —
  fragile (assumes atom 0 is local and representative) and duplicates the
  ramp/clamp logic in 8 places instead of 1.
- *Remove per-atom `lambda[]` entirely* — breaks thermo/dump/restart
  requirements already in the `fix-backmap-resolution` spec (`f_bm`,
  per-atom dump output, restart of λ state); out of scope.

### D4: Fast unit tests live outside the LAMMPS build, using CMake `FetchContent` for GoogleTest

**Decision:** Create `tests/unit/` with its own minimal `CMakeLists.txt`
that fetches GoogleTest via `FetchContent` and compiles only
`test_backmap_lambda.cpp` against `backmap_lambda.h` — no LAMMPS headers, no
LAMMPS library link. This mirrors how LAMMPS's own `unittest/` uses
GoogleTest, but avoids that suite's requirement of linking the entire LAMMPS
library (its `unittest/force-styles/` tests instantiate a full `LAMMPS`
object and run real input scripts via YAML-described configs).

**Why over alternatives:**
- *Add these tests inside LAMMPS's own `unittest/force-styles/` YAML
  harness* — correct for true end-to-end style validation (real per-atom
  forces from a running `fix backmap` + style), but requires the full
  LAMMPS+MPI build (minutes) and cannot run before/without a LAMMPS clone.
  Deferred as a follow-up (see Risks), not dropped.
- *Header-only, assert-based smoke test (no framework)* — works, but
  GoogleTest is already the LAMMPS-ecosystem convention (`unittest/`), and
  the pre-commit/CI setup for this repo already assumes C++ tooling
  parity with upstream LAMMPS practice; `FetchContent` keeps the dependency
  local to `tests/unit/` without touching `src/` or requiring `vendor/`.
- *Reuse the existing LAMMPS-tree build in `cpp-ci.yml`'s `build` job for
  these tests* — rejected: that job clones+builds all of LAMMPS
  (~5-10 min), which defeats the "fast gate before simulating" goal. The new
  `unit-tests` job needs only a C++17 compiler + CMake, no LAMMPS clone.

### D5: `LAMBDA_AT_ONSET` guard switches from per-particle to global

**Decision:** In `pair_backmap.cpp`, change:
```cpp
if (!is_cg && (li < LAMBDA_AT_ONSET || lj < LAMBDA_AT_ONSET)) continue;
```
to:
```cpp
if (!is_cg && lambda_global < LAMBDA_AT_ONSET) continue;
```
This guard is a numerical-safety measure (deferring AT-AT pair onset until
the ramp has progressed a little, to avoid LJ singularities on
badly-overlapping initial AT placements) and is orthogonal to the
same-bead/inter-bead distinction; it's simply re-expressed in terms of the
global scalar since per-particle λ is no longer read for weighting.

### D6: `.clang-format` — adopt explicit Google base matching current source style

**Decision:** Add a checked-in `.clang-format` at the repo root with
`BasedOnStyle: Google` plus the 2-space indentation and brace style already
visible in the existing `src/*.cpp`/`.h` files (matching what the LAMMPS
project itself documents in prose, since LAMMPS ships no machine-readable
`.clang-format`). This makes the existing `--fallback-style=Google` in
`.pre-commit-config.yaml`/`cpp-ci.yml` explicit rather than implicit.

**Why over alternatives:**
- *Leave it as an undeclared fallback* — works today only because nobody
  has added a stray `.clang-format` upstream in a parent directory that
  clang-format would pick up instead; fragile and surprising for
  contributors who add one locally.
- *Adopt LLVM or Mozilla base style* — would reformat the entire `src/`
  tree in one commit, a large, unrelated diff; Google style matches current
  formatting already.

## Risks / Trade-offs

- **[Risk] Removing `compute_weight()` breaks any out-of-tree code still
  calling it** → Mitigation: it's a private header-only helper inside this
  package's `src/`, not part of any published LAMMPS command syntax; no
  external callers exist. Grep confirms exactly 8 call sites, all inside
  this repo, all updated in this change.
- **[Risk] `same_bead()` requires `atom2cg` to be valid for ghost atoms
  too** (pair/bonded computes iterate over local+ghost neighbor indices) →
  Mitigation: `atom2cg` must be forward-communicated the same way `lambda[]`
  already is (`comm_forward`); verify this is covered by the existing
  MPI-correctness forward-comm path added in the `mpi-correctness`
  worktree, and add an explicit check/test if not.
- **[Risk] Fast unit tests only cover the pure formula, not real per-atom
  force correctness inside a running simulation** → Mitigation: explicitly
  scoped as tier 1 of a two-tier plan; tier 2 (LAMMPS-native integration
  test) is a named, deferred follow-up, and the existing Tier B full-system
  regression runs (Dodecane/PE/rim135/melamine/PET) remain the end-to-end
  check for this change before it's considered validated.
- **[Risk] Regression risk for the four already-validated systems** — all
  are single-bead-per-residue or otherwise mostly-uniform-λ systems, so
  `same_bead` may rarely or never evaluate false in bonded styles, meaning
  the formula change could be close to a no-op for them, or could reveal
  something previously masked by the per-particle product happening to
  equal the global value under uniform ramps → Mitigation: task list
  requires rerunning Tier B (minimum) for all four before merging, plus a
  Tier C RDF spot-check for Dodecane.

## Migration Plan

1. Implement in a dedicated worktree/branch (`fix/lambda-weighting`),
   isolated from the in-flight `mpi-correctness` and `pete-example`
   worktrees.
2. Land the fast unit tests first (or in the same PR) so CI gates the
   formula change immediately — no simulation run should be needed to catch
   an arithmetic mistake in `compute_weight3()`.
3. Build against the local LAMMPS tree, run the four validated systems'
   Tier B (+ Dodecane Tier C spot-check), then PET Tier B.
4. No rollback complexity beyond normal `git revert` — this is a pure
   internal C++ change with no data-file or on-disk format migration.

## Open Questions

- Should the deferred LAMMPS-native integration test (tier 2) be scoped as
  its own follow-up change now, or left informal until tier 1 lands and the
  formula is confirmed correct against real systems? Leaning toward: file it
  informally in `research/` as a follow-up once this change is archived,
  not as a blocking task here.
