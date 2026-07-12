## Why

The lambda-weighting formula used by every `backmap/*` style (`pair_backmap`,
`bond/angle/dihedral_backmap_{harmonic,table,ryckaert}`, `fix_backmap/pairs`)
computes weights from a **per-particle product** `λ_i × λ_j` of the two end
atoms of an interaction. This does not match the intended physics: the
correct model uses a single **global** λ (one scalar per simulation, not
per-atom), and it must distinguish three cases — CG-CG network, AT-AT
**intra-bead** (same CG bead, must be full strength once λ>0, never scaled),
and AT-AT **inter-bead** (different CG beads, scales linearly with λ). The
current code has no runtime concept of "same bead", so non-bonded AT-AT pair
interactions inside a single CG bead are incorrectly scaled down during the
ramp instead of acting at full strength. This was discovered while debugging
the PET/Dacron network Tier B instability and confirmed by the user, who also
asked to simplify the fix to always use the global λ scalar (no per-particle
product at all).

Separately, this package has zero C++ unit tests today. The only automated
C++ checks are a `clang-format --dry-run` pass and a full LAMMPS+MPI build
that greps `lmp -h` for registered style names (`.github/workflows/cpp-ci.yml`).
A change to core weighting math should not go into a full-scale simulation
(minutes-to-hours per run, across 5 validated systems) without a fast,
deterministic test of the formula itself first.

## What Changes

- **BREAKING (internal API only, no LAMMPS input-script syntax change)**:
  Remove `BackmapLambda::compute_weight(lambda_i, lambda_j, is_cg)` (the
  `λ_i × λ_j` product) from `backmap_lambda.h`. Replace it with
  `compute_weight3(same_bead, is_cg, lambda_global)` used by all call sites.
  No user-facing fix/pair/bond/angle/dihedral command syntax changes.
- Add a `lambda_global` scalar and expose it, together with the existing
  per-atom `atom2cg` bead-membership map, via `FixBackmap::extract()` so
  styles can read them without maintaining separate state.
- Update all 8 call sites (`pair_backmap.cpp`, `bond_backmap_harmonic.cpp`,
  `bond_backmap_table.cpp`, `angle_backmap_harmonic.cpp`,
  `angle_backmap_table.cpp`, `dihedral_backmap_ryckaert.cpp`,
  `dihedral_backmap_table.cpp`, `fix_backmap_pairs.cpp`) to determine
  `same_bead` from `atom2cg` and call `compute_weight3()` instead of the
  removed per-particle formula.
- Update `pair_backmap.cpp`'s `LAMBDA_AT_ONSET` safety guard from a
  per-particle check (`λ_i < 0.1 || λ_j < 0.1`) to a global check
  (`lambda_global < 0.1`).
- Add a new fast, dependency-free (no LAMMPS build required) GoogleTest
  suite validating the weighting formula and bead-membership helper in
  isolation, runnable in seconds.
- Wire the new test suite into `cpp-ci.yml` as an independent, fast CI job,
  and add a `make test-cpp` target so it runs identically on a laptop.
- Add the missing `.clang-format` file (the existing pre-commit/CI hooks
  reference `--style=file` but no such file exists in the repo today, so
  they silently fall back to Google style).

## Capabilities

### New Capabilities

- `lambda-weighting-unit-tests`: fast, standalone GoogleTest suite for the
  lambda-weighting formula and bead-membership helper, plus the CI job and
  local `make` target that run it.

### Modified Capabilities

- `fix-backmap-resolution`: `FixBackmap::extract()` gains `"atom2cg"` and
  `"lambda_global"` accessors alongside the existing `"lambda"`.
- `pair-backmap`: the AT/CG weighting formula changes from the per-particle
  product `λ_i × λ_j` / `1 − λ_i × λ_j` to a 3-way, global-λ formula that
  also distinguishes AT intra-bead from AT inter-bead pairs.
- `bonded-backmap-styles`: the shared weighting helper
  (`compute_backmap_weight` → `compute_weight3`) and all bond/angle/dihedral
  `backmap/*` styles change from the per-particle product to the same 3-way,
  global-λ formula. `fix backmap/pairs` gets the same treatment.

## Impact

- **C++ source** (`lammps-user-backmapping/src/`): `backmap_lambda.h`,
  `fix_backmap.h`/`.cpp`, `pair_backmap.cpp`, `bond_backmap_harmonic.cpp`,
  `bond_backmap_table.cpp`, `angle_backmap_harmonic.cpp`,
  `angle_backmap_table.cpp`, `dihedral_backmap_ryckaert.cpp`,
  `dihedral_backmap_table.cpp`, `fix_backmap_pairs.cpp`.
- **New C++ tests**: `tests/unit/test_backmap_lambda.cpp` +
  `tests/unit/CMakeLists.txt` (new directory, no LAMMPS dependency).
- **CI/tooling**: `.github/workflows/cpp-ci.yml` (new job), `Makefile` (new
  `test-cpp` target), `CONTRIBUTING.md` (workflow checklist), new
  `.clang-format` file.
- **Validated examples** (re-run, not modified): Dodecane, PE, rim135,
  melamine (regression check), PET/Dacron (re-check of the CG angle-table
  anomaly under the corrected model).
- **Docs**: `docs/` pages describing the weighting model; `research/`
  notebook entry correcting the earlier wrong diagnosis; `CHANGELOG.md`.
- No breaking change to LAMMPS input-script syntax — `fix`/`pair_coeff`/
  `bond_coeff`/`angle_coeff`/`dihedral_coeff` commands are unchanged. The
  break is internal (C++ helper signature only).
