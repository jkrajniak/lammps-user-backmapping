## 1. Isolate the work

- [ ] 1.1 Create a new git worktree `fix/lambda-weighting` off
  `lammps-user-backmapping` (separate from `mpi-correctness` and
  `pete-example`)

## 2. `fix_backmap` — expose lambda_global and atom2cg

- [ ] 2.1 Add `double lambda_global` member to `FixBackmap` (`fix_backmap.h`)
- [ ] 2.2 Update `lambda_global` in the same place/way the per-atom ramp
  currently updates `lambda[i] += alpha` (clamped to `[0, 1]`)
- [ ] 2.3 Extend `FixBackmap::extract()` to return `"atom2cg"` (dim=1,
  `int*`) and `"lambda_global"` (dim=0, `double*`), alongside the existing
  `"lambda"`
- [ ] 2.4 Verify `atom2cg` is forward-communicated to ghost atoms (check the
  existing `comm_forward`/`pack_forward_comm`/`unpack_forward_comm` path
  added for MPI correctness); add ghost coverage if missing
- [ ] 2.5 Add/confirm restart and `grow_arrays`/`copy_arrays` handling for
  `lambda_global` (scalar, trivial) — no per-atom growth needed beyond what
  `atom2cg` already requires

## 3. `backmap_lambda.h` — replace the weighting formula

- [ ] 3.1 Remove `compute_weight(double lambda_i, double lambda_j, bool is_cg)`
  and its doc comment entirely
- [ ] 3.2 Add `compute_weight3(bool same_bead, bool is_cg, double lambda_global)`
  implementing: `is_cg → 1 − λ_global`; `same_bead → (λ_global > 0) ? 1 : 0`;
  else `→ λ_global`
- [ ] 3.3 Add `same_bead(const int *atom2cg, int i1, int i2)` (pair/bond),
  `same_bead(const int *atom2cg, int i1, int i2, int i3)` (angle), and
  `same_bead(const int *atom2cg, int i1, int i2, int i3, int i4)` (dihedral)
  overloads (or a single variadic/array-based helper)
- [ ] 3.4 Update the file's top-of-file usage-pattern doc comment to show
  `compute_weight3()` instead of the removed `compute_weight()`

## 4. Update all 8 call sites

- [ ] 4.1 `pair_backmap.cpp`: extract `atom2cg`/`lambda_global` instead of
  per-atom `lambda[]` for weighting; replace `compute_weight(li, lj, is_cg)`
  with `compute_weight3(same_bead(atom2cg, i, j), is_cg, lambda_global)`;
  change `LAMBDA_AT_ONSET` guard to `lambda_global < LAMBDA_AT_ONSET`
- [ ] 4.2 `bond_backmap_harmonic.cpp`: same pattern for `compute()` and
  `single()`
- [ ] 4.3 `bond_backmap_table.cpp`: same pattern for `compute()` and
  `single()`
- [ ] 4.4 `angle_backmap_harmonic.cpp`: same pattern, `same_bead` over
  `i1,i2,i3`
- [ ] 4.5 `angle_backmap_table.cpp`: same pattern, `same_bead` over
  `i1,i2,i3`, both `compute()` and the energy-only path
- [ ] 4.6 `dihedral_backmap_ryckaert.cpp`: same pattern, `same_bead` over
  `i1,i2,i3,i4`
- [ ] 4.7 `dihedral_backmap_table.cpp`: same pattern, `same_bead` over
  `i1,i2,i3,i4`
- [ ] 4.8 `fix_backmap_pairs.cpp`: same pattern for its pairwise loop
- [ ] 4.9 Grep the whole `src/` tree for any remaining `compute_weight(` /
  per-atom `lam[i` weighting reads to confirm no call site was missed

## 5. Fast unit tests (must pass before any full-scale run)

- [ ] 5.1 Create `tests/unit/CMakeLists.txt`: standalone project, C++17,
  `FetchContent` for GoogleTest, one test binary linked only against
  `backmap_lambda.h` (no LAMMPS headers/library)
- [ ] 5.2 Create `tests/unit/test_backmap_lambda.cpp` covering:
  - CG weight: `w(is_cg=true, λ=0)=1`, `w(is_cg=true, λ=1)=0`, linear at
    λ=0.25/0.5/0.75, independent of `same_bead`
  - AT intra-bead: `w(same_bead=true, λ=0)=0`, `w(same_bead=true, λ=ε)=1`
    for small ε>0 and for λ=1 (step function, not scaled)
  - AT inter-bead: `w(same_bead=false, λ)=λ` exactly at λ=0, 0.3, 0.5, 1
  - `same_bead()` pair/angle/dihedral overloads: true when all indices map
    to the same bead, false for any mismatch, including the "not all same,
    not all different" mixed case for angle/dihedral
  - Clamping: `lambda_global` outside `[0,1]` is clamped before use
    (negative → treated as 0, >1 → treated as 1)
- [ ] 5.3 Build and run locally: `cmake -S tests/unit -B tests/unit/build &&
  cmake --build tests/unit/build && ctest --test-dir tests/unit/build
  --output-on-failure`; confirm all pass before proceeding to task 7

## 6. CI, local target, and pre-commit

- [ ] 6.1 Add a `unit-tests` job to `.github/workflows/cpp-ci.yml` that
  configures/builds/runs `tests/unit/` with just a C++ compiler + CMake (no
  LAMMPS clone), independent of the existing `build` job
- [ ] 6.2 Add a `test-cpp` target to the top-level `Makefile` that runs the
  same configure/build/ctest sequence as the new CI job
- [ ] 6.3 Add `.clang-format` at the repo root (`BasedOnStyle: Google` +
  existing 2-space indentation) so the pre-commit/CI `--style=file` hooks
  stop relying on an undeclared fallback
- [ ] 6.4 Update `CONTRIBUTING.md`'s "Development Workflow" checklist to
  mention `make test-cpp` alongside the existing clang-format step

## 7. Build and regression-test against real systems

- [ ] 7.1 Copy `src/*.cpp src/*.h` into
  `/Users/jakubkrajniak/Work/Science/lammps/src/BACKMAP/` and rebuild
  (`cmake --build build -j$(nproc)`)
- [ ] 7.2 Rerun Tier B (minimum) for Dodecane, PE, rim135, melamine;
  confirm no regression (no crashes, no unexpected energy/temperature
  drift vs. the last known-good run)
- [ ] 7.3 Re-check Tier C RDF for Dodecane (fastest system) to confirm no
  numerical drift from the formula change
- [ ] 7.4 Rerun PET Tier B under the corrected model; record whether the
  CG angle-table anomaly (type 36 / `table_a0.table`) persists, changes
  magnitude, or resolves (expected: persists, since that term is
  `is_cg`-tagged and evaluates to the same `1 − λ_global` at λ=0 as before)

## 8. Documentation

- [ ] 8.1 Write a new dated research notebook entry
  (`research/notebook/<date>_lambda-weighting-fix.md`) documenting the
  before/after formula, the unit-test results, and the regression results,
  superseding the incorrect conclusions in
  `research/notebook/2026-07-12_pet-angle-geometry-blowup.md`
- [ ] 8.2 Update `research/STATUS.md` for the PET row and any other
  affected component status
- [ ] 8.3 Update `CHANGELOG.md` under `[Unreleased]` → `Fixed` (weighting
  formula) and `Added` (unit tests / CI job)
- [ ] 8.4 Update any `docs/` pages describing the lambda-weighting model to
  reflect the 3-way global-λ formula
