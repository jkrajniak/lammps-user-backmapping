## ADDED Requirements

### Requirement: Standalone lambda-weighting unit test suite

The package SHALL include a fast, standalone GoogleTest suite at
`tests/unit/test_backmap_lambda.cpp` that validates
`BackmapLambda::compute_weight3()` and the `same_bead()` helper family from
`backmap_lambda.h` in isolation, without requiring a LAMMPS build or LAMMPS
headers. The suite SHALL build and run in seconds via its own
`tests/unit/CMakeLists.txt`, which fetches GoogleTest via CMake
`FetchContent` and links only against `src/backmap_lambda.h`.

This suite exists to validate the lambda-weighting physics before any
full-scale simulation run (which can take minutes to hours per system),
catching arithmetic or logic errors in the weighting formula immediately.

#### Scenario: Suite builds without a LAMMPS tree
- **WHEN** `cmake -S tests/unit -B tests/unit/build` is run in a checkout with no LAMMPS source tree present
- **THEN** configuration SHALL succeed, requiring only a C++17 compiler and CMake

#### Scenario: Suite runs in seconds
- **WHEN** `ctest --test-dir tests/unit/build` is run after building
- **THEN** all tests SHALL complete in well under 30 seconds on a typical developer machine

#### Scenario: CG weight boundary and monotonicity
- **WHEN** `compute_weight3(same_bead=false, is_cg=true, lambda_global)` is evaluated at `lambda_global` = 0, 0.25, 0.5, 0.75, 1
- **THEN** it SHALL return 1, 0.75, 0.5, 0.25, 0 respectively (linear, independent of `same_bead`)

#### Scenario: AT intra-bead step function
- **WHEN** `compute_weight3(same_bead=true, is_cg=false, lambda_global)` is evaluated at `lambda_global` = 0, 1e-6, 0.5, 1
- **THEN** it SHALL return 0, 1, 1, 1 respectively (a step function, not a scaled value, for any `lambda_global > 0`)

#### Scenario: AT inter-bead linear ramp
- **WHEN** `compute_weight3(same_bead=false, is_cg=false, lambda_global)` is evaluated at `lambda_global` = 0, 0.3, 0.5, 1
- **THEN** it SHALL return exactly 0, 0.3, 0.5, 1 respectively

#### Scenario: same_bead pair helper
- **WHEN** `same_bead(atom2cg, i, j)` is evaluated with `atom2cg[i] == atom2cg[j]` and with `atom2cg[i] != atom2cg[j]`
- **THEN** it SHALL return true in the first case and false in the second

#### Scenario: same_bead angle helper — all same
- **WHEN** `same_bead(atom2cg, i1, i2, i3)` is evaluated with all three atoms mapping to the same bead index
- **THEN** it SHALL return true

#### Scenario: same_bead angle helper — mixed membership
- **WHEN** `same_bead(atom2cg, i1, i2, i3)` is evaluated with `i1, i2` in bead A and `i3` in bead B (not all same, not all different)
- **THEN** it SHALL return false

#### Scenario: same_bead dihedral helper
- **WHEN** `same_bead(atom2cg, i1, i2, i3, i4)` is evaluated with all four atoms in the same bead, and separately with `i4` in a different bead
- **THEN** it SHALL return true in the first case and false in the second

#### Scenario: lambda_global clamping
- **WHEN** `compute_weight3()` is called with a raw `lambda_global` value of -0.1 or 1.2
- **THEN** the value SHALL be clamped to 0 or 1 respectively before use in the weight formula

#### Scenario: No per-particle formula remains
- **WHEN** the test suite is compiled against the current `backmap_lambda.h`
- **THEN** it SHALL NOT reference any function taking two separate per-particle lambda values (`lambda_i`, `lambda_j`) as the weighting inputs — `compute_weight3()` takes only `(same_bead, is_cg, lambda_global[, phase])`

### Requirement: CI integration for the unit test suite

The package's CI pipeline (`.github/workflows/cpp-ci.yml`) SHALL include a
`unit-tests` job that configures, builds, and runs `tests/unit/` using only
a C++17 compiler and CMake — it SHALL NOT clone or build the full LAMMPS
source tree. This job SHALL run independently of (and faster than) the
existing `build` job that compiles LAMMPS with the BACKMAP package.

#### Scenario: Unit-tests job does not require LAMMPS
- **WHEN** the `unit-tests` CI job runs
- **THEN** it SHALL NOT include a step that clones or builds the LAMMPS source tree

#### Scenario: Unit-tests job fails on a broken formula
- **WHEN** `compute_weight3()` is modified to return an incorrect value for any of the covered scenarios
- **THEN** the `unit-tests` CI job SHALL fail before the `build` job's result is relevant

### Requirement: Local unit test target

The top-level `Makefile` SHALL provide a `test-cpp` target that runs the
identical configure/build/test sequence as the CI `unit-tests` job, so
`make test-cpp` on a developer machine matches CI exactly.

#### Scenario: Local target matches CI
- **WHEN** a developer runs `make test-cpp` locally
- **THEN** it SHALL configure, build, and run `tests/unit/` via CMake/CTest and report the same pass/fail result as the CI `unit-tests` job for the same source tree

### Requirement: Explicit C/C++ formatting configuration

The repository SHALL include a checked-in `.clang-format` file at the repo
root (`BasedOnStyle: Google`, matching the existing 2-space indentation used
throughout `src/`), so the existing `clang-format --style=file
--fallback-style=Google` invocations in `.pre-commit-config.yaml` and
`.github/workflows/cpp-ci.yml` resolve deterministically instead of relying
on an undeclared fallback.

#### Scenario: clang-format resolves a style without falling back
- **WHEN** `clang-format --style=file` is run from the repo root
- **THEN** it SHALL find and use the checked-in `.clang-format` file rather than the `--fallback-style` default

#### Scenario: Existing source already conforms
- **WHEN** `clang-format --style=file --dry-run --Werror` is run against the existing `src/*.cpp src/*.h` files after adding `.clang-format`
- **THEN** it SHALL report no formatting violations (the new file codifies the style already in use, it does not reformat existing code)
