## MODIFIED Requirements

### Requirement: Dihedral style `backmap/ryckaert`

The `dihedral_style backmap/ryckaert` SHALL compute Ryckaert-Bellemans dihedral forces with lambda weighting:

```
E = w × Σ(n=0..5) Cn × cos^n(φ)
F_φ = w × Σ(n=1..5) n × Cn × cos^(n-1)(φ) × sin(φ)
```

where φ is the dihedral angle using the LAMMPS polymer convention (trans = 180°).

The weight `w` is computed from the lambda values of the first and last atoms of the dihedral (atoms i and l in the i-j-k-l quadruplet), using `backmap_lambda.h` utilities.

Command syntax:
```
dihedral_style backmap/ryckaert
dihedral_coeff N at/cg C0 C1 C2 C3 C4 C5
```

Where:
- `N` = dihedral type
- `at` or `cg` = weighting direction
- `C0..C5` = Ryckaert-Bellemans coefficients in kcal/mol (LAMMPS convention)

Used within `dihedral_style hybrid`:
```
dihedral_style hybrid backmap/ryckaert backmap/table
dihedral_coeff 1 backmap/ryckaert at 1.53 0.776 -1.19 -3.22 0.0 0.0
```

Force, energy, and virial contributions SHALL all be scaled by the same weight factor `w`.

The style SHALL locate `fix backmap` during initialization via `find_fix_backmap()` from `backmap_lambda.h`. If no fix is found, the style SHALL abort with: "dihedral_style backmap/ryckaert requires fix backmap".

#### Scenario: RB dihedral with AT weighting at lambda=1.0
- **WHEN** `dihedral_coeff 1 backmap/ryckaert at 1.53 0.776 -1.19 -3.22 0.0 0.0` and λ_i = λ_l = 1.0
- **THEN** the dihedral force SHALL be at full strength (weight = 1.0 × 1.0 = 1.0)

#### Scenario: RB dihedral with AT weighting at lambda=0.5
- **WHEN** `dihedral_coeff 1 backmap/ryckaert at C0..C5` and λ_i = λ_l = 0.5
- **THEN** the dihedral energy SHALL be `0.25 × Σ Cn × cos^n(φ)` (weight = 0.5 × 0.5 = 0.25)

#### Scenario: RB dihedral with CG weighting at lambda=0.5
- **WHEN** `dihedral_coeff 2 backmap/ryckaert cg C0..C5` and λ_i = λ_l = 0.5
- **THEN** the dihedral energy SHALL be `0.75 × Σ Cn × cos^n(φ)` (weight = 1 − 0.25 = 0.75)

#### Scenario: RB dihedral at lambda=0 (pure CG)
- **WHEN** λ = 0 for all atoms
- **THEN** AT dihedrals (`at`) SHALL have zero force, CG dihedrals (`cg`) SHALL have full-strength force

#### Scenario: Skip computation when weight is negligible
- **WHEN** the computed weight for a dihedral is below 1e-10
- **THEN** the style SHALL skip the force computation entirely for that dihedral

#### Scenario: Restart file support
- **WHEN** a simulation using `backmap/ryckaert` is written to a restart file
- **THEN** the dihedral style parameters (C0..C5, at/cg flag) SHALL be preserved in the restart and correctly restored on read

#### Scenario: Fix not defined
- **WHEN** a `backmap/ryckaert` dihedral is defined but no `fix backmap` exists
- **THEN** the style SHALL abort with an error message: "dihedral_style backmap/ryckaert requires fix backmap"

### Requirement: Dihedral style `backmap/table`

The `dihedral_style backmap/table` SHALL compute tabulated dihedral forces with lambda weighting:

```
E = w × E_table(φ)
F_φ = w × F_table(φ)
```

Command syntax:
```
dihedral_style backmap/table linear N
dihedral_coeff M at/cg filename keyword
```

Where:
- `M` = dihedral type
- `at` or `cg` = weighting direction
- `filename` = table file path
- `keyword` = table section name
- `N` = number of interpolation points

The table format SHALL follow LAMMPS `dihedral_style table` conventions:
- Column 1: index (integer)
- Column 2: dihedral angle in degrees (-180 to 180)
- Column 3: energy in kcal/mol
- Column 4: force (-dE/dφ) in kcal/(mol·radian)

The style SHALL support both `linear` and `spline` interpolation modes.

Force, energy, and virial contributions SHALL all be scaled by the same weight factor `w`.

The style SHALL locate `fix backmap` during initialization via `find_fix_backmap()` from `backmap_lambda.h`.

#### Scenario: CG dihedral with tabulated potential at lambda=0
- **WHEN** `dihedral_coeff 2 backmap/table cg table_d1.table ENTRY` and λ = 0
- **THEN** the dihedral force SHALL be at full strength (weight = 1 − 0 = 1.0)

#### Scenario: CG dihedral with tabulated potential at lambda=0.7
- **WHEN** `dihedral_coeff 2 backmap/table cg table_d1.table ENTRY` and λ_i = λ_l = 0.7
- **THEN** the dihedral force SHALL be `0.51 × F_table(φ)` (weight = 1 − 0.49 = 0.51)

#### Scenario: AT dihedral with tabulated potential at lambda=1
- **WHEN** `dihedral_coeff 3 backmap/table at table_d2.table ENTRY` and λ = 1
- **THEN** the dihedral force SHALL be at full strength (weight = 1.0)

#### Scenario: Skip computation when weight is negligible
- **WHEN** the computed weight for a tabulated dihedral is below 1e-10
- **THEN** the style SHALL skip the force computation entirely

#### Scenario: Restart file support
- **WHEN** a simulation using `backmap/table` dihedrals is written to a restart file
- **THEN** the table filename, keyword, interpolation mode, and at/cg flag SHALL be preserved

#### Scenario: Fix not defined
- **WHEN** a `backmap/table` dihedral is defined but no `fix backmap` exists
- **THEN** the style SHALL abort with: "dihedral_style backmap/table requires fix backmap"

#### Scenario: Invalid table file
- **WHEN** the specified table file does not exist or contains invalid data
- **THEN** the style SHALL abort with a descriptive error naming the file and the issue

### Requirement: Lambda-weighted bonded force computation

The `backmap/*` bonded styles SHALL compute bonded forces (bonds, angles, dihedrals) with lambda weighting applied at force-computation time. The weighting scheme is identical to the non-bonded pair weighting in `pair_style backmap` and matches ESPResSo++'s `FixedPairListAdressInteractionTemplate`:

- **AT weight** (keyword `at`): `w = λ_i × λ_j` — for the two atoms forming the bond, or the first and last atoms of an angle/dihedral
- **CG weight** (keyword `cg`): `w = 1 − λ_i × λ_j`

For uniform backmapping (all particles share the same λ), this reduces to `w_AT = λ²` and `w_CG = 1 − λ²`, ensuring `w_AT + w_CG = 1` always.

Each style reads per-atom λ from `fix backmap`'s per-atom array via `fix->extract()`. The fix ID is auto-detected (single `fix backmap` in the simulation) or specified via a style argument.

Force, energy, and virial contributions SHALL all be scaled by the same weight factor.

**Phase-aware weighting**: During Phase 1 of two-phase backmapping (see fix-backmap-resolution spec), CG-weighted styles (`cg` keyword) SHALL apply full strength (weight = 1.0) regardless of lambda. During Phase 2, CG styles SHALL use the standard `1 − λ_i × λ_j` weighting. AT-weighted styles SHALL always use `λ_i × λ_j` weighting in both phases.

The phase is read from `fix backmap` via `fix->extract("phase")`. The `compute_backmap_weight()` helper in `backmap_lambda.h` SHALL accept a phase parameter and return the correct weight.

#### Scenario: Cross-CG AT bond at lambda=0.5
- **WHEN** a `backmap/harmonic at` bond connects atoms i,j with λ_i = λ_j = 0.5
- **THEN** the bond force SHALL be 0.25 × F_harmonic (where F_harmonic is the unweighted harmonic force)

#### Scenario: Cross-CG CG bond at lambda=0.5
- **WHEN** a `backmap/table cg` bond connects CG atoms i,j with λ_i = λ_j = 0.5
- **THEN** the bond force SHALL be 0.75 × F_table (weight = 1 − 0.25 = 0.75)

#### Scenario: Bond force at lambda=0 (pure CG)
- **WHEN** λ = 0 for all atoms
- **THEN** AT cross bonds (`at`) SHALL have zero force, CG cross bonds (`cg`) SHALL have full-strength force

#### Scenario: Bond force at lambda=1 (pure AT)
- **WHEN** λ = 1 for all atoms
- **THEN** AT cross bonds (`at`) SHALL have full-strength force, CG cross bonds (`cg`) SHALL have zero force

#### Scenario: Skip computation when weight is negligible
- **WHEN** the computed weight is below a threshold (~1e-10)
- **THEN** the style SHALL skip the force computation entirely for efficiency (matching ESPResSo++'s `is_almost_zero()` check)

#### Scenario: Phase 1 CG forces at full strength
- **WHEN** `fix backmap` is in Phase 1 and λ = 0.5
- **THEN** CG cross bonds (`cg`) SHALL have full-strength force (weight = 1.0, ignoring lambda)
- **AND** AT cross bonds (`at`) SHALL have weight 0.25 (λ² = 0.25, phase does not affect AT weighting)

#### Scenario: Phase 2 CG forces use standard weighting
- **WHEN** `fix backmap` is in Phase 2 and λ = 0.5
- **THEN** CG cross bonds (`cg`) SHALL have weight 0.75 (1 − λ² = 0.75)
- **AND** AT cross bonds (`at`) SHALL have weight 0.25

### Requirement: Shared lambda-access helper

All `backmap/*` styles (pair, bond, angle, dihedral) SHALL share a common helper class or set of utility functions for:
- Locating `fix backmap` by scanning the fix list
- Reading per-atom λ values from the fix's per-atom array
- Computing the weight `w` given two lambda values, the `at`/`cg` flag, and the current phase
- Applying the `is_almost_zero()` check to skip negligible-weight computations

This shared code SHALL be in a single header (`backmap_lambda.h`) to avoid duplication across styles.

The `compute_backmap_weight()` function SHALL accept a phase parameter:
- Phase 1 + `is_cg=true` → return 1.0 (full strength)
- Phase 1 + `is_cg=false` → return `λ_i × λ_j`
- Phase 2 + `is_cg=true` → return `1 − λ_i × λ_j`
- Phase 2 + `is_cg=false` → return `λ_i × λ_j`

#### Scenario: Weight computation helper
- **WHEN** `compute_backmap_weight(lambda_i=0.7, lambda_j=0.7, is_cg=false, phase=2)` is called
- **THEN** it SHALL return 0.49 (0.7 × 0.7)

#### Scenario: Weight computation helper CG
- **WHEN** `compute_backmap_weight(lambda_i=0.7, lambda_j=0.7, is_cg=true, phase=2)` is called
- **THEN** it SHALL return 0.51 (1 − 0.49)

#### Scenario: Phase 1 CG weight override
- **WHEN** `compute_backmap_weight(lambda_i=0.7, lambda_j=0.7, is_cg=true, phase=1)` is called
- **THEN** it SHALL return 1.0 (full strength in Phase 1)

#### Scenario: Phase 1 AT weight unchanged
- **WHEN** `compute_backmap_weight(lambda_i=0.7, lambda_j=0.7, is_cg=false, phase=1)` is called
- **THEN** it SHALL return 0.49 (AT weighting is the same regardless of phase)
