## MODIFIED Requirements

### Requirement: Lambda-weighted bonded force computation

The `backmap/*` bonded styles SHALL compute bonded forces (bonds, angles, dihedrals) with lambda weighting applied at force-computation time, using a single **global** λ scalar (`lambda_global`, from `fix backmap`'s `extract("lambda_global")`) and per-atom CG-bead membership (`atom2cg`, from `extract("atom2cg")`) — not a per-particle product of two atoms' λ values. There are three weighting cases, matching the non-bonded weighting in `pair_style backmap`:

- **CG-tagged terms** (keyword `cg`): `w = 1 − λ_global`
- **AT-tagged, intra-bead terms** (keyword `at`, all atoms of the interaction map to the same CG bead in `atom2cg`): `w = 1` once `λ_global > 0`, else `w = 0`. Never scaled once the ramp has started.
- **AT-tagged, inter-bead terms** (keyword `at`, atoms of the interaction map to different CG beads): `w = λ_global`

Each style reads `atom2cg` and `lambda_global` from `fix backmap` via `fix->extract()`. The fix ID is auto-detected (single `fix backmap` in the simulation) or specified via a style argument.

Force, energy, and virial contributions SHALL all be scaled by the same weight factor.

**Phase-aware weighting**: During Phase 1 of two-phase backmapping (see `fix-backmap-resolution` spec), CG-weighted terms (`cg` keyword) SHALL apply full strength (weight = 1.0) regardless of `lambda_global`. During Phase 2, CG terms SHALL use the standard `1 − λ_global` weighting. AT-tagged terms (`at` keyword, both intra- and inter-bead) SHALL use the same weighting in both phases (phase does not affect AT weighting).

The phase is read from `fix backmap` via `fix->extract("phase")`. The `compute_weight3()` helper in `backmap_lambda.h` SHALL accept an optional phase parameter and return the correct weight; its default (no phase argument) behaves as Phase 2.

Reference: ESPResSo++ `FixedPairListAdressInteractionTemplate::addForces()` — lines 130-153 in `src/interaction/FixedPairListAdressInteractionTemplate.hpp` (weighting concept only; the global-λ / intra-inter-bead adaptation described here is specific to this package, confirmed with the maintainer during PET/Dacron network validation).

#### Scenario: Cross-CG AT bond (inter-bead) at lambda_global=0.5
- **WHEN** a `backmap/harmonic at` bond connects atoms i,j that map to different CG beads, with `lambda_global = 0.5`
- **THEN** the bond force SHALL be `0.5 × F_harmonic` (weight = `λ_global`)

#### Scenario: Cross-CG CG bond at lambda_global=0.5
- **WHEN** a `backmap/table cg` bond connects CG atoms i,j with `lambda_global = 0.5`
- **THEN** the bond force SHALL be `0.5 × F_table` (weight = `1 − 0.5 = 0.5`)

#### Scenario: AT intra-bead term at full strength mid-ramp
- **WHEN** an `at`-tagged interaction connects atoms that all map to the same CG bead in `atom2cg`, with `lambda_global = 0.1`
- **THEN** the force SHALL be computed at full strength (weight = 1), not scaled by 0.1 or any function of 0.1

#### Scenario: Bond force at lambda_global=0 (pure CG)
- **WHEN** `lambda_global = 0`
- **THEN** AT inter-bead cross bonds (`at`) SHALL have zero force, AT intra-bead terms SHALL have zero force (weight not yet `> 0`), CG cross bonds (`cg`) SHALL have full-strength force

#### Scenario: Bond force at lambda_global=1 (pure AT)
- **WHEN** `lambda_global = 1`
- **THEN** AT inter-bead cross bonds (`at`) SHALL have full-strength force (weight = `λ_global` = 1), AT intra-bead terms SHALL have full-strength force, CG cross bonds (`cg`) SHALL have zero force

#### Scenario: Skip computation when weight is negligible
- **WHEN** the computed weight is below a threshold (~1e-10)
- **THEN** the style SHALL skip the force computation entirely for efficiency (matching ESPResSo++'s `is_almost_zero()` check)

#### Scenario: Phase 1 CG forces at full strength
- **WHEN** `fix backmap` is in Phase 1 and `lambda_global = 0.5`
- **THEN** CG cross bonds (`cg`) SHALL have full-strength force (weight = 1.0, ignoring lambda)
- **AND** AT-tagged terms SHALL have weight 0.5 (`λ_global`, phase does not affect AT weighting)

#### Scenario: Phase 2 CG forces use standard weighting
- **WHEN** `fix backmap` is in Phase 2 and `lambda_global = 0.5`
- **THEN** CG cross bonds (`cg`) SHALL have weight 0.5 (`1 − λ_global`)
- **AND** AT-tagged terms SHALL have weight 0.5 (`λ_global`)

### Requirement: Bond style `backmap/harmonic`

The `bond_style backmap/harmonic` SHALL compute harmonic bond forces with lambda weighting:

```
F = -w × k × (r - r0)
E = w × 0.5 × k × (r - r0)²
```

where `w = BackmapLambda::compute_weight3(same_bead(atom2cg, i1, i2), is_cg, lambda_global)`: `1 − λ_global` for `cg`-tagged bonds; `1` (once `λ_global > 0`) for `at`-tagged bonds whose two atoms map to the same CG bead; `λ_global` for `at`-tagged bonds whose two atoms map to different CG beads.

Command syntax:
```
bond_style backmap/harmonic
bond_coeff N at/cg K r0
```

Where:
- `N` = bond type
- `at` or `cg` = weighting direction
- `K` = spring constant (kcal/(mol·Å²))
- `r0` = equilibrium distance (Å)

Used within `bond_style hybrid`:
```
bond_style hybrid harmonic backmap/harmonic
bond_coeff 1 harmonic 500.0 1.0          # intra-CG (static)
bond_coeff 2 backmap/harmonic at 500.0 1.0  # cross-CG AT
```

#### Scenario: Harmonic bond with AT inter-bead weighting
- **WHEN** `bond_coeff 2 backmap/harmonic at 500.0 1.0` connects atoms in different CG beads and `lambda_global = 0.8`
- **THEN** the bond force SHALL be `0.8 × (-500.0 × (r - 1.0))` (weight = `λ_global` = 0.8)

#### Scenario: Restart file support
- **WHEN** a simulation using `backmap/harmonic` is written to a restart file
- **THEN** the bond style parameters (K, r0, at/cg flag) SHALL be preserved in the restart

### Requirement: Bond style `backmap/table`

The `bond_style backmap/table` SHALL compute tabulated bond forces with lambda weighting:

```
F = w × F_table(r)
E = w × E_table(r)
```

where `w` is computed by `BackmapLambda::compute_weight3(same_bead(atom2cg, i1, i2), is_cg, lambda_global)`, following the same three-case rule as `backmap/harmonic`.

Command syntax:
```
bond_style backmap/table linear N
bond_coeff M at/cg filename keyword
```

Where:
- `M` = bond type
- `at` or `cg` = weighting direction
- `filename` = table file path
- `keyword` = table section name
- `N` = number of interpolation points

The table format SHALL follow LAMMPS `bond_style table` conventions.

#### Scenario: CG bond with tabulated potential
- **WHEN** `bond_coeff 3 backmap/table cg table_b1.table ENTRY` and `lambda_global = 0.3`
- **THEN** the bond force SHALL be `0.7 × F_table(r)` (weight = `1 − 0.3 = 0.7`)

### Requirement: Angle style `backmap/harmonic`

The `angle_style backmap/harmonic` SHALL compute harmonic angle forces with lambda weighting:

```
F = w × (-K × (θ - θ0))
E = w × 0.5 × K × (θ - θ0)²
```

The weight `w = BackmapLambda::compute_weight3(same_bead(atom2cg, i1, i2, i3), is_cg, lambda_global)`, where `same_bead` is true only if all three atoms of the angle triple map to the same CG bead in `atom2cg`.

Command syntax:
```
angle_style backmap/harmonic
angle_coeff N at/cg K theta0
```

Used within `angle_style hybrid`:
```
angle_style hybrid harmonic backmap/harmonic
angle_coeff 1 harmonic 50.0 109.47          # intra-CG (static)
angle_coeff 2 backmap/harmonic at 50.0 120.0  # cross-CG AT
```

#### Scenario: Angle with AT inter-bead weighting
- **WHEN** `angle_coeff 2 backmap/harmonic at 50.0 120.0` connects atoms spanning two different CG beads and `lambda_global = 0.6`
- **THEN** the angle force SHALL be weighted by 0.6 (`λ_global`)

#### Scenario: Angle with AT intra-bead weighting
- **WHEN** `angle_coeff 2 backmap/harmonic at 50.0 120.0` connects three atoms all mapped to the same CG bead and `lambda_global = 0.05`
- **THEN** the angle force SHALL be at full strength (weight = 1)

### Requirement: Angle style `backmap/table`

The `angle_style backmap/table` SHALL compute tabulated angle forces with lambda
weighting:

```
F = w × F_table(θ)
E = w × E_table(θ)
```

The weight `w = BackmapLambda::compute_weight3(same_bead(atom2cg, i1, i2, i3), is_cg, lambda_global)`, using the same three-case rule as `backmap/harmonic`.

Command syntax:

```
angle_style backmap/table linear N
angle_coeff M at/cg filename keyword
```

Used within `angle_style hybrid`:

```
angle_style hybrid backmap/harmonic backmap/table linear 1000
angle_coeff 2 backmap/table cg table_a1.table ENTRY
```

#### Scenario: CG tabulated angle at lambda_global=0

- **WHEN** `angle_coeff 2 backmap/table cg table_a1.table ENTRY` and `lambda_global = 0`
- **THEN** the angle force SHALL be at full tabulated strength (weight = 1.0)

#### Scenario: CG tabulated angle at lambda_global=0.5

- **WHEN** `angle_coeff 2 backmap/table cg table_a1.table ENTRY` and `lambda_global = 0.5`
- **THEN** the angle energy SHALL be scaled by 0.5 (weight = `1 − 0.5`)

#### Scenario: Missing fix backmap

- **WHEN** a `backmap/table` angle is defined but no `fix backmap` exists
- **THEN** the style SHALL abort with: "angle_style backmap/table requires fix backmap"

### Requirement: Dihedral style `backmap/ryckaert`

The `dihedral_style backmap/ryckaert` SHALL compute Ryckaert-Bellemans dihedral forces with lambda weighting:

```
E = w × Σ(n=0..5) Cn × cos^n(φ)
F_φ = w × Σ(n=1..5) n × Cn × cos^(n-1)(φ) × sin(φ)
```

where φ is the dihedral angle using the LAMMPS polymer convention (trans = 180°).

The weight `w = BackmapLambda::compute_weight3(same_bead(atom2cg, i1, i2, i3, i4), is_cg, lambda_global)`, where `same_bead` is true only if all four atoms of the dihedral quadruplet map to the same CG bead in `atom2cg`.

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

#### Scenario: RB dihedral with AT inter-bead weighting at lambda_global=1.0
- **WHEN** `dihedral_coeff 1 backmap/ryckaert at 1.53 0.776 -1.19 -3.22 0.0 0.0` connects atoms spanning two CG beads and `lambda_global = 1.0`
- **THEN** the dihedral force SHALL be at full strength (weight = `λ_global` = 1.0)

#### Scenario: RB dihedral with AT inter-bead weighting at lambda_global=0.5
- **WHEN** `dihedral_coeff 1 backmap/ryckaert at C0..C5` connects atoms spanning two CG beads and `lambda_global = 0.5`
- **THEN** the dihedral energy SHALL be `0.5 × Σ Cn × cos^n(φ)` (weight = `λ_global` = 0.5)

#### Scenario: RB dihedral with AT intra-bead weighting
- **WHEN** `dihedral_coeff 1 backmap/ryckaert at C0..C5` connects four atoms all mapped to the same CG bead and `lambda_global = 0.02`
- **THEN** the dihedral energy SHALL be computed at full strength (weight = 1)

#### Scenario: RB dihedral with CG weighting at lambda_global=0.5
- **WHEN** `dihedral_coeff 2 backmap/ryckaert cg C0..C5` and `lambda_global = 0.5`
- **THEN** the dihedral energy SHALL be `0.5 × Σ Cn × cos^n(φ)` (weight = `1 − 0.5`)

#### Scenario: RB dihedral at lambda_global=0 (pure CG)
- **WHEN** `lambda_global = 0`
- **THEN** AT dihedrals (`at`, inter- or intra-bead) SHALL have zero force, CG dihedrals (`cg`) SHALL have full-strength force

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

The weight `w = BackmapLambda::compute_weight3(same_bead(atom2cg, i1, i2, i3, i4), is_cg, lambda_global)`, using the same three-case rule as `backmap/ryckaert`.

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

#### Scenario: CG dihedral with tabulated potential at lambda_global=0
- **WHEN** `dihedral_coeff 2 backmap/table cg table_d1.table ENTRY` and `lambda_global = 0`
- **THEN** the dihedral force SHALL be at full strength (weight = `1 − 0 = 1.0`)

#### Scenario: CG dihedral with tabulated potential at lambda_global=0.7
- **WHEN** `dihedral_coeff 2 backmap/table cg table_d1.table ENTRY` and `lambda_global = 0.7`
- **THEN** the dihedral force SHALL be `0.3 × F_table(φ)` (weight = `1 − 0.7 = 0.3`)

#### Scenario: AT inter-bead dihedral with tabulated potential at lambda_global=1
- **WHEN** `dihedral_coeff 3 backmap/table at table_d2.table ENTRY` connects atoms spanning two CG beads and `lambda_global = 1`
- **THEN** the dihedral force SHALL be at full strength (weight = `λ_global` = 1.0)

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

### Requirement: Lambda access from fix

All `backmap/*` bonded styles SHALL read `atom2cg` and `lambda_global` from `fix backmap` via `fix->extract()`. They SHALL NOT maintain their own lambda state and SHALL NOT read the per-atom `"lambda"` array for weighting purposes. The fix is the single source of truth for both bead membership and the global λ scalar.

The styles SHALL locate the fix during initialization and cache a pointer to it. If no `fix backmap` is found, the style SHALL error with a descriptive message.

#### Scenario: Fix not defined
- **WHEN** a `backmap/harmonic` bond is defined but no `fix backmap` exists
- **THEN** the style SHALL abort with an error message: "bond_style backmap/harmonic requires fix backmap"

#### Scenario: Lambda updated between timesteps
- **WHEN** the fix increments `lambda_global` at timestep N
- **THEN** the bonded styles SHALL use the updated `lambda_global` value at timestep N+1

### Requirement: Shared lambda-access helper

All `backmap/*` styles (pair, bond, angle, dihedral) SHALL share a common helper header (`backmap_lambda.h`) for:
- Locating `fix backmap` by scanning the fix list
- Reading `atom2cg` and `lambda_global` from the fix's `extract()` interface
- Determining CG-bead co-membership of the atoms in an interaction (2 for pair/bond, 3 for angle, 4 for dihedral) via a `same_bead()` helper
- Computing the weight `w` given `same_bead`, the `at`/`cg` flag, `lambda_global`, and the current phase
- Applying the `is_almost_zero()` check to skip negligible-weight computations

This shared code SHALL be in a single header (`backmap_lambda.h`) to avoid duplication across styles.

The `compute_weight3(same_bead, is_cg, lambda_global, phase = 2)` function SHALL implement:
- `is_cg=true`, `phase=1` → return 1.0 (full strength, ignoring `lambda_global`)
- `is_cg=true`, `phase=2` → return `1 − lambda_global`
- `is_cg=false`, `same_bead=true` → return `(lambda_global > 0.0) ? 1.0 : 0.0` (phase does not affect this case)
- `is_cg=false`, `same_bead=false` → return `lambda_global` (phase does not affect this case)

There is no per-particle (`lambda_i`, `lambda_j`) variant of this helper; the previous `compute_weight(lambda_i, lambda_j, is_cg)` function is removed.

#### Scenario: Weight computation helper — AT inter-bead
- **WHEN** `compute_weight3(same_bead=false, is_cg=false, lambda_global=0.7)` is called
- **THEN** it SHALL return 0.7

#### Scenario: Weight computation helper — AT intra-bead
- **WHEN** `compute_weight3(same_bead=true, is_cg=false, lambda_global=0.7)` is called
- **THEN** it SHALL return 1.0

#### Scenario: Weight computation helper — CG
- **WHEN** `compute_weight3(same_bead=false, is_cg=true, lambda_global=0.7)` is called
- **THEN** it SHALL return 0.3 (`1 − 0.7`)

#### Scenario: Phase 1 CG weight override
- **WHEN** `compute_weight3(same_bead=false, is_cg=true, lambda_global=0.7, phase=1)` is called
- **THEN** it SHALL return 1.0 (full strength in Phase 1)

#### Scenario: Phase 1 AT weight unchanged
- **WHEN** `compute_weight3(same_bead=false, is_cg=false, lambda_global=0.7, phase=1)` is called
- **THEN** it SHALL return 0.7 (AT weighting is the same regardless of phase)

#### Scenario: AT intra-bead weight unaffected by lambda_global magnitude
- **WHEN** `compute_weight3(same_bead=true, is_cg=false, lambda_global=0.0001)` is called
- **THEN** it SHALL return 1.0 (step function at `lambda_global > 0`, not a scaled value)

### Requirement: Fix `backmap/pairs`

The package SHALL provide `fix backmap/pairs` for lambda-weighted explicit 1–4 LJ pairs
excluded from the neighbor list by `special_bonds`:

```
fix ID group backmap/pairs at file pairs.dat cut CUTOFF
```

The fix SHALL read `atom2cg` and `lambda_global` from `fix backmap` via `fix->extract()` (not per-atom `lambda[]`), and weight each pair using `BackmapLambda::compute_weight3(same_bead(atom2cg, i, j), is_cg, lambda_global)`, matching the same three-case rule as the other `backmap/*` styles.

#### Scenario: AT inter-bead 1–4 pair at lambda_global=0.5

- **WHEN** a pair uses keyword `at`, the two atoms map to different CG beads, and `lambda_global = 0.5`
- **THEN** the LJ force SHALL be scaled by 0.5 (weight = `lambda_global`)

#### Scenario: AT intra-bead 1–4 pair at any positive lambda_global

- **WHEN** a pair uses keyword `at`, the two atoms map to the same CG bead, and `lambda_global = 0.02`
- **THEN** the LJ force SHALL be at full strength (weight = 1)
