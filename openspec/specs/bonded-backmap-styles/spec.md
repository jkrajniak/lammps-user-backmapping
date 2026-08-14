## ADDED Requirements

### Requirement: Lambda-weighted bonded force computation

The `backmap/*` bonded styles SHALL compute bonded forces (bonds, angles, dihedrals) with lambda weighting applied at force-computation time, using the same single global lambda scalar (`lambda_global`, read from `fix backmap`) and per-atom CG-bead co-membership (`atom2cg`) as `pair_style backmap` — not a per-particle product. Each `bond_coeff`/`angle_coeff`/`dihedral_coeff` line is tagged `at` or `cg`, and the weight follows the same three-case model as the pair style (`BackmapLambda::compute_weight3`):

- All atoms of the interaction mapped to the *same* CG bead (intra-bead): weight `w = 1` always, regardless of the `at`/`cg` tag
- Otherwise, tag `cg`: weight `w = 1 − lambda_global`
- Otherwise, tag `at`: weight `w = lambda_global`

Each style reads `lambda_global` and `atom2cg` from `fix backmap` via `fix->extract()`. The fix ID is auto-detected (single `fix backmap` in the simulation) or specified via a style argument.

Force, energy, and virial contributions SHALL all be scaled by the same weight factor.

Reference: `backmap_lambda_weights.h::compute_weight3()`, shared by all `backmap/*` pair/bond/angle/dihedral styles.

#### Scenario: Cross-CG AT bond, inter-bead
- **WHEN** a `backmap/harmonic at` bond connects atoms i,j mapped to different CG beads, with `lambda_global` = 0.5
- **THEN** the bond force SHALL be weighted by 0.5 (`w = lambda_global`)

#### Scenario: Cross-CG CG bond, inter-bead
- **WHEN** a `backmap/table cg` bond connects CG atoms i,j (necessarily in different CG beads — a CG atom does not map into another CG bead), with `lambda_global` = 0.5
- **THEN** the bond force SHALL be weighted by 0.5 (`w = 1 − lambda_global`)

#### Scenario: Bond force at lambda_global=0 (pure CG)
- **WHEN** `lambda_global` = 0
- **THEN** inter-bead AT bonds (`at`) SHALL have zero force, inter-bead CG bonds (`cg`) SHALL have full-strength force

#### Scenario: Bond force at lambda_global=1 (pure AT)
- **WHEN** `lambda_global` = 1
- **THEN** inter-bead AT bonds (`at`) SHALL have full-strength force, inter-bead CG bonds (`cg`) SHALL have zero force

#### Scenario: Skip computation when weight is negligible
- **WHEN** the computed weight is below a threshold (~1e-10)
- **THEN** the style SHALL skip the force computation entirely for efficiency (matching ESPResSo++'s `is_almost_zero()` check)

### Requirement: Bond style `backmap/harmonic`

The `bond_style backmap/harmonic` SHALL compute harmonic bond forces with lambda weighting:

```
F = -w × k × (r - r0)
E = w × 0.5 × k × (r - r0)²
```

where `w` follows the three-case model above.

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

#### Scenario: Harmonic bond with AT weighting
- **WHEN** `bond_coeff 2 backmap/harmonic at 500.0 1.0`, the bond is inter-bead, and `lambda_global` = 0.8
- **THEN** the bond force SHALL be `0.8 × (-500.0 × (r - 1.0))` (weight = `lambda_global`)

#### Scenario: Restart file support
- **WHEN** a simulation using `backmap/harmonic` is written to a restart file
- **THEN** the bond style parameters (K, r0, at/cg flag) SHALL be preserved in the restart

### Requirement: Bond style `backmap/table`

The `bond_style backmap/table` SHALL compute tabulated bond forces with lambda weighting:

```
F = w × F_table(r)
E = w × E_table(r)
```

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
- **WHEN** `bond_coeff 3 backmap/table cg table_b1.table ENTRY`, the bond is inter-bead, and `lambda_global` = 0.3
- **THEN** the bond force SHALL be `0.7 × F_table(r)` (weight = `1 − lambda_global`)

### Requirement: Angle style `backmap/harmonic`

The `angle_style backmap/harmonic` SHALL compute harmonic angle forces with lambda weighting:

```
F = w × (-K × (θ - θ0))
E = w × 0.5 × K × (θ - θ0)²
```

The weight `w` follows the three-case model above, evaluated across the first and last atoms of the angle (atoms i and k in the i-j-k triple).

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

#### Scenario: Angle with AT weighting
- **WHEN** `angle_coeff 2 backmap/harmonic at 50.0 120.0`, the angle is inter-bead, and `lambda_global` = 0.6
- **THEN** the angle force SHALL be weighted by 0.6 (`w = lambda_global`)

### Requirement: Angle style `backmap/table`

The `angle_style backmap/table` SHALL compute tabulated angle forces with lambda
weighting:

```
F = w × F_table(θ)
E = w × E_table(θ)
```

The weight `w` follows the same three-case model, evaluated across the first and
last atoms of the angle (atoms i and k in the i-j-k triple).

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

- **WHEN** `angle_coeff 2 backmap/table cg table_a1.table ENTRY`, the angle is inter-bead, and `lambda_global` = 0
- **THEN** the angle force SHALL be at full tabulated strength (weight = 1.0)

#### Scenario: CG tabulated angle at lambda_global=0.5

- **WHEN** `angle_coeff 2 backmap/table cg table_a1.table ENTRY`, the angle is inter-bead, and `lambda_global` = 0.5
- **THEN** the angle energy SHALL be scaled by 0.5 (weight = `1 − lambda_global`)

#### Scenario: Missing fix backmap

- **WHEN** a `backmap/table` angle is defined but no `fix backmap` exists
- **THEN** the style SHALL abort with: "angle_style backmap/table requires fix backmap"

### Requirement: Dihedral style `ryckaert` (static)

The package SHALL provide `dihedral_style ryckaert` for intra-bead RB dihedrals without
lambda weighting:

```
dihedral_style ryckaert
dihedral_coeff N C0 C1 C2 C3 C4 C5
```

#### Scenario: Intra-bead RB dihedral

- **WHEN** all four atoms of a dihedral belong to the same CG bead
- **THEN** the generator SHALL assign a static `ryckaert` type with converted C0..C5

### Requirement: Dihedral style `backmap/ryckaert`

The `dihedral_style backmap/ryckaert` SHALL compute Ryckaert-Bellemans dihedral forces with lambda weighting:

```
E = w × Σ(n=0..5) Cn × cos^n(φ)
F_φ = w × Σ(n=1..5) n × Cn × cos^(n-1)(φ) × sin(φ)
```

where φ is the dihedral angle using the LAMMPS polymer convention (trans = 180°).

The weight `w` follows the three-case model, evaluated across all four atoms of the dihedral (atoms i,j,k,l in the i-j-k-l quadruplet) via `backmap_lambda_weights.h`'s `same_bead()` four-atom overload.

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

#### Scenario: RB dihedral with AT weighting at lambda_global=1.0
- **WHEN** `dihedral_coeff 1 backmap/ryckaert at 1.53 0.776 -1.19 -3.22 0.0 0.0`, the dihedral is inter-bead, and `lambda_global` = 1.0
- **THEN** the dihedral force SHALL be at full strength (weight = 1.0)

#### Scenario: RB dihedral with AT weighting at lambda_global=0.5
- **WHEN** `dihedral_coeff 1 backmap/ryckaert at C0..C5`, the dihedral is inter-bead, and `lambda_global` = 0.5
- **THEN** the dihedral energy SHALL be `0.5 × Σ Cn × cos^n(φ)` (weight = `lambda_global`)

#### Scenario: RB dihedral with CG weighting at lambda_global=0.5
- **WHEN** `dihedral_coeff 2 backmap/ryckaert cg C0..C5`, the dihedral is inter-bead, and `lambda_global` = 0.5
- **THEN** the dihedral energy SHALL be `0.5 × Σ Cn × cos^n(φ)` (weight = `1 − lambda_global`)

#### Scenario: RB dihedral at lambda_global=0 (pure CG)
- **WHEN** `lambda_global` = 0
- **THEN** inter-bead AT dihedrals (`at`) SHALL have zero force, inter-bead CG dihedrals (`cg`) SHALL have full-strength force

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

#### Scenario: CG dihedral with tabulated potential at lambda_global=0
- **WHEN** `dihedral_coeff 2 backmap/table cg table_d1.table ENTRY` and `lambda_global` = 0
- **THEN** the dihedral force SHALL be at full strength (weight = `1 − 0` = 1.0)

#### Scenario: CG dihedral with tabulated potential at lambda_global=0.7
- **WHEN** `dihedral_coeff 2 backmap/table cg table_d1.table ENTRY`, the dihedral is inter-bead, and `lambda_global` = 0.7
- **THEN** the dihedral force SHALL be `0.3 × F_table(φ)` (weight = `1 − lambda_global`)

#### Scenario: AT dihedral with tabulated potential at lambda_global=1
- **WHEN** `dihedral_coeff 3 backmap/table at table_d2.table ENTRY` and `lambda_global` = 1
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

### Requirement: Lambda access from fix

All `backmap/*` bonded styles SHALL read `lambda_global` and `atom2cg` from `fix backmap` via `extract()`. They SHALL NOT maintain their own lambda state. The fix is the single source of truth.

The styles SHALL locate the fix during initialization and cache a pointer to it. If no `fix backmap` is found, the style SHALL error with a descriptive message.

#### Scenario: Fix not defined
- **WHEN** a `backmap/harmonic` bond is defined but no `fix backmap` exists
- **THEN** the style SHALL abort with an error message: "bond_style backmap/harmonic requires fix backmap"

#### Scenario: Lambda updated between timesteps
- **WHEN** the fix increments `lambda_global` in `end_of_step()` at timestep N
- **THEN** the bonded styles SHALL use the updated value at timestep N+1

### Requirement: Shared lambda-access helper

All `backmap/*` styles (pair, bond, angle, dihedral) SHALL share a common helper for:
- Locating `fix backmap` by scanning the fix list (`find_fix_backmap()`)
- Reading `lambda_global` and `atom2cg` from the fix (`extract_lambda_global()`, `extract_atom2cg()`)
- Computing the weight `w` given same-bead status and the `at`/`cg` flag (`compute_weight3()`)
- Applying the `is_almost_zero()` check to skip negligible-weight computations

This shared code SHALL be split across two headers: `backmap_lambda.h` (fix lookup and `extract()` accessors, requires LAMMPS headers) and `backmap_lambda_weights.h` (pure weighting math — `compute_weight3()`, `same_bead()`, `clamp_lambda()`, `is_almost_zero()` — deliberately free of LAMMPS headers so it can be unit-tested standalone).

`compute_weight3(bool is_same_bead, bool is_cg, double lambda_global)` SHALL return:
- `is_same_bead=true` → `1.0` (full strength, regardless of `is_cg`)
- `is_same_bead=false, is_cg=true` → `1 − lambda_global`
- `is_same_bead=false, is_cg=false` → `lambda_global`

#### Scenario: Weight computation helper, AT inter-bead
- **WHEN** `compute_weight3(is_same_bead=false, is_cg=false, lambda_global=0.7)` is called
- **THEN** it SHALL return 0.7

#### Scenario: Weight computation helper, CG inter-bead
- **WHEN** `compute_weight3(is_same_bead=false, is_cg=true, lambda_global=0.7)` is called
- **THEN** it SHALL return 0.3 (`1 − 0.7`)

#### Scenario: Weight computation helper, same-bead override
- **WHEN** `compute_weight3(is_same_bead=true, is_cg=false, lambda_global=0.7)` is called
- **THEN** it SHALL return 1.0 (full strength, independent of `lambda_global`)

### Requirement: Fix `backmap/pairs`

The package SHALL provide `fix backmap/pairs` for lambda-weighted explicit 1–4 LJ pairs
excluded from the neighbor list by `special_bonds`:

```
fix ID group backmap/pairs at file pairs.dat cut CUTOFF
```

#### Scenario: AT cross 1–4 at lambda_global=0.5

- **WHEN** a pair uses keyword `at`, is inter-bead, and `lambda_global` = 0.5
- **THEN** the pair force SHALL be weighted by 0.5 (`w = lambda_global`)
