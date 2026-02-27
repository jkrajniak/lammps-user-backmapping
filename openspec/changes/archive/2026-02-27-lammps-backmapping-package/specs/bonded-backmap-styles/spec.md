## ADDED Requirements

### Requirement: Lambda-weighted bonded force computation

The `backmap/*` bonded styles SHALL compute bonded forces (bonds, angles, dihedrals) with lambda weighting applied at force-computation time. The weighting scheme is identical to the non-bonded pair weighting in `pair_style backmap` and matches ESPResSo++'s `FixedPairListAdressInteractionTemplate`:

- **AT weight** (keyword `at`): `w = λ_i × λ_j` — for the two atoms forming the bond, or the first and last atoms of an angle/dihedral
- **CG weight** (keyword `cg`): `w = 1 − λ_i × λ_j`

For uniform backmapping (all particles share the same λ), this reduces to `w_AT = λ²` and `w_CG = 1 − λ²`, ensuring `w_AT + w_CG = 1` always.

Each style reads per-atom λ from `fix backmap`'s per-atom array via `fix->extract()`. The fix ID is auto-detected (single `fix backmap` in the simulation) or specified via a style argument.

Force, energy, and virial contributions SHALL all be scaled by the same weight factor.

Reference: ESPResSo++ `FixedPairListAdressInteractionTemplate::addForces()` — lines 130-153 in `src/interaction/FixedPairListAdressInteractionTemplate.hpp`

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

### Requirement: Bond style `backmap/harmonic`

The `bond_style backmap/harmonic` SHALL compute harmonic bond forces with lambda weighting:

```
F = -w × k × (r - r0)
E = w × 0.5 × k × (r - r0)²
```

where `w` is `λ_i × λ_j` (for `at`) or `1 − λ_i × λ_j` (for `cg`).

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
- **WHEN** `bond_coeff 2 backmap/harmonic at 500.0 1.0` and λ = 0.8
- **THEN** the bond force SHALL be `0.64 × (-500.0 × (r - 1.0))` (weight = 0.8²)

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
- **WHEN** `bond_coeff 3 backmap/table cg table_b1.table ENTRY` and λ = 0.3
- **THEN** the bond force SHALL be `0.91 × F_table(r)` (weight = 1 − 0.09 = 0.91)

### Requirement: Angle style `backmap/harmonic`

The `angle_style backmap/harmonic` SHALL compute harmonic angle forces with lambda weighting:

```
F = w × (-K × (θ - θ0))
E = w × 0.5 × K × (θ - θ0)²
```

The weight `w` is computed from the lambda values of the first and last atoms of the angle (atoms i and k in the i-j-k triple).

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
- **WHEN** `angle_coeff 2 backmap/harmonic at 50.0 120.0` and λ = 0.6
- **THEN** the angle force SHALL be weighted by 0.36 (0.6²)

### Requirement: Dihedral style `backmap/ryckaert`

The `dihedral_style backmap/ryckaert` SHALL compute Ryckaert-Bellemans dihedral forces with lambda weighting:

```
E = w × Σ(n=0..5) Cn × cos^n(φ)
```

The weight `w` is computed from the lambda values of the first and last atoms of the dihedral (atoms i and l in the i-j-k-l quadruplet).

Command syntax:
```
dihedral_style backmap/ryckaert
dihedral_coeff N at/cg C0 C1 C2 C3 C4 C5
```

#### Scenario: RB dihedral with AT weighting
- **WHEN** `dihedral_coeff 1 backmap/ryckaert at ...` and λ = 1.0
- **THEN** the dihedral force SHALL be at full strength (weight = 1.0)

### Requirement: Dihedral style `backmap/table`

The `dihedral_style backmap/table` SHALL compute tabulated dihedral forces with lambda weighting.

Command syntax:
```
dihedral_style backmap/table linear N
dihedral_coeff M at/cg filename keyword
```

The table format SHALL follow LAMMPS `dihedral_style table` conventions (angle in degrees, energy in kcal/mol).

#### Scenario: CG dihedral with tabulated potential
- **WHEN** `dihedral_coeff 2 backmap/table cg table_d1.table ENTRY` and λ = 0
- **THEN** the dihedral force SHALL be at full strength (weight = 1 − 0 = 1.0)

### Requirement: Lambda access from fix

All `backmap/*` bonded styles SHALL read per-atom λ values from `fix backmap`'s per-atom array. They SHALL NOT maintain their own lambda state. The fix is the single source of truth for lambda values.

The styles SHALL locate the fix during initialization and cache a pointer to it. If no `fix backmap` is found, the style SHALL error with a descriptive message.

#### Scenario: Fix not defined
- **WHEN** a `backmap/harmonic` bond is defined but no `fix backmap` exists
- **THEN** the style SHALL abort with an error message: "bond_style backmap/harmonic requires fix backmap"

#### Scenario: Lambda updated between timesteps
- **WHEN** the fix increments λ in `end_of_step()` at timestep N
- **THEN** the bonded styles SHALL use the updated λ values at timestep N+1

### Requirement: Shared lambda-access helper

All `backmap/*` styles (pair, bond, angle, dihedral) SHALL share a common helper class or set of utility functions for:
- Locating `fix backmap` by scanning the fix list
- Reading per-atom λ values from the fix's per-atom array
- Computing the weight `w` given two lambda values and the `at`/`cg` flag
- Applying the `is_almost_zero()` check to skip negligible-weight computations

This shared code SHALL be in a single header (e.g., `backmap_lambda.h`) to avoid duplication across styles.

#### Scenario: Weight computation helper
- **WHEN** `compute_weight(lambda_i=0.7, lambda_j=0.7, is_cg=false)` is called
- **THEN** it SHALL return 0.49 (0.7 × 0.7)
- **WHEN** `compute_weight(lambda_i=0.7, lambda_j=0.7, is_cg=true)` is called
- **THEN** it SHALL return 0.51 (1 − 0.49)
