## ADDED Requirements

### Requirement: Lambda-weighted pair force computation

The pair style SHALL compute non-bonded pair forces as a weighted combination of atomistic and coarse-grained sub-styles. The weighting follows the per-pair scheme from AdResS, identical to the bonded weighting in `FixedPairListAdressInteractionTemplate` and the non-bonded weighting in `VerletListAdressInteractionTemplate`:

- AT pair force weight: `w_AT = λ_i × λ_j`
- CG pair force weight: `w_CG = 1 − λ_i × λ_j`

For uniform backmapping (all particles share the same λ), this reduces to `w_AT = λ²` and `w_CG = 1 − λ²`. This ensures `w_AT + w_CG = 1` always — a physically consistent interpolation.

Bonded forces are handled separately by `bond_style backmap/*`, `angle_style backmap/*`, and `dihedral_style backmap/*` using the same weighting scheme. See the bonded styles spec for details on the three categories of bonded interactions (intra-CG static, cross-CG AT, cross-CG CG).

Reference: Krajniak et al., JCTC 2016, DOI: 10.1021/acs.jctc.6b00595 — Section 2, Eq. 1-3

#### Scenario: Forces at lambda=0 (pure CG)
- **WHEN** λ = 0 for all particles
- **THEN** AT pair forces SHALL be zero (weight 0×0=0), CG pair forces SHALL be at full strength (weight 1−0=1)

#### Scenario: Forces at lambda=1 (pure AT)
- **WHEN** λ = 1 for all particles
- **THEN** AT pair forces SHALL be at full strength (weight 1×1=1), CG pair forces SHALL be zero (weight 1−1=0)

#### Scenario: Forces at lambda=0.5 (mid-transition)
- **WHEN** λ = 0.5 for all particles
- **THEN** AT pair forces SHALL be weighted by 0.25 (0.5×0.5), CG pair forces SHALL be weighted by 0.75 (1−0.25)

#### Scenario: Negative lambda treated as zero
- **WHEN** λ < 0 for a particle (from nonuniform initialization)
- **THEN** the effective lambda for force weighting SHALL be max(0, λ), so AT weight = 0 and CG weight = 1

### Requirement: Sub-style delegation

The pair style SHALL delegate actual force computation to two LAMMPS sub-styles:
- An atomistic sub-style for AT-AT type pairs (e.g., `lj/cut`, `lj/cut/coul/cut`)
- A coarse-grained sub-style for CG-CG type pairs (e.g., `table`)

Cross-type pairs (AT-CG) SHALL have no interaction.

The pair style SHALL support any pairwise sub-style that implements the `compute()` and `single()` methods.

#### Scenario: LJ atomistic + tabulated CG
- **WHEN** the pair style is configured with `lj/cut/coul/cut` for AT and `table` for CG
- **THEN** AT type pairs SHALL use LJ+Coulomb forces (weighted by λ²), CG type pairs SHALL use tabulated forces (weighted by 1−λ²)

#### Scenario: No cross-type interactions
- **WHEN** pair_coeff is not defined for an AT-CG type pair
- **THEN** no forces SHALL be computed between AT and CG atoms of different molecules

### Requirement: Pair style command syntax

The pair style SHALL be invoked with the following syntax:
```
pair_style backmap cutoff_at at_style at_args ... cutoff_cg cg_style cg_args ...
```

Pair coefficients SHALL specify which type pairs use which sub-style:
```
pair_coeff I J atomistic at_pair_args ...
pair_coeff I J cg cg_pair_args ...
pair_coeff I J none
```

The fix ID for `fix backmap` SHALL be passed via:
```
pair_coeff * * fix backmap_fix_id
```
or auto-detected from the defined fixes.

#### Scenario: Water2 system setup
- **WHEN** the user configures pair_style backmap for SPC water + WCG beads
- **THEN** OW-OW, OW-H, H-H pairs SHALL use the atomistic sub-style, WCG-WCG SHALL use the CG sub-style, and OW-WCG and H-WCG SHALL be set to none

### Requirement: Lambda access from fix

The pair style SHALL read per-atom λ values from `fix backmap`'s per-atom array. It SHALL NOT maintain its own lambda state. The fix is the single source of truth for lambda values.

#### Scenario: Lambda updated by fix before pair computation
- **WHEN** the fix increments λ in `end_of_step()` at timestep N
- **THEN** the pair style SHALL use the updated λ values at timestep N+1

### Requirement: Energy computation

The pair style SHALL correctly compute potential energy contributions weighted by lambda for thermodynamic output. The per-pair energy SHALL be:
```
E_pair = w_AT × E_AT + w_CG × E_CG
```

#### Scenario: Energy at lambda endpoints
- **WHEN** λ = 0
- **THEN** the pair energy SHALL equal the CG potential energy only
- **WHEN** λ = 1
- **THEN** the pair energy SHALL equal the AT potential energy only
