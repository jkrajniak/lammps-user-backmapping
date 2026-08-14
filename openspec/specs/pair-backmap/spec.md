## ADDED Requirements

### Requirement: Lambda-weighted pair force computation

The pair style SHALL compute non-bonded pair forces as a weighted combination of atomistic and coarse-grained sub-styles, using a single global lambda scalar (`lambda_global`, read from `fix backmap`) and per-atom CG-bead co-membership (`atom2cg`) — not a per-particle product. Three cases:

- CG-CG pairs: weight `w = 1 − lambda_global`
- AT-AT pairs mapped to the *same* CG bead (intra-bead): weight `w = 1` always, unconditionally — this is real intra-molecular chemistry that exists independent of the resolution ramp
- AT-AT pairs mapped to *different* CG beads (inter-bead): weight `w = lambda_global`

Bonded forces are handled separately by `bond_style backmap/*`, `angle_style backmap/*`, and `dihedral_style backmap/*` using the same three-case weighting scheme (`BackmapLambda::compute_weight3`). See the bonded styles spec for details.

Reference: Krajniak et al., JCTC 2016, DOI: 10.1021/acs.jctc.6b00595 — Section 2, Eq. 1-3 (for the original per-pair AdResS-derived scheme this three-case model supersedes)

#### Scenario: Forces at lambda_global=0 (pure CG)
- **WHEN** `lambda_global` = 0
- **THEN** CG-CG pair forces SHALL be at full strength (weight 1), AT-AT inter-bead pair forces SHALL be zero (weight 0), AT-AT intra-bead pair forces SHALL remain at full strength (weight 1)

#### Scenario: Forces at lambda_global=1 (pure AT)
- **WHEN** `lambda_global` = 1
- **THEN** CG-CG pair forces SHALL be zero (weight 0), AT-AT inter-bead pair forces SHALL be at full strength (weight 1), AT-AT intra-bead pair forces SHALL remain at full strength (weight 1)

#### Scenario: Forces at lambda_global=0.5 (mid-transition)
- **WHEN** `lambda_global` = 0.5
- **THEN** CG-CG pair forces SHALL be weighted by 0.5, AT-AT inter-bead pair forces SHALL be weighted by 0.5, AT-AT intra-bead pair forces SHALL remain at weight 1

#### Scenario: Negative lambda treated as zero
- **WHEN** `lambda_global` < 0 (e.g. from a negative `lambda0`)
- **THEN** the effective value for force weighting SHALL be `max(0, lambda_global)`

### Requirement: Sub-style delegation

The pair style SHALL delegate actual force computation to two LAMMPS sub-styles:
- An atomistic sub-style for AT-AT type pairs (e.g., `lj/cut`, `lj/cut/coul/cut`)
- A coarse-grained sub-style for CG-CG type pairs (e.g., `table`)

Cross-type pairs (AT-CG) SHALL have no interaction.

The pair style SHALL support any pairwise sub-style that implements the `compute()` and `single()` methods.

#### Scenario: LJ atomistic + tabulated CG
- **WHEN** the pair style is configured with `lj/cut/coul/cut` for AT and `table` for CG
- **THEN** AT type pairs SHALL use LJ+Coulomb forces (weighted per the three-case model above), CG type pairs SHALL use tabulated forces (weighted by `1 − lambda_global`)

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

The pair style SHALL read `lambda_global` and `atom2cg` from `fix backmap` via `extract()`. It SHALL NOT maintain its own lambda state. The fix is the single source of truth.

#### Scenario: Lambda updated by fix before pair computation
- **WHEN** the fix increments `lambda_global` in `end_of_step()` at timestep N
- **THEN** the pair style SHALL use the updated value at timestep N+1

### Requirement: Energy computation

The pair style SHALL correctly compute potential energy contributions weighted by the same three-case model as forces, for thermodynamic output. The per-pair energy SHALL be `w * E_substyle`, using the same `w` as the force computation for that pair (CG: `1 − lambda_global`; AT intra-bead: `1`; AT inter-bead: `lambda_global`).

#### Scenario: Energy at lambda_global endpoints
- **WHEN** `lambda_global` = 0
- **THEN** CG-CG pair energy SHALL be at full strength and AT-AT inter-bead pair energy SHALL be zero (AT-AT intra-bead energy is unaffected, always full strength)
- **WHEN** `lambda_global` = 1
- **THEN** CG-CG pair energy SHALL be zero and AT-AT inter-bead pair energy SHALL be at full strength
