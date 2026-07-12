## MODIFIED Requirements

### Requirement: Lambda-weighted pair force computation

The pair style SHALL compute non-bonded pair forces as a weighted combination of atomistic and coarse-grained sub-styles, using a single **global** λ scalar (`lambda_global`, from `fix backmap`'s `extract("lambda_global")`) and per-atom CG-bead membership (`atom2cg`, from `extract("atom2cg")`) — not a per-particle product of two atoms' λ values. There are three weighting cases:

- **CG-CG pairs** (`is_cg=true`): `w = 1 − λ_global`
- **AT-AT intra-bead pairs** (`is_cg=false`, both atoms map to the same CG bead in `atom2cg`): `w = 1` once `λ_global > 0`, else `w = 0`. Intra-bead AT-AT interactions are never scaled down once the ramp has started — they act at full strength.
- **AT-AT inter-bead pairs** (`is_cg=false`, atoms map to different CG beads): `w = λ_global`

This uses the shared `BackmapLambda::compute_weight3(same_bead, is_cg, lambda_global)` helper (see the `bonded-backmap-styles` spec for the same helper's use in bond/angle/dihedral styles).

Cross-type (AT-CG) pairs continue to have no interaction, per the existing "Sub-style delegation" requirement.

Reference: Krajniak et al., JCTC 2016, DOI: 10.1021/acs.jctc.6b00595 — Section 2, Eq. 1-3 (adapted here to a single global λ and explicit intra/inter-bead distinction, per user clarification during the PET/Dacron network validation).

#### Scenario: Forces at lambda_global=0 (pure CG)
- **WHEN** `lambda_global = 0`
- **THEN** CG pair forces SHALL be at full strength (weight `1 − 0 = 1`), AT-AT intra-bead forces SHALL be zero (weight 0, since `λ_global` is not `> 0`), AT-AT inter-bead forces SHALL be zero (weight `λ_global = 0`)

#### Scenario: Forces at lambda_global=1 (pure AT)
- **WHEN** `lambda_global = 1`
- **THEN** CG pair forces SHALL be zero (weight `1 − 1 = 0`), AT-AT intra-bead forces SHALL be at full strength (weight 1), AT-AT inter-bead forces SHALL be at full strength (weight `λ_global = 1`)

#### Scenario: AT-AT intra-bead pair at full strength mid-ramp
- **WHEN** `lambda_global = 0.1` and two AT atoms `i, j` map to the same CG bead in `atom2cg`
- **THEN** the pair force between `i` and `j` SHALL be computed at full strength (weight = 1), not scaled by 0.1 or any function of 0.1

#### Scenario: AT-AT inter-bead pair scales linearly with lambda_global
- **WHEN** `lambda_global = 0.5` and two AT atoms `i, j` map to different CG beads
- **THEN** the pair force between `i` and `j` SHALL be weighted by 0.5

#### Scenario: CG pair scales linearly with lambda_global
- **WHEN** `lambda_global = 0.5` and both atoms are CG type
- **THEN** the CG pair force SHALL be weighted by 0.5 (`1 − 0.5`)

### Requirement: Lambda access from fix

The pair style SHALL read `lambda_global` and `atom2cg` from `fix backmap` via `fix->extract()`. It SHALL NOT maintain its own lambda state and SHALL NOT read the per-atom `"lambda"` array for weighting purposes. The fix is the single source of truth for both the global λ scalar and bead membership.

The `LAMBDA_AT_ONSET` numerical-safety guard (deferring AT-AT pair force computation until the ramp has progressed past a small threshold, to avoid LJ singularities on badly-overlapping initial AT placements) SHALL be evaluated against `lambda_global` directly (`lambda_global < LAMBDA_AT_ONSET`), not against any per-atom value.

#### Scenario: Lambda updated by fix before pair computation
- **WHEN** the fix increments `lambda_global` at timestep N
- **THEN** the pair style SHALL use the updated `lambda_global` value at timestep N+1

#### Scenario: LAMBDA_AT_ONSET guard uses the global scalar
- **WHEN** `lambda_global = 0.05` and `LAMBDA_AT_ONSET = 0.1`
- **THEN** the pair style SHALL skip AT-AT force computation entirely for that timestep, regardless of any individual atom's per-atom `lambda[]` value

### Requirement: Energy computation

The pair style SHALL correctly compute potential energy contributions weighted by lambda for thermodynamic output, using the same 3-way weighting as forces:
```
E_pair = w × E_sub(r)
```
where `w` is the appropriate weight from `compute_weight3(same_bead, is_cg, lambda_global)` for that pair (CG, AT intra-bead, or AT inter-bead), and `E_sub` is the CG or AT sub-style's energy.

#### Scenario: Energy at lambda_global endpoints
- **WHEN** `lambda_global = 0`
- **THEN** the pair energy SHALL equal the CG potential energy only (AT intra-bead and inter-bead contributions are zero)
- **WHEN** `lambda_global = 1`
- **THEN** the pair energy SHALL equal the full-strength AT potential energy (both intra-bead and inter-bead), with zero CG contribution

#### Scenario: Intra-bead AT energy at full strength mid-ramp
- **WHEN** `lambda_global = 0.2` and two AT atoms `i, j` map to the same CG bead
- **THEN** the pair energy contribution from `i,j` SHALL be the full, unweighted AT sub-style energy (weight = 1)
