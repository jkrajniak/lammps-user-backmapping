## ADDED Requirements

### Requirement: Fix `backmap/pairs`

The package SHALL provide `fix backmap/pairs` for lambda-weighted explicit 1–4 LJ pairs
excluded from the neighbor list by `special_bonds`:

```
fix ID group backmap/pairs at file pairs.dat cut CUTOFF
```

#### Scenario: AT cross 1–4 at lambda=0.5

- **WHEN** a pair uses keyword `at` and λ_i = λ_j = 0.5
- **THEN** the LJ force SHALL be scaled by 0.25 (weight = λ_i × λ_j)
