# Design — CG cross 1–4 pairs

## Physics

Cross 1–4 pairs use the same lambda weighting as other cross-CG AT interactions:
`w = λ_i × λ_j` (keyword `at`). LJ parameters follow GROMACS func-1 lookup:
`ε = fudgeLJ × sqrt(ε_i × ε_j)`, `σ = 0.5 × (σ_i + σ_j)` for OPLS combination rule 2.

## LAMMPS integration

`fix backmap/pairs` runs in `post_force`, iterating an explicit tag-pair list from
`pairs.dat`. Minimum-image convention applies. Requires `fix backmap`.

## Python routing

| Source | Style |
|--------|-------|
| `[ cross_pairs ]` func 1 lookup | `fix backmap/pairs at` |
| `cross_interactions.pairs` YAML | same (linear builder) |
