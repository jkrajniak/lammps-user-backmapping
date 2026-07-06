# CG cross 1–4 pairs (PR2c)

## Why

Rim135 hybrid topology has **3,600** `[ cross_pairs ]` (func 1, LJ lookup). With
`special_bonds lj 0.0 0.0 0.0` these 1–4 interactions are excluded from the neighbor
list and must be re-injected with lambda weighting.

## What changes

- Add `fix backmap/pairs` for explicit lambda-weighted 1–4 LJ pairs.
- Parse `[ cross_pairs ]` from hybrid TOP; resolve LJ params (fudgeLJ + combination rule).
- Emit `pairs.dat` and wire `fix backmap/pairs at` in generated input.

## Non-goals

- Coulomb on cross 1–4 (GROMACS func-1 pairs are LJ-only lookup in rim135).
- CG-level tabulated 1–4 pairs (future).
