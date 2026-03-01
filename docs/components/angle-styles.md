# angle_style backmap/harmonic

## Syntax

```
angle_style backmap/harmonic
angle_coeff N at/cg K theta0
```

- **N** -- angle type number
- **at** or **cg** -- weighting mode
    - `at`: weight = &lambda;<sub>i</sub> &times; &lambda;<sub>k</sub>
      (fades in during backmapping)
    - `cg`: weight = 1 - &lambda;<sub>i</sub> &times; &lambda;<sub>k</sub>
      (fades out during backmapping)
- **K** -- force constant (energy/radian&sup2; units)
- **theta0** -- equilibrium angle (degrees)

## Description

A harmonic angle potential scaled by the lambda weight:

\[
E = w \times K (\theta - \theta_0)^2
\]

where \( w \) is computed from the lambda values of the **first and last**
atoms in the angle triplet (*i*-*j*-*k*), using \( \lambda_i \) and
\( \lambda_k \).

!!! note
    The equilibrium angle `theta0` is specified in **degrees** in the
    `angle_coeff` command but converted to radians internally.

## Lambda Weighting

For an angle *i*-*j*-*k*:

- **AT mode** (`at`): \( w = \lambda_i \times \lambda_k \)
- **CG mode** (`cg`): \( w = 1 - \lambda_i \times \lambda_k \)

The central atom *j*'s lambda value is not used in the weight calculation.
This is because cross-CG angles typically involve atoms at the boundaries
of adjacent beads, and the end atoms' resolution determines the relevance
of the interaction.

## Example

```
angle_style backmap/harmonic

# AT cross angle (fades in with lambda)
angle_coeff 1 at 126.67 111.0
```

In a hybrid setup:

```
angle_style hybrid harmonic backmap/harmonic

# Intra-bead AT angles (always active)
angle_coeff 1 harmonic 126.67 111.0

# Cross-bead AT angles (lambda-weighted)
angle_coeff 2 backmap/harmonic at 126.67 111.0
```

## Requirements

This angle style requires [`fix backmap`](fix-backmap.md) to be defined.
It reads per-atom lambda values from the fix at each timestep.

## Restart

Angle coefficients (K, theta0, is_cg) are written to and read from restart
files.

## Related

- [fix backmap](fix-backmap.md) -- provides per-atom lambda values
- [bond_style backmap/harmonic](bond-styles.md) -- lambda-weighted bonds
- [bond_style backmap/table](bond-styles.md#backmap-table) -- lambda-weighted
  tabulated bonds
- [pair_style backmap](pair-backmap.md) -- lambda-weighted pair interactions
