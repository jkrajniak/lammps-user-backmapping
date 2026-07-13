# angle_style backmap/harmonic

## Syntax

```
angle_style backmap/harmonic
angle_coeff N at/cg K theta0
```

- **N** -- angle type number
- **at** or **cg** -- weighting mode
    - `at`: always full strength if all three atoms map to the same CG bead
      (real intra-molecular chemistry, independent of &lambda;), otherwise
      weight = &lambda;<sub>global</sub> (fades in during backmapping)
    - `cg`: weight = 1 - &lambda;<sub>global</sub>
      (fades out during backmapping)
- **K** -- force constant (energy/radian&sup2; units)
- **theta0** -- equilibrium angle (degrees)

## Description

A harmonic angle potential scaled by the lambda weight:

\[
E = w \times K (\theta - \theta_0)^2
\]

where \( w \) is computed from the single global &lambda; value and whether
all three atoms in the angle triplet (*i*-*j*-*k*) map to the same CG bead
(see [theory: Force Weighting](../theory.md#force-weighting)).

!!! note
    The equilibrium angle `theta0` is specified in **degrees** in the
    `angle_coeff` command but converted to radians internally.

## Lambda Weighting

For an angle *i*-*j*-*k*:

- **AT mode** (`at`), all three atoms in the same CG bead: \( w = 1 \) once
  \( \lambda_\text{global} > 0 \)
- **AT mode** (`at`), atoms span different CG beads: \( w = \lambda_\text{global} \)
- **CG mode** (`cg`): \( w = 1 - \lambda_\text{global} \)

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

# angle_style backmap/table

## Syntax

```
angle_style backmap/table linear N
angle_coeff M cg filename keyword
```

- **N** -- number of interpolation points (typically 1000)
- **M** -- angle type number
- **cg** -- weighting mode (CG cross-angles use `cg`; weight = 1 − λ<sub>i</sub> × λ<sub>k</sub>)
- **filename** -- LAMMPS table file converted from GROMACS `table_a*.xvg`
- **keyword** -- table section label (usually `ENTRY`)

## Description

Tabulated angle potential for CG cross-bead angles, scaled by the same
lambda weight as [`angle_style backmap/harmonic`](#angle_style-backmapharmonic):

\[
E = w \times U(\theta)
\]

where \( U(\theta) \) comes from the table (angle in degrees, energy in
kcal/mol) and \( w = 1 - \lambda_\text{global} \) for `cg` mode.
At λ = 0 (pure CG), the full tabulated torque applies.

## Example

```
angle_style hybrid backmap/harmonic backmap/table linear 1000

angle_coeff 588 backmap/table cg table_a1.table ENTRY
angle_coeff 589 backmap/table cg table_a2.table ENTRY
```

## Requirements

Requires [`fix backmap`](fix-backmap.md). Source tables are GROMACS
`table_a*.xvg` files (angle in degrees); `backmap-prep` converts them to
LAMMPS `.table` format via `convert_tables`.

## Related

- [`angle_style backmap/harmonic`](#angle_style-backmapharmonic) -- harmonic cross angles
- [`bond_style backmap/table`](bond-styles.md#backmap-table) -- tabulated cross bonds
