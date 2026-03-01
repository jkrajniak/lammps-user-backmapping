# pair_style backmap

## Syntax

```
pair_style backmap cut_at at_style at_args ... cut_cg cg_style cg_args ...
```

- **cut_at** -- cutoff for AT sub-style (distance units)
- **at_style** -- name of the AT pair sub-style (e.g., `lj/cut/coul/cut`)
- **at_args** -- arguments passed to the AT sub-style
- **cut_cg** -- cutoff for CG sub-style (distance units)
- **cg_style** -- name of the CG pair sub-style (e.g., `table`)
- **cg_args** -- arguments passed to the CG sub-style

## Description

`pair_style backmap` is a hybrid-like pair style that delegates force
computation to two sub-styles (AT and CG) and weights the results by the
current lambda values of the interacting atoms.

For each pair of atoms *i* and *j*:

- **AT-type pairs** (tagged `atomistic` in `pair_coeff`):
  force is weighted by \( w_\text{AT} = \lambda_i \times \lambda_j \)
- **CG-type pairs** (tagged `cg` in `pair_coeff`):
  force is weighted by \( w_\text{CG} = 1 - \lambda_i \times \lambda_j \)
- **None pairs** (tagged `none`): no interaction

Both force and energy are scaled by the weight factor. Interactions with
negligible weight (< 10<sup>-10</sup>) are skipped for efficiency.

## pair_coeff

```
pair_coeff I J atomistic at_args ...
pair_coeff I J cg cg_args ...
pair_coeff I J none
```

- **I, J** -- atom type pair (1-indexed, can use `*` for ranges)
- **atomistic** -- forward coefficients to the AT sub-style, weight by
  &lambda;&lambda;
- **cg** -- forward coefficients to the CG sub-style, weight by
  1 - &lambda;&lambda;
- **none** -- no interaction between these types

The AT and CG arguments after the keyword are passed directly to the
respective sub-style's `pair_coeff` command.

## Requirements

This pair style requires [`fix backmap`](fix-backmap.md) to be defined.
It reads per-atom lambda values from the fix at each timestep.

## Example

```
# AT sub-style: LJ + Coulomb with 14 Å cutoff
# CG sub-style: tabulated with 14 Å cutoff
pair_style backmap 14.0 lj/cut/coul/cut 14.0 9.0 14.0 table linear 1000

# CG-CG interactions (fade out with lambda)
pair_coeff 1 1 cg 0.0 0.0
pair_coeff 1 2 cg 0.0 0.0
pair_coeff 2 2 cg 0.0 0.0

# Cross-resolution interactions (no interaction)
pair_coeff 1 3 none
pair_coeff 1 4 none
pair_coeff 2 3 none
pair_coeff 2 4 none

# AT-AT interactions (fade in with lambda)
pair_coeff 3 3 atomistic 0.207266 3.748000
pair_coeff 3 4 atomistic 0.156387 3.826500
pair_coeff 4 4 atomistic 0.117997 3.905000
```

In this example:

- Types 1-2 are CG beads; they interact via tabulated CG potentials
- Types 3-4 are AT atoms; they interact via LJ + Coulomb
- CG-AT cross pairs have no direct interaction (forces are distributed
  through `fix backmap`)

## Restart

Cutoff values are written to restart files. Sub-style parameters must be
re-specified after a restart.

## Related

- [fix backmap](fix-backmap.md) -- provides per-atom lambda values
- [bond_style backmap/harmonic](bond-styles.md) -- lambda-weighted bonds
