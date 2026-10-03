# pair_style backmap

## Syntax

```
pair_style backmap cut_at at_style at_args ... cut_cg cg_style cg_args ... [cg_special w12 w13 w14]
```

- **cut_at** -- cutoff for AT sub-style (distance units)
- **at_style** -- name of the AT pair sub-style (e.g., `lj/cut/coul/cut`)
- **at_args** -- arguments passed to the AT sub-style
- **cut_cg** -- cutoff for CG sub-style (distance units)
- **cg_style** -- name of the CG pair sub-style (e.g., `table`)
- **cg_args** -- arguments passed to the CG sub-style
- **cg_special** (optional) -- special-bond weights for CG-CG pairs at 1-2,
  1-3 and 1-4, replacing the `special_bonds` factors for those pairs

## Special bonds

Each listed pair is evaluated with its `special_bonds` factors (LJ and
Coulomb), as a plain pair style would; a factor below 1e-30 counts as an
exclusion. `cg_special` gives CG-CG pairs their own weights, for a CG model
that excludes fewer neighbours than the AT force field. A MARTINI model
(nrexcl = 1) with an AT force field at nrexcl = 3:

```
special_bonds lj 0.0 1.0e-100 1.0e-100 coul 0.0 1.0e-100 1.0e-100
pair_style backmap 12.0 lj/cut/coul/cut 12.0 9.0 11.0 table linear 1000 cg_special 0.0 1.0 1.0
```

The tiny nonzero factors keep 1-3 and 1-4 pairs in the neighbor list (LAMMPS
drops pairs whose factors are exactly zero); AT pairs stay excluded, CG pairs
get weight 1.

## Description

`pair_style backmap` is a hybrid-like pair style that delegates force
computation to two sub-styles (AT and CG) and weights the results by the
single global lambda value (`lambda_global`) and whether the pair's atoms
map to the same CG bead (see [theory: Force Weighting](../theory.md#force-weighting)).

For each pair of atoms *i* and *j*:

- **CG-type pairs** (tagged `cg` in `pair_coeff`):
  force is weighted by \( w_\text{CG} = 1 - \lambda_\text{global} \)
- **AT-type pairs** (tagged `atomistic` in `pair_coeff`), same CG bead:
  always full strength, unconditionally (real intra-molecular chemistry,
  independent of the resolution ramp)
- **AT-type pairs**, different CG beads:
  force is weighted by \( w_\text{AT} = \lambda_\text{global} \), evaluated
  smoothly from \( \lambda_\text{global} = 0 \) upward (see the note on
  `LAMBDA_AT_ONSET` in [theory: Force Weighting](../theory.md#force-weighting))
- **None pairs** (tagged `none`): no interaction

Both force and energy are scaled by the weight factor. Interactions with
negligible weight (< 10<sup>-10</sup>) are skipped for efficiency.

The weighted energy of each pair is tallied once, as van der Waals energy:
`evdwl` holds the whole sub-style energy (LJ and any cut-off Coulomb part), and
`ecoul` stays zero. At \( \lambda_\text{global} = 1 \) the pair energy of the
atomistic part equals that of the plain atomistic pair style.

## pair_coeff

```
pair_coeff I J atomistic at_args ...
pair_coeff I J cg cg_args ...
pair_coeff I J none
```

- **I, J** -- atom type pair (1-indexed, can use `*` for ranges)
- **atomistic** -- forward coefficients to the AT sub-style; weight is
  always full strength for pairs in the same CG bead, otherwise
  &lambda;<sub>global</sub> (see Description above)
- **cg** -- forward coefficients to the CG sub-style, weight by
  1 - &lambda;<sub>global</sub>
- **none** -- no interaction between these types

The AT and CG arguments after the keyword are passed directly to the
respective sub-style's `pair_coeff` command.

## Requirements

This pair style requires [`fix backmap`](fix-backmap.md) to be defined.
It reads the global lambda value and per-atom CG-bead membership from the
fix at each timestep.

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

- [fix backmap](fix-backmap.md) -- provides the global lambda value and
  per-atom CG-bead membership
- [bond_style backmap/harmonic](bond-styles.md) -- lambda-weighted bonds
