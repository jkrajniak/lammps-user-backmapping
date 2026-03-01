# fix backmap

## Syntax

```
fix ID group-ID backmap cg_type T alpha A lambda0 L0 [nonuniform yes/no]
```

- **ID** -- fix identifier (user-chosen)
- **group-ID** -- group of atoms the fix applies to (typically `all`)
- **cg_type** -- atom type number of CG beads (integer, 1-indexed)
- **alpha** -- lambda ramp rate per timestep (positive float)
- **lambda0** -- initial lambda value (float, default 0.0)
- **nonuniform** -- `yes` or `no`; enable staggered initial lambda (optional,
  default `no`)

## Description

`fix backmap` drives the time-dependent backmapping simulation. It manages:

1. **Lambda ramp** -- per-atom resolution parameter that increases by `alpha`
   each timestep, clamped to \[0, 1\]
2. **Molecule map** -- builds and maintains the mapping between CG beads and
   their constituent AT atoms based on molecule IDs
3. **COM tracking** -- after each timestep, updates CG bead positions to the
   center-of-mass of their AT atoms
4. **CG force distribution** -- redistributes forces on CG beads to AT atoms
   proportional to their mass fraction
5. **CG velocity zeroing** -- sets CG bead velocities to zero each step

The per-atom lambda value is stored as a per-atom vector and accessible via
`f_ID` in dump commands (e.g., `f_bm` if the fix ID is `bm`).

## Arguments

### `cg_type` (required)

The numeric atom type identifying CG beads. All atoms with this type are
treated as CG beads; all other atoms in the same molecule are AT atoms.

```
fix bm all backmap cg_type 1 alpha 0.001
```

### `alpha` (required)

Lambda increment per timestep. At each step:

\[
\lambda \leftarrow \min(\lambda + \alpha,\ 1)
\]

Smaller values produce a slower, gentler transition. Typical range:
10<sup>-4</sup> to 10<sup>-3</sup>.

### `lambda0` (optional)

Initial value of lambda for all atoms. Default: `0.0`.

### `nonuniform` (optional)

When `yes`, each atom receives a random offset to its initial lambda, creating
a staggered transition. Default: `no`.

## fix_modify Options

### `active`

Activate or deactivate the lambda ramp:

```
fix_modify bm active no    # freeze lambda (CG equilibration)
fix_modify bm active yes   # resume ramp (backmapping phase)
```

When `active no`, lambda values remain frozen at their current values. CG
force distribution and COM tracking continue operating.

## Per-Atom Data

The fix stores one per-atom value: the current lambda. Access it via:

- `f_ID` in dump commands: `dump 1 all custom 100 dump.dat id type f_bm`
- `f_ID` in thermo output
- `extract("lambda", dim)` from C++ code

## Restart

Per-atom lambda values are written to restart files and restored on restart,
allowing seamless continuation of a backmapping simulation.

## Example

```
# Define atom groups
group at_atoms type 3 4
group cg_atoms type 1 2

# Set up backmapping fix
fix bm all backmap cg_type 1 alpha 0.0001 lambda0 0.0 nonuniform no

# Phase 1: CG equilibration (lambda frozen)
fix_modify bm active no
run 10000

# Phase 2: Backmapping (lambda ramps)
fix_modify bm active yes
run 10000

# Phase 3: AT production (lambda = 1 everywhere)
run 10000
```

## Mass Validation

During initialization, the fix checks that the CG bead mass equals the sum
of AT atom masses within each molecule. A warning is issued if they differ
by more than 10<sup>-4</sup>.

## Related

- [pair_style backmap](pair-backmap.md) -- lambda-weighted non-bonded
  interactions
- [bond_style backmap/harmonic](bond-styles.md) -- lambda-weighted harmonic
  bonds
- [bond_style backmap/table](bond-styles.md#backmap-table) -- lambda-weighted
  tabulated bonds
- [angle_style backmap/harmonic](angle-styles.md) -- lambda-weighted
  harmonic angles
