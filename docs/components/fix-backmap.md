# fix backmap

## Syntax

```
fix ID group-ID backmap cg_type T1 [T2 ...] alpha A lambda0 L0 [nonuniform yes/no]
```

- **ID** -- fix identifier (user-chosen)
- **group-ID** -- group of atoms the fix applies to (typically `all`)
- **cg_type** -- one or more atom type numbers of CG beads (integers,
  1-indexed). Multiple types are supported for systems with heterogeneous CG
  bead types (e.g. end and middle beads).
- **alpha** -- lambda ramp rate per timestep (positive float)
- **lambda0** -- initial lambda value (float, default 0.0)
- **nonuniform** -- `yes` or `no`; enable staggered initial lambda (optional,
  default `no`)

## Description

`fix backmap` drives the time-dependent backmapping simulation. It manages:

1. **Lambda ramp** -- per-atom resolution parameter that increases by `alpha`
   each timestep, clamped to \[0, 1\]
2. **Bead map** -- builds and maintains the per-bead mapping between CG beads
   and their constituent AT atoms based on molecule IDs and global tag order.
   Rebuilt automatically on every neighbor list rebuild.
3. **COM tracking** -- after each timestep, updates CG bead positions to the
   center-of-mass of their AT atoms (using PBC-aware displacement)
4. **CG force distribution** -- redistributes forces on CG beads to AT atoms
   proportional to their mass fraction
5. **CG velocity zeroing** -- sets CG bead velocities to zero each step

The per-atom lambda value is stored as a per-atom vector and accessible via
`f_ID` in dump commands (e.g., `f_bm` if the fix ID is `bm`).

## Arguments

### `cg_type` (required)

One or more numeric atom types identifying CG beads. All atoms matching any
listed type are treated as CG beads; all other atoms in the same molecule are
AT atoms. Within each molecule, AT atoms are distributed evenly among CG beads
in global tag order: CG bead *i* receives AT atoms
\[i \times \text{apb},\ (i+1) \times \text{apb})\] where
`apb = num_AT / num_CG`.

```
fix bm all backmap cg_type 1 2 alpha 0.001
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

Typical **single-segment** λ ramp (CG melt equilibrated before building the
hybrid; `backmap-prep` default):

```
# Define atom groups
group at_atoms type 3 4
group cg_atoms type 1 2

fix integrate at_atoms nvt temp 298.0 298.0 100.0
fix bm all backmap cg_type 1 2 alpha 0.0001 lambda0 0.0 nonuniform no
compute at_temp at_atoms temp
thermo_modify temp at_temp

fix_modify bm active yes
run 10000
write_data system_hybrid.data
```

Optional **multi-segment** input (only if you need in-hybrid CG equilibration
or post-backmap steps in the same file):

```
fix_modify bm active no
run 10000
fix_modify bm active yes
run 10000
write_data system_hybrid.data
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
