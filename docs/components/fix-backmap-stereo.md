# fix backmap/stereo

## Syntax

```
fix ID group backmap/stereo K file
```

- **K** -- restraint constant (energy units)
- **file** -- one restraint per line: `atom1 atom2 atom3 atom4 phi_t`, atom IDs and the
  target dihedral in degrees with the IUPAC sign; `#` starts a comment

## Description

Restrains improper dihedrals with

\[
E = K \left[ 1 - \cos(\phi - \phi_t) \right]
\]

which is zero at \( \phi = \phi_t \) and \( 2K/3 \) at the mirror image of an ideal
tetrahedral center (\( \phi_t = \pm 35.26^\circ \)). `backmap-prep` writes the file
(`stereo_centers.txt`) and the fix (`stereo.lmp`) when `simulation.stereo_restraint` is set,
with two restraints per stereocenter ([Fragment orientation](../fragment-orientation.md)).

`fix restrain dihedral` has the same functional form, with the target shifted by 180 degrees,
but it adds forces to ghost atoms in `post_force`, after the reverse communication of ghost
forces. With `newton_bond on` (the default) a restraint whose atoms lie on two MPI ranks loses
part of its force, and with `newton_bond off` its energy is counted more than once. Here every
rank that owns one of the four atoms computes the restraint, applies forces to its own atoms
only, and tallies its share of the energy. Forces and energies do not depend on the
decomposition.

## Output

- `f_ID` -- total restraint energy (extensive)
- `f_ID[m]` -- energy of restraint `m`, in file order

`fix_modify ID energy yes` adds the energy to the potential energy, which minimization needs.
`stereo.lmp` sets it.

## Restrictions

Needs an atom map (molecular systems have one). Every rank that owns one of the four atoms
must also see the other three as owned or ghost atoms; otherwise the fix stops with an error.
