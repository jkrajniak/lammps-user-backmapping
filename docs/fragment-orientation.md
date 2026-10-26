# Fragment orientation

`backmap-prep` places one atomistic fragment on every CG bead. By default a fragment is only
*translated* onto its bead: it keeps the orientation it had in the source molecule. That is
enough when each fragment contains the geometry that matters, but it leaves the **relative
orientation of neighboring fragments** to chance, which matters wherever a stereocenter or a
rigid unit spans two beads.

`prep.fragment_orientation: align` rotates every fragment so that it fits its neighbors.

```yaml
prep:
  fragment_orientation: align   # default: translate
```

![Glycerol C2 of POPC placed by translation and by alignment](images/fragment-orientation.png)

*(a) The glycerol carbon C2 of a MARTINI 3 POPC lipid, seen along the C1-C2 bond. Three of its
four neighbors (C1, O21, HS) are in one fragment and the fourth (C3) in the next. Translation puts
C3 on the wrong side for this lipid; alignment puts it on the right side. (b) Signed height of C3
over the plane through C2, C1, and O21 for the 128 lipids of the bilayer, before any dynamics.
(c) Share of lipids with the handedness of the source molecule. The independent CG2AT bilayer
used as reference in the paper has 89 %.*

## What it does

For every residue, the builder

1. compares the vectors between the centers of mass of bonded fragments in the source molecule
   with the vectors between the CG beads, and takes the proper rotation (Kabsch, determinant +1,
   so a fragment is never mirrored) that makes them agree;
2. refines each rotation with the covalent bonds that join fragments: an atom bonded to a
   neighboring fragment should lie where that fragment, rotated, puts it;
3. checks every stereocenter whose neighbors lie in two fragments. If the placed geometry has
   the opposite handedness to the source molecule, the fragment gets the best-fitting
   orientation, among 2000 trial rotations, that restores it.

Fragments are rotated rigidly about their bead, so the CG bead positions and the internal
geometry of every fragment are unchanged. Beads without a bonded neighbor in the same residue
(for example a MARTINI water bead) are not rotated. Fragments cut from different source chains
are not compared.

## When to use it

Use `align` when a stereocenter, a ring, or a double bond has atoms in two beads, as in lipids
(glycerol, cis double bonds) and residues split across beads. For linear polymers whose fragments
are self-contained, `translate` gives the same structure after the ramp and keeps the output
identical to earlier releases, which is why it is the default.

## Keeping the handedness during the ramp: `simulation.stereo_restraint`

Placement is not enough. In the MARTINI POPC bilayer the fragments are placed with the right
handedness at the glycerol carbon, but a few percent of the lipids invert at the start of the ramp,
when the bonds between fragments switch on from a nearly flat arrangement, and the atomistic force
field has no term that brings them back. `simulation.stereo_restraint` (kJ/mol, robust protocol)
restrains every stereocenter of the atomistic part at the handedness it has in the hybrid
structure, from the start of the backmapping until the hybrid structure is written:

```yaml
simulation:
  protocol: robust
  stereo_restraint: 83.68   # kJ/mol (20 kcal/mol)
```

A stereocenter is an atom with four bonded atomistic neighbors that are all inequivalent (equal
atom type and equal environment count as equivalent). It is restrained on two improper
dihedrals, one of its three heavy neighbors and one with the fourth neighbor (the hydrogen when
there is one) in place of the third, each with the ideal tetrahedral value (35.26 degrees) and the
sign of its handedness in the hybrid structure. Each restraint has zero energy at that geometry
and 2K/3 at the mirror image. The heavy-atom improper alone keeps the center on its side of the
heavy neighbors but not the hydrogen: under a stiff restraint the hydrogen crosses to the side of
the heavy neighbors, the center becomes nearly planar, and it inverts as soon as the restraint is
released. The second improper rises when the hydrogen crosses. The input script includes `stereo.lmp` before the first minimization and
`stereo_release.lmp` after `<prefix>_hybrid.data` is written.

`fix restrain dihedral` has its energy minimum at the written target + 180 degrees, although its
documentation puts the minimum at the target (LAMMPS stable_22Jul2025 to develop of 2026-10-25).
`backmap-prep` writes the target shifted by 180 degrees, and the test
`python/tests/test_lammps_stereo_restraint.py` fails if LAMMPS changes this convention.

The numbers below were measured with release 1.4.0, whose restraint had its energy minimum away
from the tetrahedral geometry (fixed in the next release) and pushed every restrained center. The
restraint is soft on purpose. In the POPC bilayer (128 lipids) it reduces the inverted
centers at the end of the hybrid stage from 17 (`align` alone) to 6; a stiff restraint (about
800 kJ/mol) holds all 128 but strains the center and destabilizes the 1 fs stage, and a restraint
that decays over the stages releases the same lipids when it drops. It is off by default.

## Limits

- Handedness is taken from the **source molecule**; the CG model carries no chirality, so a
  molecule of the opposite handedness needs its own source structure.
- The orientation fixes the placement only. Whether a geometry survives the hybrid run depends
  on the transient of the ramp (see the ramp-length discussion in the paper).
- A stereocenter is detected from the source geometry (atoms closer than 0.19 nm are bonded) and
  only when its neighbors lie in two fragments of one residue.
