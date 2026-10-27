# fix backmap/ringguard

## Syntax

```
fix ID group backmap/ringguard K rc rcb file
```

- **K** -- wall constant (energy/distance^2 units)
- **rc** -- wall radius around each ring centroid (distance units)
- **rcb** -- smaller wall radius for the atoms bonded to the ring (distance units, 0 to rc)
- **file** -- one ring per line: the atom IDs of the ring (3 to 12 atoms), optionally followed
  by `|` and the IDs of the atoms bonded to the ring; `#` starts a comment

## Description

Puts a repulsive wall around the centroid of every listed ring. An atom of the fix group that is
not in the ring and lies at distance \( d < r_c \) from the centroid (the unweighted mean of the
ring atoms) adds

\[
E = K (r_c - d)^2
\]

An atom listed after the `|` of the ring adds \( K (r_{cb} - d)^2 \) for \( d < r_{cb} \) instead.

and is pushed out along the centroid-atom axis. The ring atoms share the reaction force equally,
so the total force is zero. The wall does not depend on lambda.

`backmap-prep` writes the ring list (`rings.txt`) and the fix (`ring_guard.lmp`, on `at_atoms`)
when `simulation.ring_guard` is set (robust protocol). The fix is on from the first minimization
until the hybrid structure is written, then `ring_guard_release.lmp` removes it. The rings are all
cycles of 5 to 7 atomistic atoms (the shortest cycle through each bond).

The defaults come from the epoxy example. 0.20 nm (`ring_guard`) is a suitable radius: an atom
passing through a six-membered ring comes within 0.14 nm of its centroid, while stacked rings stay
about 0.35 nm apart. The atoms bonded to a ring (its hydrogens and the first atom of each
substituent) are 0.145 nm or more from its centroid in a normal ring (axial hydrogens of
cyclohexane about 0.18 nm), so they get the smaller radius `ring_guard_bonded` (default 0.12 nm).
A substituent folded into its own ring sits at 0.07-0.11 nm and is pushed out. The default
`ring_guard_k` is 418400 kJ/mol/nm\(^2\) (1000 kcal/mol/Å\(^2\)). A ten times softer wall let atoms
through during the heat spike at the start of the ramp.

### Why

At the start of the lambda ramp the bonded terms act but the atomistic nonbonded repulsion is
still scaled down, so a substituent can fold into its own ring and a bond ends up passing through
it (ring spearing). Once the full potential is on, the defect is trapped. In the epoxy example the
CH2-NH2 arm of isophorone diamine does this for a few of the amine nitrogens, which shows up as
small peaks at short distance in the C-N radial distribution functions. The phantom particles at
ring centers of Krajniak, Zhang et al. (*J. Comput. Chem.* 2018) solve the related problem of
rings split over several CG beads, at the CG level. This fix acts on the atoms.

### Parallel evaluation

Each ring is evaluated once, by the rank that owns its first atom. The atoms near the centroid are
found with a cell grid over owned and ghost atoms rather than a neighbor list, because neighbor
lists leave out special pairs with weight 0, and an atom two bonds from a ring atom is such a
pair. Forces are accumulated on owned and ghost atoms and summed onto the owners with a reverse
communication of the fix, so forces and energies do not depend on the decomposition.

## Output

- `f_ID` -- wall energy (extensive)
- `f_ID[1]` -- wall energy
- `f_ID[2]` -- number of atom-ring contacts inside rc, a monitor of threading attempts

`fix_modify ID energy yes` adds the energy to the potential energy, which minimization needs.
`ring_guard.lmp` sets it and also writes `f_ID[1]` and `f_ID[2]` to `ring_guard.log` every 1000
steps.

## Restrictions

Needs an atom map (molecular systems have one). The rank that owns the first ring atom must see
the other ring atoms and all atoms within rc of the centroid as owned or ghost atoms. The virial
is not tallied, so the fix should not act under constant pressure (the robust protocol releases
it before any constant-pressure stage).
