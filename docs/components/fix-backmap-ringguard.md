# fix backmap/ringguard

## Syntax

```
fix ID group backmap/ringguard K rc rcb file
```

- **K** -- force constant (energy/distance^2 units)
- **rc** -- radius of the repulsive sphere around each ring centroid (distance units)
- **rcb** -- smaller radius for the atoms bonded to the ring (distance units, 0 to rc)
- **file** -- one ring per line: the atom IDs of the ring (3 to 12 atoms), optionally followed
  by `|` and the IDs of the atoms bonded to the ring; `#` starts a comment

## Description

Centres a repulsive sphere of radius \( r_c \) on the centroid of every listed ring. An atom of the fix group that is
not in the ring and lies at distance \( d < r_c \) from the centroid (the unweighted mean of the
ring atoms) adds

\[
E = K (r_c - d)^2
\]

and is pushed out along the centroid-atom axis. An atom listed after the `|` of the ring adds
\( K (r_{cb} - d)^2 \) for \( d < r_{cb} \) instead.

A stretched bond can also sweep through a ring while both of its atoms stay outside the sphere. Each
bond of the fix group with no atom in the ring therefore adds

\[
E = w(s) - w(d_\mathrm{near}), \qquad w(x) = K (r_w - x)^2 \text{ for } x < r_w
\]

when the point of the bond nearest the centroid lies between its two atoms, at distance \( s \).
Here \( d_\mathrm{near} \) is the distance of the nearer bond atom, and \( r_w = r_{cb} \) if a
bond atom is bonded to the ring, otherwise \( r_c \). The term keeps the bond out of the same sphere, adds only what the atom terms do
not already cover and vanishes when the nearest point reaches an end of the bond, so the energy is
continuous. The ring atoms share every reaction force equally, so the total force is zero. The
forces do not depend on lambda. The excluded region is a sphere, not a cylinder along the ring
axis, so atoms above the ring face (a stacked ring at about 0.35 nm) are not affected.

`backmap-prep` writes the ring list (`rings.txt`) and the fix (`ring_guard.lmp`, on `at_atoms`)
in the robust protocol whenever the system has rings: `simulation.ring_guard` defaults to 0.20 nm,
and `ring_guard: null` switches the guard off. The fix is on from the first minimization
until the hybrid structure is written, then `ring_guard_release.lmp` removes it. The rings are all
cycles of 5 to 7 atomistic atoms (the shortest cycle through each bond).

The defaults come from the epoxy example. 0.20 nm (`ring_guard`) is a suitable radius: an atom
passing through a six-membered ring comes within 0.14 nm of its centroid, while stacked rings stay
about 0.35 nm apart. The atoms bonded to a ring (its hydrogens and the first atom of each
substituent) are 0.145 nm or more from its centroid in a normal ring (axial hydrogens of
cyclohexane about 0.18 nm), so they get the smaller radius `ring_guard_bonded` (default 0.12 nm).
A substituent folded into its own ring sits at 0.07-0.11 nm and is pushed out. The default
`ring_guard_k` is 418400 kJ/mol/nm\(^2\) (1000 kcal/mol/Å\(^2\)). A ten times smaller K let atoms
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

Each atom term is evaluated by the rank that owns the atom, each bond term by the rank that owns the
atom storing the bond. The rings near an atom or a bond are found with a cell grid of ring
centroids rather than a neighbor list, because neighbor lists leave out special pairs with weight
0, and an atom two bonds from a ring atom is such a pair. A rank builds a ring around the owned
image of its first atom, or else the ghost image nearest its subdomain, and skips the ring if the
other ring atoms are not all present as minimum images of that anchor. Without this check a ghost
at the edge of the ghost shell could take the other ring atoms from across the box and place the
centroid elsewhere. Forces are accumulated on owned and ghost atoms and summed onto the owners with a reverse
communication of the fix, so forces and energies do not depend on the decomposition.

## Output

- `f_ID` -- guard energy (extensive)
- `f_ID[1]` -- guard energy
- `f_ID[2]` -- number of atom-ring and bond-ring contacts inside the sphere

`fix_modify ID energy yes` adds the energy to the potential energy, which minimization needs.
`ring_guard.lmp` sets it and also writes `f_ID[1]` and `f_ID[2]` to `ring_guard.log` every 1000
steps.

## Restrictions

Needs an atom map (molecular systems have one) and an orthogonal box. The ghost cutoff must cover
a ring plus rc around every owned atom and bond, which the pair cutoff of any atomistic force
field does. The virial
is not tallied, so the fix should not act under constant pressure (the robust protocol releases
it before any constant-pressure stage).
