# Add a ring guard against ring spearing during the ramp (`simulation.ring_guard`)

## Why

In the epoxy (RIM135) validation run, four of the 200 amine nitrogens and one oxygen end up inside
a six-membered ring of their own molecule: a bond passes through the ring (ring spearing). The
nitrogens are those of isophorone diamine, whose CH2-NH2 arm folds into its cyclohexane ring.
The defects are frozen once the atomistic potential is fully on (no change over the 1 ns window,
through the 1000 K anneal) and show up as small peaks at short distance in the C-N RDFs of the
paper's epoxy figure (`research/experiments/20261027_rim135-small-r-peaks.md`).

They are not placement errors. At insertion and through the lambda = 0 minimization and
relaxation the nitrogens sit 0.28-0.35 nm from the ring centroid, outside the ring. They move in
during the first part of the ramp (about 0.20 nm at lambda 0.09, 0.12-0.19 nm at lambda 0.19,
0.09 nm at lambda 1), when the bonded terms act but the atomistic nonbonded repulsion is still
scaled down by lambda.

The 2018 paper (Krajniak, Zhang et al., J. Comput. Chem. 39, 648) solved a related problem for the
melamine network with phantom particles at ring centers that repel the CG beads. That acts at the
CG level, between rings split over several beads. Here the rings are single beads and the defect
is atomistic and intramolecular, so the phantom has to act on atoms.

## What changes

- New LAMMPS style `fix backmap/ringguard`:
  `fix ID group backmap/ringguard K rc rcb file`. Each line of `file` lists the atom IDs of one
  ring, then `|` and the atoms bonded to the ring. The fix adds
  E = K (rc - d)^2 for d < rc,
  where d is the distance of an atom of the fix group, not in the ring, from the ring centroid (an
  unweighted mean of the ring atoms, a massless virtual site). An atom bonded to the ring adds
  K (rcb - d)^2 for d < rcb instead: its normal place is 0.145 nm or more from the centroid (axial
  hydrogens about 0.18 nm), and a substituent folded into its own ring sits at 0.07-0.11 nm. The
  reaction force is shared equally by the ring atoms, so total force is conserved. The wall is
  independent of lambda.
  - Bond walls: a stretched bond can sweep through a ring with both atoms outside rc (seen in the
    prototype). Each bond with no atom in the ring adds E = w(s) - w(d_near), w(x) = K (rw - x)^2
    for x < rw, when the point of the bond nearest the centroid (distance s) lies between its
    atoms. d_near is the distance of the nearer bond atom, rw = rcb if a bond atom is bonded to
    the ring, else rc. The term is zero when the nearest point is an end, so the energy is
    continuous, and it adds only what the atom walls miss.
  - Parallel: atom walls are evaluated by the owner of the atom, bond walls by the owner of the
    atom storing the bond, against a cell grid of ring centroids, not a neighbor list, which
    drops special pairs with weight 0 (an atom two bonds from a ring atom is one). Each rank
    builds a ring around the owned image of its first atom, else the ghost image nearest its
    subdomain, and skips it unless all ring atoms are minimum images of that anchor (a ghost at
    the edge of the shell can see the other ring atoms only across the box, which placed phantom
    centroids in the prototype). Forces are accumulated on owned and ghost atoms in a per-atom
    buffer and summed onto the owners with a reverse communication of the fix, so they do not
    depend on the decomposition (the failure mode of `fix restrain`, see `fix backmap/stereo`).
    Orthogonal boxes only.
  - Output: `f_ID` energy (with `fix_modify energy yes` part of pe, needed for minimization),
    `f_ID[1]` energy, `f_ID[2]` number of atom and bond contacts inside the wall.
  - Virial: not tallied in the prototype. The guard acts only in the backmapping stages, which
    run at constant volume; documented as a restriction.
- New `backmap-prep` options (robust protocol, default off, output unchanged when off):
  - `simulation.ring_guard`: wall radius rc in nm (proposal 0.20 nm). Physical contacts stay
    outside: stacked rings are about 0.35 nm apart and a C-H pointing at a ring puts the H about
    0.25 nm from the centroid, while an atom passing through a six-membered ring comes within
    0.14 nm.
  - `simulation.ring_guard_k`: K in kJ/mol/nm^2 (default 418400, which is 1000 kcal/mol/A^2;
    100 kcal/mol/A^2 let atoms through, see the prototype results).
  - `simulation.ring_guard_bonded`: rcb in nm for the atoms bonded to a ring (default 0.12, at
    most `ring_guard`).
  - Rings: all cycles of 5 to 7 atomistic atoms found from the bond graph (shortest cycle through
    each bond). Smaller rings (epoxide, cyclopropane) are left out: their centroid is closer than
    rc to their own substituents.
  - Files `rings.txt`, `ring_guard.lmp` (fix on `at_atoms` + `fix_modify energy yes`),
    `ring_guard_release.lmp` (`unfix`). The guard is switched on before the first minimization
    and released where the stereo restraint is released (after the staged NVT, when the hybrid
    structure is written).
- Documentation: component page `docs/components/fix-backmap-ringguard.md`, settings rows,
  changelog entry.

## Validation (prototype)

1. Unit and LAMMPS integration tests (below).
2. Epoxy example, backmapping stage only (minimization, relaxation, ramp, staged NVT), with and
   without the guard, same seeds: count of heteroatoms within 0.15 nm of a ring centroid at the
   end of the ramp (6 without the guard in the rc5 run), `f_ID[2]` along the ramp, ramp energy and
   temperature transients (`b3_ramp_energy.py`).
3. If (2) removes the defects: rerun the epoxy Tier C protocol and compare the RDFs with the
   published runs (the small-r C-N features should vanish, nothing else should change beyond
   noise). Check PET (benzene rings) and melamine network (triazine) for threaded rings with and
   without the guard.

## Prototype results (epoxy, `research/experiments/20261027_ring-guard-prototype.md`)

- No guard: 10 threaded rings at the hybrid write (6 in an earlier run).
- rc 0.20 nm, K 100 kcal/mol/A^2, all atoms: 3 threaded rings. Two N and one O passed the wall
  during the heat spike at the start of the ramp (lambda 0.09-0.29). Most wall contacts were the
  rings' own axial hydrogens, which the wall pushed outward.
- rc 0.20 nm, K 1000, atoms bonded to the ring left out: no atom two or more bonds from a ring
  passed. One defect remained, the ring amine N of an isophorone folded into its own ring
  (0.085 nm), which the exclusion let through. Hence the radius rcb for bonded atoms.
- rc 0.20 nm, K 1000, rcb 0.12 nm (the defaults): 1 threaded ring. No atom passed the walls; a
  stretched bond (0.45-0.60 nm at lambda = 0) swept through a ring with both atoms outside rc,
  was held at the wall and fell in after the release. Hence the bond walls.
- Bond walls added: the run heated to 950 K in the NVT stage and crashed. Cause: phantom ring
  centroids in parallel (above), present since the first prototype, so the guard contact counts
  of all earlier guarded runs are inflated (the threading counts come from the dumps and stand).
  After the fix, energies equal an independent numpy recomputation on 1-8 ranks.
- Defaults with bond walls and the anchor fix: running.

## Out of scope

- CG-level phantoms between beads (the 2018 remedy); no current example splits a ring over beads.
- Repair of defects in an existing backmapped structure.
- Virial contribution (needed only if the guard is kept on under constant pressure).

## Tests

- `python/tests/test_ring_guard.py`: ring detection (benzene, cyclohexane, fused naphthalene
  gives two rings, an epoxide and a chain give none, CG atoms ignored), the written files and
  units, default off.
- `python/tests/test_lammps_ring_guard.py` (integration, `BACKMAP_LMP`): zero energy and force
  for an atom outside rc; energy K(rc - d)^2 and force along the centroid-atom axis for an atom
  inside; total force zero; ring members never repel each other; forces identical for serial and
  2-rank runs with the ring split over the ranks; the contact count; bonded atoms inside and
  outside rcb; a bond through the ring with both atoms outside, a bond beside it, the bond and
  atom forces against finite differences, bond walls serial vs 2 ranks; a ring swept across the
  periodic boundary in a lattice of probe atoms with a ghost shell wider than half the box,
  serial vs 2 and 4 ranks.
