# Proposal: MPI correctness for fix backmap

## Why

`fix backmap` is validated serially but produces silently wrong results under
MPI domain decomposition. Two operations are incorrect:

1. **COM position update (`initial_integrate`)** sums each CG bead's mapped AT
   atoms from the local+ghost set on the CG-owner rank. AT atoms that live on
   non-neighbour ranks (beyond the ghost cutoff) are missed, so the COM is
   incomplete and the CG bead is driven to the wrong position.
2. **CG force distribution (`post_force`)** adds `ratio * f_cg` to AT atoms in
   the local+ghost set, including ghost AT atoms. Forces written to ghost atoms
   are never reverse-communicated to their owners, so those AT atoms miss the
   CG-force contribution. AT atoms whose CG bead is on another rank and not
   ghosted locally get nothing.

This blocks the CPC manuscript claim that the package runs in parallel, and is
the last open correctness item in `STATUS.md`.

## What changes

- Replace the CG-owner-centric loops with an atom-centric, `fix_rigid_small`-style
  communication pattern:
  - **COM**: each *local* AT atom contributes `(m, m*dx)` to a per-CG-bead
    accumulator indexed by the CG atom's local-or-ghost index; a
    `reverse_comm(4)` consolidates ghost-CG contributions to CG owners, which
    then apply the COM shift.
  - **Force**: forward-comm CG `(fx, fy, fz, mass)` plus `lambda` to ghost CG
    atoms (`comm_forward = 5`); each *local* AT atom reads its (local or ghost)
    CG bead's force and adds its share to its own `f`. All force additions land
    on local atoms, so no reverse-comm is needed.
- A per-atom `atom2cg` map (local-or-ghost CG index per AT atom) is built in
  `build_bead_map` and grown with `nmax`.
- The bead-partner atoms must still be local-or-ghost on each rank, i.e.
  `comm_modify cutoff` must be >= the maximum CG-AT distance within a bead
  (a few Å for backmapped fragments). This is the same constraint as `fix rigid`
  for bodies that span domains and is documented.

## Out of scope

- `pair_style backmap` and `bond/angle/dihedral_style backmap/*` MPI behaviour is
  unchanged (they use standard LAMMPS ghost communication and are not affected).
- Performance optimisation of the new communication (profile first).

## Evidence

- Serial regression: existing dodecane/PE/melamine Tier B runs unchanged.
- MPI test: small dodecane, serial vs 4-rank, compare final positions / energy
  within tolerance (`examples/dodecane/large/test_mpi_serial_vs_4rank.sh`).
