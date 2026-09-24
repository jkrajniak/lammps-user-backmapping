# Design: one hybrid engine

## Pipeline

```
settings.yaml ─► Settings (schema) ─► bakery XML (v2_loader)
                                          │
              LAMMPS-format inputs ──► converted to .gro/.top (new, exact units)
                                          │
                                          ▼
                              vendored bakery build
                              (fragment COM on bead, degree-aware fragments,
                               cross bonds/angles/dihedrals, exclusions)
                                          │  hyb_conf.gro + hyb_topol.top
                                          ▼
                              build_system_from_hybrid (lammps_builder)
                                          │  System
                                          ▼
                              writers: .data, in.<prefix>, tables
```

`rebuild` repositions the hybrid from an equilibrated CG frame
(`network/rebuild.py`); `cg-only` uses `build_system_from_cg`. Both already exist
for networks and become the only paths.

## Decisions

1. **Keep bakery's placement rule as the single rule.** Fragment atoms are
   shifted by the fragment's mass-weighted COM and then by the bead position.
   This is the published method (JCTC 2016, JCC 2018) and the one confirmed by
   the author. No alternative placement is kept.
2. **LAMMPS-format inputs are converted, not special-cased.** A LAMMPS CG data
   file or AT fragment (data + bounded input script) is written as temporary
   `.gro`/`.top` in GROMACS units before the bakery build. Values are printed at
   full precision so the round trip back to `real` units is exact to 1e-10
   relative. One code path means the fragment-placement and topology logic is
   the same for every input format.
3. **`prep.engine` is removed, not aliased.** Accepting and ignoring it keeps old
   settings files loading; a deprecation warning tells the user it no longer
   selects anything.
4. **Validation against GROMACS, not against either old engine.** Neither engine
   is a trustworthy reference: the linear one has the defects in the proposal,
   the network one had never been run on linear systems and writes LJ as zero
   for them. The reference is the GROMACS topology itself, evaluated by GROMACS:
   - `gmx grompp` + `gmx mdrun -rerun` on the AT-only part of `hyb_topol.top`
     (CG beads removed, cross terms at full strength), same coordinates, same
     cutoffs, no long-range corrections.
   - Compare bond, angle, dihedral and LJ energies with LAMMPS `run 0` of the
     generated input at lambda = 1 (`fix backmap lambda0 1.0`, ramp inactive).
   - Tolerance: 1e-5 relative per term (float formatting in `.gro` limits
     coordinates to 1e-3 nm; the comparison uses the same `.gro` coordinates on
     both sides).
   GROMACS 2023 is installed in user space on the compute VM
   (`~/sc/tools/gmx`, conda-forge).
5. **Placement check is a property test.** For every generated data file, every
   bead is within 0.01 A of the mass-weighted COM of its AT atoms under the bead
   map `fix backmap` uses (per molecule, tags sorted, contiguous blocks).
6. **CG bonded completeness is a property test.** The number of CG bonds, angles
   and dihedrals in the data file equals the number implied by the CG topology
   and `cross_interactions` with `cg_bonded: true`.

## Risks

- Bakery's fragment matching is degree-based. Plain `atoms:` beads map at any
  degree (`*`), which `v2_loader` already emits; the double-prefix bug hid this.
- Type numbering and style choices differ from the old linear output (e.g. all
  AT bonds in one `backmap/harmonic at` type; same-bead weight is 1, so this is
  equivalent). Hand-written inputs that hard-code type numbers must be
  regenerated, not patched.
- The paper-data parity tests (rim135, PET) must stay byte-identical: the merge
  must not change the network path's output for networks.
