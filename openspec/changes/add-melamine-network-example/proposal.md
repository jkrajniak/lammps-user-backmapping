# Add crosslinked melamine (MF) network example

## Why

`examples/melamine/` (Tier C RDF validation, 5/11 PASS best result so far) simulates
500 fully independent trimethylolmelamine monomers. Investigation this session
(`research/experiments/20260802_melamine-bakery-rerun.md`, 2026-08-06 update) found
that bakery's actual reference RDFs come from a **substantially different chemical
system**: the same 500 molecules, but with 675 real inter-molecular covalent
crosslinks (`bakery_full_unused/examples/network_backmapping/mf/backmapping/cg_topol.top`,
`; chem` bonds, e.g. `1 622` linking molecule 1 to a molecule ~200 residues away —
roughly 45% of each molecule's 3 pendant arms condensed into an ether bridge with
another molecule). Comparing our uncrosslinked monomer gas against that network
reference is very likely the dominant reason L2(g) has failed on every RDF pair,
every protocol variant, across the entire multi-week melamine campaign, independent
of any force-field tuning (three separate mechanisms -- Langevin damping, `fix bm`
residual force, ring-charge symmetrization -- were tested this session and none
moved that result).

`openspec/changes/phase3-network-backmapping/design.md` already lists **"2. MF
polymerized -- `network_backmapping/mf/` (depends on melamine base)"** as the
second item in its planned migration order (right after epoxy/RIM135, which is
done), and its `tasks.md` item 6.1 ("MF polymerized network") is an unchecked,
explicitly-flagged-as-separate-change follow-on. This change is that follow-on.

## What changes

- Keep `examples/melamine/` (uncrosslinked) exactly as-is -- no changes to that
  example, its `settings.yaml`, or its `openspec/specs/example-melamine/spec.md`.
  It stays the paper's local-structure/force-field-fidelity reference.
- Add a new `examples/melamine_network/` example: the same 500-molecule MF system,
  crosslinked to match bakery's reference network exactly (same 675-bond pattern),
  built via the already-proven Settings v2 / `network-backmap-prep` engine (same
  machinery used for epoxy/RIM135 and PET -- `atoms_by_degree` bead mapping, active
  sites, `network/v2_loader.py` + `network/lammps_builder.py`), not new bespoke code.
- Source real OPLS-AA parameters for the new C-O-C angle and C-C-O-C dihedral terms
  that appear at each crosslink site (bakery's own reference never filled these in --
  `missing_definitions.txt` / `; chem MISSING params type: A-A` throughout their
  `cg_topol.top` and `hyb_topol.top`), added to the shared `forcefield_dir` already
  used to fill melamine's missing amide/imide dihedrals.
- Run the full pipeline end-to-end on the VM: `backmap-prep build-hybrid` ->
  LAMMPS eq/ramp/production -> RDF comparison against the same
  `ref_C_C.xvg`/`ref_C_N.xvg`/`ref_O_H.xvg` files (now a fair, like-for-like
  comparison for the first time).

## Success criterion

A stable eq/ramp/production run on the real crosslinked topology, with correct
local chemistry at every site (crosslinked and unreacted): bond lengths/angles,
charges, and crosslink connectivity all matching OPLS-AA / bakery's reference.
RDF results are then computed and reported honestly, whatever they show -- this
change does not commit to matching or beating the uncrosslinked example's RDF
pass count as a target. (Options were discussed and this was the user's explicit
choice: get a structurally faithful run first, not force-field-tune toward a
number.)

## Non-goals

- Live/reactive bond formation during the simulation (matching
  `phase3-network-backmapping`'s own non-goal: static cured topology only). The
  675-bond network is imported as a fixed, pre-existing structure, exactly as
  bakery's own `cg_topol.top` already defines it -- we are not re-deriving or
  re-simulating the curing process itself.
- Standing up ESPResSo++ to reuse bakery's own backmapping engine (checked:
  not built/importable on the VM; would require a new C++/MPI/Python build for a
  one-off structure-generation step). The crosslinked AT structure is built by our
  own `network-backmap-prep` engine instead.
- Changing anything about `examples/melamine/` (the uncrosslinked example) or its
  spec.
- Chasing a specific RDF pass-count target (see Success criterion).
- General-purpose reactive-crosslink-generation tooling (we import bakery's one
  specific network; we are not building a "generate an arbitrary degree of cure"
  feature). If a different degree of cure is ever wanted later, that is a
  separate, future change.

## Impact

- **Python:** likely a small extension to `network/v2_loader.py` if MF's bakery
  source layout (single `at_topol.top` + `settings.xml` active-site/remove rules)
  doesn't already fit the input shape the loader expects from epoxy (separate
  pre-built per-degree `.itp`/`.gro` files) -- see design.md Risk section. Possibly
  a one-time conversion script to translate bakery's `cg_topol.top` crosslink list
  into the v2 settings format.
- **Examples:** new `examples/melamine_network/` directory.
- **Forcefield:** new OPLS-AA ether angle/dihedral entries added to the shared
  `forcefield_dir`.
- **Research:** new experiment note once run; `research/checkpoints.md` gets a new
  entry once a run is validated (not before).
- **Docs:** new example README; no changes to the existing melamine README.

## Dependencies

- `openspec/changes/phase3-network-backmapping` -- uses the `network-backmap-prep`
  capability (Settings v2, `atoms_by_degree`, active sites) that change
  implements. That change is still open (19/29 tasks); this one only consumes
  its already-completed, already-proven-on-epoxy-and-PET generator surface, and
  does not block on its remaining open tasks (charge-transfer-rules task 2.4 is
  not required here -- see design.md).
- Bakery reference data (read-only, VM-local):
  `bakery_full_unused/examples/network_backmapping/mf/backmapping/` (`cg_topol.top`,
  `cg_conf.gro`, `at_topol.top`, `settings.xml`, `hyb_topol.top`).
- Session finding this proposal is based on:
  `research/experiments/20260802_melamine-bakery-rerun.md` (2026-08-06 update).

## References

- Same JCC paper as the uncrosslinked example (trimethylolmelamine polymerization
  system), referenced in `phase3-network-backmapping/proposal.md`.
- `phase3-network-backmapping/design.md` migration order, item 2.
