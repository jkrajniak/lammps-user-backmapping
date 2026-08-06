# Tasks — crosslinked melamine (MF) network example

Parent change: `phase3-network-backmapping` (tasks.md item 6.1)
Design: [design.md](./design.md)

---

## 1. Understand the existing generator surface against MF's specific shape

- [ ] 1.1 Read `network/v2_loader.py` mapping-resolution code end-to-end
      (`_resolve_mapping`, `_emit_cg_bead`, `_emit_bead_extras`) against how
      epoxy's `settings.v2.yaml` expresses `degree`/`base`/`add`/`active_sites`
- [ ] 1.2 Read bakery's vendored `network/bakery/structures.py`
      (`BackmapperSettings2`) to confirm it consumes MF's exact XML shape
      (single `at_topol.top` + `degree="2"/"3"` beads + `active_site`/`remove`)
      the same way it consumes epoxy's (separate per-degree source files) --
      or identify precisely what differs
- [ ] 1.3 Record findings in this file / a short note before writing YAML --
      confirm whether a small `v2_loader.py` extension is needed (see
      design.md Risks) or whether MF fits the existing shape as-is

## 2. Import bakery's crosslink network

- [ ] 2.1 Write a one-time conversion script: parse
      `cg_topol.top`'s `; chem` bonds (675 total) into a molecule-instance +
      bead-pair list; cross-check against `hyb_topol.top`'s "AT cross bonds"
      section and `settings.xml`'s `active_site`/`remove` rules for which atoms
      are removed on which side of each crosslink
- [ ] 2.2 Verify the conversion against known facts: 500 residues, 1500 static
      + 675 chem bonds, ~45% of arms crosslinked -- script output should
      reproduce these counts exactly
- [ ] 2.3 Emit the crosslink list in whatever input shape task 1.3 determines
      (v2 YAML `active_sites`/mapping entries directly, or a supplementary
      file consumed by a small loader extension)

## 3. Force field for crosslink sites

- [ ] 3.1 Source real OPLS-AA parameters for the C-O-C angle at each ether
      bridge and the C-C-O-C dihedral(s) around it (same category of lookup
      already used for melamine's missing amide/imide dihedrals via
      `forcefield_dir`)
- [ ] 3.2 Add these to the shared `forcefield_dir`; verify explicit inline
      params (if any survive per-molecule) still take priority over table
      lookup, matching the existing fill-in behavior
- [ ] 3.3 Unit test: after generation, confirm zero remaining
      `missing_definitions.txt` entries for angle/dihedral types actually
      present in the generated system

## 4. `examples/melamine_network/` example

- [ ] 4.1 Directory layout mirroring `examples/epoxy/` (`settings.v2.yaml`,
      `README.md`, `large/` for full-scale generated files) -- assets resolve
      via `prep.data_dir` pointing at bakery's MF source, same pattern as
      epoxy's `settings.v2.yaml`
- [ ] 4.2 `settings.v2.yaml`: MF molecule, A1/A2/A3 beads, degree-1 (unreacted,
      9-atom, matches `examples/melamine/`'s existing bead definition) +
      degree-2 (reacted) mapping, active sites from task 2
- [ ] 4.3 Tier A: `backmap-prep build-hybrid` parity check against bakery's own
      `hyb_conf.gro`/`hyb_topol.top` (atom/bond/angle/dihedral counts; charge
      sum conserved within 1e-6 e)
- [ ] 4.4 Wire `backmap-prep build` -> LAMMPS `.data`/`in.*` (reuse
      `network/lammps_builder.py::build_system_from_hybrid`, already proven
      for rim135/PET)
- [ ] 4.5 New crosslink AT bonds/angles/dihedrals written as plain always-on
      `harmonic` terms (not lambda-weighted `backmap/harmonic`) -- verify in
      the generated `in.*` file, not just assumed

## 5. LAMMPS run (VM, tmux)

- [ ] 5.1 Short pilot (eq + short ramp + short production) first, matching
      this session's own pattern for the uncrosslinked example -- check
      stability (dangerous rebuilds, thermo sanity) before committing to a
      full-length run
- [ ] 5.2 Full eq/ramp/production run once the pilot is stable
- [ ] 5.3 `compare_melamine_structure.py` against the same `ref/*.xvg` files
      used by `examples/melamine/`

## 6. Validation and docs

- [ ] 6.1 `examples/melamine_network/README.md` -- system description,
      crosslink provenance (bakery's exact network, cited), run instructions
- [ ] 6.2 New `research/experiments/` note recording the run (protocol,
      stability, RDF result, honestly reported per design.md Success
      criterion)
- [ ] 6.3 Update `research/checkpoints.md` with a new entry once a run is
      validated (not before, per that file's own registry practice)
- [ ] 6.4 Check off `phase3-network-backmapping/tasks.md` item 6.1, linking to
      this change
- [ ] 6.5 Add `openspec/specs/example-melamine-network/spec.md` (this change's
      own spec delta, once implementation matches it) -- do not modify
      `openspec/specs/example-melamine/spec.md`
